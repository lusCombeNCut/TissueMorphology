#ifndef TESTVISCOELASTICSHEAREXPERIMENT3D_HPP_
#define TESTVISCOELASTICSHEAREXPERIMENT3D_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <fstream>
#include <sstream>

#include "AbstractCellBasedTestSuite.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * 3D simple shear stress-relaxation test on a pure viscoelastic ghost node
 * network (no cells).
 *
 * Analogous to TestViscoelasticShearExperiment.hpp but using
 * ViscoelasticGhostNodeEcmField<3> with a cubic lattice (18-connected).
 * Shear is applied in the x-direction with height gradient in y.
 * The z-faces are free (unconfined in depth).
 *
 * Validates the Cauchy-Born unit conversion for the 3D cubic lattice:
 *   G = k/a  (Ostoja-Starzewski, Appl. Mech. Rev. 55(1), 2002, Eq. 2.10)
 *
 * Same environment variables as the 2D test:
 *   SHEAR_G_PA       — target shear modulus in Pa (default 800)
 *   SHEAR_ARM_RATIO  — E1/E0 SLS arm ratio (default 0.2)
 *   SHEAR_TAU        — relaxation time in hours (default 1.0)
 *   SHEAR_OUTPUT_DIR — output subdirectory name (default SimpleShear_3D)
 *   SHEAR_SKIP_VTP   — set to "1" to skip VTP output (faster for sweeps)
 *
 * Output dir: $CHASTE_TEST_OUTPUT/ViscoelasticCompression/<SHEAR_OUTPUT_DIR>/
 */
class TestViscoelasticShearExperiment3d : public AbstractCellBasedTestSuite
{
private:

    static double GetEnvDouble(const char* name, double defaultVal)
    {
        const char* val = std::getenv(name);
        return val ? std::atof(val) : defaultVal;
    }

    static std::string GetEnvString(const char* name, const std::string& defaultVal)
    {
        const char* val = std::getenv(name);
        return val ? std::string(val) : defaultVal;
    }

    static bool GetEnvBool(const char* name, bool defaultVal)
    {
        const char* val = std::getenv(name);
        if (!val) return defaultVal;
        return std::string(val) == "1" || std::string(val) == "true";
    }

    /**
     * Enforce shear boundary conditions for the 3D slab:
     *   - Bottom y-face: pin x, y, z to original positions
     *   - Top y-face:    prescribe x = original_x + shearOffset; pin y, z
     */
    void EnforceShearBoundaryConditions3d(
        ViscoelasticGhostNodeEcmField<3>& rField,
        const std::vector<unsigned>& rBottomNodes,
        const std::vector<unsigned>& rTopNodes,
        const std::vector<double>& rBottomOrigX,
        const std::vector<double>& rBottomOrigY,
        const std::vector<double>& rBottomOrigZ,
        const std::vector<double>& rTopOrigX,
        const std::vector<double>& rTopOrigY,
        const std::vector<double>& rTopOrigZ,
        double shearOffset)
    {
        for (unsigned k = 0; k < rBottomNodes.size(); k++)
        {
            unsigned idx = rBottomNodes[k];
            rField.rGetNode(idx).position[0] = rBottomOrigX[k];
            rField.rGetNode(idx).position[1] = rBottomOrigY[k];
            rField.rGetNode(idx).position[2] = rBottomOrigZ[k];
        }
        for (unsigned k = 0; k < rTopNodes.size(); k++)
        {
            unsigned idx = rTopNodes[k];
            rField.rGetNode(idx).position[0] = rTopOrigX[k] + shearOffset;
            rField.rGetNode(idx).position[1] = rTopOrigY[k];
            rField.rGetNode(idx).position[2] = rTopOrigZ[k];
        }
    }

public:

    void TestSimpleShearRelaxation3d()
    {
        // ── Read parameters (env or defaults) ─────────────────────
        const double G_input_Pa  = GetEnvDouble("SHEAR_G_PA", 800.0);
        const double arm_ratio   = GetEnvDouble("SHEAR_ARM_RATIO", 0.2);
        const double tau         = GetEnvDouble("SHEAR_TAU", 1.0);
        const std::string outDir = GetEnvString("SHEAR_OUTPUT_DIR", "SimpleShear_3D");
        const bool skipVtp       = GetEnvBool("SHEAR_SKIP_VTP", false);

        // ── Unit conversion constants ─────────────────────────────
        const double kF        = 0.001;     // N/m per sim unit
        const double spacing_m = 10e-6;     // 10 µm

        // ── Cauchy-Born conversion: G = k/a → E0 = G*a/kF ────────
        // Valid for the 18-connected cubic lattice (Ostoja-Starzewski 2002, Eq. 2.10)
        const double E0 = G_input_Pa * spacing_m / kF;
        const double E1 = arm_ratio * E0;

        // ── Grid geometry ─────────────────────────────────────────
        // x: shear direction (wide)
        // y: height / shear-gradient direction
        // z: depth (free faces)
        const double spacing = 1.0;
        const double x_min = 0.0, x_max = 49.0;   // 50 nodes wide
        const double y_min = 0.0, y_max = 9.0;    // 10 nodes tall (height)
        const double z_min = 0.0, z_max = 9.0;    // 10 nodes deep

        const unsigned num_cols_x = static_cast<unsigned>((x_max - x_min) / spacing) + 1; // 50
        const unsigned num_cols_z = static_cast<unsigned>((z_max - z_min) / spacing) + 1; // 10
        const double height    = y_max - y_min;                      // 9
        const double face_area = (x_max - x_min) * (z_max - z_min); // 49 * 9 = 441

        // ── Drag: set for fast diffusion, respecting CFL ──────────
        // Cubic lattice: N_neighbours = 18 (face + edge-diagonal; body-diagonal excluded)
        // CFL = (E0+E1)*18*dt/eta < 2; solve for eta given CFL_target
        const double dt         = 0.001;
        const double cfl_target = 0.8;
        const double eta_drag   = std::max(0.01, (E0 + E1) * 18.0 * dt / cfl_target);

        // ── Shear protocol ────────────────────────────────────────
        const double gamma_max      = 0.05;
        const double shear_disp_max = gamma_max * height;
        const double shear_velocity = shear_disp_max / 1.0;
        const double t_load         = 1.0;

        // Relaxation must exceed diffusion timescale: t_diff = h^2 * eta / E0
        const double t_diff  = height * height * eta_drag / E0;
        const double t_relax = std::max(10.0 * tau, 5.0 * t_diff);

        const unsigned n_load       = static_cast<unsigned>(std::round(t_load  / dt));
        const unsigned n_relax      = static_cast<unsigned>(std::round(t_relax / dt));
        const unsigned total_steps  = n_load + n_relax;
        const unsigned vtp_interval = std::max(1u, total_steps / 100u);

        // ── Construct ghost node field ────────────────────────────
        RandomNumberGenerator::Instance()->Reseed(42);

        boost::shared_ptr<ViscoelasticGhostNodeEcmField<3>> p_field(
            new ViscoelasticGhostNodeEcmField<3>("random", spacing,
                                                  x_min, x_max,
                                                  y_min, y_max,
                                                  z_min, z_max,
                                                  "cubic"));

        p_field->SetRelaxedStiffness(E0);
        p_field->SetRelaxationModulus(E1);
        p_field->SetRelaxationTime(tau);
        p_field->SetGhostDamping(eta_drag);
        p_field->SetAnisotropyStrength(0.0);

        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
            p_field->rGetNode(i).density = 1.0;

        // ── Identify top/bottom y-face nodes ─────────────────────
        std::vector<unsigned> bottom_nodes, top_nodes;
        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            const auto& node = p_field->GetNode(i);
            if (node.position[1] < y_min + 0.5 * spacing) bottom_nodes.push_back(i);
            if (node.position[1] > y_max - 0.5 * spacing) top_nodes.push_back(i);
        }
        TS_ASSERT_EQUALS(bottom_nodes.size(), num_cols_x * num_cols_z);
        TS_ASSERT_EQUALS(top_nodes.size(),    num_cols_x * num_cols_z);

        // Store original positions of boundary nodes
        std::vector<double> bottom_orig_x(bottom_nodes.size()),
                            bottom_orig_y(bottom_nodes.size()),
                            bottom_orig_z(bottom_nodes.size());
        for (unsigned k = 0; k < bottom_nodes.size(); k++)
        {
            bottom_orig_x[k] = p_field->GetNode(bottom_nodes[k]).position[0];
            bottom_orig_y[k] = p_field->GetNode(bottom_nodes[k]).position[1];
            bottom_orig_z[k] = p_field->GetNode(bottom_nodes[k]).position[2];
        }
        std::vector<double> top_orig_x(top_nodes.size()),
                            top_orig_y(top_nodes.size()),
                            top_orig_z(top_nodes.size());
        for (unsigned k = 0; k < top_nodes.size(); k++)
        {
            top_orig_x[k] = p_field->GetNode(top_nodes[k]).position[0];
            top_orig_y[k] = p_field->GetNode(top_nodes[k]).position[1];
            top_orig_z[k] = p_field->GetNode(top_nodes[k]).position[2];
        }

        // ── Open CSV output ───────────────────────────────────────
        OutputFileHandler handler("ViscoelasticCompression/" + outDir, false);
        out_stream p_csv = handler.OpenOutputFile("shear_data.csv");
        *p_csv << "time,shear_displacement,shear_strain,shear_stress,shear_stress_Pa,total_force_x_top\n";

        std::string out_dir = handler.GetOutputDirectoryFullPath();
        out_stream p_pvd;
        unsigned vtp_counter = 0;

        if (!skipVtp)
        {
            p_pvd = handler.OpenOutputFile("ghost_ecm_results.pvd");
            *p_pvd << "<?xml version=\"1.0\"?>\n"
                   << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
                   << "  <Collection>\n";

            std::ostringstream fname;
            fname << "ghost_ecm_" << vtp_counter << ".vtp";
            p_field->WriteOutput(out_dir + fname.str(), 0.0);
            if (p_pvd)
                *p_pvd << "    <DataSet timestep=\"0\" group=\"\" part=\"0\" file=\""
                       << fname.str() << "\"/>\n";
            ++vtp_counter;
        }

        // ── Tracking ──────────────────────────────────────────────
        double peak_shear_stress    = 0.0;
        double final_shear_stress   = 0.0;
        double prev_relax_stress    = 1e10;
        bool   relaxation_monotonic = true;

        // ── Phase 1: Shear Loading ────────────────────────────────
        for (unsigned step = 0; step < n_load; step++)
        {
            double time         = (step + 1) * dt;
            double shear_offset = shear_velocity * time;

            EnforceShearBoundaryConditions3d(*p_field,
                bottom_nodes, top_nodes,
                bottom_orig_x, bottom_orig_y, bottom_orig_z,
                top_orig_x,    top_orig_y,    top_orig_z,
                shear_offset);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fx_top = 0.0;
            for (unsigned k : top_nodes)
                total_fx_top += p_field->rGetNode(k).force[0];

            EnforceShearBoundaryConditions3d(*p_field,
                bottom_nodes, top_nodes,
                bottom_orig_x, bottom_orig_y, bottom_orig_z,
                top_orig_x,    top_orig_y,    top_orig_z,
                shear_offset);

            double shear_stress    = -total_fx_top / face_area;
            double shear_stress_Pa = shear_stress * kF / spacing_m;

            *p_csv << time << "," << shear_offset << "," << (shear_offset / height)
                   << "," << shear_stress << "," << shear_stress_Pa
                   << "," << total_fx_top << "\n";

            if (shear_stress > peak_shear_stress) peak_shear_stress = shear_stress;

            if (!skipVtp && (step + 1) % vtp_interval == 0)
            {
                std::ostringstream fname;
                fname << "ghost_ecm_" << vtp_counter << ".vtp";
                p_field->WriteOutput(out_dir + fname.str(), time);
                *p_pvd << "    <DataSet timestep=\"" << time
                       << "\" group=\"\" part=\"0\" file=\"" << fname.str() << "\"/>\n";
                ++vtp_counter;
            }
        }

        // ── Phase 2: Shear Relaxation ─────────────────────────────
        for (unsigned step = 0; step < n_relax; step++)
        {
            double time = t_load + (step + 1) * dt;

            EnforceShearBoundaryConditions3d(*p_field,
                bottom_nodes, top_nodes,
                bottom_orig_x, bottom_orig_y, bottom_orig_z,
                top_orig_x,    top_orig_y,    top_orig_z,
                shear_disp_max);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fx_top = 0.0;
            for (unsigned k : top_nodes)
                total_fx_top += p_field->rGetNode(k).force[0];

            EnforceShearBoundaryConditions3d(*p_field,
                bottom_nodes, top_nodes,
                bottom_orig_x, bottom_orig_y, bottom_orig_z,
                top_orig_x,    top_orig_y,    top_orig_z,
                shear_disp_max);

            double shear_stress    = -total_fx_top / face_area;
            double shear_stress_Pa = shear_stress * kF / spacing_m;

            *p_csv << time << "," << shear_disp_max << "," << gamma_max
                   << "," << shear_stress << "," << shear_stress_Pa
                   << "," << total_fx_top << "\n";

            if (shear_stress > prev_relax_stress + 1e-4) relaxation_monotonic = false;
            prev_relax_stress = shear_stress;
            final_shear_stress = shear_stress;

            if (!skipVtp && (step + 1) % vtp_interval == 0)
            {
                std::ostringstream fname;
                fname << "ghost_ecm_" << vtp_counter << ".vtp";
                p_field->WriteOutput(out_dir + fname.str(), time);
                *p_pvd << "    <DataSet timestep=\"" << time
                       << "\" group=\"\" part=\"0\" file=\"" << fname.str() << "\"/>\n";
                ++vtp_counter;
            }
        }

        p_csv->close();

        if (!skipVtp)
        {
            *p_pvd << "  </Collection>\n</VTKFile>\n";
            p_pvd->close();
        }

        // ── Compute recovered shear modulus ────────────────────────
        double peak_shear_stress_Pa  = peak_shear_stress  * kF / spacing_m;
        double final_shear_stress_Pa = final_shear_stress * kF / spacing_m;
        double G_relaxed_meas = final_shear_stress_Pa / gamma_max;
        double G_peak_meas    = peak_shear_stress_Pa  / gamma_max;

        // ── Assertions ────────────────────────────────────────────
        TS_ASSERT_LESS_THAN(0.0, final_shear_stress);
        TS_ASSERT_LESS_THAN(0.0, peak_shear_stress + final_shear_stress);

        // ── Machine-parseable summary (one line, for sweep scripts) ─
        std::cout << "SHEAR_RESULT,"
                  << G_input_Pa << ","
                  << G_relaxed_meas << ","
                  << G_peak_meas << ","
                  << E0 << "," << E1 << ","
                  << tau << "," << eta_drag << ","
                  << relaxation_monotonic << std::endl;

        std::cout << "\n=== 3D Simple Shear Test Results ===\n"
                  << "Input G = " << G_input_Pa << " Pa  |  E0=" << E0 << "  E1=" << E1
                  << "  tau=" << tau << " h  dt=" << dt << "  eta=" << eta_drag << "\n"
                  << "Grid: 50x10x10  Lattice: cubic (18-connected)\n"
                  << "Shear strain gamma:     " << gamma_max << "\n"
                  << "Face area (sim):        " << face_area << " sim²\n"
                  << "t_diff = " << t_diff << " h  t_relax = " << t_relax << " h\n"
                  << "Peak shear stress:      " << peak_shear_stress
                  << " sim  =  " << peak_shear_stress_Pa << " Pa\n"
                  << "Final shear stress:     " << final_shear_stress
                  << " sim  =  " << final_shear_stress_Pa << " Pa\n"
                  << "Stress ratio f/p:       "
                  << (peak_shear_stress > 0 ? final_shear_stress / peak_shear_stress : 0.0) << "\n"
                  << "─── Recovered Shear Modulus ───\n"
                  << "G_peak  = peak_stress/gamma  = " << G_peak_meas  << " Pa\n"
                  << "G_relax = final_stress/gamma = " << G_relaxed_meas << " Pa\n"
                  << "Input G = " << G_input_Pa << " Pa\n"
                  << "Ratio G_relax/G_input:        " << G_relaxed_meas / G_input_Pa << "\n"
                  << "Ratio G_peak/G_input:         " << G_peak_meas / G_input_Pa << "\n"
                  << "Relaxation monotonic:   " << (relaxation_monotonic ? "yes" : "NO") << "\n"
                  << "VTP frames written:     " << vtp_counter << "\n"
                  << "Output dir: " << out_dir << "\n" << std::endl;
    }
};

#endif // TESTVISCOELASTICSHEAREXPERIMENT3D_HPP_
