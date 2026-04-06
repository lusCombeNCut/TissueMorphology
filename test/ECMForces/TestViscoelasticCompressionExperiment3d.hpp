#ifndef TESTVISCOELASTICCOMPRESSIONEXPERIMENT3D_HPP_
#define TESTVISCOELASTICCOMPRESSIONEXPERIMENT3D_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
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
 * 3D confined uniaxial compression stress-relaxation test on a pure
 * viscoelastic ghost node network (no cells).
 *
 * Analogous to the 2D test (TestViscoelasticCompressionExperiment.hpp) but
 * using ViscoelasticGhostNodeEcmField<3> with a cubic grid. Compression is
 * applied in the z-direction (top z-face moves, bottom z-face pinned, lateral
 * faces confined in x and y).
 *
 * SLS parameters and unit convention are identical to the 2D test:
 *   G = 800 Pa → E0=464, E1=23.2, tau=1.0h, eta=50 (sim units)
 *
 * Output dir: $CHASTE_TEST_OUTPUT/ViscoelasticCompression/StressRelaxation_3D/
 *
 * Reference: Fertala et al. (2025) doi:10.1101/2025.07.02.662292
 */
class TestViscoelasticCompressionExperiment3d : public AbstractCellBasedTestSuite
{
private:

    /**
     * Enforce all boundary conditions for the 3D slab:
     *   - Bottom z-face: pin x, y, z
     *   - Top z-face:    pin x, y; prescribe z
     *   - Four lateral faces (±x, ±y): pin the normal direction only
     */
    void EnforceBoundaryConditions(
        ViscoelasticGhostNodeEcmField<3>& rField,
        const std::vector<unsigned>& rBottomNodes,
        const std::vector<unsigned>& rTopNodes,
        const std::vector<unsigned>& rXMinNodes,
        const std::vector<unsigned>& rXMaxNodes,
        const std::vector<unsigned>& rYMinNodes,
        const std::vector<unsigned>& rYMaxNodes,
        const std::vector<double>& rOriginalX,
        const std::vector<double>& rOriginalY,
        const std::vector<double>& rBottomOriginalZ,
        double prescribedZTop,
        double xMin, double xMax,
        double yMin, double yMax)
    {
        for (unsigned k = 0; k < rBottomNodes.size(); k++)
        {
            unsigned idx = rBottomNodes[k];
            rField.rGetNode(idx).position[0] = rOriginalX[idx];
            rField.rGetNode(idx).position[1] = rOriginalY[idx];
            rField.rGetNode(idx).position[2] = rBottomOriginalZ[k];
        }
        for (unsigned k = 0; k < rTopNodes.size(); k++)
        {
            unsigned idx = rTopNodes[k];
            rField.rGetNode(idx).position[0] = rOriginalX[idx];
            rField.rGetNode(idx).position[1] = rOriginalY[idx];
            rField.rGetNode(idx).position[2] = prescribedZTop;
        }
        for (unsigned idx : rXMinNodes) rField.rGetNode(idx).position[0] = xMin;
        for (unsigned idx : rXMaxNodes) rField.rGetNode(idx).position[0] = xMax;
        for (unsigned idx : rYMinNodes) rField.rGetNode(idx).position[1] = yMin;
        for (unsigned idx : rYMaxNodes) rField.rGetNode(idx).position[1] = yMax;
    }

public:

    void TestUniaxialStressRelaxation3d()
    {
        // ── Grid geometry ─────────────────────────────────────────
        const double spacing = 1.0;
        const double x_min = 0.0, x_max = 9.0;   // 10 nodes → smaller for speed
        const double y_min = 0.0, y_max = 9.0;
        const double z_min = 0.0, z_max = 9.0;
        const unsigned num_cols_x = static_cast<unsigned>((x_max - x_min) / spacing) + 1; // 10
        const unsigned num_cols_y = static_cast<unsigned>((y_max - y_min) / spacing) + 1; // 10
        const double face_area = (x_max - x_min) * (y_max - y_min); // 9×9 = 81 sim²
        const double L0 = z_max - z_min;   // 9.0 initial height in z

        // ── SLS parameters (G=800 Pa → sim units) ──
        // Cauchy-Born: E0 = G * Δx_phys / kF = 800 * 10e-6 / 0.001 = 8.0
        const double E0       = 8.0;
        const double E1       = 1.6;    // SLS arm ratio 0.2
        const double tau      = 1.0;    // Relaxation time (hours)
        // CFL: (9.6)*18*0.001/1.0 = 0.17 (safe)
        const double eta_drag = 1.0;

        // ── Time-stepping ─────────────────────────────────────────
        const double dt                   = 0.001;
        const double compression_velocity = 1.0;   // sim lengths per hour
        const double t_load               = 1.0;   // loading duration (h)
        const double t_relax              = 10.0;  // relaxation duration (10*tau)
        const unsigned n_load  = static_cast<unsigned>(std::round(t_load  / dt));
        const unsigned n_relax = static_cast<unsigned>(std::round(t_relax / dt));
        const double total_displacement = compression_velocity * t_load; // 1.0
        const double final_strain       = total_displacement / L0;       // ~0.111

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
        {
            p_field->rGetNode(i).density = 1.0;
        }

        // ── Identify boundary nodes ───────────────────────────────
        std::vector<unsigned> bottom_nodes, top_nodes;
        std::vector<unsigned> xmin_nodes, xmax_nodes, ymin_nodes, ymax_nodes;
        std::vector<double>   original_x, original_y;

        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            const auto& node = p_field->GetNode(i);
            original_x.push_back(node.position[0]);
            original_y.push_back(node.position[1]);

            if (node.position[2] < z_min + 0.5 * spacing) bottom_nodes.push_back(i);
            if (node.position[2] > z_max - 0.5 * spacing) top_nodes.push_back(i);
            if (node.position[0] < x_min + 0.5 * spacing) xmin_nodes.push_back(i);
            if (node.position[0] > x_max - 0.5 * spacing) xmax_nodes.push_back(i);
            if (node.position[1] < y_min + 0.5 * spacing) ymin_nodes.push_back(i);
            if (node.position[1] > y_max - 0.5 * spacing) ymax_nodes.push_back(i);
        }

        TS_ASSERT_EQUALS(bottom_nodes.size(), num_cols_x * num_cols_y);
        TS_ASSERT_EQUALS(top_nodes.size(),    num_cols_x * num_cols_y);

        std::vector<double> bottom_original_z(bottom_nodes.size());
        for (unsigned k = 0; k < bottom_nodes.size(); k++)
        {
            bottom_original_z[k] = p_field->GetNode(bottom_nodes[k]).position[2];
        }

        // ── Open CSV output ───────────────────────────────────────
        OutputFileHandler handler("ViscoelasticCompression/StressRelaxation_3D", false);
        out_stream p_csv = handler.OpenOutputFile("compression_data.csv");
        *p_csv << "time,displacement,strain,stress,total_force_z\n";

        // ── Open PVD index for ParaView ───────────────────────────
        std::string out_dir = handler.GetOutputDirectoryFullPath();
        out_stream p_pvd = handler.OpenOutputFile("ghost_ecm_results.pvd");
        *p_pvd << "<?xml version=\"1.0\"?>\n"
               << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
               << "  <Collection>\n";

        unsigned vtp_counter = 0;
        auto WriteVtp = [&](double time)
        {
            std::ostringstream fname;
            fname << "ghost_ecm_" << vtp_counter << ".vtp";
            p_field->WriteOutput(out_dir + fname.str(), time);
            *p_pvd << "    <DataSet timestep=\"" << time
                   << "\" group=\"\" part=\"0\" file=\"" << fname.str() << "\"/>\n";
            ++vtp_counter;
        };

        WriteVtp(0.0);

        // ── Tracking variables ────────────────────────────────────
        double peak_stress       = 0.0;
        double final_stress      = 0.0;
        double prev_relax_stress = 1e10;
        bool   relaxation_monotonic = true;

        // ── Phase 1: Loading ──────────────────────────────────────
        for (unsigned step = 0; step < n_load; step++)
        {
            double time           = (step + 1) * dt;
            double displacement   = compression_velocity * time;
            double prescribed_z_top = z_max - displacement;

            EnforceBoundaryConditions(*p_field,
                bottom_nodes, top_nodes,
                xmin_nodes, xmax_nodes, ymin_nodes, ymax_nodes,
                original_x, original_y, bottom_original_z,
                prescribed_z_top, x_min, x_max, y_min, y_max);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fz_top = 0.0;
            for (unsigned k : top_nodes)
                total_fz_top += p_field->rGetNode(k).force[2];

            EnforceBoundaryConditions(*p_field,
                bottom_nodes, top_nodes,
                xmin_nodes, xmax_nodes, ymin_nodes, ymax_nodes,
                original_x, original_y, bottom_original_z,
                prescribed_z_top, x_min, x_max, y_min, y_max);

            double stress = total_fz_top / face_area;
            *p_csv << time << "," << displacement << "," << (displacement / L0)
                   << "," << stress << "," << total_fz_top << "\n";

            if (stress > peak_stress) peak_stress = stress;

            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        // ── Phase 2: Relaxation ───────────────────────────────────
        double z_top_hold = z_max - total_displacement;

        for (unsigned step = 0; step < n_relax; step++)
        {
            double time = t_load + (step + 1) * dt;

            EnforceBoundaryConditions(*p_field,
                bottom_nodes, top_nodes,
                xmin_nodes, xmax_nodes, ymin_nodes, ymax_nodes,
                original_x, original_y, bottom_original_z,
                z_top_hold, x_min, x_max, y_min, y_max);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fz_top = 0.0;
            for (unsigned k : top_nodes)
                total_fz_top += p_field->rGetNode(k).force[2];

            EnforceBoundaryConditions(*p_field,
                bottom_nodes, top_nodes,
                xmin_nodes, xmax_nodes, ymin_nodes, ymax_nodes,
                original_x, original_y, bottom_original_z,
                z_top_hold, x_min, x_max, y_min, y_max);

            double stress = total_fz_top / face_area;
            *p_csv << time << "," << total_displacement << ","
                   << final_strain << "," << stress << "," << total_fz_top << "\n";

            if (stress > prev_relax_stress + 1e-4) relaxation_monotonic = false;
            prev_relax_stress = stress;
            final_stress = stress;

            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        p_csv->close();

        *p_pvd << "  </Collection>\n</VTKFile>\n";
        p_pvd->close();

        // ── Assertions ────────────────────────────────────────────
        TS_ASSERT_LESS_THAN(final_stress, peak_stress);
        double stress_ratio = final_stress / peak_stress;
        TS_ASSERT_LESS_THAN(0.0, stress_ratio);
        TS_ASSERT_LESS_THAN(stress_ratio, 1.0);
        TS_ASSERT(relaxation_monotonic);
        TS_ASSERT_LESS_THAN(0.0, final_stress);

        std::cout << "\n=== 3D Compression Test Results (Confined) ===\n"
                  << "G = 800 Pa  |  E0=" << E0 << "  E1=" << E1
                  << "  tau=" << tau << " h  dt=" << dt << "\n"
                  << "Grid: " << num_cols_x << "x" << num_cols_y
                  << "x" << static_cast<unsigned>(L0/spacing)+1 << "\n"
                  << "Peak stress:          " << peak_stress << "\n"
                  << "Final stress:         " << final_stress << "\n"
                  << "Stress ratio f/p:     " << stress_ratio
                  << "  (step-strain theory: " << E0/(E0+E1) << ")\n"
                  << "Final strain:         " << final_strain << "\n"
                  << "Relaxation monotonic: " << (relaxation_monotonic ? "yes" : "NO") << "\n"
                  << "VTP frames written:   " << vtp_counter << "\n"
                  << "Output dir: " << out_dir << "\n" << std::endl;
    }

    /**
     * Unconfined uniaxial compression in 3D — all four lateral faces free
     * to expand (Poisson effect). Only top/bottom z-faces are driven.
     * A single bottom-face node is anchored in x and y to prevent rigid-body drift.
     * Output: ViscoelasticCompression/StressRelaxation_3D_Unconfined/
     */
    void TestUniaxialStressRelaxation3dUnconfined()
    {
        const double spacing = 1.0;
        const double x_min = 0.0, x_max = 9.0;
        const double y_min = 0.0, y_max = 9.0;
        const double z_min = 0.0, z_max = 9.0;
        const unsigned num_cols_x = static_cast<unsigned>((x_max - x_min) / spacing) + 1;
        const unsigned num_cols_y = static_cast<unsigned>((y_max - y_min) / spacing) + 1;
        const double face_area = (x_max - x_min) * (y_max - y_min);
        const double L0 = z_max - z_min;

        const double E0       = 8.0;
        const double E1       = 1.6;
        const double tau      = 1.0;
        const double eta_drag = 1.0;

        const double dt                   = 0.001;
        const double compression_velocity = 1.0;
        const double t_load               = 1.0;
        const double t_relax              = 10.0;
        const unsigned n_load  = static_cast<unsigned>(std::round(t_load  / dt));
        const unsigned n_relax = static_cast<unsigned>(std::round(t_relax / dt));
        const double total_displacement = compression_velocity * t_load;
        const double final_strain       = total_displacement / L0;

        const unsigned total_steps  = n_load + n_relax;
        const unsigned vtp_interval = std::max(1u, total_steps / 100u);

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

        std::vector<unsigned> bottom_nodes, top_nodes;
        std::vector<double>   original_x, original_y;

        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            const auto& node = p_field->GetNode(i);
            original_x.push_back(node.position[0]);
            original_y.push_back(node.position[1]);
            if (node.position[2] < z_min + 0.5 * spacing) bottom_nodes.push_back(i);
            if (node.position[2] > z_max - 0.5 * spacing) top_nodes.push_back(i);
        }
        TS_ASSERT_EQUALS(bottom_nodes.size(), num_cols_x * num_cols_y);
        TS_ASSERT_EQUALS(top_nodes.size(),    num_cols_x * num_cols_y);

        std::vector<double> bottom_original_z(bottom_nodes.size());
        for (unsigned k = 0; k < bottom_nodes.size(); k++)
            bottom_original_z[k] = p_field->GetNode(bottom_nodes[k]).position[2];

        // Anchor the centre bottom node in x and y to prevent rigid-body drift
        unsigned anchor_idx = bottom_nodes[bottom_nodes.size() / 2];

        // Lambda: enforce BCs without lateral walls
        auto EnforceBCUnconfined3d = [&](double prescribed_z_top)
        {
            for (unsigned k = 0; k < bottom_nodes.size(); k++)
                p_field->rGetNode(bottom_nodes[k]).position[2] = bottom_original_z[k];
            for (unsigned k = 0; k < top_nodes.size(); k++)
                p_field->rGetNode(top_nodes[k]).position[2] = prescribed_z_top;
            // Prevent lateral rigid-body translation
            p_field->rGetNode(anchor_idx).position[0] = original_x[anchor_idx];
            p_field->rGetNode(anchor_idx).position[1] = original_y[anchor_idx];
        };

        OutputFileHandler handler("ViscoelasticCompression/StressRelaxation_3D_Unconfined", false);
        out_stream p_csv = handler.OpenOutputFile("compression_data.csv");
        *p_csv << "time,displacement,strain,stress,total_force_z\n";

        std::string out_dir = handler.GetOutputDirectoryFullPath();
        out_stream p_pvd = handler.OpenOutputFile("ghost_ecm_results.pvd");
        *p_pvd << "<?xml version=\"1.0\"?>\n"
               << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
               << "  <Collection>\n";

        unsigned vtp_counter = 0;
        auto WriteVtp = [&](double time)
        {
            std::ostringstream fname;
            fname << "ghost_ecm_" << vtp_counter << ".vtp";
            p_field->WriteOutput(out_dir + fname.str(), time);
            *p_pvd << "    <DataSet timestep=\"" << time
                   << "\" group=\"\" part=\"0\" file=\"" << fname.str() << "\"/>\n";
            ++vtp_counter;
        };
        WriteVtp(0.0);

        double peak_stress       = 0.0;
        double final_stress      = 0.0;
        double prev_relax_stress = 1e10;
        bool   relaxation_monotonic = true;

        for (unsigned step = 0; step < n_load; step++)
        {
            double time             = (step + 1) * dt;
            double displacement     = compression_velocity * time;
            double prescribed_z_top = z_max - displacement;

            EnforceBCUnconfined3d(prescribed_z_top);
            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fz_top = 0.0;
            for (unsigned k : top_nodes)
                total_fz_top += p_field->rGetNode(k).force[2];

            EnforceBCUnconfined3d(prescribed_z_top);

            double stress = total_fz_top / face_area;
            *p_csv << time << "," << displacement << "," << (displacement / L0)
                   << "," << stress << "," << total_fz_top << "\n";
            if (stress > peak_stress) peak_stress = stress;
            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        double z_top_hold = z_max - total_displacement;
        for (unsigned step = 0; step < n_relax; step++)
        {
            double time = t_load + (step + 1) * dt;

            EnforceBCUnconfined3d(z_top_hold);
            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fz_top = 0.0;
            for (unsigned k : top_nodes)
                total_fz_top += p_field->rGetNode(k).force[2];

            EnforceBCUnconfined3d(z_top_hold);

            double stress = total_fz_top / face_area;
            *p_csv << time << "," << total_displacement << ","
                   << final_strain << "," << stress << "," << total_fz_top << "\n";
            if (stress > prev_relax_stress + 1e-4) relaxation_monotonic = false;
            prev_relax_stress = stress;
            final_stress = stress;
            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        p_csv->close();
        *p_pvd << "  </Collection>\n</VTKFile>\n";
        p_pvd->close();

        TS_ASSERT_LESS_THAN(final_stress, peak_stress);
        double stress_ratio = final_stress / peak_stress;
        TS_ASSERT_LESS_THAN(0.0, stress_ratio);
        TS_ASSERT_LESS_THAN(stress_ratio, 1.0);
        TS_ASSERT(relaxation_monotonic);
        TS_ASSERT_LESS_THAN(0.0, final_stress);

        std::cout << "\n=== 3D Compression Test Results (Unconfined) ===\n"
                  << "G = 800 Pa  |  E0=" << E0 << "  E1=" << E1
                  << "  tau=" << tau << " h  dt=" << dt << "\n"
                  << "Grid: " << num_cols_x << "x" << num_cols_y
                  << "x" << static_cast<unsigned>(L0/spacing)+1 << "\n"
                  << "Peak stress:          " << peak_stress << "\n"
                  << "Final stress:         " << final_stress << "\n"
                  << "Stress ratio f/p:     " << stress_ratio
                  << "  (step-strain theory: " << E0/(E0+E1) << ")\n"
                  << "Final strain:         " << final_strain << "\n"
                  << "Relaxation monotonic: " << (relaxation_monotonic ? "yes" : "NO") << "\n"
                  << "VTP frames written:   " << vtp_counter << "\n"
                  << "Output dir: " << out_dir << "\n" << std::endl;
    }
};

#endif // TESTVISCOELASTICCOMPRESSIONEXPERIMENT3D_HPP_
