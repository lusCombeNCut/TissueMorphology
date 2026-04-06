#ifndef TESTVISCOELASTICCOMPRESSIONEXPERIMENT_HPP_
#define TESTVISCOELASTICCOMPRESSIONEXPERIMENT_HPP_

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
 * Confined uniaxial compression stress-relaxation test on a pure
 * viscoelastic ghost node network (no cells).
 *
 * Validates that the emergent bulk mechanical response of the discrete
 * SLS spring network matches the expected Standard Linear Solid behaviour:
 *   E(t) = E_0 + E_1 * exp(-t/tau)
 *
 * Protocol (analogous to parallel plate rheometry):
 *   Phase 1 (Loading):    Compress top plate downward at constant velocity
 *   Phase 2 (Relaxation): Hold displacement, let internal stress decay
 *
 * SLS parameters derived from G=800 Pa via CryptBuddingParams::Finalise():
 *   E0_phys = G * spacing_phys  (Cauchy-Born, Ostoja-Starzewski 2002)
 *   SLS_ARM_RATIO=0.2, spacing_phys=10 µm, kF=0.001 N/m per sim unit
 *
 * Reference: Fertala et al. (2025) doi:10.1101/2025.07.02.662292
 */
class TestViscoelasticCompressionExperiment : public AbstractCellBasedTestSuite
{
private:

    /**
     * Enforce all boundary conditions: pin bottom, set top to prescribed y,
     * pin left/right walls (x only). Call before AND after each timestep.
     */
    void EnforceBoundaryConditions(
        ViscoelasticGhostNodeEcmField<2>& rField,
        const std::vector<unsigned>& rBottomNodes,
        const std::vector<unsigned>& rTopNodes,
        const std::vector<unsigned>& rLeftNodes,
        const std::vector<unsigned>& rRightNodes,
        const std::vector<double>& rOriginalX,
        const std::vector<double>& rBottomOriginalY,
        double prescribedYTop,
        double xMin, double xMax)
    {
        for (unsigned k = 0; k < rBottomNodes.size(); k++)
        {
            unsigned idx = rBottomNodes[k];
            rField.rGetNode(idx).position[0] = rOriginalX[idx];
            rField.rGetNode(idx).position[1] = rBottomOriginalY[k];
        }
        for (unsigned k = 0; k < rTopNodes.size(); k++)
        {
            unsigned idx = rTopNodes[k];
            rField.rGetNode(idx).position[0] = rOriginalX[idx];
            rField.rGetNode(idx).position[1] = prescribedYTop;
        }
        for (unsigned idx : rLeftNodes)
        {
            rField.rGetNode(idx).position[0] = xMin;
        }
        for (unsigned idx : rRightNodes)
        {
            rField.rGetNode(idx).position[0] = xMax;
        }
    }

public:

    /**
     * Enforce boundary conditions WITHOUT lateral wall constraints (unconfined).
     * Only the top and bottom plates are driven; lateral edges are free to expand.
     * A single bottom-row node is pinned in x to stop rigid-body lateral drift.
     */
    void EnforceBCUnconfined(
        ViscoelasticGhostNodeEcmField<2>& rField,
        const std::vector<unsigned>& rBottomNodes,
        const std::vector<unsigned>& rTopNodes,
        const std::vector<double>& rOriginalX,
        const std::vector<double>& rBottomOriginalY,
        double prescribedYTop,
        unsigned anchorIdx)
    {
        for (unsigned k = 0; k < rBottomNodes.size(); k++)
        {
            unsigned idx = rBottomNodes[k];
            rField.rGetNode(idx).position[1] = rBottomOriginalY[k];
        }
        for (unsigned k = 0; k < rTopNodes.size(); k++)
        {
            unsigned idx = rTopNodes[k];
            rField.rGetNode(idx).position[1] = prescribedYTop;
        }
        // Prevent net rigid-body lateral translation
        rField.rGetNode(anchorIdx).position[0] = rOriginalX[anchorIdx];
    }

    void TestUniaxialStressRelaxation()
    {
        // ── Grid geometry ────────────────────────────────────────
        const double spacing = 1.0;
        const double x_min = 0.0, x_max = 19.0;
        const double y_min = 0.0, y_max = 19.0;
        const unsigned num_cols = static_cast<unsigned>((x_max - x_min) / spacing) + 1; // 20
        const double width = x_max - x_min;   // 19.0 (for stress = force/width)
        const double L0 = y_max - y_min;       // 19.0 initial height

        // ── SLS parameters (G=800 Pa → sim units) ───────────────
        // Cauchy-Born: G = k_bond / Δx → E0 = G * Δx_phys / kF
        // G=800 Pa, Δx=10µm, kF=0.001 N/m → E0=8.0, E1=1.6 (r=0.2)
        const double E0 = 8.0;
        const double E1 = 1.6;
        const double tau = 1.0;    // Relaxation time (hours)
        // CFL: (E0+E1)*N*dt/eta < 2, N=8 for square lattice
        // (9.6)*8*0.001/1.0 = 0.077 (safe)
        const double eta_drag = 1.0;

        // ── Time-stepping ────────────────────────────────────────
        const double dt = 0.001;                   // Reduced for accuracy
        const double compression_velocity = 1.0;   // sim lengths per hour
        const double t_load = 1.0;                 // loading duration (hours)
        const double t_relax = 10.0;               // relaxation duration (10*tau)
        const unsigned n_load  = static_cast<unsigned>(std::round(t_load  / dt));
        const unsigned n_relax = static_cast<unsigned>(std::round(t_relax / dt));
        const double total_displacement = compression_velocity * t_load; // 1.0
        const double final_strain = total_displacement / L0;             // ~0.053

        // VTP output: sample every vtp_interval steps (~100 frames total)
        const unsigned total_steps = n_load + n_relax;
        const unsigned vtp_interval = std::max(1u, total_steps / 100u);

        // ── Construct ghost node field ───────────────────────────
        RandomNumberGenerator::Instance()->Reseed(42);

        boost::shared_ptr<ViscoelasticGhostNodeEcmField<2>> p_field(
            new ViscoelasticGhostNodeEcmField<2>("random", spacing,
                                                  x_min, x_max,
                                                  y_min, y_max,
                                                  "square"));

        p_field->SetRelaxedStiffness(E0);
        p_field->SetRelaxationModulus(E1);
        p_field->SetRelaxationTime(tau);
        p_field->SetGhostDamping(eta_drag);
        p_field->SetAnisotropyStrength(0.0);

        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            p_field->rGetNode(i).density    = 1.0;
        }

        // ── Identify boundary nodes ──────────────────────────────
        std::vector<unsigned> bottom_nodes, top_nodes, left_nodes, right_nodes;
        std::vector<double> original_x;

        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            const auto& node = p_field->GetNode(i);
            original_x.push_back(node.position[0]);

            if (node.position[1] < y_min + 0.5 * spacing) bottom_nodes.push_back(i);
            if (node.position[1] > y_max - 0.5 * spacing) top_nodes.push_back(i);
            if (node.position[0] < x_min + 0.5 * spacing) left_nodes.push_back(i);
            if (node.position[0] > x_max - 0.5 * spacing) right_nodes.push_back(i);
        }

        TS_ASSERT_EQUALS(bottom_nodes.size(), num_cols);
        TS_ASSERT_EQUALS(top_nodes.size(), num_cols);

        std::vector<double> bottom_original_y(bottom_nodes.size());
        for (unsigned k = 0; k < bottom_nodes.size(); k++)
        {
            bottom_original_y[k] = p_field->GetNode(bottom_nodes[k]).position[1];
        }

        // ── Open CSV output ──────────────────────────────────────
        OutputFileHandler handler("ViscoelasticCompression/StressRelaxation_2D", false);
        out_stream p_csv = handler.OpenOutputFile("compression_data.csv");
        *p_csv << "time,displacement,strain,stress,total_force_y\n";

        // ── Open PVD index for ParaView ──────────────────────────
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

        // Write initial state
        WriteVtp(0.0);

        // ── Tracking variables ───────────────────────────────────
        double peak_stress = 0.0;
        double final_stress = 0.0;
        double prev_relax_stress = 1e10;
        bool relaxation_monotonic = true;

        // ── Phase 1: Loading ─────────────────────────────────────
        for (unsigned step = 0; step < n_load; step++)
        {
            double time = (step + 1) * dt;
            double displacement = compression_velocity * time;
            double prescribed_y_top = y_max - displacement;

            EnforceBoundaryConditions(*p_field, bottom_nodes, top_nodes,
                                      left_nodes, right_nodes,
                                      original_x, bottom_original_y,
                                      prescribed_y_top, x_min, x_max);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fy_top = 0.0;
            for (unsigned k : top_nodes)
                total_fy_top += p_field->rGetNode(k).force[1];

            EnforceBoundaryConditions(*p_field, bottom_nodes, top_nodes,
                                      left_nodes, right_nodes,
                                      original_x, bottom_original_y,
                                      prescribed_y_top, x_min, x_max);

            double stress = total_fy_top / width;
            *p_csv << time << "," << displacement << "," << (displacement / L0)
                   << "," << stress << "," << total_fy_top << "\n";

            if (stress > peak_stress) peak_stress = stress;

            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        // ── Phase 2: Relaxation ──────────────────────────────────
        double y_top_hold = y_max - total_displacement;

        for (unsigned step = 0; step < n_relax; step++)
        {
            double time = t_load + (step + 1) * dt;

            EnforceBoundaryConditions(*p_field, bottom_nodes, top_nodes,
                                      left_nodes, right_nodes,
                                      original_x, bottom_original_y,
                                      y_top_hold, x_min, x_max);

            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fy_top = 0.0;
            for (unsigned k : top_nodes)
                total_fy_top += p_field->rGetNode(k).force[1];

            EnforceBoundaryConditions(*p_field, bottom_nodes, top_nodes,
                                      left_nodes, right_nodes,
                                      original_x, bottom_original_y,
                                      y_top_hold, x_min, x_max);

            double stress = total_fy_top / width;
            *p_csv << time << "," << total_displacement << ","
                   << final_strain << "," << stress << "," << total_fy_top << "\n";

            if (stress > prev_relax_stress + 1e-4) relaxation_monotonic = false;
            prev_relax_stress = stress;
            final_stress = stress;

            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        p_csv->close();

        // Close PVD
        *p_pvd << "  </Collection>\n</VTKFile>\n";
        p_pvd->close();

        // ── Assertions ───────────────────────────────────────────
        TS_ASSERT_LESS_THAN(final_stress, peak_stress);  // relaxation occurred

        double stress_ratio = final_stress / peak_stress;
        TS_ASSERT_LESS_THAN(0.0, stress_ratio);   // solid-like (E0 arm remains)
        TS_ASSERT_LESS_THAN(stress_ratio, 1.0);    // viscous arm relaxed

        TS_ASSERT(relaxation_monotonic);
        TS_ASSERT_LESS_THAN(0.0, final_stress);    // still under compression

        std::cout << "\n=== Compression Test Results (Confined) ===\n"
                  << "G = 800 Pa  |  E0=" << E0 << "  E1=" << E1
                  << "  tau=" << tau << " h  dt=" << dt << "\n"
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
     * Unconfined uniaxial compression — lateral edges free to expand (Poisson effect).
     * Only top/bottom plates are driven. A single bottom-row node is x-anchored
     * to prevent net rigid-body lateral drift without constraining Poisson expansion.
     * Output: ViscoelasticCompression/StressRelaxation_2D_Unconfined/
     */
    void TestUniaxialStressRelaxationUnconfined()
    {
        const double spacing = 1.0;
        const double x_min = 0.0, x_max = 19.0;
        const double y_min = 0.0, y_max = 19.0;
        const double width = x_max - x_min;
        const double L0    = y_max - y_min;

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
        boost::shared_ptr<ViscoelasticGhostNodeEcmField<2>> p_field(
            new ViscoelasticGhostNodeEcmField<2>("random", spacing,
                                                  x_min, x_max,
                                                  y_min, y_max,
                                                  "square"));
        p_field->SetRelaxedStiffness(E0);
        p_field->SetRelaxationModulus(E1);
        p_field->SetRelaxationTime(tau);
        p_field->SetGhostDamping(eta_drag);
        p_field->SetAnisotropyStrength(0.0);
        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
            p_field->rGetNode(i).density = 1.0;

        std::vector<unsigned> bottom_nodes, top_nodes;
        std::vector<double>   original_x;
        for (unsigned i = 0; i < p_field->GetNumNodes(); i++)
        {
            const auto& node = p_field->GetNode(i);
            original_x.push_back(node.position[0]);
            if (node.position[1] < y_min + 0.5 * spacing) bottom_nodes.push_back(i);
            if (node.position[1] > y_max - 0.5 * spacing) top_nodes.push_back(i);
        }

        std::vector<double> bottom_original_y(bottom_nodes.size());
        for (unsigned k = 0; k < bottom_nodes.size(); k++)
            bottom_original_y[k] = p_field->GetNode(bottom_nodes[k]).position[1];

        // Anchor the middle bottom node in x to prevent rigid-body drift
        unsigned anchor_idx = bottom_nodes[bottom_nodes.size() / 2];

        OutputFileHandler handler("ViscoelasticCompression/StressRelaxation_2D_Unconfined", false);
        out_stream p_csv = handler.OpenOutputFile("compression_data.csv");
        *p_csv << "time,displacement,strain,stress,total_force_y\n";

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
            double prescribed_y_top = y_max - displacement;

            EnforceBCUnconfined(*p_field, bottom_nodes, top_nodes,
                                original_x, bottom_original_y,
                                prescribed_y_top, anchor_idx);
            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fy_top = 0.0;
            for (unsigned k : top_nodes)
                total_fy_top += p_field->rGetNode(k).force[1];

            EnforceBCUnconfined(*p_field, bottom_nodes, top_nodes,
                                original_x, bottom_original_y,
                                prescribed_y_top, anchor_idx);

            double stress = total_fy_top / width;
            *p_csv << time << "," << displacement << "," << (displacement / L0)
                   << "," << stress << "," << total_fy_top << "\n";
            if (stress > peak_stress) peak_stress = stress;
            if ((step + 1) % vtp_interval == 0) WriteVtp(time);
        }

        double y_top_hold = y_max - total_displacement;
        for (unsigned step = 0; step < n_relax; step++)
        {
            double time = t_load + (step + 1) * dt;

            EnforceBCUnconfined(*p_field, bottom_nodes, top_nodes,
                                original_x, bottom_original_y,
                                y_top_hold, anchor_idx);
            p_field->ClearGhostForces();
            p_field->UpdateGhostNodePositions(dt);

            double total_fy_top = 0.0;
            for (unsigned k : top_nodes)
                total_fy_top += p_field->rGetNode(k).force[1];

            EnforceBCUnconfined(*p_field, bottom_nodes, top_nodes,
                                original_x, bottom_original_y,
                                y_top_hold, anchor_idx);

            double stress = total_fy_top / width;
            *p_csv << time << "," << total_displacement << ","
                   << final_strain << "," << stress << "," << total_fy_top << "\n";
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

        std::cout << "\n=== Compression Test Results (Unconfined) ===\n"
                  << "G = 800 Pa  |  E0=" << E0 << "  E1=" << E1
                  << "  tau=" << tau << " h  dt=" << dt << "\n"
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

#endif // TESTVISCOELASTICCOMPRESSIONEXPERIMENT_HPP_
