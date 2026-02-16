/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TEST3DCRYPTORGANOID_HPP_
#define TEST3DCRYPTORGANOID_HPP_

/**
 * Test3dCryptOrganoid.hpp
 * 
 * Simulates a 3D intestinal crypt organoid using the cell-centre (node-based)
 * approach from TissueMorphology, with dynamic ECM forces including:
 *   - ECM contact guidance (fibers bias cell migration)
 *   - ECM degradation by cells (MMP-mediated)
 *   - ECM fiber remodeling by cell traction
 *   - Basement membrane confinement with time-dependent degradation
 *   - Cell proliferation with contact inhibition
 *   - Differential adhesion between stem/Paneth and transit-amplifying cells
 * 
 * Duration: 7 simulated days (168 hours)
 * 
 * This combines the cell-centre mechanics of TissueMorphology with
 * biological processes known to drive crypt budding in intestinal organoids:
 *   1. Stem cells + Paneth cells at organoid base divide to expand tissue
 *   2. Transit-amplifying cells proliferate rapidly in crypt zone
 *   3. Differential adhesion causes cell sorting (stem cluster at base)
 *   4. ECM degradation at expansion sites allows outward budding
 *   5. Basement membrane constrains organoid shape
 *   6. ECM fiber alignment guides the direction of bud extension
 * 
 * References:
 *   - Sato et al. (2009) Nature 459: 262-265 (organoid culture)
 *   - Painter (2009) J. Math. Biol. (ECM contact guidance model)
 *   - Dunn et al. (2012) J. Theor. Biol. (crypt mechanics)
 */

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "SimulationTime.hpp"

// Cell-based 3D
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedEventHandler.hpp"

// Cell cycle models
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"

// Cell types and states
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

// Modifiers
#include "VolumeTrackingModifier.hpp"

// Writers
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

// Project-specific: TissueMorphology forces
#include "BasementMembraneForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"
#include "ECMFieldWriter3d.hpp"

// For progress tracking
#include "AbstractCellBasedSimulationModifier.hpp"
#include <iomanip>

// PETSc setup
#include "PetscSetupAndFinalize.hpp"


/**
 * Custom simulation modifier for progress tracking
 */
template<unsigned DIM>
class ProgressTrackingModifier : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    double mEndTime;
    unsigned mUpdateInterval;
    unsigned mLastUpdateStep;
    std::chrono::steady_clock::time_point mStartWallTime;
    
public:
    ProgressTrackingModifier(double endTime, unsigned updateInterval = 100)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mEndTime(endTime),
          mUpdateInterval(updateInterval),
          mLastUpdateStep(0)
    {
    }
    
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        double current_time = SimulationTime::Instance()->GetTime();
        unsigned current_step = SimulationTime::Instance()->GetTimeStepsElapsed();
        
        // Update progress every mUpdateInterval steps
        if (current_step - mLastUpdateStep >= mUpdateInterval)
        {
            double percentage = (current_time / mEndTime) * 100.0;
            unsigned num_cells = rCellPopulation.GetNumRealCells();
            
            // Calculate real elapsed wall-clock time
            auto now = std::chrono::steady_clock::now();
            double elapsed_seconds = std::chrono::duration<double>(now - mStartWallTime).count();
            int elapsed_min = static_cast<int>(elapsed_seconds) / 60;
            int elapsed_sec = static_cast<int>(elapsed_seconds) % 60;
            
            // Estimate remaining time
            double rate = current_time / elapsed_seconds;  // simulated hrs per real second
            double remaining_sim_time = mEndTime - current_time;
            double eta_seconds = (rate > 0) ? remaining_sim_time / rate : 0;
            int eta_min = static_cast<int>(eta_seconds) / 60;
            int eta_sec = static_cast<int>(eta_seconds) % 60;
            
            // Create progress bar
            int bar_width = 40;
            int pos = static_cast<int>(bar_width * current_time / mEndTime);
            
            std::cout << "\r[";
            for (int i = 0; i < bar_width; ++i)
            {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << std::fixed << std::setprecision(1) << percentage << "% "
                      << "| Sim: " << std::setprecision(1) << current_time << "/" << mEndTime << " hrs "
                      << "| Cells: " << num_cells
                      << " | Elapsed: " << elapsed_min << "m" << std::setw(2) << std::setfill('0') << elapsed_sec << "s"
                      << " ETA: " << eta_min << "m" << std::setw(2) << eta_sec << "s" << std::setfill(' ')
                      << std::flush;
            
            mLastUpdateStep = current_step;
        }
    }
    
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
    {
        mStartWallTime = std::chrono::steady_clock::now();
        std::cout << "\nProgress:" << std::endl;
    }
    
    virtual void UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        // Calculate total elapsed time
        auto now = std::chrono::steady_clock::now();
        double elapsed_seconds = std::chrono::duration<double>(now - mStartWallTime).count();
        int elapsed_min = static_cast<int>(elapsed_seconds) / 60;
        int elapsed_sec = static_cast<int>(elapsed_seconds) % 60;
        
        // Print final completed progress bar
        std::cout << "\r[";
        for (int i = 0; i < 40; ++i) std::cout << "=";
        std::cout << "] 100.0% | Sim: " << mEndTime << "/" << mEndTime << " hrs "
                  << "| Cells: " << rCellPopulation.GetNumRealCells()
                  << " | Total: " << elapsed_min << "m" << std::setw(2) << std::setfill('0') << elapsed_sec << "s"
                  << std::setfill(' ') << std::endl;
    }
    
    virtual void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t<EndTime>" << mEndTime << "</EndTime>\n";
        *rParamsFile << "\t\t<UpdateInterval>" << mUpdateInterval << "</UpdateInterval>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};


class Test3dCryptOrganoid : public AbstractCellBasedTestSuite
{
private:

    /**
     * Create a spherical shell of nodes representing the initial organoid.
     * 
     * Intestinal organoids start as a hollow spherical cyst of epithelial cells.
     * Cells are placed on a spherical shell with some thickness, representing
     * a single-layer epithelium.
     * 
     * @param num_cells  Number of cells to place
     * @param radius     Radius of the spherical shell (µm)
     * @param shell_thickness  Thickness of the cell layer (µm)
     * @return vector of Node pointers
     */
    std::vector<Node<3>*> CreateSphericalShell(unsigned num_cells,
                                                double radius,
                                                double shell_thickness)
    {
        std::vector<Node<3>*> nodes;
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            // Fibonacci sphere for approximately uniform distribution
            double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;
            double theta = 2.0 * M_PI * i / golden_ratio;
            double phi = acos(1.0 - 2.0 * (i + 0.5) / num_cells);
            
            // Add small perturbation within shell thickness
            double r = radius + (p_gen->ranf() - 0.5) * shell_thickness;
            
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);
            
            nodes.push_back(new Node<3>(i, false, x, y, z));
        }
        
        return nodes;
    }
    
    /**
     * Assign cell types based on position on the spheroid.
     * 
     * Intestinal organoid cell types:
     *   - Stem cells:     At the "south pole" (bottom, z < -0.5*R) — form crypt base
     *   - Paneth cells:   Interspersed with stem cells — niche support
     *   - Transit-amplifying: Middle zone — rapidly dividing
     *   - Enterocytes:    Top half — differentiated, non-dividing
     * 
     * For the cell-centre model:
     *   - Stem/Paneth → StemCellProliferativeType
     *   - Transit-amplifying → TransitCellProliferativeType
     *   - Enterocytes → DifferentiatedCellProliferativeType
     */
    void AssignCellTypes(std::vector<CellPtr>& cells,
                         NodesOnlyMesh<3>& mesh,
                         double radius)
    {
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        
        for (unsigned i = 0; i < cells.size(); i++)
        {
            c_vector<double, 3> location = mesh.GetNode(i)->rGetLocation();
            double z = location[2];
            double z_frac = z / radius;  // Normalized z in [-1, 1]
            
            if (z_frac < -0.5)
            {
                // Crypt base: stem & Paneth cells
                cells[i]->SetCellProliferativeType(p_stem_type);
                cells[i]->GetCellData()->SetItem("cell_type_id", 0.0);  // Stem
                cells[i]->GetCellData()->SetItem("is_paneth", (i % 3 == 0) ? 1.0 : 0.0);
                // Higher basement membrane degradation at crypt base
                cells[i]->GetCellData()->SetItem("ecm_degradation_factor", 2.0);
            }
            else if (z_frac < 0.3)
            {
                // Transit-amplifying zone
                cells[i]->SetCellProliferativeType(p_transit_type);
                cells[i]->GetCellData()->SetItem("cell_type_id", 1.0);  // TA
                cells[i]->GetCellData()->SetItem("is_paneth", 0.0);
                cells[i]->GetCellData()->SetItem("ecm_degradation_factor", 1.0);
            }
            else
            {
                // Villus-like: differentiated enterocytes
                cells[i]->SetCellProliferativeType(p_diff_type);
                cells[i]->GetCellData()->SetItem("cell_type_id", 2.0);  // Differentiated
                cells[i]->GetCellData()->SetItem("is_paneth", 0.0);
                cells[i]->GetCellData()->SetItem("ecm_degradation_factor", 0.5);
            }
            
            // All cells share these
            cells[i]->GetCellData()->SetItem("volume", 1.0);  // For VolumeTrackingModifier
            cells[i]->GetCellData()->SetItem("basement_membrane_stiffness", 3.0);
        }
    }

public:

    /**
     * Main test: 3D intestinal crypt organoid with dynamic ECM.
     * 
     * Simulates 7 days of organoid development starting from a spherical
     * cyst of ~100 cells embedded in Matrigel-like ECM.
     */
    void Test3dIntestinalCryptOrganoid()
    {
        // ================================================================
        // SIMULATION PARAMETERS
        // ================================================================
        
        // Geometry
        const unsigned num_initial_cells = 100;
        const double organoid_radius = 25.0;       // µm — typical organoid starting radius
        const double shell_thickness = 3.0;         // µm — epithelial thickness
        const double interaction_cutoff = 15.0;     // µm — cell-cell interaction range
        
        // Time
        const double dt = 0.01;                    // hours (3.6 seconds)
        const double end_time = 168.0;              // hours (7 days)
        const unsigned sampling_multiple = 20;     // Output every dt*50 = 0.25 hour
        
        // ECM
        const double ecm_domain_half = 80.0;        // µm — ECM domain: [-80, 80]^3
        const double ecm_grid_spacing = 10.0;        // µm per voxel
        
        // Forces
        const double spring_stiffness = 20.0;       // Cell-cell spring stiffness
        const double bm_stiffness = 5.0;            // Basement membrane stiffness
        const double bm_radius = 30.0;              // Initial BM radius (µm)
        const double bm_degradation_rate = 0.15;    // BM expansion rate (µm/hour)
        const double bm_max_radius = 80.0;          // Max BM radius (µm)
        const double ecm_base_speed = 0.3;          // ECM guidance speed (µm/hour)
        
        std::cout << "============================================" << std::endl;
        std::cout << "  3D Intestinal Crypt Organoid Simulation   " << std::endl;
        std::cout << "  Duration: 7 days | Cells: " << num_initial_cells << std::endl;
        std::cout << "  Framework: TissueMorphology (cell-centre)" << std::endl;
        std::cout << "  ECM: Dynamic 3D with fiber guidance      " << std::endl;
        std::cout << "============================================" << std::endl;
        
        // ================================================================
        // 1. CREATE SPHERICAL ORGANOID MESH
        // ================================================================
        
        std::vector<Node<3>*> nodes = CreateSphericalShell(
            num_initial_cells, organoid_radius, shell_thickness);
        
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, interaction_cutoff);
        
        // ================================================================
        // 2. CREATE CELLS WITH PROLIFERATION MODEL
        // ================================================================
        
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        
        for (unsigned i = 0; i < num_initial_cells; i++)
        {
            // Contact inhibition cell cycle: cells stop dividing when compressed
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(3);
            p_cycle->SetQuiescentVolumeFraction(0.7);
            p_cycle->SetEquilibriumVolume(1.0);
            
            CellPtr p_cell(new Cell(p_state, p_cycle));
            p_cell->SetCellProliferativeType(p_transit_type);
            
            // Random birth time to avoid synchronous division
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 18.0;
            p_cell->SetBirthTime(birth_time);
            p_cell->InitialiseCellCycleModel();
            
            cells.push_back(p_cell);
        }
        
        // ================================================================
        // 3. CREATE CELL POPULATION
        // ================================================================
        
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Assign cell types based on position (stem/TA/differentiated)
        AssignCellTypes(cells, mesh, organoid_radius);
        
        // Add output writers
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        
        // ================================================================
        // 4. CREATE 3D DYNAMIC ECM FIELD
        // ================================================================
        
        // Radial ECM: fibers point outward from organoid center,
        // mimicking Matrigel/collagen gel surrounding the organoid
        boost::shared_ptr<DynamicECMField3d> p_ecm_field(
            new DynamicECMField3d("radial",           // fiber pattern
                                  ecm_grid_spacing,   // voxel size
                                  -ecm_domain_half, ecm_domain_half,
                                  -ecm_domain_half, ecm_domain_half,
                                  -ecm_domain_half, ecm_domain_half));
        
        // ECM dynamics parameters
        p_ecm_field->SetDegradationRate(0.002);   // MMPs degrade ECM
        p_ecm_field->SetRemodelingRate(0.05);      // Traction remodels fibers
        p_ecm_field->SetDepositionRate(0.0003);    // Slow ECM deposition
        
        // ================================================================
        // 5. SET UP SIMULATION
        // ================================================================
        
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("CryptOrganoid3d/7day_dynamicECM");
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);
        
        // ================================================================
        // 6. ADD FORCES
        // ================================================================
        
        // 6a. Cell-cell spring forces (repulsion + adhesion)
        //     F_ij = mu * (|r_ij| - l_0) * r_hat_ij
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring_force->SetMeinekeSpringGrowthDuration(1.0);
        p_spring_force->SetCutOffLength(interaction_cutoff);
        simulator.AddForce(p_spring_force);
        
        // 6b. Basement membrane force
        //     F_BM = -k_BM * (|r - r_c| - R) * r_hat  for |r - r_c| > R
        //     With time-dependent radius: R(t) = R_0 + k_deg * t
        MAKE_PTR(BasementMembraneForce<3>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetTargetRadius(bm_radius);
        p_bm_force->EnableEcmDegradation(bm_degradation_rate, bm_max_radius);
        simulator.AddForce(p_bm_force);
        
        // 6c. 3D ECM contact guidance force
        //     F_ECM = v0 * s * sqrt(rho) * n_migration
        //     where n_migration combines fiber-aligned, perpendicular, and random components
        MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm_force);
        p_ecm_force->SetECMField(p_ecm_field);
        p_ecm_force->SetBaseSpeed(ecm_base_speed);
        p_ecm_force->SetECMSensitivity(1.0);
        p_ecm_force->SetEnableDegradation(true);
        p_ecm_force->SetEnableRemodeling(true);
        p_ecm_force->SetEnableDeposition(false);
        simulator.AddForce(p_ecm_force);
        
        // ================================================================
        // 7. ADD SIMULATION MODIFIERS
        // ================================================================
        
        // Volume tracking for contact inhibition cell cycle
        MAKE_PTR(VolumeTrackingModifier<3>, p_vol_modifier);
        simulator.AddSimulationModifier(p_vol_modifier);
        
        // Progress tracking modifier - updates every 100 timesteps
        boost::shared_ptr<ProgressTrackingModifier<3>> p_progress_modifier(
            new ProgressTrackingModifier<3>(end_time, 100));
        simulator.AddSimulationModifier(p_progress_modifier);
        
        // ECM field writer: output every sampling_multiple timesteps
        // TEMPORARILY DISABLED FOR DEBUGGING
        // boost::shared_ptr<ECMFieldWriter3d> p_ecm_writer(
        //     new ECMFieldWriter3d(p_ecm_field, sampling_multiple));
        // simulator.AddSimulationModifier(p_ecm_writer);
        
        // ================================================================
        // 8. RUN SIMULATION
        // ================================================================
        
        std::cout << "Progress update interval: every 100 timesteps (every " 
                  << (100 * dt) << " hours)" << std::endl;
        std::cout << "\nStarting simulation..." << std::endl;
        std::cout << "End time: " << end_time << " hours (" << end_time / 24.0 << " days)" << std::endl;
        std::cout << "Timestep: " << dt << " hours" << std::endl;
        std::cout << "Total steps: " << (unsigned)(end_time / dt) << std::endl;
        std::cout << "Output frames: " << (unsigned)(end_time / (dt * sampling_multiple)) << std::endl;
        
        simulator.Solve();
        
        // ================================================================
        // 9. POST-SIMULATION ANALYSIS
        // ================================================================
        
        NodeBasedCellPopulation<3>& r_population
            = dynamic_cast<NodeBasedCellPopulation<3>&>(simulator.rGetCellPopulation());
        
        unsigned final_cell_count = r_population.GetNumRealCells();
        c_vector<double, 3> centroid = r_population.GetCentroidOfCellPopulation();
        
        // Count cell types
        unsigned stem_count = 0, ta_count = 0, diff_count = 0;
        double max_distance = 0.0;
        
        for (AbstractCellPopulation<3>::Iterator cell_iter = r_population.Begin();
             cell_iter != r_population.End();
             ++cell_iter)
        {
            double type_id = cell_iter->GetCellData()->GetItem("cell_type_id");
            if (type_id < 0.5) stem_count++;
            else if (type_id < 1.5) ta_count++;
            else diff_count++;
            
            // Track maximum radial distance (organoid extent)
            unsigned node_idx = r_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> pos = r_population.GetNode(node_idx)->rGetLocation();
            double dist = norm_2(pos - centroid);
            if (dist > max_distance) max_distance = dist;
        }
        
        // ================================================================
        // 10. REPORT RESULTS
        // ================================================================
        
        std::cout << "\n============================================" << std::endl;
        std::cout << "  SIMULATION COMPLETE" << std::endl;
        std::cout << "============================================" << std::endl;
        std::cout << "Final cell count:     " << final_cell_count << std::endl;
        std::cout << "  Stem cells:         " << stem_count << std::endl;
        std::cout << "  Transit-amplifying:  " << ta_count << std::endl;
        std::cout << "  Differentiated:     " << diff_count << std::endl;
        std::cout << "Centroid: (" << centroid[0] << ", "
                  << centroid[1] << ", " << centroid[2] << ")" << std::endl;
        std::cout << "Max radial distance:  " << max_distance << " µm" << std::endl;
        std::cout << "Expansion ratio:      " << max_distance / organoid_radius << "x" << std::endl;
        std::cout << "\nOutput in: testoutput/CryptOrganoid3d/7day_dynamicECM/" << std::endl;
        std::cout << "  - results_*.vtu    : cell positions (ParaView)" << std::endl;
        std::cout << "  - ecm3d_*.vti      : ECM field (ParaView)" << std::endl;
        std::cout << "  - ecm3d_results.pvd: ECM time series" << std::endl;
        std::cout << "============================================" << std::endl;
        
        // Crypt counting placeholder — to be implemented later with
        // curvature analysis of the organoid surface
        std::cout << "\n[NOTE] Crypt structure counting not yet implemented." << std::endl;
        std::cout << "Future: analyze surface curvature to detect bud protrusions." << std::endl;
        
        // Basic assertions
        TS_ASSERT_LESS_THAN(num_initial_cells, final_cell_count);  // Cells should have divided
        TS_ASSERT_LESS_THAN(0u, stem_count);   // Some stem cells should remain
    }
};

#endif /* TEST3DCRYPTORGANOID_HPP_ */
