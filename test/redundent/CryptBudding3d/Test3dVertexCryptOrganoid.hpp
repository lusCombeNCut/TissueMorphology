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

#ifndef TEST3DVERTEXCRYPTORGANOID_HPP_
#define TEST3DVERTEXCRYPTORGANOID_HPP_

/**
 * Test3dVertexCryptOrganoid.hpp
 * 
 * Simulates a 3D intestinal crypt organoid using a VERTEX MODEL from
 * the OrganoidChaste project combined with TissueMorphology ECM forces.
 * 
 * Unlike the node-based version (Test3dCryptOrganoid), this uses a
 * finite-thickness monolayer vertex mesh where each cell is a 3D polyhedron
 * with distinct apical, basal, and lateral faces. This enables:
 *   - Epithelial buckling and folding (crypt invagination)
 *   - Apical constriction driving bud formation
 *   - T1 transitions (cell neighbour rearrangements)
 *   - Realistic cell shape changes during morphogenesis
 *   - Lumen pressure effects
 * 
 * Forces:
 *   - SurfaceTensionForce: apical/basal/lateral surface tensions (OrganoidChaste)
 *   - GeometricalTargetVolumeModifier: volume homeostasis (OrganoidChaste)
 *   - BasementMembraneForce: BM confinement with degradation (TissueMorphology)
 *   - DynamicECMContactGuidanceForce3d: ECM fiber guidance (TissueMorphology)
 * 
 * The vertex model naturally captures buckling instabilities that arise when
 * differential growth in a confined epithelium exceeds a critical threshold,
 * which is the primary mechanism for crypt formation in intestinal organoids.
 * 
 * References:
 *   - Drozdowski & Schwarz (2025) OrganoidChaste (finite-thickness vertex model)
 *   - Hannezo et al. (2011) PRL (buckling instability in epithelia)
 *   - Sato et al. (2009) Nature 459: 262-265 (organoid culture)
 *   - Painter (2009) J. Math. Biol. (ECM contact guidance model)
 */

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <string>
#include <vector>
#include <ctime>
#include <sstream>
#include <iomanip>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "SimulationTime.hpp"

// OrganoidChaste: 3D vertex model infrastructure
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"

// OrganoidChaste: vertex forces and modifiers
#include "SurfaceTensionForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"

// OrganoidChaste: writers
#include "CellThicknessWriter.hpp"
#include "FaceTypeWriter.hpp"

// Cell cycle models
#include "NoCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"

// Cell types and states
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellsGenerator.hpp"

// Standard Chaste writers
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"

// TissueMorphology: ECM forces
#include "BasementMembraneForce.hpp"
#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"
#include "ECMFieldWriter3d.hpp"

// PETSc setup
#include "PetscSetupAndFinalize.hpp"


class Test3dVertexCryptOrganoid : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test: 3D vertex-based organoid — relaxation only (no growth).
     * 
     * Creates a spherical monolayer vertex mesh and lets it relax under
     * surface tension + basement membrane forces. This verifies the
     * combined force framework works before adding cell division.
     */
    void TestVertexOrganoidRelaxation()
    {
        // ================================================================
        // PARAMETERS
        // ================================================================
        const unsigned num_cells = 100;
        const double sphere_radius = 10.0;  // Inner radius of sphere
        
        // Surface tension parameters (dimensionless, OrganoidChaste convention)
        const double gamma_apical = 0.85;   // Apical surface tension
        const double gamma_basal = 0.85;    // Basal surface tension  
        const double gamma_lateral = 0.7;   // Lateral surface tension
        
        // Basement membrane
        const double bm_stiffness = 2.0;
        const double bm_radius = 12.0;      // Slightly larger than sphere
        
        // Time
        const double dt = 0.001;
        const double relax_time = 5.0;      // Short relaxation
        const unsigned sampling_multiple = 50;
        
        std::cout << "============================================" << std::endl;
        std::cout << "  3D Vertex Organoid — Relaxation Test      " << std::endl;
        std::cout << "  Cells: " << num_cells << " | Radius: " << sphere_radius << std::endl;
        std::cout << "  Model: Finite-thickness vertex monolayer  " << std::endl;
        std::cout << "============================================" << std::endl;
        
        // ================================================================
        // 1. CREATE SPHERICAL VERTEX MESH
        // ================================================================
        
        // Compute cell height from surface tension parameters
        // (OrganoidChaste convention: h ~ f(gamma_a, gamma_b))
        double height = 2.0 / 3.0 / sqrt(3.0)
                        * cbrt((9.0 / 2.0) * (9.0 / 2.0))
                        * cbrt((gamma_apical + gamma_basal) / gamma_lateral
                               * (gamma_apical + gamma_basal) / gamma_lateral)
                        * 1.0;
        
        // T1 transition threshold length
        double t1_length = 0.66 / cbrt(3.0 * 3.0 * (1.0 + gamma_lateral) / gamma_lateral);
        
        std::cout << "Cell height: " << height << std::endl;
        std::cout << "T1 threshold: " << t1_length << std::endl;
        
        FiniteThicknessRandomizedSphereMeshGenerator generator(
            num_cells, t1_length, 0.001, height, sphere_radius);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        
        // Enable protorosette handling for T1 transitions
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells);
        
        // ================================================================
        // 2. CREATE CELLS (no proliferation for relaxation)
        // ================================================================
        
        std::vector<CellPtr> cells;
        CellsGenerator<NoCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
        
        // ================================================================
        // 3. CREATE VERTEX-BASED CELL POPULATION
        // ================================================================
        
        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);
        
        // Add writers
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();
        
        // ================================================================
        // 4. SET UP SIMULATION (FiniteThicknessSimulation3d)
        // ================================================================
        
        // Generate timestamp for output
        std::time_t now = std::time(nullptr);
        std::tm* ltm = std::localtime(&now);
        std::stringstream timestamp;
        timestamp << std::setfill('0')
                  << std::setw(4) << (1900 + ltm->tm_year)
                  << std::setw(2) << (1 + ltm->tm_mon)
                  << std::setw(2) << ltm->tm_mday << "_"
                  << std::setw(2) << ltm->tm_hour
                  << std::setw(2) << ltm->tm_min
                  << std::setw(2) << ltm->tm_sec;
        
        std::string output_dir = "VertexCryptOrganoid/relax_" + timestamp.str();
        
        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetEndTime(relax_time);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetDt(dt);
        
        std::cout << "\nOutput: testoutput/" << output_dir << "/" << std::endl;
        
        // ================================================================
        // 5. ADD FORCES
        // ================================================================
        
        // 5a. Surface tension force (vertex model mechanics)
        //     Drives cell shape via apical/basal/lateral tensions
        MAKE_PTR(SurfaceTensionForce<3>, p_tension_force);
        p_tension_force->CreateSurfaceTensionParametersForCells(
            gamma_apical, gamma_basal, gamma_lateral, p_mesh);
        p_tension_force->SetSimulatedAnnealingParameters(0.0, 190.0, 10.0);
        p_tension_force->SetSimulationInstance(&simulator);
        simulator.AddForce(p_tension_force);
        
        // 5b. Basement membrane force (TissueMorphology)
        //     Confines cells to spherical region
        MAKE_PTR(BasementMembraneForce<3>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetTargetRadius(bm_radius);
        simulator.AddForce(p_bm_force);
        
        // ================================================================
        // 6. ADD MODIFIERS
        // ================================================================
        
        // Geometrical target volume: maintains cell volume homeostasis
        double outer_radius = sphere_radius + height;
        double avg_cell_volume = 4.0 / 3.0 * M_PI
                                 * (outer_radius * outer_radius * outer_radius
                                    - sphere_radius * sphere_radius * sphere_radius)
                                 / num_cells;
        
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_vol_modifier,
                      (&cell_population));
        p_vol_modifier->SetGrowthDuration(0.0);
        p_vol_modifier->SetT1AdaptationDuration(0.100);
        p_vol_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_vol_modifier);
        
        // ================================================================
        // 7. RUN RELAXATION
        // ================================================================
        
        std::cout << "\nRelaxing vertex mesh..." << std::endl;
        simulator.Solve();
        
        std::cout << "\n============================================" << std::endl;
        std::cout << "  RELAXATION COMPLETE" << std::endl;
        std::cout << "  Final cells: " << cell_population.GetNumRealCells() << std::endl;
        std::cout << "============================================" << std::endl;
        
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_cells);
    }
    
    /**
     * Main test: 3D vertex organoid with growth, ECM, and buckling.
     * 
     * Uses the OrganoidChaste finite-thickness vertex model for cell
     * mechanics, combined with TissueMorphology ECM forces. Cell division
     * in a confined organoid drives buckling → crypt formation.
     * 
     * Protocol:
     *   Phase 1: Relax mesh geometry (no growth, no ECM)
     *   Phase 2: Enable growth + T1 transitions + ECM guidance
     */
    void TestVertexOrganoidWithGrowthAndECM()
    {
        // ================================================================
        // PARAMETERS
        // ================================================================
        
        // Geometry
        const unsigned num_cells = 200;
        const double sphere_radius = 10.0;
        
        // Surface tension (dimensionless)
        const double gamma_apical = 0.85;
        const double gamma_basal = 0.85;
        const double gamma_lateral = 0.7;
        
        // Basement membrane (TissueMorphology)
        const double bm_stiffness = 3.0;
        const double bm_degradation_rate = 0.002;  // Slow expansion
        const double bm_max_radius = 25.0;
        
        // ECM field (TissueMorphology)
        const double ecm_domain_half = 40.0;
        const double ecm_grid_spacing = 5.0;
        const double ecm_base_speed = 0.1;
        
        // Time
        const double dt_relax = 0.001;
        const double dt_grow = 0.006;
        const double relax_time = 5.0;
        const double grow_time = 100.0;         // Total simulation time
        const unsigned sampling_relax = 20;
        const unsigned sampling_grow = 100;
        
        // Cell cycle
        const double stem_g1_duration = 25.0;
        const double transit_g1_duration = 15.0;
        
        std::cout << "============================================" << std::endl;
        std::cout << "  3D Vertex Crypt Organoid Simulation       " << std::endl;
        std::cout << "  Cells: " << num_cells << " | Model: vertex monolayer" << std::endl;
        std::cout << "  Phase 1: Relaxation (" << relax_time << " time units)" << std::endl;
        std::cout << "  Phase 2: Growth (" << grow_time << " time units)" << std::endl;
        std::cout << "  ECM: Dynamic 3D with fiber guidance       " << std::endl;
        std::cout << "  BM: Degrading (rate=" << bm_degradation_rate << ")" << std::endl;
        std::cout << "============================================" << std::endl;
        
        // ================================================================
        // 1. CREATE SPHERICAL VERTEX MESH
        // ================================================================
        
        double height = 2.0 / 3.0 / sqrt(3.0)
                        * cbrt((9.0 / 2.0) * (9.0 / 2.0))
                        * cbrt((gamma_apical + gamma_basal) / gamma_lateral
                               * (gamma_apical + gamma_basal) / gamma_lateral)
                        * 1.0;
        
        double t1_length = 0.66 / cbrt(3.0 * 3.0 * (1.0 + gamma_lateral) / gamma_lateral);
        
        std::cout << "\nCell height: " << height << std::endl;
        std::cout << "Inner radius: " << sphere_radius << std::endl;
        std::cout << "Outer radius: " << sphere_radius + height << std::endl;
        
        FiniteThicknessRandomizedSphereMeshGenerator generator(
            num_cells, t1_length, 0.001, height, sphere_radius);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        
        p_mesh->SetProtorosetteFormationProbability(1.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        
        TS_ASSERT_EQUALS(p_mesh->GetNumElements(), num_cells);
        
        // ================================================================
        // 2. CREATE CELLS WITH PROLIFERATION
        // ================================================================
        
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel;
            p_model->SetDimension(3);
            
            CellPtr p_cell(new Cell(p_state, p_model));
            
            // Assign cell types based on position on sphere
            // Use the element centroid z-coordinate
            c_vector<double, 3> centroid = p_mesh->GetCentroidOfElement(i);
            double z_frac = centroid[2] / (sphere_radius + height / 2.0);
            
            if (z_frac < -0.5)
            {
                // Crypt base: stem cells (slow cycle)
                p_cell->SetCellProliferativeType(p_stem_type);
                p_model->SetStemCellG1Duration(stem_g1_duration);
                p_model->SetTransitCellG1Duration(transit_g1_duration);
                p_model->SetMaxTransitGenerations(2);
            }
            else if (z_frac < 0.3)
            {
                // Transit-amplifying zone (faster cycle)
                p_cell->SetCellProliferativeType(p_transit_type);
                p_model->SetStemCellG1Duration(stem_g1_duration);
                p_model->SetTransitCellG1Duration(transit_g1_duration);
                p_model->SetMaxTransitGenerations(2);
            }
            else
            {
                // Differentiated (no division)
                p_cell->SetCellProliferativeType(p_diff_type);
                p_model->SetStemCellG1Duration(1e6);  // Effectively no division
                p_model->SetTransitCellG1Duration(1e6);
                p_model->SetMaxTransitGenerations(0);
            }
            
            p_model->SetSDuration(25.0);
            p_model->SetG2Duration(0.001);
            p_model->SetMDuration(0.001);
            
            // Random birth time to desynchronize
            double birth_time = -RandomNumberGenerator::Instance()->ranf() * 10.0;
            p_cell->SetBirthTime(birth_time);
            
            cells.push_back(p_cell);
        }
        
        // Count initial cell types
        unsigned initial_stem = 0, initial_ta = 0, initial_diff = 0;
        for (unsigned i = 0; i < cells.size(); i++)
        {
            if (cells[i]->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                initial_stem++;
            else if (cells[i]->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                initial_ta++;
            else
                initial_diff++;
        }
        std::cout << "\nInitial cell types:" << std::endl;
        std::cout << "  Stem: " << initial_stem << std::endl;
        std::cout << "  Transit-amplifying: " << initial_ta << std::endl;
        std::cout << "  Differentiated: " << initial_diff << std::endl;
        
        // ================================================================
        // 3. CREATE CELL POPULATION
        // ================================================================
        
        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellRearrangementLocations(false);
        cell_population.SetRestrictVertexMovementBoolean(false);
        cell_population.SetDoInitialVolumeRelaxation(true);
        
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellThicknessWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddFaceWriter<FaceTypeWriter>();
        
        // ================================================================
        // 4. CREATE ECM FIELD (TissueMorphology)
        // ================================================================
        
        boost::shared_ptr<DynamicECMField3d> p_ecm_field(
            new DynamicECMField3d("radial",
                                  ecm_grid_spacing,
                                  -ecm_domain_half, ecm_domain_half,
                                  -ecm_domain_half, ecm_domain_half,
                                  -ecm_domain_half, ecm_domain_half));
        
        p_ecm_field->SetDegradationRate(0.002);
        p_ecm_field->SetRemodelingRate(0.05);
        p_ecm_field->SetDepositionRate(0.0003);
        
        // ================================================================
        // 5. PHASE 1: RELAXATION (no growth, no ECM)
        // ================================================================
        
        std::time_t now = std::time(nullptr);
        std::tm* ltm = std::localtime(&now);
        std::stringstream timestamp;
        timestamp << std::setfill('0')
                  << std::setw(4) << (1900 + ltm->tm_year)
                  << std::setw(2) << (1 + ltm->tm_mon)
                  << std::setw(2) << ltm->tm_mday << "_"
                  << std::setw(2) << ltm->tm_hour
                  << std::setw(2) << ltm->tm_min
                  << std::setw(2) << ltm->tm_sec;
        
        std::string output_dir = "VertexCryptOrganoid/growth_" + timestamp.str();
        
        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetEndTime(relax_time);
        simulator.SetSamplingTimestepMultiple(sampling_relax);
        simulator.SetDt(dt_relax);
        
        std::cout << "\nOutput: testoutput/" << output_dir << "/" << std::endl;
        
        // Surface tension force
        MAKE_PTR(SurfaceTensionForce<3>, p_tension_force);
        p_tension_force->CreateSurfaceTensionParametersForCells(
            gamma_apical, gamma_basal, gamma_lateral, p_mesh);
        p_tension_force->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension_force->SetSimulationInstance(&simulator);
        simulator.AddForce(p_tension_force);
        
        // Basement membrane (TissueMorphology) — no degradation during relaxation
        MAKE_PTR(BasementMembraneForce<3>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetTargetRadius(sphere_radius + height + 1.0);
        simulator.AddForce(p_bm_force);
        
        // Target volume modifier
        double outer_radius = sphere_radius + height;
        double avg_cell_volume = 4.0 / 3.0 * M_PI
                                 * (outer_radius * outer_radius * outer_radius
                                    - sphere_radius * sphere_radius * sphere_radius)
                                 / num_cells;
        
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_vol_modifier,
                      (&cell_population));
        p_vol_modifier->SetGrowthDuration(0.0);
        p_vol_modifier->SetT1AdaptationDuration(0.100);
        p_vol_modifier->SetReferenceTargetVolume(avg_cell_volume);
        simulator.AddSimulationModifier(p_vol_modifier);
        
        std::cout << "\n--- PHASE 1: Relaxation ---" << std::endl;
        simulator.Solve();
        std::cout << "Relaxation complete. Cells: "
                  << cell_population.GetNumRealCells() << std::endl;
        
        // ================================================================
        // 6. PHASE 2: GROWTH + ECM + BUCKLING
        // ================================================================
        
        std::cout << "\n--- PHASE 2: Growth + ECM + Buckling ---" << std::endl;
        
        simulator.SetEndTime(grow_time);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetSamplingTimestepMultiple(sampling_grow);
        simulator.SetDt(dt_grow);
        
        // Enable active T1 transitions (enables buckling)
        p_tension_force->SetSimulatedAnnealingParameters(0.003, 1900000.0, 1.0);
        p_tension_force->SetSimulationInstance(&simulator);
        p_tension_force->SetPerformActiveT1Swaps(true);
        p_tension_force->SetT1TransitionParameters(2.0, false);
        
        // Enable BM degradation
        p_bm_force->EnableEcmDegradation(bm_degradation_rate, bm_max_radius);
        
        // Add ECM guidance force (TissueMorphology)
        MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm_force);
        p_ecm_force->SetECMField(p_ecm_field);
        p_ecm_force->SetBaseSpeed(ecm_base_speed);
        p_ecm_force->SetECMSensitivity(1.0);
        p_ecm_force->SetEnableDegradation(true);
        p_ecm_force->SetEnableRemodeling(true);
        p_ecm_force->SetEnableDeposition(false);
        simulator.AddForce(p_ecm_force);
        
        // Add ECM field writer
        boost::shared_ptr<ECMFieldWriter3d> p_ecm_writer(
            new ECMFieldWriter3d(p_ecm_field, sampling_grow));
        simulator.AddSimulationModifier(p_ecm_writer);
        
        // Run growth phase
        simulator.Solve();
        
        // ================================================================
        // 7. POST-SIMULATION ANALYSIS
        // ================================================================
        
        unsigned final_cell_count = cell_population.GetNumRealCells();
        
        // Count cell types
        unsigned stem_count = 0, ta_count = 0, diff_count = 0;
        for (AbstractCellPopulation<3>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_iter->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                stem_count++;
            else if (cell_iter->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                ta_count++;
            else
                diff_count++;
        }
        
        std::cout << "\n============================================" << std::endl;
        std::cout << "  SIMULATION COMPLETE" << std::endl;
        std::cout << "============================================" << std::endl;
        std::cout << "Final cell count:     " << final_cell_count << std::endl;
        std::cout << "  Stem cells:         " << stem_count << std::endl;
        std::cout << "  Transit-amplifying: " << ta_count << std::endl;
        std::cout << "  Differentiated:     " << diff_count << std::endl;
        std::cout << "\nOutput in: testoutput/" << output_dir << "/" << std::endl;
        std::cout << "  - results_*.vtu     : vertex mesh (ParaView)" << std::endl;
        std::cout << "  - ecm3d_*.vti       : ECM field (ParaView)" << std::endl;
        std::cout << "============================================" << std::endl;
        
        // Assertions
        TS_ASSERT_LESS_THAN(num_cells, final_cell_count);  // Cells should have divided
    }
};

#endif /* TEST3DVERTEXCRYPTORGANOID_HPP_ */
