#ifndef TEST3DVERTEXECMINTERACTION_HPP_
#define TEST3DVERTEXECMINTERACTION_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <string>
#include <vector>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"
#include "SimulationTime.hpp"

// OrganoidChaste: 3D vertex model infrastructure
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "FiniteThicknessHoneycombVertexMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"

// OrganoidChaste: vertex forces and modifiers
#include "SurfaceTensionForce.hpp"
#include "GeometricalTargetVolumeModifier.hpp"

// Cell cycle models
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellsGenerator.hpp"

// TissueMorphology: ECM forces
#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"
#include "ECMFieldWriter3d.hpp"

// PETSc setup
#include "PetscSetupAndFinalize.hpp"

class Test3dVertexECMInteraction : public AbstractCellBasedTestSuite
{
public:
    void TestVertexECMAlignment()
    {
        // ================================================================
        // PARAMETERS
        // ================================================================
        const unsigned num_across = 3;
        const unsigned num_up = 3;
        const double height = 1.0;
        
        // Surface tension
        const double gamma_apical = 0.1;
        const double gamma_basal = 0.1;
        const double gamma_lateral = 0.1;
        
        // ECM
        const double ecm_grid_spacing = 1.0;
        const double ecm_domain_size = 10.0;
        const double ecm_base_speed = 1.0;
        
        // Time
        const double dt = 0.01;
        const double end_time = 0.5; // Short run
        
        // ================================================================
        // 1. CREATE MESH
        // ================================================================
        FiniteThicknessHoneycombVertexMeshGenerator generator(num_across, num_up, false, 0.01, 0.001, height);
        MutableMonolayerVertexMesh<3, 3>* p_mesh = generator.GetMesh();
        
        // ================================================================
        // 2. CREATE CELLS
        // ================================================================
        std::vector<CellPtr> cells;
        // Use FixedG1GenerationalCellCycleModel to satisfy GeometricalTargetVolumeModifier
        CellsGenerator<FixedG1GenerationalCellCycleModel, 3> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh->GetNumElements());
        
        // Set long G1 duration to prevent division
        for (auto& p_cell : cells)
        {
            static_cast<FixedG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel())->SetStemCellG1Duration(1000.0);
            static_cast<FixedG1GenerationalCellCycleModel*>(p_cell->GetCellCycleModel())->SetTransitCellG1Duration(1000.0);
        }
        
        // ================================================================
        // 3. CREATE POPULATION
        // ================================================================
        MonolayerVertexBasedCellPopulation<3> cell_population(*p_mesh, cells);
        cell_population.SetDoInitialVolumeRelaxation(true);
        
        // ================================================================
        // 4. CREATE ECM FIELD
        // ================================================================
        // "radial" pattern is well-defined in DynamicECMField3d
        boost::shared_ptr<DynamicECMField3d> p_ecm_field(
            new DynamicECMField3d("radial",
                                  ecm_grid_spacing,
                                  -ecm_domain_size, ecm_domain_size,
                                  -ecm_domain_size, ecm_domain_size,
                                  -ecm_domain_size, ecm_domain_size));
        
        // ================================================================
        // 5. SET UP SIMULATION
        // ================================================================
        FiniteThicknessSimulation3d simulator(cell_population);
        simulator.SetOutputDirectory("Test3dVertexECMInteraction");
        simulator.SetEndTime(end_time);
        simulator.SetDt(dt);
        
        // ================================================================
        // 6. ADD FORCES
        // ================================================================
        MAKE_PTR(SurfaceTensionForce<3>, p_tension_force);
        p_tension_force->CreateSurfaceTensionParametersForCells(
            gamma_apical, gamma_basal, gamma_lateral, p_mesh);
        p_tension_force->SetSimulationInstance(&simulator);
        simulator.AddForce(p_tension_force);
        
        MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm_force);
        p_ecm_force->SetECMField(p_ecm_field);
        p_ecm_force->SetBaseSpeed(ecm_base_speed);
        p_ecm_force->SetECMSensitivity(1.0);
        p_ecm_force->SetEnableDegradation(false);
        p_ecm_force->SetEnableRemodeling(false);
        simulator.AddForce(p_ecm_force);
        
        // ================================================================
        // 7. ADD MODIFIERS
        // ================================================================
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_vol_modifier, (&cell_population));
        p_vol_modifier->SetGrowthDuration(0.0); // Fixes the crash
        p_vol_modifier->SetReferenceTargetVolume(0.5 * sqrt(3.0) * height);
        simulator.AddSimulationModifier(p_vol_modifier);
        
        // ================================================================
        // 8. RUN
        // ================================================================
        // Initial positions
        std::vector<c_vector<double, 3>> initial_positions;
        for (unsigned i = 0; i < cell_population.GetNumRealCells(); ++i)
        {
            initial_positions.push_back(cell_population.GetLocationOfCellCentre(cell_population.GetCellUsingLocationIndex(i)));
        }

        simulator.Solve();
        
        // ================================================================
        // 9. VERIFY
        // ================================================================
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), num_across * num_up);
        
        // Check that cells moved (due to ECM force)
        bool any_movement = false;
        for (unsigned i = 0; i < cell_population.GetNumRealCells(); ++i)
        {
            c_vector<double, 3> final_pos = cell_population.GetLocationOfCellCentre(cell_population.GetCellUsingLocationIndex(i));
            if (norm_2(final_pos - initial_positions[i]) > 1e-6)
            {
                any_movement = true;
                break;
            }
        }
        TS_ASSERT(any_movement);

        // Check that CellData was updated by the force
        CellPtr p_first_cell = *cell_population.Begin();
        TS_ASSERT_THROWS_NOTHING(p_first_cell->GetCellData()->GetItem("ecm_density"));
        TS_ASSERT_THROWS_NOTHING(p_first_cell->GetCellData()->GetItem("ecm_fiber_x"));
    }
};

#endif /* TEST3DVERTEXECMINTERACTION_HPP_ */
