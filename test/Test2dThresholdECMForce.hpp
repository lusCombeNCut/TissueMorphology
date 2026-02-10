#ifndef TEST2DTHRESHOLDECMFORCE_HPP_
#define TEST2DTHRESHOLDECMFORCE_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <string>

// Core Chaste includes
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"

// Cell-based includes for 2D
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "NodesOnlyMesh.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "SimulationTime.hpp"

// VTK output writers
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Project-specific includes
#include "OrganoidCellFactory.hpp"
#include "BasementMembraneForce.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Test suite for 2D organoid development with threshold ECM force
 * Studies basement membrane degradation and organoid expansion in 2D
 */
class Test2dThresholdECMForce : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test 2D organoid development with threshold-based ECM degradation.
     * Studies how basement membrane stiffness changes affect organoid expansion in 2D.
     */
    void Test2dThresholdBasedECMDegradation()
    {
        // Reset simulation time before test
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Create circular organoid for 2D studies
        std::vector<Node<2>*> nodes;
        unsigned num_cells = 50;  // Moderate number for 2D
        double organoid_radius = 5.0;  // Smaller than basement membrane (6μm)
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            double theta = 2.0 * M_PI * p_gen->ranf();
            double r = organoid_radius * sqrt(p_gen->ranf());  // Uniform in disk
            
            double x = r * cos(theta);
            double y = r * sin(theta);
            
            nodes.push_back(new Node<2>(i, false, x, y));
        }
        
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        
        // Create heterogeneous cell population
        OrganoidCellFactory<2,2> cell_factory(num_cells, 0.25);
        std::vector<CellPtr> cells;
        cell_factory.CreateCells(cells, num_cells);
        
        // Add cell type heterogeneity for realistic behavior
        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Inner cells (core) vs outer cells (surface) have different properties
            c_vector<double,2> cell_location = mesh.GetNode(i)->rGetLocation();
            double distance_from_center = norm_2(cell_location);
            
            if (distance_from_center < 3.0)  // Inner cells
            {
                cells[i]->GetCellData()->SetItem("cell_type", 0.0);  // Core cells
                cells[i]->GetCellData()->SetItem("basement_membrane_stiffness", 1.0);
            }
            else  // Outer cells
            {
                cells[i]->GetCellData()->SetItem("cell_type", 1.0);  // Surface cells
                cells[i]->GetCellData()->SetItem("basement_membrane_stiffness", 5.0);
            }
        }
        
        // Set up 2D simulation
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        
        // Update cell population to initialize neighbor relationships
        cell_population.Update();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Organoid2d/ThresholdECMForce");
        simulator.SetDt(0.005);  // Small timestep for stability
        simulator.SetSamplingTimestepMultiple(40);   // Output every 0.2 time units
        simulator.SetEndTime(20.0); // Long simulation (100 outputs)
        
        // Increase movement threshold to allow initial settling
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Add differential adhesion forces
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(35.0);
        simulator.AddForce(p_spring_force);
        
        // Add basement membrane force with ECM degradation for realistic expansion
        MAKE_PTR(BasementMembraneForce<2>, p_basement_force);
        p_basement_force->SetBasementMembraneParameter(2.0);  // Softer initially
        p_basement_force->SetTargetRadius(6.0);  // Start smaller
        p_basement_force->EnableEcmDegradation(0.4, 20.0);  // Allow significant expansion
        simulator.AddForce(p_basement_force);
        
        // Run simulation
        simulator.Solve();
        
        std::cout << "=== 2D Threshold ECM Force Test Complete ===" << std::endl;
        std::cout << "Check Organoid2d/ThresholdECMForce for VTK output" << std::endl;
        std::cout << "Final cell count: " << cell_population.GetNumRealCells() << std::endl;
    }
};

#endif /*TEST2DTHRESHOLDECMFORCE_HPP_*/
