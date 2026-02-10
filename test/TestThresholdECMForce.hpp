#ifndef TESTTHRESHOLDECMFORCE_HPP_
#define TESTTHRESHOLDECMFORCE_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <string>

// Core Chaste includes
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"

// Cell-based includes for 3D
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
#include "VtkMeshWriter.hpp"

// Project-specific includes
#include "OrganoidCellFactory.hpp"
#include "BasementMembraneForce.hpp"

// Enable parallel execution
#include "PetscSetupAndFinalize.hpp"

/**
 * Test suite for 3D organoid development with threshold ECM force
 * Studies basement membrane degradation and organoid expansion
 */
class TestThresholdECMForce : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test 3D organoid development with threshold-based ECM degradation.
     * Studies how basement membrane stiffness changes affect organoid expansion.
     */
    void TestThresholdBasedECMDegradation()
    {
        // Reset simulation time before test
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Create spherical organoid for long-term studies
        std::vector<Node<3>*> nodes;
        unsigned num_cells = 50;  // More cells for interesting dynamics
        double organoid_radius = 18.0;  // Smaller than basement membrane (20μm)
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            double theta = 2.0 * M_PI * p_gen->ranf();
            double phi = acos(1.0 - 2.0 * p_gen->ranf());
            double r = organoid_radius * pow(p_gen->ranf(), 1.0/3.0);
            
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);
            
            nodes.push_back(new Node<3>(i, false, x, y, z));
        }
        
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5.0);
        
        // Create heterogeneous cell population
        OrganoidCellFactory<3,3> cell_factory(num_cells, 0.25);
        std::vector<CellPtr> cells;
        cell_factory.CreateCells(cells, num_cells);
        
        // Add cell type heterogeneity for realistic folding
        for (unsigned i = 0; i < cells.size(); i++)
        {
            // Inner cells (core) vs outer cells (surface) have different properties
            c_vector<double,3> cell_location = mesh.GetNode(i)->rGetLocation();
            double distance_from_center = norm_2(cell_location);
            
            if (distance_from_center < 12.0)  // Inner cells
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
        
        // Set up long-term simulation
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        
        // Update cell population to initialize neighbor relationships
        cell_population.Update();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Organoid3d/ThresholdECMForce");
        simulator.SetDt(0.005);  // Larger timestep for faster execution
        simulator.SetSamplingTimestepMultiple(40);   // Output every 0.2 time units
        simulator.SetEndTime(50.0); // Very long simulation (100 outputs)
        
        // Increase movement threshold to allow initial settling
        cell_population.SetAbsoluteMovementThreshold(500.0);
        
        // Add differential adhesion forces
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(35.0);
        simulator.AddForce(p_spring_force);
        
        // Add basement membrane force with ECM degradation for realistic expansion
        MAKE_PTR(BasementMembraneForce<3>, p_basement_force);
        p_basement_force->SetBasementMembraneParameter(2.0);  // Softer initially
        p_basement_force->SetTargetRadius(20.0);  // Start smaller
        p_basement_force->EnableEcmDegradation(0.4, 80.0);  // Allow significant expansion
        simulator.AddForce(p_basement_force);
        
        // Run simulation
        simulator.Solve();
        
        std::cout << "=== Threshold ECM Force Test Complete ===" << std::endl;
        std::cout << "Check Organoid3d/ThresholdECMForce for VTK output" << std::endl;
        std::cout << "Final cell count: " << cell_population.GetNumRealCells() << std::endl;
    }
};

#endif /*TESTTHRESHOLDECMFORCE_HPP_*/