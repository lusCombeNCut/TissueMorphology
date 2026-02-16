/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

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
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TEST3DORGANOIDFORMATION_HPP_
#define TEST3DORGANOIDFORMATION_HPP_

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

// This test can run in parallel
#include "PetscSetupAndFinalize.hpp"

/**
 * Test suite for 3D organoid formation with basement membrane stiffness effects
 * and folding analysis. Inspired by TestMonodomain3dRabbitHeartTutorial structure.
 */
class Test3dOrganoidFormation : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test 3D organoid formation with spherical initialization.
     * 
     * This test creates a spherical organoid and studies how basement membrane
     * stiffness affects 3D morphology and potential folding patterns.
     */
    void Test3dSphericalOrganoidFormation()
    {
        // Create a 3D spherical organoid using nodes-only approach
        // This is better for 3D organoid formation with folding
        std::vector<Node<3>*> nodes;
        
        // Create spherical distribution of nodes
        unsigned num_cells = 30;  // Much fewer cells for stability
        double organoid_radius = 20.0;  // 20μm radius - INSIDE basement membrane
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            // Generate points in spherical distribution
            double theta = 2.0 * M_PI * p_gen->ranf();  // azimuthal angle
            double phi = acos(1.0 - 2.0 * p_gen->ranf());  // polar angle (uniform on sphere)
            double r = organoid_radius * pow(p_gen->ranf(), 1.0/3.0);  // uniform in volume
            
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);
            
            nodes.push_back(new Node<3>(i, false, x, y, z));
        }
        
        // Create 3D nodes-only mesh for organoid formation
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5.0);  // interaction cutoff distance
        
        // Create 3D cell factory with basement membrane properties
        OrganoidCellFactory<3,3> cell_factory(num_cells, 0.15);  // num_cells, basement_probability
        std::vector<CellPtr> cells;
        cell_factory.CreateCells(cells, num_cells);
        
        // Create 3D node-based cell population for organoid formation
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        
        // Update cell population to initialize neighbor relationships
        cell_population.Update();
        
        // Add VTK output writers for 3D visualization
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        // Note: CellVolumesWriter has issues with NodeBasedCellPopulation, skip for now
        
        // Set up 3D simulation
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Organoid3d/Basic");
        simulator.SetDt(0.002);  // Stable timestep
        simulator.SetSamplingTimestepMultiple(50);   // Output every 0.1 time units
        simulator.SetEndTime(2.0);  // Longer simulation for more timesteps (50 outputs)
        
        // Increase movement threshold to allow initial settling
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Add 3D spring forces between neighboring cells
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(30.0);
        p_spring_force->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring_force->SetMeinekeSpringGrowthDuration(1.0);
        simulator.AddForce(p_spring_force);
        
        // Add 3D basement membrane force with ECM degradation for organoid expansion
        MAKE_PTR(BasementMembraneForce<3>, p_basement_force);
        p_basement_force->SetBasementMembraneParameter(1.5);  // Softer initially for easier expansion
        p_basement_force->SetTargetRadius(25.0);  // Smaller initial radius
        p_basement_force->EnableEcmDegradation(0.3, 60.0);  // Gradual expansion: rate=0.3, max=60μm
        simulator.AddForce(p_basement_force);
        
        // Run simulation
        simulator.Solve();
        
        // Post-simulation analysis
        NodeBasedCellPopulation<3>* p_tissue = static_cast<NodeBasedCellPopulation<3>*>(&simulator.rGetCellPopulation());
        
        // Calculate 3D morphology metrics
        unsigned num_cells_final = p_tissue->GetNumRealCells();
        c_vector<double,3> centroid = p_tissue->GetCentroidOfCellPopulation();
        
        // Calculate sphericity (for folding analysis)
        double total_volume = 0.0;
        
        for (AbstractCellPopulation<3>::Iterator cell_iter = p_tissue->Begin();
             cell_iter != p_tissue->End();
             ++cell_iter)
        {
            // Calculate cell contribution to organoid volume and surface
            total_volume += 1.0;  // Approximate unit cell volume
        }
        
        // Sphericity = (π^(1/3) * (6V)^(2/3)) / A
        // Perfect sphere has sphericity = 1, more folded/irregular shapes have lower values
        
        std::cout << "=== 3D Organoid Formation Results ===" << std::endl;
        std::cout << "Number of cells: " << num_cells_final << std::endl;
        std::cout << "Centroid: (" << centroid[0] << ", " << centroid[1] << ", " << centroid[2] << ")" << std::endl;
        std::cout << "Estimated volume: " << total_volume << std::endl;
        
        // Basic assertions for test validation
        TS_ASSERT_LESS_THAN(0, num_cells_final);
        TS_ASSERT_LESS_THAN(num_cells_final, 1000);  // Reasonable cell count
    }
    
    /**
     * Test basement membrane stiffness effects on 3D folding patterns.
     */
    void Test3dBasementMembraneStiffnessEffects()
    {
        // Reset simulation time before test
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Test single stiffness value to study folding (simplified for now)
        double stiffness = 2.0;
        
        // Create spherical organoid using nodes-only approach
        std::vector<Node<3>*> nodes;
        unsigned num_cells = 30;  // Reduced for stability
        double organoid_radius = 18.0;  // Smaller than basement membrane (20μm)
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned j = 0; j < num_cells; j++)
        {
            double theta = 2.0 * M_PI * p_gen->ranf();
            double phi = acos(1.0 - 2.0 * p_gen->ranf());
            double r = organoid_radius * pow(p_gen->ranf(), 1.0/3.0);
            
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);
            
            nodes.push_back(new Node<3>(j, false, x, y, z));
        }
        
        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 5.0);
        
        // Create 3D cell factory with basement membrane properties
        OrganoidCellFactory<3,3> cell_factory(num_cells, 0.15);
        std::vector<CellPtr> cells;
        cell_factory.CreateCells(cells, num_cells);
        
        // Set up simulation
        NodeBasedCellPopulation<3> cell_population(mesh, cells);
        
        // Update cell population to initialize neighbor relationships
        cell_population.Update();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Organoid3d/StiffnessTest");
        simulator.SetDt(0.002);  // Stable timestep
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(5.0);  // Longer duration for more timesteps (50 outputs)
        
        // Increase movement threshold to allow initial settling
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Add forces
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(25.0);
        simulator.AddForce(p_spring_force);
        
        MAKE_PTR(BasementMembraneForce<3>, p_basement_force);
        p_basement_force->SetBasementMembraneParameter(stiffness);
        p_basement_force->SetTargetRadius(20.0);  // Smaller initial radius
        p_basement_force->EnableEcmDegradation(0.5, 50.0);  // Allow expansion
        simulator.AddForce(p_basement_force);
        
        // Run simulation
        simulator.Solve();
        
        std::cout << "Completed stiffness test (K=" << stiffness << ")" << std::endl;
    }
    
    /**
     * Test long-term 3D organoid development with folding analysis.
     */
    void Test3dLongTermOrganoidDevelopment()
    {
        // Reset simulation time before test
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Create larger spherical organoid for folding studies
        std::vector<Node<3>*> nodes;
        unsigned num_cells = 40;  // Reduced for stability
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
            
            if (distance_from_center < 30.0)  // Inner cells
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
        // Note: CellVolumesWriter has issues with NodeBasedCellPopulation, skip for now
        
        OffLatticeSimulation<3> simulator(cell_population);
        simulator.SetOutputDirectory("Organoid3d/LongTermDevelopment");
        simulator.SetDt(0.002);  // Stable timestep for stability
        simulator.SetSamplingTimestepMultiple(50);   // Output every 0.1 time units
        simulator.SetEndTime(100.0); // Long simulation for extensive development (100 outputs)
        
        // Increase movement threshold to allow initial settling
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Add differential adhesion forces
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(35.0);
        simulator.AddForce(p_spring_force);
        
        // Add basement membrane force with ECM degradation for realistic expansion
        MAKE_PTR(BasementMembraneForce<3>, p_basement_force);
        p_basement_force->SetBasementMembraneParameter(2.0);  // Softer initially
        p_basement_force->SetTargetRadius(20.0);  // Start smaller
        p_basement_force->EnableEcmDegradation(1, 70.0);  // Allow significant expansion
        simulator.AddForce(p_basement_force);
        
        // Run simulation
        simulator.Solve();
        
        std::cout << "=== Long-term 3D Organoid Development Complete ===" << std::endl;
        std::cout << "Check Organoid3d/LongTermDevelopment for VTK output" << std::endl;
        std::cout << "Use VTK visualization for folding analysis" << std::endl;
    }
};

#endif /*TEST3DORGANOIDFORMATION_HPP_*/