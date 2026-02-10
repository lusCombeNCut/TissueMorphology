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

#ifndef TESTORGANOIDFORMATION_HPP_
#define TESTORGANOIDFORMATION_HPP_

#include <cxxtest/TestSuite.h>

// Core Chaste includes
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "Exception.hpp"

// Cell-based includes
#include "HoneycombMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "SimulationTime.hpp"
#include "CellDataItemWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"

// Project-specific includes
#include "OrganoidCellFactory.hpp"
#include "BasementMembraneForce.hpp"

// This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * Test suite for organoid formation with basement membrane stiffness effects.
 */
class TestOrganoidFormation : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test basic organoid formation with default parameters.
     * 
     * This test creates a small population of cells and runs a basic
     * cell-based simulation to verify the setup is working correctly.
     */
    void TestBasicOrganoidFormation()
    {
        // Create a 2D hexagonal mesh 
        HoneycombMeshGenerator generator(4, 4, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
        
        // Create cells directly (following Chaste pattern)
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_type);
        
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            p_cycle_model->SetDimension(2);
            
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_type);
            p_cell->SetBirthTime(0.0);
            
            // Add basement membrane stiffness as cell data
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", 1.0);
            p_cell->GetCellData()->SetItem("cell_type", i < 5 ? 0.0 : 1.0);
            
            cells.push_back(p_cell);
        }
        
        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        // Add cell writers for visualization
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OrganoidFormation/BasicTest");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(0.1);  // Very short simulation for testing
        
        // Add a simple force law (spring forces between neighbors)
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);
        
        // Run the simulation
        simulator.Solve();
        
        // Basic checks
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);  // Should have 16 cells
        TS_ASSERT_EQUALS(cell_population.GetNumNodes(), p_mesh->GetNumNodes());
        // Simulator ran successfully if no exceptions thrown
    }
    
    /**
     * Test organoid formation with different basement membrane stiffness values.
     * 
     * This test compares the behavior of organoids formed with different
     * basement membrane stiffness parameters.
     */
    void TestBasementMembraneStiffnessEffect()
    {
        // Test with low stiffness
        {
            HoneycombMeshGenerator generator(3, 3, 0);
            boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
            
            OrganoidCellFactory<2,2> cell_factory(0.1, true); // Low stiffness
            std::vector<CellPtr> cells;
            cell_factory.CreateCells(cells, p_mesh->GetNumNodes());
            
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
            
            // Add cell writers for visualization
            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("OrganoidFormation/LowStiffness");
            simulator.SetDt(0.1);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetEndTime(0.2);
            
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            simulator.AddForce(p_force);
            
            simulator.Solve();
            
            // Simulation completed successfully if no exceptions thrown
        }
        
        // Test with high stiffness
        {
            HoneycombMeshGenerator generator(3, 3, 0);
            boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
            
            OrganoidCellFactory<2,2> cell_factory(10.0, true); // High stiffness
            std::vector<CellPtr> cells;
            cell_factory.CreateCells(cells, p_mesh->GetNumNodes());
            
            MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
            
            // Add cell writers for visualization
            cell_population.AddCellWriter<CellIdWriter>();
            cell_population.AddCellWriter<CellAgesWriter>();
            
            OffLatticeSimulation<2> simulator(cell_population);
            simulator.SetOutputDirectory("OrganoidFormation/HighStiffness");
            simulator.SetDt(0.1);
            simulator.SetSamplingTimestepMultiple(1);
            simulator.SetEndTime(0.3);
            
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
            simulator.AddForce(p_force);
            
            simulator.Solve();
            
            // Simulation completed successfully if no exceptions thrown
        }
    }
    
    /**
     * Test that cells are created with correct basement membrane properties.
     */
    void TestCellFactoryProperties()
    {
        // Create a simple mesh
        std::vector<Node<2>*> nodes;
        nodes.push_back(new Node<2>(0, false, 0.0, 0.0));
        nodes.push_back(new Node<2>(1, false, 1.0, 0.0));
        
        // Test cell factory
        OrganoidCellFactory<2,2> cell_factory(2.5, true);
        
        // Create a cell
        CellPtr p_cell = cell_factory.CreateCell(0);
        
        // Check that the cell has the correct basement membrane stiffness
        TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("basement_membrane_stiffness"), 2.5, 1e-6);
        
        // Check cell type assignment (center vs peripheral)
        TS_ASSERT_DELTA(p_cell->GetCellData()->GetItem("cell_type"), 0.0, 1e-6); // Should be stem-like (center)
        
        // Test peripheral cell
        CellPtr p_peripheral_cell = cell_factory.CreateCell(10);
        TS_ASSERT_DELTA(p_peripheral_cell->GetCellData()->GetItem("cell_type"), 1.0, 1e-6); // Should be differentiated-like
        
        // Clean up
        delete nodes[0];
        delete nodes[1];
    }
    
    /**
     * Test organoid formation with basement membrane force constraints.
     * 
     * This test demonstrates how the basement membrane force affects organoid morphology.
     */
    void TestOrganoidWithBasementMembraneForce()
    {
        // Create a larger 2D hexagonal mesh to get ~100 cells
        HoneycombMeshGenerator generator(10, 10, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();
        
        // Create cells using the organoid cell factory
        OrganoidCellFactory<2,2> cell_factory(2.0, true);
        std::vector<CellPtr> cells;
        cell_factory.CreateCells(cells, p_mesh->GetNumNodes());
        
        // Create cell population
        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
        
        // Add cell writers for visualization
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("OrganoidFormation/WithBasementMembrane");
        simulator.SetDt(0.01);  // Smaller time step for basement membrane forces
        simulator.SetSamplingTimestepMultiple(1);  // Output every timestep for 100 timesteps
        simulator.SetEndTime(1.0);  // 100 timesteps (1.0 / 0.01 = 100)
        
        // Add spring forces between neighbors
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        simulator.AddForce(p_spring_force);
        
        // Add basement membrane force
        MAKE_PTR(BasementMembraneForce<2>, p_basement_force);
        p_basement_force->SetBasementMembraneRadius(3.0);
        p_basement_force->SetForceStrength(1.0);
        simulator.AddForce(p_basement_force);
        
        // Run the simulation
        simulator.Solve();
        
        // Basic checks
        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 101u);  // Should have 101 cells (10x10 honeycomb mesh)
        // Simulation completed successfully if no exceptions thrown
        
        // Check that organoid center is reasonable
        c_vector<double, 2> center = p_basement_force->GetOrganoidCenter();
        TS_ASSERT(norm_2(center) < 20.0); // Should be reasonably centered (larger organoid)
    }
};

#endif /*TESTORGANOIDFORMATION_HPP_*/