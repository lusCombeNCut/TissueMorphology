#ifndef TEST2DAPICALCONSTRICTION_HPP_
#define TEST2DAPICALCONSTRICTION_HPP_

#include <cxxtest/TestSuite.h>

// Core Chaste includes
#include "AbstractCellBasedTestSuite.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

// Cell-based includes
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"

// Cell cycle and properties
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellData.hpp"

// Output writers
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"

// Project-specific includes
#include "ApicalConstrictionForce.hpp"

// This test runs sequentially
#include "FakePetscSetup.hpp"

/**
 * Test suite for apical constriction force in 2D.
 * Validates the basic mechanics of apical wedging and tissue invagination.
 */
class Test2dApicalConstriction : public AbstractCellBasedTestSuite
{
public:

    /**
     * Test basic apical constriction in a 2D tissue sheet.
     * 
     * Creates a monolayer with apical cells at the center,
     * applies constriction force, and observes invagination.
     */
    void TestBasicApicalConstriction()
    {
        // Reset simulation time
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Create a 2D tissue sheet (grid of cells)
        std::vector<Node<2>*> nodes;
        unsigned grid_width = 10;
        unsigned grid_height = 10;
        double cell_spacing = 1.0;
        
        for (unsigned i = 0; i < grid_width; i++)
        {
            for (unsigned j = 0; j < grid_height; j++)
            {
                unsigned index = i * grid_height + j;
                double x = i * cell_spacing;
                double y = j * cell_spacing;
                nodes.push_back(new Node<2>(index, false, x, y));
            }
        }
        
        // Create mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);
        
        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        
        for (unsigned i = 0; i < nodes.size(); i++)
        {
            UniformG1GenerationalCellCycleModel* p_model = new UniformG1GenerationalCellCycleModel();
            p_model->SetDimension(2);
            
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-10.0);  // Mature cells
            p_cell->InitialiseCellCycleModel();
            
            // Mark central cells as apical
            double x = nodes[i]->rGetLocation()[0];
            double y = nodes[i]->rGetLocation()[1];
            double center_x = (grid_width - 1) * cell_spacing / 2.0;
            double center_y = (grid_height - 1) * cell_spacing / 2.0;
            double distance_from_center = sqrt((x - center_x)*(x - center_x) + (y - center_y)*(y - center_y));
            
            if (distance_from_center < 2.0)  // Central 2-unit radius
            {
                p_cell->GetCellData()->SetItem("is_apical", 1.0);
            }
            else
            {
                p_cell->GetCellData()->SetItem("is_apical", 0.0);
            }
            
            cells.push_back(p_cell);
        }
        
        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.Update();
        
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ApicalConstriction2d/Basic");
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10.0);
        
        // Add spring force for cell-cell adhesion
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(30.0);
        simulator.AddForce(p_spring_force);
        
        // Add apical constriction force
        MAKE_PTR(ApicalConstrictionForce<2>, p_constriction_force);
        p_constriction_force->SetConstrictionStrength(10.0);
        p_constriction_force->SetTargetReduction(0.6);  // 60% area reduction
        simulator.AddForce(p_constriction_force);
        
        // Run simulation
        simulator.Solve();
        
        std::cout << "=== Apical Constriction Test Complete ===" << std::endl;
        std::cout << "Check ApicalConstriction2d/Basic for output" << std::endl;
        std::cout << "Expected: Central cells should form invagination pit" << std::endl;
    }
};

#endif /*TEST2DAPICALCONSTRICTION_HPP_*/
