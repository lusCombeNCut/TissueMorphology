/*
 * Test2dInvasiveFront.hpp
 * 
 * Reproduces Section 3.1 "Invasive Cellular Front Pushing into ECM" from Metzcar et al. 2025
 * Based on Painter 2009's computational study of tumor invasion into structured ECM
 * 
 * Uses the three-factor morphogenesis model:
 * 1. BasementMembraneForce - ECM constraint with orientation-dependent resistance
 * 2. DifferentialAdhesionForce - Cell-cell adhesion (uniform for invasion front)
 * 3. ApicalConstrictionForce - Not used in this scenario (all cells non-apical)
 * 
 * Tests four ECM orientation scenarios:
 * 1. Random - baseline diffusive invasion
 * 2. Parallel - ECM aligned with domain edge (restrictive)
 * 3. Perpendicular - ECM oriented perpendicular to edge (fast invasion)
 * 4. Mixed - combination of parallel and perpendicular regions
 */

#ifndef TEST2DINVASIVEFRONT_HPP_
#define TEST2DINVASIVEFRONT_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

// Include our three-factor model forces
#include "ECMContactGuidanceForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "ApicalConstrictionForce.hpp"

#include "PetscSetupAndFinalize.hpp"

class Test2dInvasiveFront : public AbstractCellBasedTestSuite
{
public:
    
    /**
     * Helper method to run invasion simulation with specified ECM orientation
     * 
     * Uses three-factor model:
     * 1. ECMContactGuidanceForce - structured ECM guides cell migration
     * 2. DifferentialAdhesionForce - cell-cell adhesion (uniform for this scenario)
     * 3. ApicalConstrictionForce - not used (no apical cells in invasion front)
     * 
     * Parameters from Metzcar 2025:
     * - Domain: 600 µm (x) × 1000 µm (y)
     * - Cell influx: 30 cells every 180 minutes at bottom edge  
     * - Simulation: 5 days (but we'll do 1 day for testing)
     * - Cell-cell adhesion: 0.4
     * - Repulsion: 25.0 µm
     * - Base speed: 1.25 µm/min
     * - ECM sensitivity: s = 1.0
     * - Anisotropy: a = 1.0
     */
    void RunInvasionSimulation(std::string ecmType)
    {
        // Create initial cell population at bottom edge
        std::vector<Node<2>*> nodes;
        unsigned num_initial_cells = 30;
        
        // Place initial cells along bottom edge (y ≈ 0-50 µm)
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        for (unsigned i = 0; i < num_initial_cells; i++)
        {
            double x = 50.0 + p_gen->ranf() * 500.0;  // Spread across domain (leaving margins)
            double y = p_gen->ranf() * 50.0;          // Near bottom edge
            nodes.push_back(new Node<2>(i, false, x, y));
        }
        
        // Create mesh with appropriate interaction distance
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 50.0);  // 50 µm interaction cutoff
        
        // Create cells - all uniform "invasive front" type
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        
        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            // Cell cycle model with no division (just migration)
            UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetMaxTransitGenerations(0);  // No proliferation
            
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(0.0);
            
            // Mark all cells as "invasive_front" type (for differential adhesion)
            p_cell->GetCellData()->SetItem("cell_type", 0.0);  // All same type
            
            cells.push_back(p_cell);
        }
        
        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("InvasiveFront2d/" + ecmType);
        simulator.SetSamplingTimestepMultiple(60);  // Output every 60 minutes
        simulator.SetDt(1.0);  // dt = 1 minute
        
        // Run for 1 day initially (1440 minutes) - can extend to 5 days later
        simulator.SetEndTime(1440.0);
        
        // Force 1: Cell-cell adhesion (uniform for all cells in invasion front)
        // Using GeneralisedLinearSpringForce since all cells are the same type
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_adhesion_force);
        p_adhesion_force->SetCutOffLength(50.0);
        
        // From paper: cell-cell adhesion = 0.4, repulsion = 25.0
        p_adhesion_force->SetMeinekeSpringStiffness(0.4);
        simulator.AddForce(p_adhesion_force);
        
        // Force 2: ECM Contact Guidance (ECMContactGuidanceForce)
        // This is the key force - structured ECM guides migration
        MAKE_PTR(ECMContactGuidanceForce<2>, p_ecm_force);
        p_ecm_force->SetECMOrientationType(ecmType);
        p_ecm_force->SetBaseSpeed(1.25);       // µm/min
        p_ecm_force->SetECMSensitivity(1.0);   // Full ECM response
        p_ecm_force->SetAnisotropy(1.0);       // Highly anisotropic ECM
        p_ecm_force->SetDomainWidth(600.0);    // For mixed pattern
        simulator.AddForce(p_ecm_force);
        
        // Force 3: ApicalConstrictionForce not used in this scenario
        // (No apical cells in invasive front)
        
        // TODO: Add cell birth at boundary every 180 minutes
        // This requires implementing a custom AbstractCellBasedSimulationModifier
        // For now, we start with 30 cells and watch them migrate
        
        // Run simulation
        simulator.Solve();
        
        // Clean up
        for (unsigned i = 0; i < nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    
    void TestRandomECMInvasion()
    {
        // Random ECM: baseline diffusive invasion pattern
        RunInvasionSimulation("random");
    }
    
    void TestParallelECMInvasion()
    {
        // Parallel ECM: slowest invasion (fibers block perpendicular movement)
        RunInvasionSimulation("parallel");
    }
    
    void TestPerpendicularECMInvasion()
    {
        // Perpendicular ECM: fastest invasion (fibers facilitate penetration)
        RunInvasionSimulation("perpendicular");
    }
    
    void TestMixedECMInvasion()
    {
        // Mixed ECM: heterogeneous pattern (fast in center, slow on sides)
        RunInvasionSimulation("mixed");
    }
};

#endif /* TEST2DINVASIVEFRONT_HPP_ */
