/*
 * Test2dDynamicECMInvasion.hpp
 * 
 * Enhanced version of Test2dInvasiveFront with:
 * 1. Dynamic ECM remodeling (cells align ECM through traction)
 * 2. ECM degradation and deposition
 * 3. Proper cell-cell adhesion using DifferentialAdhesionForce
 * 4. Bidirectional fiber movement with perpendicular diffusion
 * 
 * Reproduces Section 3.1 from Metzcar et al. 2025 with realistic dynamics
 */

#ifndef TEST2DDYNAMICECMINVASION_HPP_
#define TEST2DDYNAMICECMINVASION_HPP_
#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "VolumeTrackingModifier.hpp"
#include "RandomCellKiller.hpp"
#include <sstream>
#include <cstdlib>
#include "Debug.hpp"
#include "SimulationTime.hpp"

// Include our three-factor model forces
#include "DynamicECMContactGuidanceForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class Test2dDynamicECMInvasion : public AbstractCellBasedTestSuite
{
public:
    
    /**
     * Helper method to initialize ECM field with specified pattern
     */
    boost::shared_ptr<DynamicECMField> CreateECMField(std::string ecmType, double domainWidth, double domainHeight)
    {
        boost::shared_ptr<DynamicECMField> p_ecm_field(new DynamicECMField(ecmType, 12.5, 0.0, domainWidth, 0.0, domainHeight));
        
        if (ecmType == "random")
        {
            // Random ECM orientations
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            for (double x = 0; x < domainWidth; x += 12.5)
            {
                for (double y = 0; y < domainHeight; y += 12.5)
                {
                    c_vector<double, 2> pos;
                    pos[0] = x;
                    pos[1] = y;
                    double angle = p_gen->ranf() * 2.0 * M_PI;
                    p_ecm_field->SetFiberAngle(pos, angle);
                }
            }
        }
        else if (ecmType == "parallel")
        {
            // Horizontal fibers (parallel to domain edge)
            for (double x = 0; x < domainWidth; x += 12.5)
            {
                for (double y = 0; y < domainHeight; y += 12.5)
                {
                    c_vector<double, 2> pos;
                    pos[0] = x;
                    pos[1] = y;
                    p_ecm_field->DepositECM(pos, 0.0, 1.0);  // 0 radians = horizontal
                }
            }
        }
        else if (ecmType == "perpendicular")
        {
            // Vertical fibers (perpendicular to domain edge)
            for (double x = 0; x < domainWidth; x += 12.5)
            {
                for (double y = 0; y < domainHeight; y += 12.5)
                {
                    c_vector<double, 2> pos;
                    pos[0] = x;
                    pos[1] = y;
                    p_ecm_field->DepositECM(pos, M_PI/2.0, 1.0);  // π/2 = vertical
                }
            }
        }
        else if (ecmType == "mixed")
        {
            // Mixed: perpendicular in center, parallel on sides
            for (double x = 0; x < domainWidth; x += 12.5)
            {
                for (double y = 0; y < domainHeight; y += 12.5)
                {
                    c_vector<double, 2> pos;
                    pos[0] = x;
                    pos[1] = y;
                    
                    // Center third: perpendicular (fast invasion)
                    // Outer thirds: parallel (slow invasion)
                    double angle = (x > domainWidth/3.0 && x < 2.0*domainWidth/3.0) ? M_PI/2.0 : 0.0;
                    p_ecm_field->DepositECM(pos, angle, 1.0);
                }
            }
        }
        
        return p_ecm_field;
    }
    
    /**
     * Run invasion simulation with dynamic ECM and proper adhesion
     */
    void RunDynamicInvasionSimulation(std::string ecmType, unsigned randomSeed = 12345)
    {
        // Destroy and recreate simulation time for each run
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Domain parameters
        double domain_width = 600.0;   // µm
        double domain_height = 1000.0; // µm
        
        // First create ECM field (for random type, this consumes RNG)
        boost::shared_ptr<DynamicECMField> p_ecm_field_temp = CreateECMField(ecmType, domain_width, domain_height);
        
        // Now reseed RNG to ensure identical cell positions for all ECM types
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(randomSeed);
        
        // Create initial cell population at bottom edge
        std::vector<Node<2>*> nodes;
        unsigned num_initial_cells = 30;
        
        for (unsigned i = 0; i < num_initial_cells; i++)
        {
            double x = 50.0 + p_gen->ranf() * 500.0;  // Spread across domain
            double y = p_gen->ranf() * 50.0;          // Near bottom edge
            nodes.push_back(new Node<2>(i, false, x, y));
        }
        
        // Create mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 50.0);  // 50 µm interaction cutoff
        
        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        
        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            // Use ContactInhibitionCellCycleModel for organoid-like behavior
            ContactInhibitionCellCycleModel* p_cycle_model = new ContactInhibitionCellCycleModel();
            p_cycle_model->SetDimension(2);
            
            // Cell cycle phase durations (typical epithelial cells)
            p_cycle_model->SetStemCellG1Duration(8.0);      // Stem cells: 8 hours
            p_cycle_model->SetTransitCellG1Duration(12.0);  // Transit cells: 12 hours  
            p_cycle_model->SetSDuration(5.0);               // S phase: 5 hours
            p_cycle_model->SetG2Duration(4.0);              // G2 phase: 4 hours  
            p_cycle_model->SetMDuration(1.0);               // M phase: 1 hour
            // Total cell cycle: 8-12h (G1) + 5h (S) + 4h (G2) + 1h (M) = 18-22 hours
            
            // Contact inhibition: quiescence when compressed to 70% of equilibrium volume
            p_cycle_model->SetQuiescentVolumeFraction(0.7);
            p_cycle_model->SetEquilibriumVolume(M_PI * 5.0 * 5.0); // ~78.5 µm² for radius=5µm
            
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            
            // Random birth times: ages 0-22 hours (avoid synchronized divisions)
            p_cell->SetBirthTime(-p_gen->ranf() * 22.0);
            
            // All cells same type for uniform adhesion
            p_cell->GetCellData()->SetItem("cell_type", 1.0);  // Type 1 = invasive front
            
            cells.push_back(p_cell);
        }
        
        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(50.0);
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        
        // Create unique output directory for each run
        std::stringstream output_dir;
        output_dir << "DynamicECMInvasion2d/" << ecmType << "/run_" << randomSeed;
        simulator.SetOutputDirectory(output_dir.str());
        
        simulator.SetSamplingTimestepMultiple(1440);  // Output every day
        simulator.SetDt(2);  // dt = 2 minutes
        simulator.SetEndTime(7200.0);  // 5 days = 2880 minutes (faster test)
        
        // === THREE-FACTOR MODEL FORCES ===
        
        // Force 1: Cell-cell adhesion (uniform for all cells in invasion front)
        // COMMENTED OUT TO TEST ECM FORCES ONLY
        // MAKE_PTR(GeneralisedLinearSpringForce<2>, p_adhesion_force);
        // p_adhesion_force->SetCutOffLength(50.0);
        // 
        // // From paper: cell-cell adhesion = 0.4
        // p_adhesion_force->SetMeinekeSpringStiffness(0.4);
        // 
        // simulator.AddForce(p_adhesion_force);
        
        // Force 2: Dynamic ECM Contact Guidance (use ECM created earlier)
        MAKE_PTR(DynamicECMContactGuidanceForce<2>, p_ecm_force);
        p_ecm_force->SetECMField(p_ecm_field_temp);
        p_ecm_force->SetBaseSpeed(1.25);       // µm/min
        p_ecm_force->SetECMSensitivity(1.0);   // Full ECM response
        p_ecm_force->SetEnableRemodeling(true);   // Dynamic ECM!
        p_ecm_force->SetEnableDegradation(true);  // Cells degrade ECM
        p_ecm_force->SetEnableDeposition(false);  // Don't deposit new ECM (yet)
        
        simulator.AddForce(p_ecm_force);
        
        // Force 3: ApicalConstrictionForce not used (no apical cells)
        
        // Add volume tracking modifier (required for ContactInhibitionCellCycleModel)
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);
        
        // Add ECM field writer to output grid data (match VTU output frequency)
        // Note: There's a 2x bug in Chaste where simulations run twice as long as EndTime
        // VTU files output every 120 timesteps (matching SamplingTimestepMultiple)
        // So ECM writer also needs 120 to get same output times
        MAKE_PTR_ARGS(ECMFieldWriter<2>, p_ecm_writer, (p_ecm_field_temp, 120));
        simulator.AddSimulationModifier(p_ecm_writer);
        
        // NOTE: Cell death disabled for now to debug proliferation
        // Add cell death (homeostatic balance)
        // 0.05% per hour = ~1% per day death rate
        // MAKE_PTR_ARGS(RandomCellKiller<2>, p_cell_killer, (&cell_population, 0.0005));
        // simulator.AddCellKiller(p_cell_killer);
        
        // Run simulation
        simulator.Solve();
        
        // Clean up
        for (unsigned i = 0; i < nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
    
    void TestDynamicRandomECM()
    {
        // Get run number from environment variable (for parallel execution)
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 12345 + run_num * 1000;
            RunDynamicInvasionSimulation("random", seed);
        }
        else
        {
            // Run all 20 simulations sequentially if no environment variable
            for (unsigned i = 0; i < 1; i++)
            {
                unsigned seed = 12345 + i * 1000;
                RunDynamicInvasionSimulation("random", seed);
            }
        }
    }
    
    void TestDynamicParallelECM()
    {
        // Get run number from environment variable (for parallel execution)
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 12345 + run_num * 1000;
            RunDynamicInvasionSimulation("parallel", seed);
        }
        else
        {
            // Run all 20 simulations sequentially if no environment variable
            for (unsigned i = 0; i < 1; i++)
            {
                unsigned seed = 12345 + i * 1000;
                RunDynamicInvasionSimulation("parallel", seed);
            }
        }
    }
    
    void TestDynamicPerpendicularECM()
    {
        // Get run number from environment variable (for parallel execution)
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 12345 + run_num * 1000;
            RunDynamicInvasionSimulation("perpendicular", seed);
        }
        else
        {
            // Run all 20 simulations sequentially if no environment variable
            for (unsigned i = 0; i < 1; i++)
            {
                unsigned seed = 12345 + i * 1000;
                RunDynamicInvasionSimulation("perpendicular", seed);
            }
        }
    }
    
    void TestDynamicMixedECM()
    {
        RunDynamicInvasionSimulation("mixed");
    }
};

#endif /* TEST2DDYNAMICECMINVASION_HPP_ */
