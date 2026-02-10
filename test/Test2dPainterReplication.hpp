/*
 * Test2dPainterReplication.hpp
 * 
 * Reproduces aspects of Painter 2009's computational study:
 * - Constant influx of cells at domain boundary (30 cells every 180 minutes)
 * - Static ECM (no remodeling)
 * - Strong ECM contact guidance
 * - 5 day simulation duration
 * - Analyze 95th percentile of cell invasion
 * 
 * Reference: Painter 2009, Metzcar et al. 2025
 */

#ifndef TEST2DPAINTERREPLICATION_HPP_
#define TEST2DPAINTERREPLICATION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "RandomNumberGenerator.hpp"
#include <sstream>
#include <cstdlib>
#include "SimulationTime.hpp"

// Include our forces
#include "DynamicECMContactGuidanceForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class Test2dPainterReplication : public AbstractCellBasedTestSuite
{
private:
    /**
     * Create ECM field with specified orientation pattern
     */
    boost::shared_ptr<DynamicECMField> CreateStaticECMField(std::string ecmType, 
                                                             double width, double height)
    {
        // ECM grid spacing
        double grid_spacing = 25.0;  // 25 µm grid spacing
        
        // Create ECM field using same constructor as Test2dDynamicECMInvasion
        boost::shared_ptr<DynamicECMField> p_ecm_field(new DynamicECMField(ecmType, grid_spacing, 0.0, width, 0.0, height));
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        if (ecmType == "random")
        {
            // Random ECM orientations (re-seed to create different pattern)
            for (double x = 0; x < width; x += grid_spacing)
            {
                for (double y = 0; y < height; y += grid_spacing)
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
            for (double x = 0; x < width; x += grid_spacing)
            {
                for (double y = 0; y < height; y += grid_spacing)
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
            for (double x = 0; x < width; x += grid_spacing)
            {
                for (double y = 0; y < height; y += grid_spacing)
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
            for (double x = 0; x < width; x += grid_spacing)
            {
                for (double y = 0; y < height; y += grid_spacing)
                {
                    c_vector<double, 2> pos;
                    pos[0] = x;
                    pos[1] = y;
                    
                    // Center third: perpendicular (fast invasion)
                    // Outer thirds: parallel (slow invasion)
                    double angle = (x > width/3.0 && x < 2.0*width/3.0) ? M_PI/2.0 : 0.0;
                    p_ecm_field->DepositECM(pos, angle, 1.0);
                }
            }
        }
        
        return p_ecm_field;
    }
    
    /**
     * Create initial batch of cells at bottom edge
     */
    std::vector<Node<2>*> CreateInitialCells(unsigned num_cells, double domain_width, 
                                             unsigned& next_node_index)
    {
        std::vector<Node<2>*> nodes;
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (unsigned i = 0; i < num_cells; i++)
        {
            double x = 50.0 + p_gen->ranf() * 500.0;  // Spread across domain width
            double y = p_gen->ranf() * 50.0;          // Near bottom edge
            nodes.push_back(new Node<2>(next_node_index++, false, x, y));
        }
        
        return nodes;
    }
    
    /**
     * Custom cell killer that adds cells at domain boundary
     * (Simulates constant influx)
     */
    class CellInfluxModifier : public AbstractCellBasedSimulationModifier<2>
    {
    private:
        double mInfluxInterval;  // Minutes between influx events
        double mLastInfluxTime;
        unsigned mCellsPerInflux;
        double mDomainWidth;
        unsigned mNextNodeIndex;
        
    public:
        CellInfluxModifier()
            : AbstractCellBasedSimulationModifier<2>(),
              mInfluxInterval(180.0),
              mLastInfluxTime(0.0),
              mCellsPerInflux(30),
              mDomainWidth(600.0),
              mNextNodeIndex(1000000)
        {
        }
        
        CellInfluxModifier(double influxInterval, unsigned cellsPerInflux, double domainWidth)
            : AbstractCellBasedSimulationModifier<2>(),
              mInfluxInterval(influxInterval),
              mLastInfluxTime(0.0),
              mCellsPerInflux(cellsPerInflux),
              mDomainWidth(domainWidth),
              mNextNodeIndex(1000000)  // Start from high index to avoid conflicts
        {
        }
        
        void SetInfluxInterval(double interval) { mInfluxInterval = interval; }
        void SetCellsPerInflux(unsigned cells) { mCellsPerInflux = cells; }
        void SetNextNodeIndex(unsigned index) { mNextNodeIndex = index; }
        
        void UpdateAtEndOfTimeStep(AbstractCellPopulation<2,2>& rCellPopulation)
        {
            double current_time = SimulationTime::Instance()->GetTime();
            
            // Check if it's time to add new cells
            if (current_time - mLastInfluxTime >= mInfluxInterval)
            {
                mLastInfluxTime = current_time;
                
                // Add new cells at bottom edge
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                NodeBasedCellPopulation<2>* p_node_population = static_cast<NodeBasedCellPopulation<2>*>(&rCellPopulation);
                
                for (unsigned i = 0; i < mCellsPerInflux; i++)
                {
                    double x = 50.0 + p_gen->ranf() * 500.0;
                    double y = p_gen->ranf() * 50.0;
                    
                    // Create new node
                    Node<2>* p_new_node = new Node<2>(mNextNodeIndex, false, x, y);
                    p_new_node->ClearAppliedForce();  // Important for velocity output
                    
                    // Add node to mesh - this returns the index where it was added
                    unsigned new_node_index = p_node_population->rGetMesh().AddNode(p_new_node);
                    
                    // Create new cell with uniform cell cycle (no division)
                    MAKE_PTR(WildTypeCellMutationState, p_state);
                    MAKE_PTR(TransitCellProliferativeType, p_transit_type);
                    
                    UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
                    p_cycle_model->SetDimension(2);
                    p_cycle_model->SetMinCellCycleDuration(1e6);  // Effectively no division
                    p_cycle_model->SetMaxCellCycleDuration(1e6);
                    
                    CellPtr p_new_cell(new Cell(p_state, p_cycle_model));
                    p_new_cell->SetCellProliferativeType(p_transit_type);
                    p_new_cell->SetBirthTime(current_time);
                    p_new_cell->GetCellData()->SetItem("cell_type", 1.0);
                    p_new_cell->GetCellData()->SetItem("ecm_angle", 0.0);
                    p_new_cell->GetCellData()->SetItem("ecm_density", 0.0);
                    p_new_cell->GetCellData()->SetItem("ecm_orientation_x", 1.0);
                    p_new_cell->GetCellData()->SetItem("ecm_orientation_y", 0.0);
                    p_new_cell->GetCellData()->SetItem("migration_direction_x", 0.0);
                    p_new_cell->GetCellData()->SetItem("migration_direction_y", 1.0);
                    
                    // Manually add cell to population (mimicking AddCell but without division)
                    p_node_population->rGetCells().push_back(p_new_cell);
                    p_node_population->SetCellUsingLocationIndex(new_node_index, p_new_cell);
                    
                    mNextNodeIndex++;
                }
                
                // Update cell population after adding cells
                p_node_population->Update(false);  // Don't call CheckForStepSizeException
            }
        }
        
        void SetupSolve(AbstractCellPopulation<2,2>& rCellPopulation, std::string outputDirectory)
        {
            mLastInfluxTime = SimulationTime::Instance()->GetTime();
        }
        
        void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<2,2>& rCellPopulation) {}
        void OutputSimulationModifierParameters(out_stream& rParamsFile) {}
    };
    
public:
    void RunPainterSimulation(std::string ecmType, unsigned randomSeed = 12345)
    {
        // Destroy and recreate simulation time
        SimulationTime::Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
        
        // Domain parameters (from Metzcar 2025)
        double domain_width = 600.0;   // µm
        double domain_height = 1000.0; // µm
        
        // Create static ECM field
        boost::shared_ptr<DynamicECMField> p_ecm_field = CreateStaticECMField(ecmType, domain_width, domain_height);
        
        // Seed RNG
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        p_gen->Reseed(randomSeed);
        
        // Create initial cell population (30 cells at bottom edge)
        unsigned next_node_index = 0;
        std::vector<Node<2>*> initial_nodes = CreateInitialCells(30, domain_width, next_node_index);
        
        // Create mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(initial_nodes, 50.0);
        
        // Create cells with uniform cell cycle (no proliferation)
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);
        
        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetMinCellCycleDuration(1e6);  // No division
            p_cycle_model->SetMaxCellCycleDuration(1e6);
            
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);  // Random ages
            p_cell->GetCellData()->SetItem("cell_type", 1.0);
            
            cells.push_back(p_cell);
        }
        
        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(500.0);
        
        // Set up simulation
        OffLatticeSimulation<2> simulator(cell_population);
        
        std::stringstream output_dir;
        output_dir << "PainterReplication2d/" << ecmType << "/run_" << randomSeed;
        simulator.SetOutputDirectory(output_dir.str());
        
        simulator.SetSamplingTimestepMultiple(7200);  // Output every 60 minutes
        simulator.SetDt(5);  // in minutes
        simulator.SetEndTime(7200);  // 5 days = 7200 minutes
        // 5 days in mintues = 5 * 24 * 60 = 7200 minutes
        // Force 1: Cell-cell adhesion (COMMENTED OUT - using only ECM forces like Test2dDynamicECMInvasion)
        // MAKE_PTR(DifferentialAdhesionForce<2>, p_adhesion_force);
        // p_adhesion_force->SetCutOffLength(50.0);
        // p_adhesion_force->SetMeinekeSpringStiffness(0.4);  // Cell-cell adhesion = 0.4
        // p_adhesion_force->SetMeinekeSpringGrowthDuration(1.0);
        // p_adhesion_force->SetMeinekeDivisionRestingSpringLength(0.5);
        // 
        // simulator.AddForce(p_adhesion_force);
        
        // Force: Dynamic ECM Contact Guidance (with remodeling like Test2dDynamicECMInvasion)
        MAKE_PTR(DynamicECMContactGuidanceForce<2>, p_ecm_force);
        p_ecm_force->SetECMField(p_ecm_field);
        p_ecm_force->SetBaseSpeed(2.5);       // Base speed = 1.25 μm/min
        p_ecm_force->SetECMSensitivity(1.0);   // Full ECM sensitivity
        p_ecm_force->SetEnableRemodeling(true);   // Dynamic ECM!
        p_ecm_force->SetEnableDegradation(true);  // Cells degrade ECM
        p_ecm_force->SetEnableDeposition(false);  // Don't deposit new ECM (yet)
        
        simulator.AddForce(p_ecm_force);
        
        // Add cell influx modifier (30 cells every 180 minutes)
        MAKE_PTR(CellInfluxModifier, p_influx_modifier);
        p_influx_modifier->SetInfluxInterval(180.0);
        p_influx_modifier->SetCellsPerInflux(30);
        p_influx_modifier->SetNextNodeIndex(initial_nodes.size());
        simulator.AddSimulationModifier(p_influx_modifier);
        
        // Add ECM field writer
        MAKE_PTR_ARGS(ECMFieldWriter<2>, p_ecm_writer, (p_ecm_field, 120));
        simulator.AddSimulationModifier(p_ecm_writer);
        
        // Run simulation
        simulator.Solve();
        
        // Clean up
        for (unsigned i = 0; i < initial_nodes.size(); i++)
        {
            delete initial_nodes[i];
        }
    }
    
    void TestPainterRandomECM()
    {
        // Get run number from environment variable (for parallel execution)
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 50000 + run_num * 1000;
            RunPainterSimulation("random", seed);
        }
        else
        {
            // Run single simulation
            RunPainterSimulation("random", 50000);
        }
    }
    
    void TestPainterParallelECM()
    {
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 50000 + run_num * 1000;
            RunPainterSimulation("parallel", seed);
        }
        else
        {
            RunPainterSimulation("parallel", 50000);
        }
    }
    
    void TestPainterPerpendicularECM()
    {
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 50000 + run_num * 1000;
            RunPainterSimulation("perpendicular", seed);
        }
        else
        {
            RunPainterSimulation("perpendicular", 50000);
        }
    }
    
    void TestPainterMixedECM()
    {
        const char* run_num_env = std::getenv("RUN_NUMBER");
        if (run_num_env != nullptr)
        {
            unsigned run_num = std::atoi(run_num_env);
            unsigned seed = 50000 + run_num * 1000;
            RunPainterSimulation("mixed", seed);
        }
        else
        {
            RunPainterSimulation("mixed", 50000);
        }
    }
};

#endif /* TEST2DPAINTERREPLICATION_HPP_ */
