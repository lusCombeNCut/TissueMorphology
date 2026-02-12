/*
 * Test2dCryptBuddingNodeBased.hpp
 *
 * 2D node-based simulations of intestinal crypt budding at varying
 * ECM/substrate stiffness levels. Designed to answer:
 *   "Do changes in stiffness of in-silico hydrogel simulations have the
 *    same effect as changes in in-vivo ECM stiffness?"
 *
 * Approach:
 *   - Start from a flat monolayer of cells on a honeycomb mesh
 *   - Stem cells at the base proliferate; differentiated cells slough at top
 *   - A radial basement membrane force confines cells, with stiffness as the
 *     key independent variable
 *   - Higher BM stiffness → more confinement → fewer/shallower buds
 *   - Lower BM stiffness → easier budding → more crypts
 *
 * The simulation sweeps over ECM stiffness values by reading the
 * environment variables ECM_STIFFNESS and RUN_NUMBER.
 * Each (stiffness, run) pair uses a distinct random seed for reproducibility.
 *
 * Post-processing: count number of crypt buds from VTU output using the
 * companion Python script (scripts/analyse_crypt_budding.py).
 *
 * No ECM fibre modelling or degradation — baseline mechanical model.
 *
 * References:
 *   - Dunn et al. (2012) J. Theor. Biol. (crypt mechanics)
 *   - Osborne et al. (2010) Phil. Trans. R. Soc. A (crypt geometry)
 *   - Painter (2009) J. Math. Biol. (ECM contact guidance)
 */

#ifndef TEST2DCRYPTBUDDINGNODEBASED_HPP_
#define TEST2DCRYPTBUDDINGNODEBASED_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

// Mesh & population
#include "HoneycombMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

// Simulation
#include "OffLatticeSimulation.hpp"

// Forces
#include "GeneralisedLinearSpringForce.hpp"

// Cell cycle
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"

// Cell types & states
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

// Cell killers
#include "PlaneBasedCellKiller.hpp"

// Modifiers & writers
#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

// Project-specific
#include "BasementMembraneForce.hpp"

// Boundary conditions
#include "PlaneBoundaryCondition.hpp"

#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"


/**
 * Simulation modifier that writes a simple summary CSV every output step.
 * Records: time, num_cells, centroid_y, max_y, min_y, y_range, stiffness
 */
template<unsigned DIM>
class CryptBuddingSummaryModifier : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    std::string mOutputDir;
    double mStiffness;
    bool mHeaderWritten;
    unsigned mSamplingMultiple;
    unsigned mLastOutputStep;

public:
    CryptBuddingSummaryModifier(double stiffness, unsigned samplingMultiple)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mStiffness(stiffness),
          mHeaderWritten(false),
          mSamplingMultiple(samplingMultiple),
          mLastOutputStep(0)
    {
    }

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                    std::string outputDirectory)
    {
        mOutputDir = outputDirectory;
    }

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        unsigned current_step = SimulationTime::Instance()->GetTimeStepsElapsed();
        if (current_step - mLastOutputStep < mSamplingMultiple && current_step > 0)
        {
            return;
        }
        mLastOutputStep = current_step;

        double current_time = SimulationTime::Instance()->GetTime();
        unsigned num_cells = rCellPopulation.GetNumRealCells();

        // Compute y-statistics for crypt depth analysis
        double sum_y = 0.0;
        double max_y = -1e10;
        double min_y = 1e10;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetNode(idx)->rGetLocation();
            double y = pos[1];
            sum_y += y;
            if (y > max_y) max_y = y;
            if (y < min_y) min_y = y;
        }

        double mean_y = (num_cells > 0) ? sum_y / num_cells : 0.0;
        double y_range = max_y - min_y;

        // Write to summary CSV
        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory()
                               + mOutputDir + "/crypt_summary.csv";

        std::ofstream file;
        if (!mHeaderWritten)
        {
            file.open(filename.c_str());
            file << "time,num_cells,mean_y,max_y,min_y,y_range,stiffness" << std::endl;
            mHeaderWritten = true;
        }
        else
        {
            file.open(filename.c_str(), std::ios::app);
        }

        file << std::fixed << std::setprecision(4)
             << current_time << ","
             << num_cells << ","
             << mean_y << ","
             << max_y << ","
             << min_y << ","
             << y_range << ","
             << mStiffness
             << std::endl;
        file.close();
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t<Stiffness>" << mStiffness << "</Stiffness>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};


class Test2dCryptBuddingNodeBased : public AbstractCellBasedTestSuite
{
private:

    /**
     * Run a single crypt budding simulation at the given BM stiffness.
     *
     * @param bmStiffness  basement membrane stiffness parameter
     * @param randomSeed   seed for random number generator (reproducibility)
     * @param outputSubDir sub-directory label (stiffness + run ID)
     */
    void RunCryptBuddingSimulation(double bmStiffness, unsigned randomSeed,
                                   const std::string& outputSubDir)
    {
        // ================================================================
        // PARAMETERS
        // ================================================================

        // Domain geometry — flat monolayer that can buckle downward
        const unsigned cells_across = 20;   // width of initial monolayer
        const unsigned cells_up = 6;        // initial height (rows)
        const double interaction_cutoff = 3.0;  // cell interaction radius

        // Time parameters (hours)
        const double dt = 0.005;            // timestep
        const double end_time = 200.0;      // ~8.3 days
        const unsigned sampling_multiple = 200; // output every dt*200 = 1 hour

        // Basement membrane / substrate
        const double bm_radius = 15.0;      // radius of confinement region
        // bmStiffness passed as argument — the key independent variable

        // Spring force
        const double spring_stiffness = 30.0;
        const double spring_cutoff = 1.5;   // natural spring rest length

        // Cell cycle
        const double stem_g1_min = 12.0;    // hours
        const double stem_g1_max = 14.0;
        const double ta_g1_min = 4.0;
        const double ta_g1_max = 6.0;
        const double quiescent_fraction = 0.8;  // contact inhibition threshold

        // Sloughing
        const double slough_height = 20.0;  // cells above this y are killed

        // ================================================================
        // SETUP
        // ================================================================

        // Seed RNG
        RandomNumberGenerator::Instance()->Reseed(randomSeed);

        // Create honeycomb mesh
        HoneycombMeshGenerator generator(cells_across, cells_up, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();

        // Convert to nodes-only mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, interaction_cutoff);

        // ================================================================
        // CREATE CELLS
        // ================================================================

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            // Use contact inhibition cell cycle
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(quiescent_fraction);
            p_cycle->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle));

            // Assign cell type based on y-position
            c_vector<double, 2> pos = mesh.GetNode(i)->rGetLocation();
            double y = pos[1];
            double y_max = (double)(cells_up - 1) * sqrt(3.0) / 2.0;
            double y_frac = y / y_max;

            if (y_frac < 0.2)
            {
                // Bottom 20%: stem cells
                p_cell->SetCellProliferativeType(p_stem_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 0.0);
            }
            else if (y_frac < 0.6)
            {
                // Middle: transit-amplifying
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 1.0);
            }
            else
            {
                // Top: differentiated (non-dividing — will slough)
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 2.0);
            }

            // Desynchronize birth times
            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();

            // Volume tracking data
            p_cell->GetCellData()->SetItem("volume", 1.0);

            // Store stiffness in cell data (used by BasementMembraneForce)
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", bmStiffness);

            cells.push_back(p_cell);
        }

        // ================================================================
        // CREATE CELL POPULATION
        // ================================================================

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(50.0);

        // Writers
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // ================================================================
        // SIMULATION
        // ================================================================

        OffLatticeSimulation<2> simulator(cell_population);

        std::string output_dir = "CryptBudding2d_NodeBased/" + outputSubDir;
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);

        // ================================================================
        // FORCES
        // ================================================================

        // 1. Cell-cell spring force (repulsion + adhesion)
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(spring_cutoff);
        p_spring_force->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring_force->SetMeinekeSpringGrowthDuration(1.0);
        simulator.AddForce(p_spring_force);

        // 2. Basement membrane / substrate stiffness force
        //    This is the KEY parameter — models ECM/hydrogel stiffness
        //    Higher stiffness → stronger confinement → fewer buds
        MAKE_PTR(BasementMembraneForce<2>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bmStiffness);
        p_bm_force->SetBasementMembraneRadius(bm_radius);
        // No degradation in this baseline experiment
        p_bm_force->DisableEcmDegradation();
        simulator.AddForce(p_bm_force);

        // ================================================================
        // BOUNDARY CONDITIONS
        // ================================================================

        // Pin cells at the bottom (y = 0 plane)
        c_vector<double, 2> bottom_point = zero_vector<double>(2);
        c_vector<double, 2> bottom_normal = zero_vector<double>(2);
        bottom_normal[1] = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc,
                      (&cell_population, bottom_point, bottom_normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // ================================================================
        // CELL KILLERS
        // ================================================================

        // Slough cells at the top to mimic villus shedding
        c_vector<double, 2> top_point = zero_vector<double>(2);
        top_point[1] = slough_height;
        c_vector<double, 2> top_normal = zero_vector<double>(2);
        top_normal[1] = 1.0;
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer,
                      (&cell_population, top_point, top_normal));
        simulator.AddCellKiller(p_killer);

        // ================================================================
        // MODIFIERS
        // ================================================================

        // Volume tracking for contact inhibition
        MAKE_PTR(VolumeTrackingModifier<2>, p_vol_modifier);
        simulator.AddSimulationModifier(p_vol_modifier);

        // Custom summary writer
        boost::shared_ptr<CryptBuddingSummaryModifier<2>> p_summary(
            new CryptBuddingSummaryModifier<2>(bmStiffness, sampling_multiple));
        simulator.AddSimulationModifier(p_summary);

        // ================================================================
        // RUN
        // ================================================================

        std::cout << "============================================" << std::endl;
        std::cout << "  2D Node-Based Crypt Budding Simulation" << std::endl;
        std::cout << "  BM Stiffness:  " << bmStiffness << std::endl;
        std::cout << "  Random Seed:   " << randomSeed << std::endl;
        std::cout << "  End Time:      " << end_time << " hours" << std::endl;
        std::cout << "  Output:        testoutput/" << output_dir << std::endl;
        std::cout << "============================================" << std::endl;

        simulator.Solve();

        // ================================================================
        // POST-SIMULATION SUMMARY
        // ================================================================

        unsigned final_cells = cell_population.GetNumRealCells();

        // Count cell types
        unsigned stem_count = 0, ta_count = 0, diff_count = 0;
        double max_y = -1e10, min_y = 1e10;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double type_id = cell_iter->GetCellData()->GetItem("cell_type_id");
            if (type_id < 0.5) stem_count++;
            else if (type_id < 1.5) ta_count++;
            else diff_count++;

            unsigned idx = cell_population.GetLocationIndexUsingCell(*cell_iter);
            double y = cell_population.GetNode(idx)->rGetLocation()[1];
            if (y > max_y) max_y = y;
            if (y < min_y) min_y = y;
        }

        std::cout << "\n============================================" << std::endl;
        std::cout << "  SIMULATION COMPLETE" << std::endl;
        std::cout << "  Final cells:   " << final_cells << std::endl;
        std::cout << "    Stem:        " << stem_count << std::endl;
        std::cout << "    TA:          " << ta_count << std::endl;
        std::cout << "    Diff:        " << diff_count << std::endl;
        std::cout << "  Y range:       [" << min_y << ", " << max_y << "]" << std::endl;
        std::cout << "  BM stiffness:  " << bmStiffness << std::endl;
        std::cout << "============================================" << std::endl;

        // Basic sanity assertions
        TS_ASSERT_LESS_THAN(0u, final_cells);
    }

public:

    /**
     * Parameterised test: reads ECM_STIFFNESS and RUN_NUMBER from
     * environment variables (set by HPC array job script).
     *
     * Default: stiffness = 5.0, run = 0 (for local testing).
     */
    void TestCryptBuddingStiffnessSweep()
    {
        // Read stiffness from environment
        double bm_stiffness = 5.0;
        const char* stiffness_env = std::getenv("ECM_STIFFNESS");
        if (stiffness_env != nullptr)
        {
            bm_stiffness = std::atof(stiffness_env);
        }

        // Read run number from environment
        unsigned run_number = 0;
        const char* run_env = std::getenv("RUN_NUMBER");
        if (run_env != nullptr)
        {
            run_number = std::atoi(run_env);
        }

        // Deterministic seed from stiffness and run number
        unsigned seed = static_cast<unsigned>(bm_stiffness * 10000) + run_number * 137;

        // Output directory encodes stiffness and run
        std::stringstream subdir;
        subdir << "stiffness_" << std::fixed << std::setprecision(1) << bm_stiffness
               << "/run_" << run_number;

        RunCryptBuddingSimulation(bm_stiffness, seed, subdir.str());
    }
};

#endif /* TEST2DCRYPTBUDDINGNODEBASED_HPP_ */
