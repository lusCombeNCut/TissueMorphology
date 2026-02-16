/*
 * Test2dCryptBuddingNodeBased.hpp
 *
 * 2D node-based simulations of intestinal crypt budding from an organoid.
 *
 * Initial condition: cells arranged in an annular ring representing an
 * early organoid cross-section with a central lumen. The apical surface
 * faces inward (toward the lumen) and the basal surface faces outward
 * (toward the ECM/basement membrane).
 *
 * Key independent variable: ECM_STIFFNESS (environment variable)
 *   - Controls basement membrane stiffness (radial confinement)
 *   - Higher stiffness -> stronger confinement -> fewer/shallower buds
 *   - Lower stiffness -> easier budding -> more crypts
 *
 * Physics:
 *   - DifferentialAdhesionForce: cell-type-specific spring constants
 *   - BasementMembraneForce: inward radial restoring force (ECM)
 *   - LumenPressureForce: outward radial pressure maintaining lumen
 *   - ApicalConstrictionForce: actomyosin contraction on lumen-facing cells
 *   - ContactInhibitionCellCycleModel: density-dependent proliferation
 *
 * Crypt detection:
 *   - Post-processing in polar coordinates (analyse_crypt_budding.py)
 *   - Outward radial protrusions from the organoid surface = crypt buds
 *
 * References:
 *   - Dunn et al. (2012) J. Theor. Biol. (crypt mechanics)
 *   - Hannezo et al. (2011) PRL (epithelial buckling)
 *   - Serra et al. (2019) Nature (organoid symmetry breaking)
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
#include <map>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

// Mesh & population
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

// Simulation
#include "OffLatticeSimulation.hpp"

// Cell cycle
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
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

// Project-specific forces
#include "DifferentialAdhesionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"

#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"


/**
 * Summary modifier -- tracks radial statistics relevant to crypt budding
 * from a circular organoid: mean/max/min radius, number of cells, etc.
 */
template<unsigned DIM>
class OrganoidBuddingSummaryModifier : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    std::string mOutputDir;
    double mStiffness;
    double mEndTime;
    bool mHeaderWritten;
    unsigned mSamplingMultiple;
    unsigned mLogInterval;
    unsigned mLastOutputStep;
    unsigned mLastLogStep;

public:
    OrganoidBuddingSummaryModifier(double stiffness, unsigned samplingMultiple,
                                   double endTime = 168.0)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mStiffness(stiffness),
          mEndTime(endTime),
          mHeaderWritten(false),
          mSamplingMultiple(samplingMultiple),
          mLogInterval(samplingMultiple / 6),
          mLastOutputStep(0),
          mLastLogStep(0)
    {
        if (mLogInterval == 0) mLogInterval = 1;
    }

    void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
                    std::string outputDirectory)
    {
        mOutputDir = outputDirectory;
    }

    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        unsigned current_step = SimulationTime::Instance()->GetTimeStepsElapsed();
        double current_time = SimulationTime::Instance()->GetTime();
        unsigned num_cells = rCellPopulation.GetNumRealCells();

        // Progress log
        if (current_step == 0 || current_step - mLastLogStep >= mLogInterval)
        {
            mLastLogStep = current_step;
            double pct = (mEndTime > 0.0) ? (current_time / mEndTime) * 100.0 : 0.0;
            std::cout << "[Progress] t=" << std::fixed << std::setprecision(1)
                      << current_time << "h / " << mEndTime << "h  ("
                      << std::setprecision(1) << pct << "%)  cells="
                      << num_cells << std::endl;
        }

        // CSV at sampling interval
        if (current_step - mLastOutputStep < mSamplingMultiple && current_step > 0)
        {
            return;
        }
        mLastOutputStep = current_step;

        // Compute centroid
        c_vector<double, DIM> centroid = zero_vector<double>(DIM);
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            centroid += rCellPopulation.GetNode(idx)->rGetLocation();
        }
        if (num_cells > 0) centroid /= static_cast<double>(num_cells);

        // Compute radial statistics from centroid
        double sum_r = 0.0, sum_r2 = 0.0;
        double max_r = -1e10, min_r = 1e10;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetNode(idx)->rGetLocation();
            double r = norm_2(pos - centroid);
            sum_r += r;
            sum_r2 += r * r;
            if (r > max_r) max_r = r;
            if (r < min_r) min_r = r;
        }
        double mean_r = (num_cells > 0) ? sum_r / num_cells : 0.0;
        double var_r = (num_cells > 1) ? (sum_r2 / num_cells - mean_r * mean_r) : 0.0;

        // Write CSV
        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory()
                               + mOutputDir + "/crypt_summary.csv";
        std::ofstream file;
        if (!mHeaderWritten)
        {
            file.open(filename.c_str());
            file << "time,num_cells,mean_r,var_r,max_r,min_r,r_range,centroid_x,centroid_y,stiffness"
                 << std::endl;
            mHeaderWritten = true;
        }
        else
        {
            file.open(filename.c_str(), std::ios::app);
        }

        file << std::fixed << std::setprecision(4)
             << current_time << ","
             << num_cells << ","
             << mean_r << ","
             << var_r << ","
             << max_r << ","
             << min_r << ","
             << (max_r - min_r) << ","
             << centroid[0] << ","
             << centroid[1] << ","
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
     * Run a single node-based organoid crypt budding simulation.
     *
     * @param bmStiffness  basement membrane / ECM stiffness parameter
     * @param randomSeed   RNG seed for reproducibility
     * @param outputSubDir output sub-directory label
     */
    void RunCryptBuddingSimulation(double bmStiffness, unsigned randomSeed,
                                   const std::string& outputSubDir)
    {
        // ================================================================
        // PARAMETERS
        // ================================================================

        // Organoid geometry -- annular ring (circular cross-section)
        const double organoid_radius = 8.0;
        const unsigned num_cells_in_ring = 80;
        const double interaction_cutoff = 2.5;

        // Center of the organoid
        const double center_x = 0.0;
        const double center_y = 0.0;

        // Time parameters (hours)
        const double dt = 0.005;
        const double relaxation_time = 10.0;       // let geometry settle before growth
        const double end_time = 200.0;             // ~8 days (after relaxation)
        const unsigned sampling_multiple = 200;    // output every dt*200 = 1 hour

        // Basement membrane / ECM (outer confinement)
        const double bm_radius = organoid_radius + 2.0;
        // bmStiffness is the key variable

        // ECM degradation
        const double ecm_degradation_rate = 0.02;
        const double ecm_max_radius = organoid_radius * 4.0;

        // Lumen pressure (inner expansion)
        const double lumen_pressure = 2.0;
        const double lumen_eq_radius = organoid_radius + 1.0;  // above initial r so outward pressure exists at startup

        // Apical constriction
        const double apical_constriction_strength = 3.0;

        // Spring forces (differential adhesion)
        const double spring_stiffness = 30.0;
        const double spring_cutoff = 1.5;
        const double apical_apical_adhesion = 1.2;
        const double basal_basal_adhesion = 1.0;
        const double apical_basal_adhesion = 0.5;

        // Cell cycle
        const double quiescent_fraction = 0.7;

        // Radial sloughing boundary
        const double max_radius_for_slough = organoid_radius * 5.0;

        // ================================================================
        // SETUP
        // ================================================================

        RandomNumberGenerator::Instance()->Reseed(randomSeed);

        // Create nodes in an annular ring
        std::vector<Node<2>*> nodes;
        for (unsigned i = 0; i < num_cells_in_ring; i++)
        {
            double theta = 2.0 * M_PI * static_cast<double>(i) / num_cells_in_ring;

            // Slight radial noise to break perfect symmetry
            double r_noise = organoid_radius +
                (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.3;

            double x = center_x + r_noise * cos(theta);
            double y = center_y + r_noise * sin(theta);

            nodes.push_back(new Node<2>(i, false, x, y));
        }

        // Create nodes-only mesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, interaction_cutoff);

        // Clean up node pointers (mesh takes ownership via copy)
        for (unsigned i = 0; i < nodes.size(); i++)
        {
            delete nodes[i];
        }

        // ================================================================
        // CREATE CELLS with apical-basal polarity
        // ================================================================

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(quiescent_fraction);
            p_cycle->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle));

            // Assign cell type based on angular position
            // Bottom of ring (theta ~ -pi/2) = crypt base = stem cells
            // Top (theta ~ pi/2) = villus-like = differentiated
            c_vector<double, 2> pos = mesh.GetNode(i)->rGetLocation();
            double theta = atan2(pos[1] - center_y, pos[0] - center_x);

            // Normalize: bottom = 0, going clockwise
            double angle_from_bottom = theta + M_PI / 2.0;
            if (angle_from_bottom < 0.0) angle_from_bottom += 2.0 * M_PI;
            double frac = angle_from_bottom / (2.0 * M_PI);

            if (frac < 0.2 || frac > 0.8)
            {
                // Bottom 40% arc: stem cells (crypt base)
                p_cell->SetCellProliferativeType(p_stem_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 0.0);
            }
            else if (frac < 0.35 || frac > 0.65)
            {
                // Flanks: transit-amplifying
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 1.0);
            }
            else
            {
                // Top: differentiated
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 2.0);
            }

            // Random birth time desynchronization
            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();

            // Volume tracking
            p_cell->GetCellData()->SetItem("volume", 1.0);

            // BM stiffness per cell
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", bmStiffness);

            // Apical marker (for DifferentialAdhesionForce & ApicalConstrictionForce)
            p_cell->GetCellData()->SetItem("is_apical", 1.0);

            cells.push_back(p_cell);
        }

        // ================================================================
        // CELL POPULATION
        // ================================================================

        NodeBasedCellPopulation<2> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(50.0);

        // Writers
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // ================================================================
        // SIMULATION
        // ================================================================

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(outputSubDir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);

        // ================================================================
        // FORCES
        // ================================================================

        // 1. Cell-cell interactions with differential adhesion
        MAKE_PTR(DifferentialAdhesionForce<2>, p_spring_force);
        p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
        p_spring_force->SetCutOffLength(spring_cutoff);
        p_spring_force->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring_force->SetMeinekeSpringGrowthDuration(1.0);
        p_spring_force->SetApicalApicalAdhesion(apical_apical_adhesion);
        p_spring_force->SetBasalBasalAdhesion(basal_basal_adhesion);
        p_spring_force->SetApicalBasalAdhesion(apical_basal_adhesion);
        simulator.AddForce(p_spring_force);

        // 2. Basement membrane -- radial confinement (KEY parameter)
        MAKE_PTR(BasementMembraneForce<2>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bmStiffness);
        p_bm_force->SetBasementMembraneRadius(bm_radius);
        c_vector<double, 2> organoid_center;
        organoid_center[0] = center_x;
        organoid_center[1] = center_y;
        p_bm_force->SetOrganoidCenter(organoid_center);
        p_bm_force->EnableEcmDegradation(ecm_degradation_rate, ecm_max_radius);
        simulator.AddForce(p_bm_force);

        // 3. Lumen pressure -- outward radial force maintaining lumen
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetPressureStrength(lumen_pressure);
        p_lumen_force->SetLumenEquilibriumRadius(lumen_eq_radius);
        p_lumen_force->SetTrackCenter(true);
        simulator.AddForce(p_lumen_force);

        // 4. Apical constriction -- drives invagination
        MAKE_PTR(ApicalConstrictionForce<2>, p_apical_force);
        p_apical_force->SetConstrictionStrength(apical_constriction_strength);
        simulator.AddForce(p_apical_force);

        // ================================================================
        // CELL KILLERS -- bounding box sloughing
        // ================================================================
        {
            c_vector<double, 2> point, normal;

            // Top
            point = zero_vector<double>(2); point[1] = max_radius_for_slough;
            normal = zero_vector<double>(2); normal[1] = 1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_k1, (&cell_population, point, normal));
            simulator.AddCellKiller(p_k1);

            // Bottom
            point = zero_vector<double>(2); point[1] = -max_radius_for_slough;
            normal = zero_vector<double>(2); normal[1] = -1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_k2, (&cell_population, point, normal));
            simulator.AddCellKiller(p_k2);

            // Right
            point = zero_vector<double>(2); point[0] = max_radius_for_slough;
            normal = zero_vector<double>(2); normal[0] = 1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_k3, (&cell_population, point, normal));
            simulator.AddCellKiller(p_k3);

            // Left
            point = zero_vector<double>(2); point[0] = -max_radius_for_slough;
            normal = zero_vector<double>(2); normal[0] = -1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_k4, (&cell_population, point, normal));
            simulator.AddCellKiller(p_k4);
        }

        // ================================================================
        // MODIFIERS
        // ================================================================

        MAKE_PTR(VolumeTrackingModifier<2>, p_vol_modifier);
        simulator.AddSimulationModifier(p_vol_modifier);

        boost::shared_ptr<OrganoidBuddingSummaryModifier<2>> p_summary(
            new OrganoidBuddingSummaryModifier<2>(bmStiffness, sampling_multiple, end_time));
        simulator.AddSimulationModifier(p_summary);

        // ================================================================
        // RUN
        // ================================================================

        std::cout << "============================================" << std::endl;
        std::cout << "  2D Node-Based Organoid Crypt Budding" << std::endl;
        std::cout << "  BM Stiffness:     " << bmStiffness << std::endl;
        std::cout << "  BM Radius:        " << bm_radius << std::endl;
        std::cout << "  Lumen Pressure:   " << lumen_pressure << std::endl;
        std::cout << "  Lumen Eq Radius:  " << lumen_eq_radius << std::endl;
        std::cout << "  Apical Constr.:   " << apical_constriction_strength << std::endl;
        std::cout << "  Spring Stiffness: " << spring_stiffness << std::endl;
        std::cout << "  Organoid Radius:  " << organoid_radius << std::endl;
        std::cout << "  Cells in Ring:    " << num_cells_in_ring << std::endl;
        std::cout << "  Random Seed:      " << randomSeed << std::endl;
        std::cout << "  dt:               " << dt << std::endl;
        std::cout << "  Relaxation:       " << relaxation_time << " hours" << std::endl;
        std::cout << "  ECM Degrad Rate:  " << ecm_degradation_rate << std::endl;
        std::cout << "  End Time:         " << end_time << " hours" << std::endl;
        std::cout << "  Output:           testoutput/" << outputSubDir << std::endl;
        std::cout << "============================================" << std::endl;

        // ================================================================
        // PHASE 1: RELAXATION — let geometry settle, no proliferation
        // ================================================================

        // Temporarily set all cells to differentiated (non-dividing)
        std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> original_types;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            original_types[*cell_iter] = cell_iter->GetCellProliferativeType();
            cell_iter->SetCellProliferativeType(p_diff_type);
        }

        simulator.SetEndTime(relaxation_time);

        std::cout << "\n--- Phase 1: Relaxation (" << relaxation_time << "h, no proliferation) ---" << std::endl;
        simulator.Solve();
        std::cout << "Relaxation complete. Cells: " << cell_population.GetNumRealCells() << std::endl;

        // ================================================================
        // PHASE 2: GROWTH — restore cell types and run main simulation
        // ================================================================

        // Restore original proliferative types
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            if (original_types.count(*cell_iter))
            {
                cell_iter->SetCellProliferativeType(original_types[*cell_iter]);
            }
        }

        simulator.SetEndTime(relaxation_time + end_time);

        std::cout << "\n--- Phase 2: Growth (" << end_time << "h with proliferation) ---" << std::endl;
        simulator.Solve();

        // ================================================================
        // POST-SIMULATION SUMMARY
        // ================================================================

        unsigned final_cells = cell_population.GetNumRealCells();

        unsigned stem_count = 0, ta_count = 0, diff_count = 0;
        double max_r = 0.0, min_r = 1e10;
        c_vector<double, 2> final_centroid = zero_vector<double>(2);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            unsigned idx = cell_population.GetLocationIndexUsingCell(*cell_iter);
            final_centroid += cell_population.GetNode(idx)->rGetLocation();
        }
        if (final_cells > 0) final_centroid /= static_cast<double>(final_cells);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double type_id = cell_iter->GetCellData()->GetItem("cell_type_id");
            if (type_id < 0.5) stem_count++;
            else if (type_id < 1.5) ta_count++;
            else diff_count++;

            unsigned idx = cell_population.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 2> pos = cell_population.GetNode(idx)->rGetLocation();
            double r = norm_2(pos - final_centroid);
            if (r > max_r) max_r = r;
            if (r < min_r) min_r = r;
        }

        std::cout << "\n============================================" << std::endl;
        std::cout << "  SIMULATION COMPLETE" << std::endl;
        std::cout << "  Final cells:      " << final_cells << std::endl;
        std::cout << "    Stem:           " << stem_count << std::endl;
        std::cout << "    TA:             " << ta_count << std::endl;
        std::cout << "    Diff:           " << diff_count << std::endl;
        std::cout << "  Radius range:     [" << min_r << ", " << max_r << "]" << std::endl;
        std::cout << "  BM stiffness:     " << bmStiffness << std::endl;
        std::cout << "============================================" << std::endl;

        TS_ASSERT_LESS_THAN(0u, final_cells);
    }

public:

    /**
     * Parameterised test: reads ECM_STIFFNESS and RUN_NUMBER from
     * environment variables (set by HPC array job script).
     */
    void TestCryptBuddingStiffnessSweep()
    {
        double bm_stiffness = 5.0;
        const char* stiffness_env = std::getenv("ECM_STIFFNESS");
        if (stiffness_env != nullptr)
        {
            bm_stiffness = std::atof(stiffness_env);
        }

        unsigned run_number = 0;
        const char* run_env = std::getenv("RUN_NUMBER");
        if (run_env != nullptr)
        {
            run_number = std::atoi(run_env);
        }

        unsigned seed = static_cast<unsigned>(bm_stiffness * 10000) + run_number * 137;

        std::stringstream subdir;
        subdir << "CryptBudding2d_NodeBased/stiffness_" << std::fixed
               << std::setprecision(1) << bm_stiffness << "/run_" << run_number;

        RunCryptBuddingSimulation(bm_stiffness, seed, subdir.str());
    }
};

#endif /* TEST2DCRYPTBUDDINGNODEBASED_HPP_ */
