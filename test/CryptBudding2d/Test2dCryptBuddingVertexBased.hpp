/*
 * Test2dCryptBuddingVertexBased.hpp
 *
 * 2D vertex-based simulations of intestinal crypt budding at varying
 * ECM/substrate stiffness levels. Companion to Test2dCryptBuddingNodeBased.
 *
 * The vertex model naturally captures epithelial buckling instabilities:
 * when proliferative pressure in a confined monolayer exceeds a critical
 * threshold (dependent on substrate stiffness), the tissue buckles to
 * form crypt-like invaginations.
 *
 * Key independent variable: ECM_STIFFNESS (environment variable)
 *   - Modelled as the area elasticity parameter in the Nagai-Honda force
 *     (resistance to cell deformation ≈ substrate confinement)
 *   - Also influences the basement membrane restoring force strength
 *
 * Approach:
 *   - Start from a flat honeycomb vertex mesh (epithelial monolayer)
 *   - Stem cells at the base proliferate with contact inhibition
 *   - Nagai-Honda force governs cell mechanics (area + perimeter energy)
 *   - Basement membrane force provides radial confinement
 *   - Cells sloughed at top boundary
 *   - No ECM degradation or fibre modelling (baseline)
 *
 * Output:
 *   - VTU files for ParaView visualization
 *   - crypt_summary.csv with time-series statistics
 *   - Final cell positions for crypt counting (Python post-processing)
 *
 * References:
 *   - Nagai & Honda (2001) Phil. Mag. B (vertex model)
 *   - Hannezo et al. (2011) PRL (epithelial buckling)
 *   - Dunn et al. (2012) J. Theor. Biol. (crypt mechanics)
 */

#ifndef TEST2DCRYPTBUDDINGVERTEXBASED_HPP_
#define TEST2DCRYPTBUDDINGVERTEXBASED_HPP_

#include <cxxtest/TestSuite.h>
#include <cmath>
#include <cstdlib>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iomanip>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

// Vertex mesh & population
#include "HoneycombVertexMeshGenerator.hpp"
#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"

// Simulation
#include "OffLatticeSimulation.hpp"

// Forces
#include "NagaiHondaForce.hpp"
#include "FarhadifarForce.hpp"

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
#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellVolumesWriter.hpp"

// Project-specific
#include "BasementMembraneForce.hpp"

// Boundary conditions
#include "PlaneBoundaryCondition.hpp"

#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"


/**
 * Vertex-specific summary modifier — tracks cell shape statistics
 * relevant to buckling detection: area variance, elongation, etc.
 */
template<unsigned DIM>
class VertexCryptSummaryModifier : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    std::string mOutputDir;
    double mStiffness;
    bool mHeaderWritten;
    unsigned mSamplingMultiple;
    unsigned mLastOutputStep;

public:
    VertexCryptSummaryModifier(double stiffness, unsigned samplingMultiple)
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

        // Compute centroid y-statistics
        double sum_y = 0.0, sum_y2 = 0.0;
        double max_y = -1e10, min_y = 1e10;
        double sum_area = 0.0;

        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned elem_idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> centroid = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double y = centroid[1];
            sum_y += y;
            sum_y2 += y * y;
            if (y > max_y) max_y = y;
            if (y < min_y) min_y = y;

            // Cell area from vertex mesh
            if (p_vertex_pop)
            {
                double area = p_vertex_pop->rGetMesh().GetVolumeOfElement(elem_idx);
                sum_area += area;
            }
        }

        double mean_y = (num_cells > 0) ? sum_y / num_cells : 0.0;
        double var_y = (num_cells > 1) ? (sum_y2 / num_cells - mean_y * mean_y) : 0.0;
        double mean_area = (num_cells > 0) ? sum_area / num_cells : 0.0;

        // Write CSV
        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory()
                               + mOutputDir + "/crypt_summary.csv";

        std::ofstream file;
        if (!mHeaderWritten)
        {
            file.open(filename.c_str());
            file << "time,num_cells,mean_y,var_y,max_y,min_y,y_range,mean_area,stiffness" << std::endl;
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
             << var_y << ","
             << max_y << ","
             << min_y << ","
             << (max_y - min_y) << ","
             << mean_area << ","
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


class Test2dCryptBuddingVertexBased : public AbstractCellBasedTestSuite
{
private:

    /**
     * Run a single vertex-based crypt budding simulation.
     *
     * @param ecmStiffness      ECM stiffness parameter (varies BM + area elasticity)
     * @param randomSeed        RNG seed for reproducibility
     * @param outputSubDir      sub-directory for output
     */
    void RunVertexCryptBuddingSimulation(double ecmStiffness, unsigned randomSeed,
                                          const std::string& outputSubDir)
    {
        // ================================================================
        // PARAMETERS
        // ================================================================

        // Mesh geometry — flat honeycomb monolayer
        const unsigned cells_across = 28;   // width (doubled)
        const unsigned cells_up = 16;       // height (doubled)
        const bool flat_bottom = true;

        // Time (hours)
        const double dt = 0.002;
        const double end_time = 200.0;      // ~8.3 days
        const unsigned sampling_multiple = 500; // output every dt*500 = 1 hour

        // Nagai-Honda force parameters
        // Area elasticity scales with ECM stiffness — stiffer substrate
        // resists cell deformation more, suppressing buckling
        const double nagai_honda_deformation = ecmStiffness;  // KEY parameter
        const double nagai_honda_membrane = 10.0;
        const double nagai_honda_adhesion = 1.0;
        const double nagai_honda_boundary_adhesion = 2.0;

        // Basement membrane
        const double bm_stiffness = ecmStiffness * 0.5;  // scales with ECM
        const double bm_radius = 24.0;

        // ECM degradation — allows organoid to grow beyond initial BM radius
        const double ecm_degradation_rate = 0.05;  // radius units per hour
        const double ecm_max_radius = 50.0;         // biological ceiling

        // Target area
        const double target_area = 1.0;

        // Sloughing
        const double slough_height = 28.0;

        // ================================================================
        // SETUP
        // ================================================================

        RandomNumberGenerator::Instance()->Reseed(randomSeed);

        // Create vertex mesh
        HoneycombVertexMeshGenerator generator(cells_across, cells_up,
                                                flat_bottom);
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh = generator.GetMesh();

        // Set T1 and T2 thresholds
        p_mesh->SetCellRearrangementThreshold(0.1);
        p_mesh->SetT2Threshold(0.01);
        p_mesh->SetCellRearrangementRatio(1.5);

        // ================================================================
        // CREATE CELLS
        // ================================================================

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        // Get mesh bounds for cell type assignment
        double max_y_mesh = 0.0;
        for (unsigned i = 0; i < p_mesh->GetNumElements(); i++)
        {
            c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(i);
            if (centroid[1] > max_y_mesh) max_y_mesh = centroid[1];
        }

        for (unsigned i = 0; i < p_mesh->GetNumElements(); i++)
        {
            // Contact inhibition cell cycle
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(0.6);
            p_cycle->SetEquilibriumVolume(target_area);

            CellPtr p_cell(new Cell(p_state, p_cycle));

            // Assign type based on y-position of element centroid
            c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(i);
            double y_frac = centroid[1] / max_y_mesh;

            if (y_frac < 0.25)
            {
                p_cell->SetCellProliferativeType(p_stem_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 0.0);
            }
            else if (y_frac < 0.6)
            {
                p_cell->SetCellProliferativeType(p_transit_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 1.0);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff_type);
                p_cell->GetCellData()->SetItem("cell_type_id", 2.0);
            }

            // Random birth time
            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();

            // Target area and volume data
            p_cell->GetCellData()->SetItem("target area", target_area);
            p_cell->GetCellData()->SetItem("volume", target_area);

            // BM stiffness in cell data
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", bm_stiffness);

            cells.push_back(p_cell);
        }

        // ================================================================
        // CREATE VERTEX-BASED CELL POPULATION
        // ================================================================

        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Writers
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // ================================================================
        // SIMULATION
        // ================================================================

        OffLatticeSimulation<2> simulator(cell_population);

        std::string output_dir = "CryptBudding2d_VertexBased/" + outputSubDir;
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);

        // ================================================================
        // FORCES
        // ================================================================

        // 1. Nagai-Honda vertex force — cell mechanics
        //    Area elasticity parameter = ECM stiffness effect
        MAKE_PTR(NagaiHondaForce<2>, p_nh_force);
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation);
        p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane);
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_adhesion);
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_boundary_adhesion);
        simulator.AddForce(p_nh_force);

        // 2. Basement membrane restoring force
        //    Stiffness co-varies with ECM stiffness
        MAKE_PTR(BasementMembraneForce<2>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetBasementMembraneRadius(bm_radius);
        // Enable ECM degradation so the organoid can keep growing
        p_bm_force->EnableEcmDegradation(ecm_degradation_rate, ecm_max_radius);
        simulator.AddForce(p_bm_force);

        // ================================================================
        // BOUNDARY CONDITIONS
        // ================================================================

        // Pin bottom boundary
        c_vector<double, 2> bottom_point = zero_vector<double>(2);
        c_vector<double, 2> bottom_normal = zero_vector<double>(2);
        bottom_normal[1] = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc,
                      (&cell_population, bottom_point, bottom_normal));
        p_bc->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // ================================================================
        // CELL KILLERS
        // ================================================================

        // Slough at top
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

        // Target area modifier — maintains cell area homeostasis
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_area_modifier);
        p_area_modifier->SetReferenceTargetArea(target_area);
        simulator.AddSimulationModifier(p_area_modifier);

        // Volume tracking for contact inhibition
        MAKE_PTR(VolumeTrackingModifier<2>, p_vol_modifier);
        simulator.AddSimulationModifier(p_vol_modifier);

        // Custom summary writer
        boost::shared_ptr<VertexCryptSummaryModifier<2>> p_summary(
            new VertexCryptSummaryModifier<2>(ecmStiffness, sampling_multiple));
        simulator.AddSimulationModifier(p_summary);

        // ================================================================
        // RUN
        // ================================================================

        std::cout << "============================================" << std::endl;
        std::cout << "  2D Vertex-Based Crypt Budding Simulation" << std::endl;
        std::cout << "  ECM Stiffness:   " << ecmStiffness << std::endl;
        std::cout << "    → NH deformation: " << nagai_honda_deformation << std::endl;
        std::cout << "    → BM stiffness:   " << bm_stiffness << std::endl;
        std::cout << "  Random Seed:     " << randomSeed << std::endl;
        std::cout << "  Mesh:            " << cells_across << " x " << cells_up << std::endl;
        std::cout << "  End Time:        " << end_time << " hours" << std::endl;
        std::cout << "  Output:          testoutput/" << output_dir << std::endl;
        std::cout << "============================================" << std::endl;

        simulator.Solve();

        // ================================================================
        // POST-SIMULATION
        // ================================================================

        unsigned final_cells = cell_population.GetNumRealCells();

        unsigned stem_count = 0, ta_count = 0, diff_count = 0;
        double max_y = -1e10, min_y = 1e10;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double type_id = cell_iter->GetCellData()->GetItem("cell_type_id");
            if (type_id < 0.5) stem_count++;
            else if (type_id < 1.5) ta_count++;
            else diff_count++;

            c_vector<double, 2> centroid =
                cell_population.GetLocationOfCellCentre(*cell_iter);
            if (centroid[1] > max_y) max_y = centroid[1];
            if (centroid[1] < min_y) min_y = centroid[1];
        }

        std::cout << "\n============================================" << std::endl;
        std::cout << "  SIMULATION COMPLETE" << std::endl;
        std::cout << "  Final cells:     " << final_cells << std::endl;
        std::cout << "    Stem:          " << stem_count << std::endl;
        std::cout << "    TA:            " << ta_count << std::endl;
        std::cout << "    Diff:          " << diff_count << std::endl;
        std::cout << "  Y range:         [" << min_y << ", " << max_y << "]" << std::endl;
        std::cout << "  ECM stiffness:   " << ecmStiffness << std::endl;
        std::cout << "============================================" << std::endl;

        TS_ASSERT_LESS_THAN(0u, final_cells);
    }

public:

    /**
     * Parameterised test: reads ECM_STIFFNESS and RUN_NUMBER from
     * environment variables (set by HPC array job script).
     *
     * Default: stiffness = 5.0, run = 0.
     */
    void TestVertexCryptBuddingStiffnessSweep()
    {
        double ecm_stiffness = 5.0;
        const char* stiffness_env = std::getenv("ECM_STIFFNESS");
        if (stiffness_env != nullptr)
        {
            ecm_stiffness = std::atof(stiffness_env);
        }

        unsigned run_number = 0;
        const char* run_env = std::getenv("RUN_NUMBER");
        if (run_env != nullptr)
        {
            run_number = std::atoi(run_env);
        }

        unsigned seed = static_cast<unsigned>(ecm_stiffness * 10000) + run_number * 137;

        std::stringstream subdir;
        subdir << "stiffness_" << std::fixed << std::setprecision(1) << ecm_stiffness
               << "/run_" << run_number;

        RunVertexCryptBuddingSimulation(ecm_stiffness, seed, subdir.str());
    }
};

#endif /* TEST2DCRYPTBUDDINGVERTEXBASED_HPP_ */
