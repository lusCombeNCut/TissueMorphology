/*
 * Test2dCryptBuddingVertexBased.hpp
 *
 * 2D vertex-based simulations of crypt budding from an organoid.
 *
 * Initial condition: cells arranged as an annular ring of wedge-shaped
 * vertex elements, representing an early organoid cross-section with a
 * central lumen. The apical surface faces inward (lumen) and the basal
 * surface faces outward (ECM/basement membrane).
 *
 * Key independent variable: ECM_STIFFNESS (environment variable)
 *   - Controls basement membrane stiffness & Nagai-Honda area elasticity
 *   - Higher stiffness -> stronger confinement -> fewer/shallower buds
 *   - Lower stiffness -> easier budding -> more crypts
 *
 * Physics:
 *   - NagaiHondaForce: cell mechanics (area + perimeter energy)
 *   - BasementMembraneForce: inward radial restoring force (ECM)
 *   - LumenPressureForce: outward radial pressure maintaining lumen
 *   - ApicalConstrictionForce: actomyosin contraction on lumen-facing cells
 *   - ContactInhibitionCellCycleModel: density-dependent proliferation
 *
 * References:
 *   - Nagai & Honda (2001) Phil. Mag. B (vertex model)
 *   - Hannezo et al. (2011) PRL (epithelial buckling)
 *   - Serra et al. (2019) Nature (organoid symmetry breaking)
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
#include <numeric>
#include <iomanip>

// Core Chaste
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"

// Vertex mesh & population
#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"

// Simulation
#include "OffLatticeSimulation.hpp"

// Forces
#include "NagaiHondaForce.hpp"

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
#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellVolumesWriter.hpp"

// Project-specific forces
#include "BasementMembraneForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"

#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"
#include "PetscSetupAndFinalize.hpp"


/**
 * Summary modifier for vertex organoid — tracks radial statistics:
 * mean/max/min radius from centroid, mean cell area, etc.
 */
template<unsigned DIM>
class VertexOrganoidSummaryModifier : public AbstractCellBasedSimulationModifier<DIM>
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
    VertexOrganoidSummaryModifier(double stiffness, unsigned samplingMultiple,
                                  double endTime = 200.0)
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
            centroid += rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }
        if (num_cells > 0) centroid /= static_cast<double>(num_cells);

        // Compute radial statistics and mean area
        double sum_r = 0.0, sum_r2 = 0.0;
        double max_r = -1e10, min_r = 1e10;
        double sum_area = 0.0;

        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double r = norm_2(pos - centroid);
            sum_r += r;
            sum_r2 += r * r;
            if (r > max_r) max_r = r;
            if (r < min_r) min_r = r;

            if (p_vertex_pop)
            {
                unsigned elem_idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                double area = p_vertex_pop->rGetMesh().GetVolumeOfElement(elem_idx);
                sum_area += area;
            }
        }

        double mean_r = (num_cells > 0) ? sum_r / num_cells : 0.0;
        double var_r = (num_cells > 1) ? (sum_r2 / num_cells - mean_r * mean_r) : 0.0;
        double mean_area = (num_cells > 0) ? sum_area / num_cells : 0.0;

        // Write CSV
        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory()
                               + mOutputDir + "/crypt_summary.csv";
        std::ofstream file;
        if (!mHeaderWritten)
        {
            file.open(filename.c_str());
            file << "time,num_cells,mean_r,var_r,max_r,min_r,r_range,mean_area,centroid_x,centroid_y,stiffness"
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
             << mean_area << ","
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


class Test2dCryptBuddingVertexBased : public AbstractCellBasedTestSuite
{
private:

    /**
     * Build a MutableVertexMesh<2,2> as an annular ring of wedge-shaped
     * elements. Each element is a quadrilateral with inner (apical) and
     * outer (basal) edges.
     *
     * @param numElements  number of cells in the ring
     * @param innerRadius  apical (lumen) radius
     * @param outerRadius  basal (ECM) radius
     * @param cellRearrangementThreshold  T1 swap threshold
     * @param t2Threshold  T2 swap threshold
     * @return  shared pointer to the constructed mesh
     */
    boost::shared_ptr<MutableVertexMesh<2,2> > MakeAnnularVertexMesh(
        unsigned numElements,
        double innerRadius,
        double outerRadius,
        double cellRearrangementThreshold = 0.1,
        double t2Threshold = 0.01)
    {
        // Create vertex nodes: 2 rings of nodes (inner + outer)
        // plus shared corners between adjacent wedges
        std::vector<Node<2>*> nodes;

        double dtheta = 2.0 * M_PI / static_cast<double>(numElements);

        // Inner ring nodes: indices 0 .. numElements-1
        for (unsigned i = 0; i < numElements; i++)
        {
            double theta = i * dtheta;
            double x = innerRadius * cos(theta);
            double y = innerRadius * sin(theta);
            nodes.push_back(new Node<2>(i, true, x, y));  // boundary
        }

        // Outer ring nodes: indices numElements .. 2*numElements-1
        for (unsigned i = 0; i < numElements; i++)
        {
            double theta = i * dtheta;
            double x = outerRadius * cos(theta);
            double y = outerRadius * sin(theta);
            nodes.push_back(new Node<2>(numElements + i, true, x, y));  // boundary
        }

        // Create elements: each is a quadrilateral
        // Vertices: inner[i], outer[i], outer[i+1], inner[i+1] (CCW)
        std::vector<VertexElement<2,2>*> elements;

        for (unsigned i = 0; i < numElements; i++)
        {
            unsigned next_i = (i + 1) % numElements;

            std::vector<Node<2>*> element_nodes;
            element_nodes.push_back(nodes[i]);                    // inner left
            element_nodes.push_back(nodes[numElements + i]);      // outer left
            element_nodes.push_back(nodes[numElements + next_i]); // outer right
            element_nodes.push_back(nodes[next_i]);               // inner right

            elements.push_back(new VertexElement<2,2>(i, element_nodes));
        }

        // Create the mutable vertex mesh
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh(
            new MutableVertexMesh<2,2>(nodes, elements,
                                       cellRearrangementThreshold,
                                       t2Threshold));

        p_mesh->SetCellRearrangementRatio(1.5);
        p_mesh->SetProtorosetteFormationProbability(0.0);
        p_mesh->SetProtorosetteResolutionProbabilityPerTimestep(1.0);
        p_mesh->SetCheckForInternalIntersections(true);

        return p_mesh;
    }


    /**
     * Run a single vertex-based organoid crypt budding simulation.
     *
     * @param ecmStiffness  ECM stiffness (varies BM + area elasticity)
     * @param randomSeed    RNG seed for reproducibility
     * @param outputSubDir  sub-directory for output
     */
    void RunVertexCryptBuddingSimulation(double ecmStiffness, unsigned randomSeed,
                                          const std::string& outputSubDir)
    {
        // ================================================================
        // PARAMETERS
        // ================================================================

        // Organoid geometry — annular ring (circular cross-section)
        const unsigned num_cells = 40;
        const double inner_radius = 6.0;    // apical / lumen
        const double outer_radius = 8.0;    // basal / ECM
        const double organoid_radius = 0.5 * (inner_radius + outer_radius); // mid-shell

        // Time (hours)
        // Adaptive dt for stability — annular vertex mesh needs smaller steps
        const double dt = (ecmStiffness < 1.0) ? 0.0002
                        : (ecmStiffness < 2.0) ? 0.0005
                        : (ecmStiffness < 5.0) ? 0.0005
                        :                        0.0005;
        const double end_time = 168.0;      // ~8 days
        const unsigned sampling_multiple = static_cast<unsigned>(1.0 / dt);

        // Nagai-Honda force parameters
        const double nagai_honda_deformation = ecmStiffness;  // KEY parameter
        const double nagai_honda_membrane = 10.0;
        const double nagai_honda_adhesion = 1.0;
        const double nagai_honda_boundary_adhesion = 2.0;

        // Basement membrane (radial confinement)
        const double bm_stiffness = ecmStiffness * 0.5;
        const double bm_radius = outer_radius + 2.0;

        // ECM degradation
        const double ecm_degradation_rate = 0.02;
        const double ecm_max_radius = outer_radius * 4.0;

        // Lumen pressure
        const double lumen_pressure = 2.0;
        const double lumen_eq_radius = outer_radius + 1.0;  // above outer_radius so outward pressure at startup

        // Apical constriction
        const double apical_constriction_strength = 3.0;

        // Target area — approximate area of a wedge element
        const double dtheta = 2.0 * M_PI / static_cast<double>(num_cells);
        const double target_area = 0.5 * dtheta * (outer_radius * outer_radius - inner_radius * inner_radius);

        // Cell cycle
        const double quiescent_fraction = 0.6;

        // Sloughing bounding box
        const double max_radius_for_slough = outer_radius * 5.0;

        // T1/T2 thresholds — larger values improve mesh stability
        // by triggering rearrangements before degenerate configurations form
        double t1_threshold = (ecmStiffness < 2.0) ? 0.2 : 0.15;
        double t2_threshold = 0.05;

        // ================================================================
        // SETUP
        // ================================================================

        RandomNumberGenerator::Instance()->Reseed(randomSeed);

        // Build annular vertex mesh
        boost::shared_ptr<MutableVertexMesh<2,2> > p_mesh =
            MakeAnnularVertexMesh(num_cells, inner_radius, outer_radius,
                                  t1_threshold, t2_threshold);

        // Add slight radial noise to break perfect symmetry
        for (unsigned i = 0; i < p_mesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> pos = p_mesh->GetNode(i)->rGetLocation();
            double r = norm_2(pos);
            if (r > 1e-6)
            {
                double noise = (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.2;
                c_vector<double, 2> unit_r = pos / r;
                pos += noise * unit_r;
                ChastePoint<2> new_point(pos[0], pos[1]);
                p_mesh->GetNode(i)->SetPoint(new_point);
            }
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

        for (unsigned i = 0; i < p_mesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(quiescent_fraction);
            p_cycle->SetEquilibriumVolume(target_area);

            CellPtr p_cell(new Cell(p_state, p_cycle));

            // Assign cell type based on angular position of element centroid
            c_vector<double, 2> centroid = p_mesh->GetCentroidOfElement(i);
            double theta = atan2(centroid[1], centroid[0]);

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

            // Random birth time
            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();

            // Target area & volume data (for contact inhibition)
            p_cell->GetCellData()->SetItem("target area", target_area);
            p_cell->GetCellData()->SetItem("volume", target_area);

            // BM stiffness per cell
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", bm_stiffness);

            // All cells are apical (face lumen)
            p_cell->GetCellData()->SetItem("is_apical", 1.0);

            cells.push_back(p_cell);
        }

        // ================================================================
        // CELL POPULATION
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
        simulator.SetOutputDirectory(outputSubDir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_multiple);
        simulator.SetEndTime(end_time);

        // ================================================================
        // FORCES
        // ================================================================

        // 1. Nagai-Honda vertex force — cell mechanics
        MAKE_PTR(NagaiHondaForce<2>, p_nh_force);
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(nagai_honda_deformation);
        p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(nagai_honda_membrane);
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(nagai_honda_adhesion);
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(nagai_honda_boundary_adhesion);
        simulator.AddForce(p_nh_force);

        // 2. Basement membrane — radial confinement (KEY parameter)
        MAKE_PTR(BasementMembraneForce<2>, p_bm_force);
        p_bm_force->SetBasementMembraneParameter(bm_stiffness);
        p_bm_force->SetBasementMembraneRadius(bm_radius);
        c_vector<double, 2> center_vec;
        center_vec[0] = 0.0;
        center_vec[1] = 0.0;
        p_bm_force->SetOrganoidCenter(center_vec);
        p_bm_force->EnableEcmDegradation(ecm_degradation_rate, ecm_max_radius);
        simulator.AddForce(p_bm_force);

        // 3. Lumen pressure — outward radial force maintaining lumen
        MAKE_PTR(LumenPressureForce<2>, p_lumen_force);
        p_lumen_force->SetPressureStrength(lumen_pressure);
        p_lumen_force->SetLumenEquilibriumRadius(lumen_eq_radius);
        p_lumen_force->SetTrackCenter(true);
        simulator.AddForce(p_lumen_force);

        // 4. Apical constriction — drives invagination
        MAKE_PTR(ApicalConstrictionForce<2>, p_apical_force);
        p_apical_force->SetConstrictionStrength(apical_constriction_strength);
        simulator.AddForce(p_apical_force);

        // ================================================================
        // CELL KILLERS — bounding box sloughing
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

        // Target area modifier — maintains cell area homeostasis
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_area_modifier);
        p_area_modifier->SetReferenceTargetArea(target_area);
        simulator.AddSimulationModifier(p_area_modifier);

        // Volume tracking for contact inhibition
        MAKE_PTR(VolumeTrackingModifier<2>, p_vol_modifier);
        simulator.AddSimulationModifier(p_vol_modifier);

        // Custom summary writer
        boost::shared_ptr<VertexOrganoidSummaryModifier<2>> p_summary(
            new VertexOrganoidSummaryModifier<2>(ecmStiffness, sampling_multiple, end_time));
        simulator.AddSimulationModifier(p_summary);

        // ================================================================
        // RUN
        // ================================================================

        std::cout << "============================================" << std::endl;
        std::cout << "  2D Vertex-Based Organoid Crypt Budding" << std::endl;
        std::cout << "  ECM Stiffness:    " << ecmStiffness << std::endl;
        std::cout << "    -> NH deform:   " << nagai_honda_deformation << std::endl;
        std::cout << "    -> BM stiffness:" << bm_stiffness << std::endl;
        std::cout << "  Lumen Pressure:   " << lumen_pressure << std::endl;
        std::cout << "  Apical Constr.:   " << apical_constriction_strength << std::endl;
        std::cout << "  Inner/Outer R:    " << inner_radius << " / " << outer_radius << std::endl;
        std::cout << "  Cells in Ring:    " << num_cells << std::endl;
        std::cout << "  Target Area:      " << target_area << std::endl;
        std::cout << "  Random Seed:      " << randomSeed << std::endl;
        std::cout << "  dt:               " << dt << std::endl;
        std::cout << "  T1 threshold:     " << t1_threshold << std::endl;
        std::cout << "  End Time:         " << end_time << " hours" << std::endl;
        std::cout << "  Output:           testoutput/" << outputSubDir << std::endl;
        std::cout << "============================================" << std::endl;

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
            final_centroid += cell_population.GetLocationOfCellCentre(*cell_iter);
        }
        if (final_cells > 0) final_centroid /= static_cast<double>(final_cells);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End(); ++cell_iter)
        {
            double type_id = cell_iter->GetCellData()->GetItem("cell_type_id");
            if (type_id < 0.5) stem_count++;
            else if (type_id < 1.5) ta_count++;
            else diff_count++;

            c_vector<double, 2> pos = cell_population.GetLocationOfCellCentre(*cell_iter);
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
        std::cout << "  ECM stiffness:    " << ecmStiffness << std::endl;
        std::cout << "============================================" << std::endl;

        TS_ASSERT_LESS_THAN(0u, final_cells);
    }

public:

    /**
     * Parameterised test: reads ECM_STIFFNESS and RUN_NUMBER from
     * environment variables (set by HPC array job script).
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
        subdir << "CryptBudding2d_VertexBased/stiffness_" << std::fixed
               << std::setprecision(1) << ecm_stiffness << "/run_" << run_number;

        RunVertexCryptBuddingSimulation(ecm_stiffness, seed, subdir.str());
    }
};

#endif /* TEST2DCRYPTBUDDINGVERTEXBASED_HPP_ */
