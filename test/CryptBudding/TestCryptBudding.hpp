/*
 * TestCryptBudding.hpp
 *
 * Unified crypt budding / organoid simulation test.
 *
 * Selects model type and optional force modules via environment variables
 * so that a single test executable drives all four model configurations:
 *
 *   MODEL_TYPE        = node2d | vertex2d | node3d | vertex3d
 *   ECM_STIFFNESS     = <double>  (default 5.0)
 *   RUN_NUMBER        = <int>     (default 0, used for random seed)
 *
 * Feature toggles (0 = off, 1 = on; defaults shown):
 *   ENABLE_LUMEN_PRESSURE        = 1
 *   ENABLE_APICAL_CONSTRICTION   = 1
 *   ENABLE_ECM_GUIDANCE          = 0  (3D only)
 *   ENABLE_RELAXATION            = 1
 *   ENABLE_SLOUGHING             = 1  (2D only)
 *   ENABLE_DIFFERENTIAL_ADHESION = 1  (node-based only)
 *
 * All models share:
 *   - ContactInhibitionCellCycleModel (density-dependent proliferation)
 *   - Three-tier cell type assignment (Stem / TA / Differentiated)
 *   - BasementMembraneForce with ECM degradation
 *   - VolumeTrackingModifier
 *   - Relaxation phase (geometry settling before proliferation)
 *   - Summary CSV writer and progress logging
 *
 * Model-specific cell–cell force:
 *   - node2d  : DifferentialAdhesionForce (linear spring, 2D)
 *   - vertex2d: NagaiHondaForce (area + perimeter energy, 2D)
 *   - node3d  : DifferentialAdhesionForce (linear spring, 3D)
 *   - vertex3d: SurfaceTensionForce (OrganoidChaste finite-thickness, 3D)
 *
 * NagaiHondaForce vs SurfaceTensionForce:
 *   NagaiHondaForce is Chaste's built-in 2D vertex model (Nagai & Honda 2001).
 *   SurfaceTensionForce is OrganoidChaste's 3D finite-thickness monolayer
 *   vertex model (Drozdowski & Schwarz 2025) with per-face surface tensions
 *   (apical, basal, lateral) and volume elasticity. They cannot substitute
 *   for each other: NagaiHondaForce operates on MutableVertexMesh<2,2> /
 *   VertexBasedCellPopulation<2>, while SurfaceTensionForce requires
 *   MutableMonolayerVertexMesh<3,3> / MonolayerVertexBasedCellPopulation<3>.
 *   Both are energy-based vertex models, but SurfaceTensionForce is the
 *   natural 3D extension for finite-thickness epithelia.
 */

#ifndef TESTCRYPTBUDDING_HPP_
#define TESTCRYPTBUDDING_HPP_

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
#include <chrono>

// ===================================================================
// Core Chaste
// ===================================================================
#include "AbstractCellBasedTestSuite.hpp"
#include "CheckpointArchiveTypes.hpp"
#include "SmartPointers.hpp"
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"

// ===================================================================
// Mesh & population — 2D node
// ===================================================================
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"

// ===================================================================
// Mesh & population — 2D vertex
// ===================================================================
#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"

// ===================================================================
// Mesh & population — 3D vertex (OrganoidChaste)
// ===================================================================
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "MutableMonolayerVertexMesh.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"

// ===================================================================
// Simulation
// ===================================================================
#include "OffLatticeSimulation.hpp"

// ===================================================================
// Forces — vertex
// ===================================================================
#include "NagaiHondaForce.hpp"
#include "SurfaceTensionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"

// ===================================================================
// Cell cycle
// ===================================================================
#include "ContactInhibitionCellCycleModel.hpp"
#include "NoCellCycleModel.hpp"

// ===================================================================
// Cell types & states
// ===================================================================
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

// ===================================================================
// Cell killers
// ===================================================================
#include "PlaneBasedCellKiller.hpp"

// ===================================================================
// Modifiers & writers
// ===================================================================
#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "FaceTypeWriter.hpp"

// ===================================================================
// TissueMorphology project forces
// ===================================================================
#include "DifferentialAdhesionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"

// ===================================================================
// TissueMorphology & OrganoidChaste — 3D ECM (optional)
// ===================================================================
#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"
#include "ECMFieldWriter3d.hpp"

// OrganoidChaste sub-force for vertex3d lumen pressure
#include "LumenPressureSubForce.hpp"

#include "PetscSetupAndFinalize.hpp"


// ===================================================================
// Shared parameter block — single source of truth
// ===================================================================

/**
 * Parameters used across all four model types.
 * Values are set once in the public test method and passed through.
 */
struct CryptBuddingParams
{
    // --- Environment-driven ---
    std::string modelType;
    double ecmStiffness;
    unsigned runNumber;
    unsigned randomSeed;

    // --- Feature toggles ---
    bool enableLumenPressure;
    bool enableApicalConstriction;
    bool enableEcmGuidance;
    bool enableRelaxation;
    bool enableSloughing;
    bool enableDifferentialAdhesion;

    // --- Time ---
    double dt;
    double relaxationTime;
    double endTime;
    unsigned samplingMultiple;

    // --- Geometry: 2D ring ---
    double organoidRadius2d;
    unsigned numCells2dNode;
    unsigned numCells2dVertex;
    double innerRadius2d;
    double outerRadius2d;
    double interactionCutoff2d;

    // --- Geometry: 3D sphere ---
    unsigned numCells3dNode;
    unsigned numCells3dVertex;
    double organoidRadius3d;
    double shellThickness3d;
    double interactionCutoff3d;
    double sphereRadius3dVertex;

    // --- Basement membrane ---
    double bmStiffnessNode;
    double bmStiffnessVertex;
    double bmRadius2d;
    double bmRadius3d;
    double ecmDegradationRate;
    double ecmMaxRadius2d;
    double ecmMaxRadius3d;

    // --- Lumen pressure ---
    double lumenPressure;
    double lumenEqRadius2d;

    // --- Apical constriction ---
    double apicalConstrictionStrength;

    // --- Cell–cell springs (node-based) ---
    double springStiffness;
    double springCutoff;
    double apicalApicalAdhesion;
    double basalBasalAdhesion;
    double apicalBasalAdhesion;

    // --- Nagai–Honda (vertex2d) ---
    double nhMembraneSurface;
    double nhCellCellAdhesion;
    double nhBoundaryAdhesion;

    // --- Surface tension (vertex3d) ---
    double gammaApical;
    double gammaBasal;
    double gammaLateral;

    // --- Cell cycle ---
    double quiescentFraction;

    // --- Sloughing ---
    double sloughRadiusFactor;

    // --- ECM guidance (3D) ---
    double ecmDomainHalf;
    double ecmGridSpacing;
    double ecmBaseSpeed;

    // --- T1/T2 thresholds (2D vertex) ---
    double t1Threshold2d;
    double t2Threshold2d;

    /**
     * Populate with defaults, then override from env vars.
     */
    void InitFromEnvironment()
    {
        // --- Read env vars ---
        const char* e;

        e = std::getenv("MODEL_TYPE");
        modelType = (e != nullptr) ? std::string(e) : "node2d";

        e = std::getenv("ECM_STIFFNESS");
        ecmStiffness = (e != nullptr) ? std::atof(e) : 5.0;

        e = std::getenv("RUN_NUMBER");
        runNumber = (e != nullptr) ? static_cast<unsigned>(std::atoi(e)) : 0u;

        e = std::getenv("ENABLE_LUMEN_PRESSURE");
        enableLumenPressure = (e != nullptr) ? (std::atoi(e) != 0) : true;

        e = std::getenv("ENABLE_APICAL_CONSTRICTION");
        enableApicalConstriction = (e != nullptr) ? (std::atoi(e) != 0) : true;

        e = std::getenv("ENABLE_ECM_GUIDANCE");
        enableEcmGuidance = (e != nullptr) ? (std::atoi(e) != 0) : false;

        e = std::getenv("ENABLE_RELAXATION");
        enableRelaxation = (e != nullptr) ? (std::atoi(e) != 0) : true;

        e = std::getenv("ENABLE_SLOUGHING");
        enableSloughing = (e != nullptr) ? (std::atoi(e) != 0) : true;

        e = std::getenv("ENABLE_DIFFERENTIAL_ADHESION");
        enableDifferentialAdhesion = (e != nullptr) ? (std::atoi(e) != 0) : true;

        // --- Seed ---
        randomSeed = static_cast<unsigned>(ecmStiffness * 10000) + runNumber * 137;

        // --- Time ---
        relaxationTime = 10.0;
        endTime        = 168.0;

        // --- 2D ring geometry ---
        organoidRadius2d    = 8.0;
        numCells2dNode      = 80;
        numCells2dVertex    = 40;
        innerRadius2d       = 6.0;
        outerRadius2d       = 8.0;
        interactionCutoff2d = 2.5;

        // --- 3D sphere geometry ---
        numCells3dNode       = 100;
        numCells3dVertex     = 200;
        organoidRadius3d     = 25.0;
        shellThickness3d     = 3.0;
        interactionCutoff3d  = 15.0;
        sphereRadius3dVertex = 10.0;

        // --- Basement membrane ---
        bmStiffnessNode   = ecmStiffness;
        bmStiffnessVertex = ecmStiffness * 0.5;
        bmRadius2d        = organoidRadius2d + 2.0;   // node2d uses outerRadius2d + 2 for vertex2d
        bmRadius3d        = 30.0;
        ecmDegradationRate = 0.02;
        ecmMaxRadius2d     = organoidRadius2d * 4.0;
        ecmMaxRadius3d     = 80.0;

        // --- Lumen pressure ---
        lumenPressure   = 2.0;
        lumenEqRadius2d = organoidRadius2d + 1.0;

        // --- Apical constriction ---
        apicalConstrictionStrength = 3.0;

        // --- Spring forces (node-based) ---
        springStiffness       = 30.0;
        springCutoff          = 1.5;
        apicalApicalAdhesion  = 1.2;
        basalBasalAdhesion    = 1.0;
        apicalBasalAdhesion   = 0.5;

        // --- Nagai–Honda (vertex2d) ---
        nhMembraneSurface = 10.0;
        nhCellCellAdhesion = 1.0;
        nhBoundaryAdhesion = 2.0;

        // --- Surface tension (vertex3d) ---
        gammaApical  = 0.85;
        gammaBasal   = 0.85;
        gammaLateral = 0.7;

        // --- Cell cycle ---
        quiescentFraction = 0.7;

        // --- Sloughing ---
        sloughRadiusFactor = 5.0;

        // --- ECM guidance (3D) ---
        ecmDomainHalf  = 80.0;
        ecmGridSpacing = 10.0;
        ecmBaseSpeed   = 0.3;

        // --- T1/T2 (vertex2d) ---
        t1Threshold2d = (ecmStiffness < 2.0) ? 0.2 : 0.15;
        t2Threshold2d = 0.05;

        // --- Model-specific dt & sampling ---
        if (modelType == "node2d")
        {
            dt = 0.005;
            samplingMultiple = 200;
        }
        else if (modelType == "vertex2d")
        {
            dt = (ecmStiffness < 1.0) ? 0.0002
               : (ecmStiffness < 5.0) ? 0.0005
               :                        0.0005;
            samplingMultiple = static_cast<unsigned>(1.0 / dt);
            endTime = 168.0;
        }
        else if (modelType == "node3d")
        {
            dt = 0.01;
            samplingMultiple = 20;
            endTime = 168.0;
        }
        else if (modelType == "vertex3d")
        {
            dt = 0.001;              // relaxation dt (growth phase bumps to 0.006)
            samplingMultiple = 20;
            endTime = 100.0;
        }
        else
        {
            dt = 0.005;
            samplingMultiple = 200;
        }
    }
};


// ===================================================================
// Summary modifier — writes radial statistics to CSV
// ===================================================================

template<unsigned DIM>
class CryptBuddingSummaryModifier : public AbstractCellBasedSimulationModifier<DIM>
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
    CryptBuddingSummaryModifier(double stiffness, unsigned samplingMultiple,
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
        double current_time  = SimulationTime::Instance()->GetTime();
        unsigned num_cells   = rCellPopulation.GetNumRealCells();

        // Progress log — use std::cerr so output is visible in real time
        // (CTest captures and buffers stdout, only showing it after completion)
        if (current_step == 0 || current_step - mLastLogStep >= mLogInterval)
        {
            mLastLogStep = current_step;
            double pct = (mEndTime > 0.0) ? (current_time / mEndTime) * 100.0 : 0.0;
            std::cerr << "[Progress] t=" << std::fixed << std::setprecision(1)
                      << current_time << "h / " << mEndTime << "h  ("
                      << std::setprecision(1) << pct << "%)  cells="
                      << num_cells << std::endl;
        }

        // CSV at sampling interval
        if (current_step - mLastOutputStep < mSamplingMultiple && current_step > 0)
            return;
        mLastOutputStep = current_step;

        // Compute centroid
        c_vector<double, DIM> centroid = zero_vector<double>(DIM);
        for (typename AbstractCellPopulation<DIM>::Iterator it = rCellPopulation.Begin();
             it != rCellPopulation.End(); ++it)
        {
            centroid += rCellPopulation.GetLocationOfCellCentre(*it);
        }
        if (num_cells > 0) centroid /= static_cast<double>(num_cells);

        double sum_r = 0.0, sum_r2 = 0.0;
        double max_r = -1e10, min_r = 1e10;

        for (typename AbstractCellPopulation<DIM>::Iterator it = rCellPopulation.Begin();
             it != rCellPopulation.End(); ++it)
        {
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*it);
            double r = norm_2(pos - centroid);
            sum_r  += r;
            sum_r2 += r * r;
            if (r > max_r) max_r = r;
            if (r < min_r) min_r = r;
        }

        double mean_r = (num_cells > 0) ? sum_r / num_cells : 0.0;
        double var_r  = (num_cells > 1) ? (sum_r2 / num_cells - mean_r * mean_r) : 0.0;

        // Write CSV
        std::string filename = OutputFileHandler::GetChasteTestOutputDirectory()
                               + mOutputDir + "/crypt_summary.csv";
        std::ofstream file;
        if (!mHeaderWritten)
        {
            file.open(filename.c_str());
            file << "time,num_cells,mean_r,var_r,max_r,min_r,r_range,stiffness" << std::endl;
            mHeaderWritten = true;
        }
        else
        {
            file.open(filename.c_str(), std::ios::app);
        }

        file << std::fixed << std::setprecision(4)
             << current_time << "," << num_cells << ","
             << mean_r << "," << var_r << ","
             << max_r << "," << min_r << ","
             << (max_r - min_r) << "," << mStiffness
             << std::endl;
        file.close();
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t<Stiffness>" << mStiffness << "</Stiffness>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};


// ===================================================================
// Utility: assign cell type by fractional position
//   frac ∈ [0,1], 0 = crypt base
//   Bottom 40% → Stem, flanks 30% → TA, top 30% → Diff
// ===================================================================

inline void AssignCellTypeByFraction(
    CellPtr pCell, double frac,
    boost::shared_ptr<AbstractCellProliferativeType> pStemType,
    boost::shared_ptr<AbstractCellProliferativeType> pTransitType,
    boost::shared_ptr<AbstractCellProliferativeType> pDiffType)
{
    if (frac < 0.2 || frac > 0.8)
    {
        pCell->SetCellProliferativeType(pStemType);
        pCell->GetCellData()->SetItem("cell_type_id", 0.0);
    }
    else if (frac < 0.35 || frac > 0.65)
    {
        pCell->SetCellProliferativeType(pTransitType);
        pCell->GetCellData()->SetItem("cell_type_id", 1.0);
    }
    else
    {
        pCell->SetCellProliferativeType(pDiffType);
        pCell->GetCellData()->SetItem("cell_type_id", 2.0);
    }
}


// ===================================================================
// Utility: angular fraction for a 2D position on the ring
//   Bottom of ring (θ ≈ -π/2) → 0, top (θ ≈ π/2) → 0.5
// ===================================================================

inline double AngularFraction2d(double x, double y,
                                double cx = 0.0, double cy = 0.0)
{
    double theta = atan2(y - cy, x - cx);
    double angle_from_bottom = theta + M_PI / 2.0;
    if (angle_from_bottom < 0.0) angle_from_bottom += 2.0 * M_PI;
    return angle_from_bottom / (2.0 * M_PI);
}


// ===================================================================
// Utility: z-fraction for a 3D position on the sphere
//   z_frac ∈ [-1, 1];  map to [0, 1] for AssignCellTypeByFraction
// ===================================================================

inline double ZFractionToRingFraction(double z, double R)
{
    // z/R ∈ [-1, 1].  Map: bottom (-1) → 0, top (+1) → 0.5
    // Use same thresholds as 2D: bottom 40% = stem, flanks 30% = TA, top 30% = diff
    double z_frac = z / R;   // [-1, 1]

    if (z_frac < -0.5)      return 0.0;   // bottom → stem
    else if (z_frac < 0.3)  return 0.3;   // mid → TA
    else                     return 0.5;   // top → diff
}


// ===================================================================
// Utility: add plane-based cell killers as bounding box
// ===================================================================

template<unsigned DIM>
void AddBoundingBoxKillers(OffLatticeSimulation<DIM>& rSimulator,
                           AbstractCellPopulation<DIM>& rPopulation,
                           double halfWidth)
{
    for (unsigned d = 0; d < DIM; d++)
    {
        // Positive plane
        {
            c_vector<double, DIM> point  = zero_vector<double>(DIM);
            c_vector<double, DIM> normal = zero_vector<double>(DIM);
            point[d]  =  halfWidth;
            normal[d] =  1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<DIM>, p_killer,
                          (&rPopulation, point, normal));
            rSimulator.AddCellKiller(p_killer);
        }
        // Negative plane
        {
            c_vector<double, DIM> point  = zero_vector<double>(DIM);
            c_vector<double, DIM> normal = zero_vector<double>(DIM);
            point[d]  = -halfWidth;
            normal[d] = -1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<DIM>, p_killer,
                          (&rPopulation, point, normal));
            rSimulator.AddCellKiller(p_killer);
        }
    }
}


// ===================================================================
// Utility: print banner
// ===================================================================

inline void PrintBanner(const CryptBuddingParams& p)
{
    std::cerr << "\n============================================" << std::endl;
    std::cerr << "  Crypt Budding — Unified Test" << std::endl;
    std::cerr << "  Model:          " << p.modelType << std::endl;
    std::cerr << "  ECM Stiffness:  " << p.ecmStiffness << std::endl;
    std::cerr << "  Run Number:     " << p.runNumber << std::endl;
    std::cerr << "  Seed:           " << p.randomSeed << std::endl;
    std::cerr << "  dt:             " << p.dt << std::endl;
    std::cerr << "  Relaxation:     " << (p.enableRelaxation ? "ON" : "OFF")
              << " (" << p.relaxationTime << "h)" << std::endl;
    std::cerr << "  End Time:       " << p.endTime << "h" << std::endl;
    std::cerr << "  Features:" << std::endl;
    std::cerr << "    Lumen Pressure:        " << (p.enableLumenPressure ? "ON" : "OFF") << std::endl;
    std::cerr << "    Apical Constriction:    " << (p.enableApicalConstriction ? "ON" : "OFF") << std::endl;
    std::cerr << "    ECM Guidance (3D):      " << (p.enableEcmGuidance ? "ON" : "OFF") << std::endl;
    std::cerr << "    Sloughing:              " << (p.enableSloughing ? "ON" : "OFF") << std::endl;
    std::cerr << "    Differential Adhesion:  " << (p.enableDifferentialAdhesion ? "ON" : "OFF") << std::endl;
    std::cerr << "============================================\n" << std::endl;
}


// ===================================================================
// Main test class
// ===================================================================

class TestCryptBudding : public AbstractCellBasedTestSuite
{
private:

    // ===============================================================
    // MODEL 1: 2D Node-Based
    // ===============================================================

    void RunNode2d(const CryptBuddingParams& p, const std::string& outputDir)
    {
        RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

        // --- Create annular ring of nodes ---
        std::vector<Node<2>*> nodes;
        for (unsigned i = 0; i < p.numCells2dNode; i++)
        {
            double theta = 2.0 * M_PI * i / p.numCells2dNode;
            double r_noise = p.organoidRadius2d +
                (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.3;
            double x = r_noise * cos(theta);
            double y = r_noise * sin(theta);
            nodes.push_back(new Node<2>(i, false, x, y));
        }

        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, p.interactionCutoff2d);
        for (unsigned i = 0; i < nodes.size(); i++) delete nodes[i];

        // --- Create cells ---
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem);
        MAKE_PTR(TransitCellProliferativeType, p_ta);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle));
            c_vector<double, 2> pos = mesh.GetNode(i)->rGetLocation();
            double frac = AngularFraction2d(pos[0], pos[1]);
            AssignCellTypeByFraction(p_cell, frac, p_stem, p_ta, p_diff);

            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("volume", 1.0);
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessNode);
            p_cell->GetCellData()->SetItem("is_apical", 1.0);
            cells.push_back(p_cell);
        }

        // --- Population ---
        NodeBasedCellPopulation<2> population(mesh, cells);
        population.SetAbsoluteMovementThreshold(50.0);
        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CellAgesWriter>();
        population.AddCellWriter<CellVolumesWriter>();
        population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // --- Simulation ---
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory(outputDir);
        simulator.SetDt(p.dt);
        simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
        simulator.SetEndTime(p.endTime);

        // --- Forces ---
        // Cell–cell
        if (p.enableDifferentialAdhesion)
        {
            MAKE_PTR(DifferentialAdhesionForce<2>, p_spring);
            p_spring->SetMeinekeSpringStiffness(p.springStiffness);
            p_spring->SetCutOffLength(p.springCutoff);
            p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
            p_spring->SetMeinekeSpringGrowthDuration(1.0);
            p_spring->SetApicalApicalAdhesion(p.apicalApicalAdhesion);
            p_spring->SetBasalBasalAdhesion(p.basalBasalAdhesion);
            p_spring->SetApicalBasalAdhesion(p.apicalBasalAdhesion);
            simulator.AddForce(p_spring);
        }
        else
        {
            MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring);
            p_spring->SetMeinekeSpringStiffness(p.springStiffness);
            p_spring->SetCutOffLength(p.springCutoff);
            simulator.AddForce(p_spring);
        }

        // Basement membrane
        MAKE_PTR(BasementMembraneForce<2>, p_bm);
        p_bm->SetBasementMembraneParameter(p.bmStiffnessNode);
        p_bm->SetBasementMembraneRadius(p.bmRadius2d);
        c_vector<double, 2> center2d = zero_vector<double>(2);
        p_bm->SetOrganoidCenter(center2d);
        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius2d);
        simulator.AddForce(p_bm);

        // Lumen pressure
        if (p.enableLumenPressure)
        {
            MAKE_PTR(LumenPressureForce<2>, p_lumen);
            p_lumen->SetPressureStrength(p.lumenPressure);
            p_lumen->SetLumenEquilibriumRadius(p.lumenEqRadius2d);
            p_lumen->SetTrackCenter(true);
            simulator.AddForce(p_lumen);
        }

        // Apical constriction
        if (p.enableApicalConstriction)
        {
            MAKE_PTR(ApicalConstrictionForce<2>, p_ac);
            p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
            simulator.AddForce(p_ac);
        }

        // Sloughing
        if (p.enableSloughing)
        {
            AddBoundingBoxKillers<2>(simulator, population,
                                    p.organoidRadius2d * p.sloughRadiusFactor);
        }

        // Modifiers
        MAKE_PTR(VolumeTrackingModifier<2>, p_vol);
        simulator.AddSimulationModifier(p_vol);

        boost::shared_ptr<CryptBuddingSummaryModifier<2>> p_summary(
            new CryptBuddingSummaryModifier<2>(p.ecmStiffness, p.samplingMultiple,
                                               p.relaxationTime + p.endTime));
        simulator.AddSimulationModifier(p_summary);

        // --- Run with optional relaxation ---
        if (p.enableRelaxation)
        {
            // Save original types, set all to differentiated
            std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
            for (AbstractCellPopulation<2>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                origTypes[*it] = it->GetCellProliferativeType();
                it->SetCellProliferativeType(p_diff);
            }

            simulator.SetEndTime(p.relaxationTime);
            std::cerr << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
            simulator.Solve();

            // Restore types
            for (AbstractCellPopulation<2>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
            }

            simulator.SetEndTime(p.relaxationTime + p.endTime);
            std::cerr << "--- Phase 2: Growth (" << p.endTime << "h) ---" << std::endl;
            simulator.Solve();
        }
        else
        {
            simulator.Solve();
        }

        // Post-simulation
        unsigned final_cells = population.GetNumRealCells();
        std::cerr << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
        TS_ASSERT_LESS_THAN(0u, final_cells);
    }


    // ===============================================================
    // MODEL 2: 2D Vertex-Based  (NagaiHondaForce)
    // ===============================================================

    /**
     * Build closed annular vertex mesh.
     */
    boost::shared_ptr<MutableVertexMesh<2,2>> MakeAnnularVertexMesh(
        unsigned numElements, double innerR, double outerR,
        double t1, double t2)
    {
        std::vector<Node<2>*> vNodes;
        double dtheta = 2.0 * M_PI / static_cast<double>(numElements);

        for (unsigned i = 0; i < numElements; i++)
        {
            double theta = i * dtheta;
            vNodes.push_back(new Node<2>(i, true,
                                         innerR * cos(theta), innerR * sin(theta)));
        }
        for (unsigned i = 0; i < numElements; i++)
        {
            double theta = i * dtheta;
            vNodes.push_back(new Node<2>(numElements + i, true,
                                         outerR * cos(theta), outerR * sin(theta)));
        }

        std::vector<VertexElement<2,2>*> elements;
        for (unsigned i = 0; i < numElements; i++)
        {
            unsigned ni = (i + 1) % numElements;
            std::vector<Node<2>*> en;
            en.push_back(vNodes[i]);
            en.push_back(vNodes[numElements + i]);
            en.push_back(vNodes[numElements + ni]);
            en.push_back(vNodes[ni]);
            elements.push_back(new VertexElement<2,2>(i, en));
        }

        boost::shared_ptr<MutableVertexMesh<2,2>> pMesh(
            new MutableVertexMesh<2,2>(vNodes, elements, t1, t2));
        pMesh->SetCellRearrangementRatio(1.5);
        pMesh->SetProtorosetteFormationProbability(0.0);
        pMesh->SetProtorosetteResolutionProbabilityPerTimestep(1.0);
        pMesh->SetCheckForInternalIntersections(true);
        return pMesh;
    }

    void RunVertex2d(const CryptBuddingParams& p, const std::string& outputDir)
    {
        RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

        // --- Build annular mesh ---
        boost::shared_ptr<MutableVertexMesh<2,2>> pMesh =
            MakeAnnularVertexMesh(p.numCells2dVertex, p.innerRadius2d, p.outerRadius2d,
                                  p.t1Threshold2d, p.t2Threshold2d);

        // Add noise
        for (unsigned i = 0; i < pMesh->GetNumNodes(); i++)
        {
            c_vector<double, 2> pos = pMesh->GetNode(i)->rGetLocation();
            double r = norm_2(pos);
            if (r > 1e-6)
            {
                double noise = (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.2;
                pos += noise * (pos / r);
                ChastePoint<2> pt(pos[0], pos[1]);
                pMesh->GetNode(i)->SetPoint(pt);
            }
        }

        // Target area
        double dtheta = 2.0 * M_PI / static_cast<double>(p.numCells2dVertex);
        double target_area = 0.5 * dtheta
            * (p.outerRadius2d * p.outerRadius2d - p.innerRadius2d * p.innerRadius2d);

        // --- Create cells ---
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem);
        MAKE_PTR(TransitCellProliferativeType, p_ta);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i = 0; i < pMesh->GetNumElements(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(target_area);

            CellPtr p_cell(new Cell(p_state, p_cycle));
            c_vector<double, 2> centroid = pMesh->GetCentroidOfElement(i);
            double frac = AngularFraction2d(centroid[0], centroid[1]);
            AssignCellTypeByFraction(p_cell, frac, p_stem, p_ta, p_diff);

            p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", target_area);
            p_cell->GetCellData()->SetItem("volume", target_area);
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessVertex);
            p_cell->GetCellData()->SetItem("is_apical", 1.0);
            cells.push_back(p_cell);
        }

        // --- Population ---
        VertexBasedCellPopulation<2> population(*pMesh, cells);
        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CellAgesWriter>();
        population.AddCellWriter<CellVolumesWriter>();
        population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // --- Simulation ---
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory(outputDir);
        simulator.SetDt(p.dt);
        simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
        simulator.SetEndTime(p.endTime);

        // --- Forces ---
        // Nagai–Honda
        MAKE_PTR(NagaiHondaForce<2>, p_nh);
        p_nh->SetNagaiHondaDeformationEnergyParameter(p.ecmStiffness);
        p_nh->SetNagaiHondaMembraneSurfaceEnergyParameter(p.nhMembraneSurface);
        p_nh->SetNagaiHondaCellCellAdhesionEnergyParameter(p.nhCellCellAdhesion);
        p_nh->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(p.nhBoundaryAdhesion);
        simulator.AddForce(p_nh);

        // Basement membrane
        MAKE_PTR(BasementMembraneForce<2>, p_bm);
        p_bm->SetBasementMembraneParameter(p.bmStiffnessVertex);
        p_bm->SetBasementMembraneRadius(p.outerRadius2d + 2.0);
        c_vector<double, 2> center2d = zero_vector<double>(2);
        p_bm->SetOrganoidCenter(center2d);
        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.outerRadius2d * 4.0);
        simulator.AddForce(p_bm);

        // Lumen pressure
        if (p.enableLumenPressure)
        {
            MAKE_PTR(LumenPressureForce<2>, p_lumen);
            p_lumen->SetPressureStrength(p.lumenPressure);
            p_lumen->SetLumenEquilibriumRadius(p.outerRadius2d + 1.0);
            p_lumen->SetTrackCenter(true);
            simulator.AddForce(p_lumen);
        }

        // Apical constriction
        if (p.enableApicalConstriction)
        {
            MAKE_PTR(ApicalConstrictionForce<2>, p_ac);
            p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
            simulator.AddForce(p_ac);
        }

        // Sloughing
        if (p.enableSloughing)
        {
            AddBoundingBoxKillers<2>(simulator, population,
                                    p.outerRadius2d * p.sloughRadiusFactor);
        }

        // Modifiers
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_area);
        p_area->SetReferenceTargetArea(target_area);
        simulator.AddSimulationModifier(p_area);

        MAKE_PTR(VolumeTrackingModifier<2>, p_vol);
        simulator.AddSimulationModifier(p_vol);

        boost::shared_ptr<CryptBuddingSummaryModifier<2>> p_summary(
            new CryptBuddingSummaryModifier<2>(p.ecmStiffness, p.samplingMultiple, p.endTime));
        simulator.AddSimulationModifier(p_summary);

        // --- Run with optional relaxation ---
        if (p.enableRelaxation)
        {
            std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
            for (AbstractCellPopulation<2>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                origTypes[*it] = it->GetCellProliferativeType();
                it->SetCellProliferativeType(p_diff);
            }

            simulator.SetEndTime(p.relaxationTime);
            std::cerr << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
            simulator.Solve();

            for (AbstractCellPopulation<2>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
            }

            simulator.SetEndTime(p.relaxationTime + p.endTime);
            std::cerr << "--- Phase 2: Growth (" << p.endTime << "h) ---" << std::endl;
            simulator.Solve();
        }
        else
        {
            simulator.Solve();
        }

        unsigned final_cells = population.GetNumRealCells();
        std::cerr << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
        TS_ASSERT_LESS_THAN(0u, final_cells);
    }


    // ===============================================================
    // MODEL 3: 3D Node-Based
    // ===============================================================

    void RunNode3d(const CryptBuddingParams& p, const std::string& outputDir)
    {
        RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

        // --- Create spherical shell (Fibonacci lattice) ---
        std::vector<Node<3>*> nodes;
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double golden = (1.0 + sqrt(5.0)) / 2.0;

        for (unsigned i = 0; i < p.numCells3dNode; i++)
        {
            double theta = 2.0 * M_PI * i / golden;
            double phi = acos(1.0 - 2.0 * (i + 0.5) / p.numCells3dNode);
            double r = p.organoidRadius3d + (p_gen->ranf() - 0.5) * p.shellThickness3d;
            double x = r * sin(phi) * cos(theta);
            double y = r * sin(phi) * sin(theta);
            double z = r * cos(phi);
            nodes.push_back(new Node<3>(i, false, x, y, z));
        }

        NodesOnlyMesh<3> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, p.interactionCutoff3d);
        for (unsigned i = 0; i < nodes.size(); i++) delete nodes[i];

        // --- Create cells ---
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem);
        MAKE_PTR(TransitCellProliferativeType, p_ta);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

        for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(3);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_cycle));
            c_vector<double, 3> pos = mesh.GetNode(i)->rGetLocation();
            double frac = ZFractionToRingFraction(pos[2], p.organoidRadius3d);
            AssignCellTypeByFraction(p_cell, frac, p_stem, p_ta, p_diff);

            p_cell->SetBirthTime(-p_gen->ranf() * 18.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("volume", 1.0);
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessNode);
            p_cell->GetCellData()->SetItem("is_apical", 1.0);
            cells.push_back(p_cell);
        }

        // --- Population ---
        NodeBasedCellPopulation<3> population(mesh, cells);
        population.SetAbsoluteMovementThreshold(50.0);
        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CellAgesWriter>();
        population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

        // --- Simulation ---
        OffLatticeSimulation<3> simulator(population);
        simulator.SetOutputDirectory(outputDir);
        simulator.SetDt(p.dt);
        simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
        simulator.SetEndTime(p.endTime);

        // --- Forces ---
        // Cell–cell
        if (p.enableDifferentialAdhesion)
        {
            MAKE_PTR(DifferentialAdhesionForce<3>, p_spring);
            p_spring->SetMeinekeSpringStiffness(p.springStiffness);
            p_spring->SetCutOffLength(p.interactionCutoff3d);
            p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
            p_spring->SetMeinekeSpringGrowthDuration(1.0);
            p_spring->SetApicalApicalAdhesion(p.apicalApicalAdhesion);
            p_spring->SetBasalBasalAdhesion(p.basalBasalAdhesion);
            p_spring->SetApicalBasalAdhesion(p.apicalBasalAdhesion);
            simulator.AddForce(p_spring);
        }
        else
        {
            MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring);
            p_spring->SetMeinekeSpringStiffness(p.springStiffness);
            p_spring->SetCutOffLength(p.interactionCutoff3d);
            simulator.AddForce(p_spring);
        }

        // Basement membrane
        MAKE_PTR(BasementMembraneForce<3>, p_bm);
        p_bm->SetBasementMembraneParameter(p.bmStiffnessNode);
        p_bm->SetTargetRadius(p.bmRadius3d);
        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);
        simulator.AddForce(p_bm);

        // Lumen pressure
        if (p.enableLumenPressure)
        {
            MAKE_PTR(LumenPressureForce<3>, p_lumen);
            p_lumen->SetPressureStrength(p.lumenPressure);
            p_lumen->SetLumenEquilibriumRadius(p.organoidRadius3d + 1.0);
            p_lumen->SetTrackCenter(true);
            simulator.AddForce(p_lumen);
        }

        // Apical constriction
        if (p.enableApicalConstriction)
        {
            MAKE_PTR(ApicalConstrictionForce<3>, p_ac);
            p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
            simulator.AddForce(p_ac);
        }

        // ECM guidance (3D only)
        boost::shared_ptr<DynamicECMField3d> pEcmField;
        if (p.enableEcmGuidance)
        {
            pEcmField.reset(new DynamicECMField3d(
                "radial", p.ecmGridSpacing,
                -p.ecmDomainHalf, p.ecmDomainHalf,
                -p.ecmDomainHalf, p.ecmDomainHalf,
                -p.ecmDomainHalf, p.ecmDomainHalf));
            pEcmField->SetDegradationRate(0.002);
            pEcmField->SetRemodelingRate(0.05);
            pEcmField->SetDepositionRate(0.0003);

            MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm);
            p_ecm->SetECMField(pEcmField);
            p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
            p_ecm->SetECMSensitivity(1.0);
            p_ecm->SetEnableDegradation(true);
            p_ecm->SetEnableRemodeling(true);
            p_ecm->SetEnableDeposition(false);
            simulator.AddForce(p_ecm);
        }

        // Sloughing
        if (p.enableSloughing)
        {
            AddBoundingBoxKillers<3>(simulator, population,
                                    p.organoidRadius3d * p.sloughRadiusFactor);
        }

        // Modifiers
        MAKE_PTR(VolumeTrackingModifier<3>, p_vol);
        simulator.AddSimulationModifier(p_vol);

        boost::shared_ptr<CryptBuddingSummaryModifier<3>> p_summary(
            new CryptBuddingSummaryModifier<3>(p.ecmStiffness, p.samplingMultiple,
                                               p.relaxationTime + p.endTime));
        simulator.AddSimulationModifier(p_summary);

        // --- Run with optional relaxation ---
        if (p.enableRelaxation)
        {
            std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
            for (AbstractCellPopulation<3>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                origTypes[*it] = it->GetCellProliferativeType();
                it->SetCellProliferativeType(p_diff);
            }

            simulator.SetEndTime(p.relaxationTime);
            std::cerr << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
            simulator.Solve();

            for (AbstractCellPopulation<3>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
            }

            simulator.SetEndTime(p.relaxationTime + p.endTime);
            std::cerr << "--- Phase 2: Growth (" << p.endTime << "h) ---" << std::endl;
            simulator.Solve();
        }
        else
        {
            simulator.Solve();
        }

        unsigned final_cells = population.GetNumRealCells();
        std::cerr << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
        TS_ASSERT_LESS_THAN(0u, final_cells);
    }


    // ===============================================================
    // MODEL 4: 3D Vertex-Based (OrganoidChaste)
    //
    // Uses SurfaceTensionForce (not NagaiHondaForce) because the 3D
    // monolayer mesh has distinct apical/basal/lateral faces that need
    // face-specific tensions, plus volume elasticity. This is the
    // OrganoidChaste standard approach (Drozdowski & Schwarz 2025).
    //
    // Also uses:
    //   - FiniteThicknessSimulation3d (extends OffLatticeSimulation<3>)
    //     to handle target-volume updates after remeshing
    //   - GeometricalTargetVolumeModifier for volume homeostasis
    //   - LumenPressureSubForce (OrganoidChaste) for lumen pressure
    //     (extends AbstractForce, acts on apical faces of monolayer)
    //   - ContactInhibitionCellCycleModel for consistency with other
    //     models (replacing the original FixedG1GenerationalCellCycleModel)
    // ===============================================================

    void RunVertex3d(const CryptBuddingParams& p, const std::string& outputDir)
    {
        RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

        // --- Compute mesh parameters (OrganoidChaste convention) ---
        double height = 2.0 / 3.0 / sqrt(3.0)
                        * cbrt((9.0 / 2.0) * (9.0 / 2.0))
                        * cbrt((p.gammaApical + p.gammaBasal) / p.gammaLateral
                               * (p.gammaApical + p.gammaBasal) / p.gammaLateral)
                        * 1.0;
        double t1_length = 0.66 / cbrt(3.0 * 3.0
                           * (1.0 + p.gammaLateral) / p.gammaLateral);

        std::cerr << "  Cell height:    " << height << std::endl;
        std::cerr << "  T1 threshold:   " << t1_length << std::endl;

        // --- Generate sphere mesh ---
        FiniteThicknessRandomizedSphereMeshGenerator generator(
            p.numCells3dVertex, t1_length, 0.001, height, p.sphereRadius3dVertex);
        MutableMonolayerVertexMesh<3, 3>* pMesh = generator.GetMesh();
        pMesh->SetProtorosetteFormationProbability(1.0);
        pMesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);
        TS_ASSERT_EQUALS(pMesh->GetNumElements(), p.numCells3dVertex);

        // Compute average cell volume for equilibrium volume
        double outerR = p.sphereRadius3dVertex + height;
        double avgVol = 4.0 / 3.0 * M_PI
                        * (outerR * outerR * outerR
                           - p.sphereRadius3dVertex * p.sphereRadius3dVertex * p.sphereRadius3dVertex)
                        / p.numCells3dVertex;

        // --- Create cells with ContactInhibitionCellCycleModel ---
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem);
        MAKE_PTR(TransitCellProliferativeType, p_ta);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

        for (unsigned i = 0; i < p.numCells3dVertex; i++)
        {
            ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
            p_cycle->SetDimension(3);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(avgVol);

            CellPtr p_cell(new Cell(p_state, p_cycle));

            c_vector<double, 3> centroid = pMesh->GetCentroidOfElement(i);
            double z_frac = centroid[2] / (p.sphereRadius3dVertex + height / 2.0);

            if (z_frac < -0.5)
            {
                p_cell->SetCellProliferativeType(p_stem);
                p_cell->GetCellData()->SetItem("cell_type_id", 0.0);
            }
            else if (z_frac < 0.3)
            {
                p_cell->SetCellProliferativeType(p_ta);
                p_cell->GetCellData()->SetItem("cell_type_id", 1.0);
            }
            else
            {
                p_cell->SetCellProliferativeType(p_diff);
                p_cell->GetCellData()->SetItem("cell_type_id", 2.0);
            }

            p_cell->SetBirthTime(-RandomNumberGenerator::Instance()->ranf() * 10.0);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("volume", avgVol);
            p_cell->GetCellData()->SetItem("is_apical", 1.0);
            cells.push_back(p_cell);
        }

        // --- Population (OrganoidChaste monolayer) ---
        MonolayerVertexBasedCellPopulation<3> population(*pMesh, cells);
        population.SetOutputCellRearrangementLocations(false);
        population.SetRestrictVertexMovementBoolean(false);
        population.SetDoInitialVolumeRelaxation(true);

        population.AddCellWriter<CellVolumesWriter>();
        population.AddCellWriter<CellThicknessWriter>();
        population.AddCellWriter<CellProliferativeTypesWriter>();
        population.AddCellWriter<CellIdWriter>();
        population.AddCellWriter<CellAgesWriter>();
        population.AddFaceWriter<FaceTypeWriter>();

        // --- Simulation (FiniteThicknessSimulation3d) ---
        FiniteThicknessSimulation3d simulator(population);
        simulator.SetOutputDirectory(outputDir);
        simulator.SetDt(p.dt);   // relaxation dt
        simulator.SetSamplingTimestepMultiple(p.samplingMultiple);

        // --- Forces ---
        // Surface tension (OrganoidChaste vertex force)
        MAKE_PTR(SurfaceTensionForce<3>, p_tension);
        p_tension->CreateSurfaceTensionParametersForCells(
            p.gammaApical, p.gammaBasal, p.gammaLateral, pMesh);
        p_tension->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
        p_tension->SetSimulationInstance(&simulator);

        simulator.AddForce(p_tension);

        // NOTE: LumenPressureSubForce is added in Phase 2 only.
        // Adding it during relaxation destabilises the mesh → DIVERGED_ITS.

        // Basement membrane (TissueMorphology)
        MAKE_PTR(BasementMembraneForce<3>, p_bm);
        p_bm->SetBasementMembraneParameter(p.bmStiffnessVertex);
        p_bm->SetTargetRadius(p.sphereRadius3dVertex + height + 1.0);
        simulator.AddForce(p_bm);

        // ECM guidance (3D only, optional)
        boost::shared_ptr<DynamicECMField3d> pEcmField;
        if (p.enableEcmGuidance)
        {
            pEcmField.reset(new DynamicECMField3d(
                "radial", p.ecmGridSpacing,
                -p.ecmDomainHalf, p.ecmDomainHalf,
                -p.ecmDomainHalf, p.ecmDomainHalf,
                -p.ecmDomainHalf, p.ecmDomainHalf));
            pEcmField->SetDegradationRate(0.002);
            pEcmField->SetRemodelingRate(0.05);
            pEcmField->SetDepositionRate(0.0003);
        }

        // --- Modifiers ---
        // Volume tracking (needed for ContactInhibitionCellCycleModel)
        MAKE_PTR(VolumeTrackingModifier<3>, p_vol);
        simulator.AddSimulationModifier(p_vol);

        // GeometricalTargetVolumeModifier (OrganoidChaste: volume homeostasis)
        MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_gvol, (&population));
        p_gvol->SetGrowthDuration(0.0);
        p_gvol->SetT1AdaptationDuration(0.100);
        p_gvol->SetReferenceTargetVolume(avgVol);
        simulator.AddSimulationModifier(p_gvol);

        // Note: no sloughing killers for vertex3d because
        // MutableMonolayerVertexMesh does not support DeleteElementPriorToReMesh in 3D.

        // ================================================================
        // PHASE 1: RELAXATION (no growth)
        // ================================================================
        if (p.enableRelaxation)
        {
            // Temporarily suppress proliferation
            std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
            for (AbstractCellPopulation<3>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                origTypes[*it] = it->GetCellProliferativeType();
                it->SetCellProliferativeType(p_diff);
            }

            simulator.SetEndTime(p.relaxationTime);
            std::cerr << "--- Phase 1: Relaxation (" << p.relaxationTime << " hours) ---" << std::endl;
            simulator.Solve();
            std::cerr << "Relaxation complete. Cells: " << population.GetNumRealCells() << std::endl;

            // Restore types
            for (AbstractCellPopulation<3>::Iterator it = population.Begin();
                 it != population.End(); ++it)
            {
                if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
            }

            // ================================================================
            // PHASE 2: GROWTH + BUCKLING
            // ================================================================
            double dt_grow = 0.006;
            unsigned sampling_grow = 100;

            simulator.SetEndTime(p.relaxationTime + p.endTime);
            simulator.SetSamplingTimestepMultiple(sampling_grow);
            simulator.SetDt(dt_grow);

            // Enable active T1 transitions (allows buckling)
            p_tension->SetSimulatedAnnealingParameters(0.003, 1900000.0, 1.0);
            p_tension->SetSimulationInstance(&simulator);
            p_tension->SetPerformActiveT1Swaps(true);
            p_tension->SetT1TransitionParameters(2.0, false);

            // Enable BM degradation
            p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);

            // Add lumen pressure in Phase 2 (not during relaxation)
            if (p.enableLumenPressure)
            {
                MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_sub, (p.lumenPressure));
                simulator.AddForce(p_lumen_sub);
            }

            // Add ECM guidance if enabled
            if (p.enableEcmGuidance && pEcmField)
            {
                MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm);
                p_ecm->SetECMField(pEcmField);
                p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
                p_ecm->SetECMSensitivity(1.0);
                p_ecm->SetEnableDegradation(true);
                p_ecm->SetEnableRemodeling(true);
                p_ecm->SetEnableDeposition(false);
                simulator.AddForce(p_ecm);
            }

            std::cerr << "--- Phase 2: Growth (" << p.endTime << " hours) ---" << std::endl;
            simulator.Solve();
        }
        else
        {
            simulator.SetEndTime(p.endTime);
            simulator.Solve();
        }

        unsigned final_cells = population.GetNumRealCells();
        std::cerr << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
        TS_ASSERT_LESS_THAN(0u, final_cells);
    }


public:

    /**
     * Single entry point — reads MODEL_TYPE and dispatches.
     *
     * Usage:
     *   MODEL_TYPE=node2d ECM_STIFFNESS=5.0 RUN_NUMBER=0 \
     *     ctest -R TestCryptBudding --output-on-failure
     *
     *   MODEL_TYPE=vertex3d ECM_STIFFNESS=2.0 ENABLE_ECM_GUIDANCE=1 \
     *     ctest -R TestCryptBudding --output-on-failure
     */
    void TestCryptBuddingSweep()
    {
        CryptBuddingParams params;
        params.InitFromEnvironment();

        // Build output directory
        std::stringstream subdir;
        subdir << "CryptBudding/" << params.modelType
               << "/stiffness_" << std::fixed << std::setprecision(1) << params.ecmStiffness
               << "/run_" << params.runNumber;

        PrintBanner(params);

        if (params.modelType == "node2d")
        {
            RunNode2d(params, subdir.str());
        }
        else if (params.modelType == "vertex2d")
        {
            RunVertex2d(params, subdir.str());
        }
        else if (params.modelType == "node3d")
        {
            RunNode3d(params, subdir.str());
        }
        else if (params.modelType == "vertex3d")
        {
            RunVertex3d(params, subdir.str());
        }
        else
        {
            EXCEPTION("Unknown MODEL_TYPE: " + params.modelType
                      + ". Options: node2d, vertex2d, node3d, vertex3d");
        }
    }
};

#endif /* TESTCRYPTBUDDING_HPP_ */
