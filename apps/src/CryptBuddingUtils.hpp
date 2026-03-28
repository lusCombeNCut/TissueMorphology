/*
 * CryptBuddingUtils.hpp
 *
 * Utility functions shared across all crypt budding model runners:
 *   - AssignCellTypeByFraction  (uniform random cell type assignment)
 *   - AngularFraction2d / ZFractionToRingFraction
 *   - AddBoundingBoxKillers
 *   - MakeAnnularVertexMesh
 *   - PrintBanner / PrintUsage / ParseArguments
 */
#ifndef CRYPTBUDDINGUTILS_HPP_
#define CRYPTBUDDINGUTILS_HPP_

#include <cmath>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "PlaneBasedCellKiller.hpp"
#include "MutableVertexMesh.hpp"
#include "ChastePoint.hpp"
#include "CryptBuddingParams.hpp"

// ===================================================================
// Cell type assignment — uniform random
// ===================================================================

inline void AssignCellTypeByFraction(
    CellPtr pCell, double /*frac*/,
    boost::shared_ptr<AbstractCellProliferativeType> pStemType,
    boost::shared_ptr<AbstractCellProliferativeType> pTransitType,
    boost::shared_ptr<AbstractCellProliferativeType> pDiffType,
    double stemFraction, double transitFraction)
{
    // Uniformly random cell type assignment: each cell independently draws
    // from the target proportions, giving a spatially uniform distribution
    // across the organoid surface (rather than z-band segregation).
    // NOTE: stemFraction and transitFraction must always be provided from config file
    double u = RandomNumberGenerator::Instance()->ranf();
    if (u < stemFraction)
    {
        pCell->SetCellProliferativeType(pStemType);
        pCell->GetCellData()->SetItem("cell_type_id", 0.0);
    }
    else if (u < stemFraction + transitFraction)
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
// Geometric helpers
// ===================================================================

inline double AngularFraction2d(double x, double y,
                                double cx = 0.0, double cy = 0.0)
{
    double theta = atan2(y - cy, x - cx);
    double angle_from_bottom = theta + M_PI / 2.0;
    if (angle_from_bottom < 0.0) angle_from_bottom += 2.0 * M_PI;
    return angle_from_bottom / (2.0 * M_PI);
}

inline double ZFractionToRingFraction(double z, double R)
{
    double z_frac = z / R;
    if (z_frac < -0.5)      return 0.0;
    else if (z_frac < 0.3)  return 0.3;
    else                     return 0.5;
}

// ===================================================================
// Bounding-box cell killers (all DIM axes)
// ===================================================================

template<unsigned DIM>
void AddBoundingBoxKillers(OffLatticeSimulation<DIM>& rSimulator,
                           AbstractCellPopulation<DIM>& rPopulation,
                           double halfWidth)
{
    for (unsigned d = 0; d < DIM; d++)
    {
        {
            c_vector<double, DIM> point  = zero_vector<double>(DIM);
            c_vector<double, DIM> normal = zero_vector<double>(DIM);
            point[d]  =  halfWidth;
            normal[d] =  1.0;
            MAKE_PTR_ARGS(PlaneBasedCellKiller<DIM>, p_killer,
                          (&rPopulation, point, normal));
            rSimulator.AddCellKiller(p_killer);
        }
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
// Build annular vertex mesh for 2D vertex model
// ===================================================================

inline boost::shared_ptr<MutableVertexMesh<2,2>> MakeAnnularVertexMesh(
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
    // Note: SetCheckForInternalIntersections(true) causes PerformIntersectionSwap 
    // assertion failures - disabled to avoid crashes
    pMesh->SetCheckForInternalIntersections(false);
    // Note: T3 swaps have unimplemented edge cases that cause NEVER_REACHED crashes
    pMesh->SetCheckForT3Swaps(false);
    return pMesh;
}

// ===================================================================
// Banner & usage
// ===================================================================

inline void PrintBanner(const CryptBuddingParams& p)
{
#ifndef TM_GIT_HASH
#define TM_GIT_HASH "unknown"
#endif
    const bool is2d  = (p.modelType == "node2d"   || p.modelType == "vertex2d");
    const bool isNode = (p.modelType == "node2d"  || p.modelType == "node3d");
    const bool is3d  = (p.modelType == "node3d"   || p.modelType == "vertex3d");

    auto onoff = [](bool v) { return v ? "ON " : "OFF"; };

    std::cout << "\n============================================================\n"
              << "  Crypt Budding Simulation\n"
              << "  Git:      " << TM_GIT_HASH << "\n"
              << "  Model:    " << p.modelType
              << "   Run: " << p.runNumber
              << "   Seed: " << p.randomSeed << "\n"
              << "  Time:     end=" << p.endTime << "h"
              << "  dt=" << p.dt;
    if (p.enableRelaxation)
        std::cout << "  relax=" << p.relaxationTime << "h";
    std::cout << "\n";
    std::cout << "  Springs:  k=" << p.springStiffness
              << "  cutoff=" << p.springCutoff;
    // Warn if cutoff ≤ rest length (1.0 by default) — eliminates all attraction
    if (p.springCutoff <= 1.0)
        std::cout << "  *** WARNING: cutoff <= rest length, no attractive spring force! ***";
    std::cout << "\n";
    std::cout << "  ECM:      stiffness=" << p.ecmStiffness;
    if (p.enableGhostNodeECM)
        std::cout << " (ghost node)";
    std::cout << "\n";
    std::cout << "  Cells:    stem=" << p.stemFraction
              << "  TA=" << p.transitFraction
              << "  diff=" << (1.0 - p.stemFraction - p.transitFraction) << "\n";

    std::cout << "  Features:\n"
              << "    Lumen Pressure:        " << onoff(p.enableLumenPressure)  << "\n"
              << "    Apical Constriction:   " << onoff(p.enableApicalConstriction) << "\n"
              << "    ECM Confinement:       " << onoff(p.enableEcmConfinement) << "\n"
              << "    Ghost Node ECM:        " << onoff(p.enableGhostNodeECM)   << "\n"
              << "    Relaxation:            " << onoff(p.enableRelaxation)     << "\n"
              << "    Differential Adhesion: " << onoff(p.enableDifferentialAdhesion) << "\n";
    if (is2d)
        std::cout << "    Sloughing:             " << onoff(p.enableSloughing)  << "\n";
    if (isNode)
    {
        std::cout << "    Topology Springs:      " << onoff(p.useTopologyBasedSprings) << "\n"
                  << "    Curvature Bending:     " << onoff(p.enableCurvatureBending)  << "\n"
                  << "    Cell Polarity:         " << onoff(p.enableCellPolarity)      << "\n";
    }
    if (is3d)
        std::cout << "    ECM Guidance:          " << onoff(p.enableEcmGuidance) << "\n";
    if (p.enableContinuousPvd)
        std::cout << "    Continuous PVD:        ON\n";
    std::cout << "============================================================\n" << std::endl;
}

inline void PrintUsage()
{
    std::cout << "Usage: CryptBuddingApp -model <type> [options]\n"
              << "\nRequired:\n"
              << "  -model <type>  where type is one of:\n"
              << "      node2d   - Overlapping spheres (distance-based neighbors)\n"
              << "      mesh2d   - Delaunay mesh (topological neighbors) [recommended]\n"
              << "      vertex2d - Vertex model with polygonal cells\n"
              << "      node3d   - 3D overlapping spheres\n"
              << "      vertex3d - 3D vertex model\n"
              << "\nConfig file (loads parameters without recompiling):\n"
              << "  -config <file.json>   Load parameters from JSON file (recommended)\n"
              << "  -config <file.ini>    Load parameters from legacy INI file\n"
              << "  -saveconfig <file>    Save current parameters (.json or .ini) and exit\n"
              << "\nOptions (override config file values):\n"
              << "  -stiffness <double>   ECM stiffness (default: 5.0)\n"
              << "  -run <int>            run/replicate number (default: 0)\n"
              << "  -lumen <0|1>          lumen pressure (default: 1)\n"
              << "  -apical <0|1>         apical constriction (default: 1)\n"
              << "  -ecm <0|1>            ECM guidance, 3D only (default: 0)\n"
              << "  -relax <0|1>          relaxation phase (default: 1)\n"
              << "  -slough <0|1>         sloughing, 2D only (default: 1)\n"
              << "  -diffadh <0|1>        differential adhesion (default: 1)\n"
              << "  -bending <0|1>        curvature bending force (default: 1)\n"
              << "  -endtime <double>     override end time\n"
              << "  -dt <double>          override timestep (relaxation in vertex3d)\n"
              << "  -dtgrow <double>      growth phase timestep, vertex3d (default: 0.002)\n"
              << "  -continuous-pvd       keep .pvd files valid during simulation\n"
              << "  -v, -verbose          verbose output (verbosity level 2)\n"
              << "  -q, -quiet            suppress non-error output (verbosity level 0)\n"
              << "\nVertex model parameters (Nagai-Honda):\n"
              << "  -adhesion <double>    cell-cell adhesion strength (default: 1.0)\n"
              << "  -membrane <double>    membrane surface tension (default: 10.0)\n"
              << "  -pressure <double>    lumen pressure strength (default: 1.0)\n"
              << "  -help                 print this message\n"
              << "\nExample:\n"
              << "  ./CryptBuddingApp -model node2d -config params-Node2d.json\n"
              << "  ./CryptBuddingApp -model node2d -saveconfig default.json\n"
              << std::endl;
}

// ===================================================================
// CLI argument parsing
// ===================================================================

inline CryptBuddingParams ParseArguments(int argc, char* argv[])
{
    CryptBuddingParams p;
    p.SetDefaults();  // Initialises sentinels/internal flags only; physics values overridden below

    // ── Step 1: Auto-load system defaults from default_params.json ────────────────
    // This is the canonical source of default physics parameters (not the C++ SetDefaults()).
    // Search order:
    //   a) $CRYPTBUDDING_DEFAULTS environment variable
    //   b) <dir-of -config arg>/config/default_params.json
    {
        // Peek at the -config argument so we can derive the defaults path from it
        std::string configArg;
        for (int i = 1; i + 1 < argc; ++i)
            if (std::string(argv[i]) == "-config") { configArg = argv[i + 1]; break; }

        std::string defaultsPath;

        // a) Env var
        if (const char* env = std::getenv("CRYPTBUDDING_DEFAULTS"))
            if (env[0]) defaultsPath = env;

        // b) Relative to the -config file: <config_dir>/config/default_params.json
        if (defaultsPath.empty() && !configArg.empty())
        {
            size_t sep = configArg.find_last_of("/\\");
            std::string dir = (sep != std::string::npos) ? configArg.substr(0, sep) : ".";
            std::string candidate = dir + "/config/default_params.json";
            if (std::ifstream(candidate).good()) defaultsPath = candidate;
        }

        if (!defaultsPath.empty())
        {
            p.LoadFromJson(defaultsPath);
            // Reset override-tracking flags: the defaults JSON must not lock out
            // Finalise() auto-tuning (e.g. dt, endTime) for the user's config.
            p.endTimeOverridden    = false;
            p.dtOverridden         = false;
            p.t1ThresholdOverridden = false;
            p.t2ThresholdOverridden = false;
        }
        else
        {
            std::cerr << "Note: default_params.json not found — using compiled-in defaults.\n"
                      << "      Set $CRYPTBUDDING_DEFAULTS or ensure -config is in the apps/ directory.\n";
        }
    }

    // ── Step 2: Load user config file (overrides defaults) ───────────────────────
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-config" && i + 1 < argc)
        {
            std::string configFile = argv[++i];
            bool ok = false;
            // Dispatch on file extension: .json → JSON loader, anything else → INI loader
            if (configFile.size() > 5 &&
                configFile.substr(configFile.size() - 5) == ".json")
                ok = p.LoadFromJson(configFile);
            else
                ok = p.LoadFromFile(configFile);
            if (!ok)
                std::cerr << "Warning: Failed to load config file: " << configFile << std::endl;
            break;  // Only process first -config
        }
    }

    // Second pass: process all arguments (CLI overrides config file)
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-help" || arg == "--help")
        {
            PrintUsage();
            exit(0);
        }
        else if (arg == "-config" && i + 1 < argc)
        {
            ++i;  // Skip, already processed
        }
        else if (arg == "-saveconfig" && i + 1 < argc)
        {
            // Need model type first for proper defaults
            for (int j = 1; j < argc; j++)
            {
                if (std::string(argv[j]) == "-model" && j + 1 < argc)
                {
                    p.modelType = argv[j + 1];
                    break;
                }
            }
            p.Finalise();
            {
                std::string outFile = argv[++i];
                if (outFile.size() > 5 &&
                    outFile.substr(outFile.size() - 5) == ".json")
                    p.SaveToJson(outFile);
                else
                    p.SaveToFile(outFile);
            }
            std::cout << "Config saved. Exiting." << std::endl;
            exit(0);
        }
        else if (arg == "-model" && i + 1 < argc)
        {
            p.modelType = argv[++i];
        }
        else if (arg == "-stiffness" && i + 1 < argc)
        {
            p.ecmStiffness = std::stod(argv[++i]);
        }
        else if (arg == "-run" && i + 1 < argc)
        {
            p.runNumber = static_cast<unsigned>(std::stoi(argv[++i]));
        }
        else if (arg == "-lumen" && i + 1 < argc)
        {
            p.enableLumenPressure = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-apical" && i + 1 < argc)
        {
            p.enableApicalConstriction = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-ecm" && i + 1 < argc)
        {
            p.enableEcmGuidance = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-relax" && i + 1 < argc)
        {
            p.enableRelaxation = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-slough" && i + 1 < argc)
        {
            p.enableSloughing = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-diffadh" && i + 1 < argc)
        {
            p.enableDifferentialAdhesion = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-bending" && i + 1 < argc)
        {
            p.enableCurvatureBending = (std::stoi(argv[++i]) != 0);
        }
        else if (arg == "-endtime" && i + 1 < argc)
        {
            p.endTime = std::stod(argv[++i]);
            p.endTimeOverridden = true;
        }
        else if (arg == "-dt" && i + 1 < argc)
        {
            p.dt = std::stod(argv[++i]);
            p.dtOverridden = true;
        }
        else if (arg == "-continuous-pvd")
        {
            p.enableContinuousPvd = true;
        }
        else if (arg == "-v" || arg == "-verbose")
        {
            p.verbosity = 2;
        }
        else if (arg == "-q" || arg == "-quiet")
        {
            p.verbosity = 0;
        }
        else if (arg == "-adhesion" && i + 1 < argc)
        {
            p.nhCellCellAdhesion = std::stod(argv[++i]);
        }
        else if (arg == "-membrane" && i + 1 < argc)
        {
            p.nhMembraneSurface = std::stod(argv[++i]);
        }
        else if (arg == "-pressure" && i + 1 < argc)
        {
            p.lumenPressure = std::stod(argv[++i]);
        }
        else
        {
            std::cerr << "WARNING: Unknown argument '" << arg << "'" << std::endl;
        }
    }

    return p;
}

#endif // CRYPTBUDDINGUTILS_HPP_
