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
#include <string>
#include <iostream>
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
    std::cout << "\n============================================" << std::endl;
    std::cout << "  Crypt Budding Simulation" << std::endl;
    std::cout << "  Git Commit:     " << TM_GIT_HASH << std::endl;
    std::cout << "  Model:          " << p.modelType << std::endl;
    std::cout << "  ECM Stiffness:  " << p.ecmStiffness << std::endl;
    std::cout << "  Run Number:     " << p.runNumber << std::endl;
    std::cout << "  Seed:           " << p.randomSeed << std::endl;
    std::cout << "  dt:             " << p.dt << std::endl;
    if (p.modelType == "vertex3d")
        std::cout << "  dt (growth):    " << p.dtGrow << std::endl;
    std::cout << "  Relaxation:     " << (p.enableRelaxation ? "ON" : "OFF")
              << " (" << p.relaxationTime << "h)" << std::endl;
    std::cout << "  End Time:       " << p.endTime << "h" << std::endl;
    std::cout << "  Cell fractions: stem=" << p.stemFraction
              << " TA=" << p.transitFraction
              << " diff=" << (1.0 - p.stemFraction - p.transitFraction) << std::endl;
    if (p.modelType == "vertex3d")
    {
        std::cout << "  Gamma scaling:  stem=" << p.gammaStemScale
                  << " TA=" << p.gammaTransitScale
                  << " diff/Paneth=" << p.gammaDiffScale << std::endl;
    }
    std::cout << "  Features:" << std::endl;
    std::cout << "    Lumen Pressure:         " << (p.enableLumenPressure ? "ON" : "OFF") << std::endl;
    std::cout << "    Apical Constriction:     " << (p.enableApicalConstriction ? "ON" : "OFF") << std::endl;
    std::cout << "    ECM Confinement:         " << (p.enableEcmConfinement ? "ON" : "OFF") << std::endl;
    if (p.modelType == "node3d" || p.modelType == "vertex3d")
        std::cout << "    ECM Guidance (3D):       " << (p.enableEcmGuidance ? "ON" : "OFF") << std::endl;
    if (p.modelType == "node2d" || p.modelType == "vertex2d")
        std::cout << "    Sloughing:               " << (p.enableSloughing ? "ON" : "OFF") << std::endl;
    if (p.enableDifferentialAdhesion)
        std::cout << "    Differential Adhesion:   ON" << std::endl;
    if (p.modelType == "node3d")
    {
        std::cout << "    Curvature Bending:       " << (p.enableCurvatureBending ? "ON" : "OFF") << std::endl;
        std::cout << "    Cell Polarity:           " << (p.enableCellPolarity ? "ON" : "OFF") << std::endl;
        std::cout << "    Topology-Based Springs:  " << (p.useTopologyBasedSprings ? "ON" : "OFF") << std::endl;
    }
    if (p.modelType == "node2d")
        std::cout << "    Topology-Based Springs:  " << (p.useTopologyBasedSprings ? "ON" : "OFF") << std::endl;
    if (p.enableEcmConfinement)
        std::cout << "    ECM Grid Type:           " << p.ecmGridType << std::endl;
    if (p.enableContinuousPvd)
        std::cout << "    Continuous PVD:          ON" << std::endl;
    std::cout << "============================================\n" << std::endl;
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
              << "  -config <file.ini>    Load parameters from INI file\n"
              << "  -saveconfig <file>    Save current parameters to INI file and exit\n"
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
              << "\nVertex model parameters (Nagai-Honda):\n"
              << "  -adhesion <double>    cell-cell adhesion strength (default: 1.0)\n"
              << "  -membrane <double>    membrane surface tension (default: 10.0)\n"
              << "  -pressure <double>    lumen pressure strength (default: 1.0)\n"
              << "  -help                 print this message\n"
              << "\nExample:\n"
              << "  ./CryptBuddingApp -model mesh2d -config params.ini\n"
              << "  ./CryptBuddingApp -model mesh2d -saveconfig default.ini\n"
              << std::endl;
}

// ===================================================================
// CLI argument parsing
// ===================================================================

inline CryptBuddingParams ParseArguments(int argc, char* argv[])
{
    CryptBuddingParams p;
    p.SetDefaults();

    // First pass: look for -config to load base parameters
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-config" && i + 1 < argc)
        {
            std::string configFile = argv[++i];
            if (!p.LoadFromFile(configFile))
            {
                std::cerr << "Warning: Failed to load config file: " << configFile << std::endl;
            }
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
            p.SaveToFile(argv[++i]);
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
        else if (arg == "-dtgrow" && i + 1 < argc)
        {
            p.dtGrow = std::stod(argv[++i]);
        }
        else if (arg == "-continuous-pvd")
        {
            p.enableContinuousPvd = true;
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
