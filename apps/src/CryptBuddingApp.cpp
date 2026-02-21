/*
 * CryptBuddingApp.cpp
 *
 * Standalone Chaste app for crypt budding / organoid simulations.
 * Equivalent to TestCryptBudding.hpp but as a proper executable with
 * CLI arguments instead of environment variables.
 *
 * To rebuild: make -j8 -C projects/TissueMorphology/apps CryptBuddingApp
 * 
 * Usage:
 *    cd /<Project_Dir>/TissueMorphology/apps
 *   ./CryptBuddingApp -model node2d -stiffness 5.0 -run 0
 *   ./CryptBuddingApp -model vertex3d -stiffness 2.0 -run 3 -lumen 0 -ecm 1
 *
 * Arguments:
 *   -model <node2d|vertex2d|node3d|vertex3d>   (required)
 *   -stiffness <double>                         ECM stiffness  (default: 5.0)
 *   -run <int>                                  run/replicate  (default: 0)
 *   -lumen <0|1>                                lumen pressure (default: 1)
 *   -apical <0|1>                               apical constriction (default: 1)
 *   -ecm <0|1>                                  ECM guidance 3D (default: 0)
 *   -relax <0|1>                                relaxation phase (default: 1)
 *   -slough <0|1>                               sloughing 2D (default: 1)
 *   -diffadh <0|1>                              differential adhesion (default: 1)
 *   -endtime <double>                           override end time
 *   -dt <double>                                override timestep
 *   -dtgrow <double>                            growth phase timestep vertex3d
 *   -help                                       print usage
 */

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "PetscTools.hpp"
#include "SimulationTime.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingUtils.hpp"
#include "PvdFixUtils.hpp"
#include "RunNode2d.hpp"
#include "RunVertex2d.hpp"
#include "RunNode3d.hpp"
#include "RunVertex3d.hpp"

#include <sstream>
#include <iomanip>
#include <iostream>

int main(int argc, char* argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;
    std::string outputSubdir;  // Stored here so catch block can access it

    try
    {
        if (argc < 2)
        {
            PrintUsage();
            ExecutableSupport::PrintError("No arguments provided. Use -model <type>.", true);
            exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        }
        else
        {
            SimulationTime::Instance()->SetStartTime(0.0);

            CryptBuddingParams params = ParseArguments(argc, argv);

            if (params.modelType.empty())
            {
                PrintUsage();
                EXCEPTION("Required argument -model not provided");
            }

            params.Finalise();
            PrintBanner(params);

            // Build output directory (includes git commit hash for traceability)
#ifndef TM_GIT_HASH
#define TM_GIT_HASH "unknown"
#endif
            std::stringstream subdir;
            subdir << "CryptBudding/" << TM_GIT_HASH
                   << "/" << params.modelType
                   << "/stiffness_" << std::fixed << std::setprecision(1) << params.ecmStiffness
                   << "/run_" << params.runNumber;
            outputSubdir = subdir.str();

            if (params.modelType == "node2d")
            {
                RunNode2d(params, outputSubdir);
            }
            else if (params.modelType == "vertex2d")
            {
                RunVertex2d(params, outputSubdir);
            }
            else if (params.modelType == "node3d")
            {
                RunNode3d(params, outputSubdir);
            }
            else if (params.modelType == "vertex3d")
            {
                RunVertex3d(params, outputSubdir);
            }
            else
            {
                EXCEPTION("Unknown -model: " + params.modelType
                          + ". Options: node2d, vertex2d, node3d, vertex3d");
            }
        }
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;

        // Failsafe: close any incomplete PVD files so that ParaView can
        if (!outputSubdir.empty())
        {
            std::cout << "\n  === PVD Failsafe: checking for incomplete PVD files ===" << std::endl;
            FixAllPvdFiles(outputSubdir);
        }
    }

    ExecutableSupport::WriteMachineInfoFile("machine_info");
    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
