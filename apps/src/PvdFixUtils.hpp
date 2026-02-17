/*
 * PvdFixUtils.hpp — Utilities for fixing incomplete PVD files
 *
 * When a Chaste simulation terminates abnormally (e.g. numerical
 * instability / DIVERGED_ITS), CloseWritersFiles() is never called
 * and the results.pvd file is left without its XML closing tags:
 *
 *     </Collection>
 *   </VTKFile>
 *
 * The functions here detect and repair such files so that ParaView
 * can still open partial results.
 */
#ifndef PVDFIXUTILS_HPP_
#define PVDFIXUTILS_HPP_

#include "OutputFileHandler.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <dirent.h>
#include <sys/stat.h>

/**
 * Fix a single PVD file that is missing its XML closing tags.
 * If the file already contains </VTKFile> it is left untouched.
 */
inline void FixPvdFileEnding(const std::string& pvdPath)
{
    std::ifstream inFile(pvdPath.c_str());
    if (!inFile.is_open())
    {
        return;
    }

    std::string content;
    std::string line;
    while (std::getline(inFile, line))
    {
        content += line + "\n";
    }
    inFile.close();

    if (content.empty())
    {
        return;
    }

    // Already properly closed
    if (content.find("</VTKFile>") != std::string::npos)
    {
        return;
    }

    std::cout << "  Fixing incomplete PVD file: " << pvdPath << std::endl;

    std::ofstream outFile(pvdPath.c_str(), std::ios::app);
    if (!outFile.is_open())
    {
        std::cerr << "  WARNING: Could not open PVD file for fixing: " << pvdPath << std::endl;
        return;
    }

    if (content.find("</Collection>") == std::string::npos)
    {
        outFile << "    </Collection>\n";
    }
    outFile << "</VTKFile>\n";
    outFile.close();
}

/**
 * Scan an output directory for all results.pvd files and fix any that
 * are missing proper XML closing tags.  Searches the directory itself
 * and all results_from_time_* subdirectories.
 *
 * @param outputSubdir  The relative output subdirectory (as passed to
 *                      OffLatticeSimulation::SetOutputDirectory).
 */
inline void FixAllPvdFiles(const std::string& outputSubdir)
{
    std::string baseDir = OutputFileHandler::GetChasteTestOutputDirectory() + outputSubdir;

    // Check for PVD directly in the output dir
    std::string pvd = baseDir + "/results.pvd";
    struct stat st;
    if (stat(pvd.c_str(), &st) == 0)
    {
        FixPvdFileEnding(pvd);
    }

    // Scan for results_from_time_* subdirectories
    DIR* dir = opendir(baseDir.c_str());
    if (dir == nullptr)
    {
        return;
    }

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr)
    {
        std::string name(entry->d_name);
        if (name.find("results_from_time_") == 0)
        {
            std::string subPvd = baseDir + "/" + name + "/results.pvd";
            if (stat(subPvd.c_str(), &st) == 0)
            {
                FixPvdFileEnding(subPvd);
            }
        }
    }
    closedir(dir);
}

#endif // PVDFIXUTILS_HPP_
