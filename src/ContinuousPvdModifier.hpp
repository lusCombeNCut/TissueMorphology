/*
 * ContinuousPvdModifier.hpp
 *
 * A simulation modifier that keeps all .pvd files in the output directory
 * continuously valid for ParaView by maintaining shadow copies with proper
 * XML closing tags.
 *
 * Without this, .pvd files are only properly closed when CloseWritersFiles()
 * is called at the end of a successful simulation. If the simulation crashes
 * or is interrupted (Ctrl-C, numerical blowup, etc.), the PVD files lack
 * their closing </Collection></VTKFile> tags and ParaView cannot open them.
 *
 * This modifier, enabled via the CLI flag -continuous-pvd, runs at each
 * output timestep and for each .pvd file creates a valid shadow copy
 * (e.g. results.pvd → results_valid.pvd) that contains all entries written
 * so far plus proper closing tags. ParaView should open the *_valid.pvd
 * files for live monitoring of running simulations.
 *
 * At the end of a successful simulation, the original .pvd files are
 * properly closed by Chaste's CloseWritersFiles(), so the shadow copies
 * become redundant.
 *
 * Usage:
 *   Add as a simulation modifier AFTER all writers have been registered.
 *   Enable via -continuous-pvd CLI flag.
 */
#ifndef CONTINUOUSPVDMODIFIER_HPP_
#define CONTINUOUSPVDMODIFIER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "SimulationTime.hpp"
#include "OutputFileHandler.hpp"
#include "SimProfiler.hpp"

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <dirent.h>

template<unsigned DIM>
class ContinuousPvdModifier : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:

    /** Absolute path to the results directory */
    std::string mOutputDirectory;

    /** How often to update shadow PVD files (same as simulation sampling multiple) */
    unsigned mSamplingTimestepMultiple;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<
            AbstractCellBasedSimulationModifier<DIM, DIM>>(*this);
        archive & mOutputDirectory;
        archive & mSamplingTimestepMultiple;
    }

    /**
     * Find all .pvd files in a directory (non-recursive).
     */
    std::vector<std::string> FindPvdFiles(const std::string& directory) const
    {
        std::vector<std::string> pvdFiles;
        DIR* dir = opendir(directory.c_str());
        if (dir == nullptr)
        {
            return pvdFiles;
        }

        struct dirent* entry;
        while ((entry = readdir(dir)) != nullptr)
        {
            std::string name(entry->d_name);
            // Match *.pvd but exclude *_valid.pvd (our own shadow files)
            if (name.size() > 4
                && name.substr(name.size() - 4) == ".pvd"
                && name.find("_valid.pvd") == std::string::npos)
            {
                pvdFiles.push_back(directory + name);
            }
        }
        closedir(dir);
        return pvdFiles;
    }

    /**
     * Create a shadow copy of a PVD file with proper XML closing tags.
     *
     * Given "results.pvd", creates "results_valid.pvd" containing all
     * the content of the original plus </Collection></VTKFile> closing
     * tags. The original file is not modified.
     *
     * The shadow file is atomically renamed from a .tmp so that ParaView
     * never sees a half-written file.
     */
    void CreateValidShadow(const std::string& pvdPath) const
    {
        // Read current content of the original PVD
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

        // Build shadow filename: results.pvd → results_valid.pvd
        std::string shadowPath = pvdPath.substr(0, pvdPath.size() - 4) + "_valid.pvd";
        std::string tmpPath = shadowPath + ".tmp";

        // Write shadow file
        std::ofstream outFile(tmpPath.c_str());
        if (!outFile.is_open())
        {
            return;
        }

        outFile << content;

        // Add closing tags if missing
        if (content.find("</VTKFile>") == std::string::npos)
        {
            if (content.find("</Collection>") == std::string::npos)
            {
                // Match indent style from original
                if (content.find("  <Collection>") != std::string::npos)
                {
                    outFile << "  </Collection>\n";
                }
                else
                {
                    outFile << "    </Collection>\n";
                }
            }
            outFile << "</VTKFile>\n";
        }

        outFile.close();

        // Atomic rename (POSIX): ParaView never sees a partial file
        std::rename(tmpPath.c_str(), shadowPath.c_str());
    }

public:

    ContinuousPvdModifier(unsigned samplingMultiple = 1)
        : AbstractCellBasedSimulationModifier<DIM, DIM>(),
          mSamplingTimestepMultiple(samplingMultiple)
    {
    }

    virtual ~ContinuousPvdModifier() {}

    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                    std::string outputDirectory) override
    {
        OutputFileHandler handler(outputDirectory, false);
        mOutputDirectory = handler.GetOutputDirectoryFullPath();
    }

    /**
     * At each output timestep, create valid shadow copies of all PVD files
     * in the results directory so ParaView can open them at any time.
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer _prof("ContinuousPvd");
        if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingTimestepMultiple != 0)
        {
            return;  // Not an output timestep
        }

        // Process PVD files in the main output directory
        std::vector<std::string> pvdFiles = FindPvdFiles(mOutputDirectory);
        for (const std::string& pvdFile : pvdFiles)
        {
            CreateValidShadow(pvdFile);
        }

        // Also check results_from_time_* subdirectories
        DIR* dir = opendir(mOutputDirectory.c_str());
        if (dir != nullptr)
        {
            struct dirent* entry;
            while ((entry = readdir(dir)) != nullptr)
            {
                std::string name(entry->d_name);
                if (name.find("results_from_time_") == 0)
                {
                    std::string subdir = mOutputDirectory + name + "/";
                    std::vector<std::string> subPvds = FindPvdFiles(subdir);
                    for (const std::string& pvdFile : subPvds)
                    {
                        CreateValidShadow(pvdFile);
                    }
                }
            }
            closedir(dir);
        }
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<ContinuousPvdModifier/>\n";
        AbstractCellBasedSimulationModifier<DIM, DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};

// Explicit instantiation
template class ContinuousPvdModifier<2>;
template class ContinuousPvdModifier<3>;

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ContinuousPvdModifier)

#endif // CONTINUOUSPVDMODIFIER_HPP_
