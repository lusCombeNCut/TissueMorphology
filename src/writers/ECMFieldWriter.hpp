/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

*/

#ifndef ECMFIELDWRITER_HPP_
#define ECMFIELDWRITER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "DynamicECMField.hpp"
#include "SimProfiler.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include <boost/shared_ptr.hpp>

/**
 * A modifier that writes the ECM field grid to VTK files at each output timestep.
 */
template<unsigned DIM>
class ECMFieldWriter : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    /** Pointer to the ECM field to write */
    boost::shared_ptr<DynamicECMField> mpECMField;
    
    /** Output directory */
    std::string mOutputDirectory;
    
    /** Sampling timestep multiple (how often to write) */
    unsigned mSamplingTimestepMultiple;
    
    /** Counter for output files */
    unsigned mOutputCounter;
    
    /** PVD file stream */
    out_stream mpPvdFile;

public:
    /**
     * Constructor.
     */
    ECMFieldWriter(boost::shared_ptr<DynamicECMField> pECMField,
                   unsigned samplingTimestepMultiple = 12)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mpECMField(pECMField),
          mOutputDirectory(""),
          mSamplingTimestepMultiple(samplingTimestepMultiple),
          mOutputCounter(0)
    {
    }
    
    /**
     * Destructor.
     */
    virtual ~ECMFieldWriter()
    {
        // Close PVD file if open
        if (mpPvdFile)
        {
            *mpPvdFile << "  </Collection>\n";
            *mpPvdFile << "</VTKFile>\n";
            mpPvdFile->close();
        }
    }
    
    /**
     * Update at end of timestep - write ECM field if needed.
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        ScopedTimer _prof("ECMFieldWriter");
        // Check if we should output this timestep
        if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingTimestepMultiple == 0)
        {
            unsigned current_timestep = SimulationTime::Instance()->GetTimeStepsElapsed();
            double current_time = SimulationTime::Instance()->GetTime();
            
            // File extension depends on grid type (vti for square, vtu for hex)
            std::string ext = mpECMField->GetOutputExtension();
            
            // Create filename using timestep (to match VTU file naming convention)
            std::stringstream fname_stream;
            fname_stream << "ecm_grid_" << current_timestep << ext;
            std::string ecm_filename = fname_stream.str();
            
            std::string full_path = mOutputDirectory + "/" + ecm_filename;
            
            // Write ECM field (picks VTI or VTU based on grid type)
            mpECMField->WriteOutput(full_path, current_time);
            
            // Add entry to PVD file (use timestep number, not real time,
            // so playback synchronizes with results.pvd)
            if (mpPvdFile)
            {
                *mpPvdFile << "    <DataSet timestep=\"" << current_timestep 
                          << "\" group=\"\" part=\"0\" file=\"" << ecm_filename << "\"/>\n";
            }
            
            mOutputCounter++;
        }
    }
    
    /**
     * Setup solver - write initial ECM grid and create PVD file.
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
    {
        // outputDirectory is relative to CHASTE_TEST_OUTPUT
        // Use OutputFileHandler to get the absolute path
        OutputFileHandler output_file_handler(outputDirectory, false);
        std::string full_output_dir = output_file_handler.GetOutputDirectoryFullPath();
        
        mOutputDirectory = full_output_dir;
        
        // Create PVD file
        std::string pvd_filename = full_output_dir + "ecm_results.pvd";
        mpPvdFile = output_file_handler.OpenOutputFile("ecm_results.pvd");
        
        // Write PVD header
        *mpPvdFile << "<?xml version=\"1.0\"?>\n";
        *mpPvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        *mpPvdFile << "  <Collection>\n";
        
        // Write initial ECM field (timestep 0)
        std::string ext = mpECMField->GetOutputExtension();
        std::string ecm_filename = "ecm_grid_0" + ext;
        std::string full_path = full_output_dir + ecm_filename;
        mpECMField->WriteOutput(full_path, 0.0);
        
        // Add initial timestep to PVD (timestep 0)
        *mpPvdFile << "    <DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"" << ecm_filename << "\"/>\n";
    }
    
    /**
     * Output parameters to file.
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
        
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};

#endif /* ECMFIELDWRITER_HPP_ */
