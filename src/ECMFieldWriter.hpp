/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

*/

#ifndef ECMFIELDWRITER_HPP_
#define ECMFIELDWRITER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "DynamicECMField.hpp"
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
        // Check if we should output this timestep
        if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingTimestepMultiple == 0)
        {
            unsigned current_timestep = SimulationTime::Instance()->GetTimeStepsElapsed();
            double current_time = SimulationTime::Instance()->GetTime();
            
            // Create VTI filename using timestep (to match VTU file naming convention)
            std::stringstream vti_filename_stream;
            vti_filename_stream << "ecm_grid_" << current_timestep << ".vti";
            std::string vti_filename = vti_filename_stream.str();
            
            std::string full_vti_path = mOutputDirectory + "/" + vti_filename;
            
            // Write ECM field to VTI file (XML ImageData format)
            mpECMField->WriteToVTI(full_vti_path, current_time);
            
            // Add entry to PVD file
            if (mpPvdFile)
            {
                *mpPvdFile << "    <DataSet timestep=\"" << current_time 
                          << "\" group=\"\" part=\"0\" file=\"" << vti_filename << "\"/>\n";
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
        
        // Write initial ECM field (timestep 0) as VTI
        std::string vti_filename = "ecm_grid_0.vti";
        std::string full_vti_path = full_output_dir + vti_filename;
        mpECMField->WriteToVTI(full_vti_path, 0.0);
        
        // Add initial timestep to PVD
        *mpPvdFile << "    <DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"" << vti_filename << "\"/>\n";
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
