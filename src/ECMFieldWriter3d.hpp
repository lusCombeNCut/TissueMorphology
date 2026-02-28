/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef ECMFIELDWRITER3D_HPP_
#define ECMFIELDWRITER3D_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "DynamicECMField3d.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include <boost/shared_ptr.hpp>

/**
 * A simulation modifier that writes the 3D ECM field to VTI files
 * at regular intervals, with a PVD index for time-series playback in ParaView.
 */
class ECMFieldWriter3d : public AbstractCellBasedSimulationModifier<3>
{
private:
    /** Pointer to the 3D ECM field */
    boost::shared_ptr<DynamicECMField3d> mpECMField;
    
    /** Output directory (absolute path) */
    std::string mOutputDirectory;
    
    /** How often to write (every N timesteps) */
    unsigned mSamplingTimestepMultiple;
    
    /** Counter for output files */
    unsigned mOutputCounter;
    
    /** PVD file stream for time-series index */
    out_stream mpPvdFile;

    /** If true, write only a z-center slice instead of the full 3D volume */
    bool mSliceOnly;

public:
    /**
     * Constructor.
     * 
     * @param pECMField Pointer to the 3D ECM field to visualize
     * @param samplingTimestepMultiple Write every N timesteps
     * @param sliceOnly If true, save only the center z-slice (2D) instead of the full 3D volume
     */
    ECMFieldWriter3d(boost::shared_ptr<DynamicECMField3d> pECMField,
                     unsigned samplingTimestepMultiple = 12,
                     bool sliceOnly = true)
        : AbstractCellBasedSimulationModifier<3>(),
          mpECMField(pECMField),
          mOutputDirectory(""),
          mSamplingTimestepMultiple(samplingTimestepMultiple),
          mOutputCounter(0),
          mSliceOnly(sliceOnly)
    {
    }
    
    /**
     * Destructor — close PVD file.
     */
    virtual ~ECMFieldWriter3d()
    {
        if (mpPvdFile)
        {
            *mpPvdFile << "  </Collection>\n";
            *mpPvdFile << "</VTKFile>\n";
            mpPvdFile->close();
        }
    }
    
    /**
     * Called at the end of each timestep. Writes VTI file if the sampling
     * multiple is reached.
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<3, 3>& rCellPopulation)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingTimestepMultiple == 0)
        {
            unsigned current_timestep = SimulationTime::Instance()->GetTimeStepsElapsed();
            double current_time = SimulationTime::Instance()->GetTime();
            
            std::stringstream vti_filename_stream;
            vti_filename_stream << "ecm3d_" << current_timestep << ".vti";
            std::string vti_filename = vti_filename_stream.str();
            
            std::string full_vti_path = mOutputDirectory + "/" + vti_filename;
            
            if (mSliceOnly)
            {
                mpECMField->WriteSliceToVTI(full_vti_path, current_time);
            }
            else
            {
                mpECMField->WriteToVTI(full_vti_path, current_time);
            }
            
            if (mpPvdFile)
            {
                *mpPvdFile << "    <DataSet timestep=\"" << current_time
                           << "\" group=\"\" part=\"0\" file=\"" << vti_filename << "\"/>\n";
            }
            
            mOutputCounter++;
        }
    }
    
    /**
     * Called at simulation start. Creates output directory and PVD file,
     * writes the initial ECM state.
     */
    virtual void SetupSolve(AbstractCellPopulation<3, 3>& rCellPopulation,
                             std::string outputDirectory)
    {
        OutputFileHandler output_file_handler(outputDirectory, false);
        mOutputDirectory = output_file_handler.GetOutputDirectoryFullPath();
        
        // Create PVD time-series index file
        mpPvdFile = output_file_handler.OpenOutputFile("ecm3d_results.pvd");
        *mpPvdFile << "<?xml version=\"1.0\"?>\n";
        *mpPvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        *mpPvdFile << "  <Collection>\n";
        
        // Write initial state
        std::string initial_vti = mOutputDirectory + "/ecm3d_0.vti";
        if (mSliceOnly)
        {
            mpECMField->WriteSliceToVTI(initial_vti, 0.0);
        }
        else
        {
            mpECMField->WriteToVTI(initial_vti, 0.0);
        }
        
        if (mpPvdFile)
        {
            *mpPvdFile << "    <DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"ecm3d_0.vti\"/>\n";
        }
        
        mOutputCounter = 1;
    }
    
    /**
     * Required override — no action needed at end of solve beyond destructor.
     */
    virtual void UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<3, 3>& rCellPopulation)
    {
    }
    
    /**
     * Required override for simulation modifier.
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<SamplingTimestepMultiple>"
                     << mSamplingTimestepMultiple
                     << "</SamplingTimestepMultiple>\n";
        AbstractCellBasedSimulationModifier<3>::OutputSimulationModifierParameters(rParamsFile);
    }
};

#endif /* ECMFIELDWRITER3D_HPP_ */
