/*
 * ViscoelasticGhostNodeEcmWriter.hpp
 *
 * Simulation modifier that writes the viscoelastic ghost node ECM field to VTP
 * files (PolyData) at each sampling interval, and maintains a PVD collection
 * file for time-series playback in ParaView.
 *
 * Identical to GhostNodeEcmWriter but operates on ViscoelasticGhostNodeEcmField.
 */
#ifndef VISCOELASTICGHOSTNODEECMWRITER_HPP_
#define VISCOELASTICGHOSTNODEECMWRITER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "SimProfiler.hpp"
#include "OutputFileHandler.hpp"
#include "SimulationTime.hpp"
#include <boost/shared_ptr.hpp>

template<unsigned DIM>
class ViscoelasticGhostNodeEcmWriter : public AbstractCellBasedSimulationModifier<DIM>
{
private:
    boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> mpGhostField;
    std::string mOutputDirectory;
    unsigned mSamplingTimestepMultiple;
    unsigned mOutputCounter;
    out_stream mpPvdFile;

public:
    ViscoelasticGhostNodeEcmWriter(boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> pField,
                                   unsigned samplingTimestepMultiple = 12)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mpGhostField(pField),
          mOutputDirectory(""),
          mSamplingTimestepMultiple(samplingTimestepMultiple),
          mOutputCounter(0)
    {
    }

    virtual ~ViscoelasticGhostNodeEcmWriter()
    {
        if (mpPvdFile)
        {
            *mpPvdFile << "  </Collection>\n";
            *mpPvdFile << "</VTKFile>\n";
            mpPvdFile->close();
        }
    }

    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
    {
        ScopedTimer _prof("ViscoelasticGhostNodeEcmWriter");

        if (SimulationTime::Instance()->GetTimeStepsElapsed() % mSamplingTimestepMultiple == 0)
        {
            unsigned current_timestep = SimulationTime::Instance()->GetTimeStepsElapsed();
            double current_time = SimulationTime::Instance()->GetTime();

            std::stringstream fname_stream;
            fname_stream << "ghost_ecm_" << current_timestep << ".vtp";
            std::string filename = fname_stream.str();

            std::string full_path = mOutputDirectory + "/" + filename;
            mpGhostField->WriteOutput(full_path, current_time);

            if (mpPvdFile)
            {
                *mpPvdFile << "    <DataSet timestep=\"" << current_timestep
                           << "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n";
            }

            mOutputCounter++;
        }
    }

    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
    {
        OutputFileHandler output_file_handler(outputDirectory, false);
        mOutputDirectory = output_file_handler.GetOutputDirectoryFullPath();

        mpPvdFile = output_file_handler.OpenOutputFile("ghost_ecm_results.pvd");
        *mpPvdFile << "<?xml version=\"1.0\"?>\n";
        *mpPvdFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        *mpPvdFile << "  <Collection>\n";

        std::string filename = "ghost_ecm_0.vtp";
        std::string full_path = mOutputDirectory + "/" + filename;
        mpGhostField->WriteOutput(full_path, 0.0);

        *mpPvdFile << "    <DataSet timestep=\"0\" group=\"\" part=\"0\" file=\"" << filename << "\"/>\n";
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<SamplingTimestepMultiple>" << mSamplingTimestepMultiple << "</SamplingTimestepMultiple>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};

#endif // VISCOELASTICGHOSTNODEECMWRITER_HPP_
