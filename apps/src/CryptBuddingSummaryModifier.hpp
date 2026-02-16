/*
 * CryptBuddingSummaryModifier.hpp
 *
 * Simulation modifier that writes radial statistics to CSV and monitors
 * numerical stability by tracking cell volume ratios.
 */
#ifndef CRYPTBUDDINGSUMMARYMODIFIER_HPP_
#define CRYPTBUDDINGSUMMARYMODIFIER_HPP_

#include <cmath>
#include <string>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "AbstractCellBasedSimulationModifier.hpp"
#include "SimulationTime.hpp"
#include "OutputFileHandler.hpp"
#include "Exception.hpp"

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
    double mReferenceVolume;      // expected average cell volume
    unsigned mInstabilityCount;   // consecutive instability warnings

public:
    CryptBuddingSummaryModifier(double stiffness, unsigned samplingMultiple,
                                double endTime = 200.0, double referenceVolume = 1.0)
        : AbstractCellBasedSimulationModifier<DIM>(),
          mStiffness(stiffness),
          mEndTime(endTime),
          mHeaderWritten(false),
          mSamplingMultiple(samplingMultiple),
          mLogInterval(samplingMultiple / 6),
          mLastOutputStep(0),
          mLastLogStep(0),
          mReferenceVolume(referenceVolume),
          mInstabilityCount(0)
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

        // Progress log — stdout is real-time in app mode (no CTest buffering)
        if (current_step == 0 || current_step - mLastLogStep >= mLogInterval)
        {
            mLastLogStep = current_step;
            double pct = (mEndTime > 0.0) ? (current_time / mEndTime) * 100.0 : 0.0;
            std::cout << "\r[Progress] t=" << std::fixed << std::setprecision(1)
                      << current_time << "h / " << mEndTime << "h  ("
                      << std::setprecision(1) << pct << "%)  cells="
                      << num_cells << "    " << std::flush;
        }

        // ---- Numerical stability monitoring ----
        double max_vol = 0.0;
        double min_vol = 1e30;
        for (typename AbstractCellPopulation<DIM>::Iterator it = rCellPopulation.Begin();
             it != rCellPopulation.End(); ++it)
        {
            double vol = it->GetCellData()->GetItem("volume");
            if (vol > max_vol) max_vol = vol;
            if (vol < min_vol) min_vol = vol;
        }

        double vol_ratio = (mReferenceVolume > 0.0) ? max_vol / mReferenceVolume : 0.0;

        if (vol_ratio > 10.0)
        {
            mInstabilityCount++;
            std::cerr << "[INSTABILITY WARNING] t=" << std::fixed << std::setprecision(3)
                      << current_time << "  step=" << current_step
                      << "  max_volume=" << std::setprecision(2) << max_vol
                      << " (" << vol_ratio << "x reference)"
                      << "  min_volume=" << std::setprecision(4) << min_vol
                      << "  (consecutive warnings: " << mInstabilityCount << ")"
                      << std::endl;

            if (vol_ratio > 1000.0)
            {
                std::cerr << "[FATAL INSTABILITY] Volume ratio " << vol_ratio
                          << "x exceeds catastrophic threshold (1000x reference)."
                          << " Simulation is numerically diverging. Consider reducing dt."
                          << std::endl;
                EXCEPTION("Numerical instability detected: cell volume "
                          + std::to_string(max_vol)
                          + " exceeds 1000x reference volume "
                          + std::to_string(mReferenceVolume)
                          + " at t=" + std::to_string(current_time));
            }
        }
        else
        {
            mInstabilityCount = 0;  // reset if volumes return to normal
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
            file << "time,num_cells,mean_r,var_r,max_r,min_r,r_range,max_vol,min_vol,stiffness" << std::endl;
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
             << (max_r - min_r) << ","
             << max_vol << "," << min_vol << ","
             << mStiffness
             << std::endl;
        file.close();
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t<Stiffness>" << mStiffness << "</Stiffness>\n";
        AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};

#endif // CRYPTBUDDINGSUMMARYMODIFIER_HPP_
