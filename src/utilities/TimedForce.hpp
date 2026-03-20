/*
 * TimedForce.hpp
 *
 * Decorator that wraps any AbstractForce and accumulates wall-clock
 * time spent in AddForceContribution().  Call PrintTimingSummary()
 * at the end of the simulation to see the breakdown.
 *
 * Usage:
 *   auto p_spring = boost::make_shared<SurfaceSpringForce<3>>();
 *   // ... configure p_spring ...
 *   auto p_timed = boost::make_shared<TimedForce<3>>(p_spring, "SurfaceSpring");
 *   simulator.AddForce(p_timed);
 *
 *   // After Solve():
 *   TimedForce<DIM>::PrintAllTimings();
 */
#ifndef TIMEDFORCE_HPP_
#define TIMEDFORCE_HPP_

#include "AbstractForce.hpp"
#include "SimProfiler.hpp"
#include <boost/shared_ptr.hpp>
#include <chrono>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

template <unsigned DIM>
class TimedForce : public AbstractForce<DIM>
{
private:
    boost::shared_ptr<AbstractForce<DIM>> mpWrappedForce;
    std::string mLabel;
    double mTotalSeconds;
    unsigned long mCallCount;

    // Global registry for all TimedForce instances (for PrintAllTimings)
    static std::vector<TimedForce<DIM>*>& GetRegistry()
    {
        static std::vector<TimedForce<DIM>*> registry;
        return registry;
    }

public:
    TimedForce(boost::shared_ptr<AbstractForce<DIM>> pForce,
               const std::string& label)
        : AbstractForce<DIM>(),
          mpWrappedForce(pForce),
          mLabel(label),
          mTotalSeconds(0.0),
          mCallCount(0)
    {
        GetRegistry().push_back(this);
    }

    virtual ~TimedForce()
    {
        auto& reg = GetRegistry();
        reg.erase(std::remove(reg.begin(), reg.end(), this), reg.end());
    }

    void AddForceContribution(
        AbstractCellPopulation<DIM>& rCellPopulation) override
    {
        ScopedTimer timer(mLabel);
        auto t0 = std::chrono::steady_clock::now();
        mpWrappedForce->AddForceContribution(rCellPopulation);
        auto t1 = std::chrono::steady_clock::now();
        mTotalSeconds += std::chrono::duration<double>(t1 - t0).count();
        mCallCount++;
    }

    const std::string& GetLabel() const { return mLabel; }
    double GetTotalSeconds() const { return mTotalSeconds; }
    unsigned long GetCallCount() const { return mCallCount; }
    double GetAvgMicroseconds() const
    {
        return mCallCount > 0 ? (mTotalSeconds * 1e6 / mCallCount) : 0.0;
    }

    void ResetTiming()
    {
        mTotalSeconds = 0.0;
        mCallCount = 0;
    }

    /** Access the wrapped force for configuration */
    boost::shared_ptr<AbstractForce<DIM>> GetWrappedForce() { return mpWrappedForce; }

    /**
     * Print a timing summary for ALL TimedForce instances, sorted by
     * total time descending.
     */
    static void PrintAllTimings()
    {
        auto& reg = GetRegistry();
        if (reg.empty()) return;

        // Sort by total time descending
        std::vector<TimedForce<DIM>*> sorted(reg.begin(), reg.end());
        std::sort(sorted.begin(), sorted.end(),
                  [](const TimedForce* a, const TimedForce* b)
                  { return a->mTotalSeconds > b->mTotalSeconds; });

        double grandTotal = 0.0;
        for (auto* tf : sorted) grandTotal += tf->mTotalSeconds;

        std::cout << "\n============================================\n";
        std::cout << "  FORCE PROFILING SUMMARY\n";
        std::cout << "============================================\n";
        std::cout << std::left << std::setw(28) << "Force"
                  << std::right << std::setw(10) << "Total(s)"
                  << std::setw(8) << "%"
                  << std::setw(12) << "Calls"
                  << std::setw(12) << "Avg(µs)"
                  << "\n";
        std::cout << std::string(70, '-') << "\n";

        for (auto* tf : sorted)
        {
            double pct = grandTotal > 0
                             ? 100.0 * tf->mTotalSeconds / grandTotal
                             : 0.0;
            std::cout << std::left << std::setw(28) << tf->mLabel
                      << std::right << std::fixed << std::setprecision(2)
                      << std::setw(10) << tf->mTotalSeconds
                      << std::setw(7) << pct << "%"
                      << std::setw(12) << tf->mCallCount
                      << std::setw(12) << std::setprecision(1)
                      << tf->GetAvgMicroseconds()
                      << "\n";
        }

        std::cout << std::string(70, '-') << "\n";
        std::cout << std::left << std::setw(28) << "TOTAL FORCES"
                  << std::right << std::fixed << std::setprecision(2)
                  << std::setw(10) << grandTotal << "\n";
        std::cout << "============================================\n\n";
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<Label>" << mLabel << "</Label>\n";
        mpWrappedForce->OutputForceParameters(rParamsFile);
    }
};

#endif // TIMEDFORCE_HPP_
