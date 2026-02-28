/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

*/

#ifndef SIMPROFILER_HPP_
#define SIMPROFILER_HPP_

#include <chrono>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <fstream>

/**
 * Lightweight profiling timer for Chaste simulations.
 *
 * Usage:
 *   SimProfiler::Instance().Start("ECMConfinement");
 *   ... expensive work...
 *   SimProfiler::Instance().Stop("ECMConfinement");
 *
 *   SimProfiler::Instance().PrintSummary();
 *
 * All times are wall-clock microseconds.  The profiler is a process-wide
 * singleton so it can be accessed from any force / modifier without passing
 * pointers around.
 */
class SimProfiler
{
public:
    /** Access the singleton instance. */
    static SimProfiler& Instance()
    {
        static SimProfiler instance;
        return instance;
    }

    /** Start timing a named section. */
    void Start(const std::string& name)
    {
        mRunning[name] = Clock::now();
    }

    /** Stop timing a named section and accumulate elapsed time. */
    void Stop(const std::string& name)
    {
        auto it = mRunning.find(name);
        if (it == mRunning.end()) return;

        auto elapsed = std::chrono::duration_cast<Microseconds>(
            Clock::now() - it->second).count();

        auto& entry = mSections[name];
        entry.totalMicroseconds += elapsed;
        entry.callCount++;
        if (elapsed > entry.maxMicroseconds) entry.maxMicroseconds = elapsed;

        mRunning.erase(it);
    }

    /** Mark the beginning of a new timestep (for per-step statistics). */
    void BeginTimestep()
    {
        mTimestepStart = Clock::now();
        mTimestepCount++;
    }

    /** Mark the end of the current timestep. */
    void EndTimestep()
    {
        auto elapsed = std::chrono::duration_cast<Microseconds>(
            Clock::now() - mTimestepStart).count();
        mTimestepTotalUs += elapsed;
        if (elapsed > mTimestepMaxUs) mTimestepMaxUs = elapsed;
    }

    /** Print a human-readable summary table to stdout. */
    void PrintSummary() const
    {
        if (mSections.empty())
        {
            std::cout << "[SimProfiler] No sections recorded.\n";
            return;
        }

        // Sort by total time descending
        std::vector<std::pair<std::string, SectionData>> sorted(
            mSections.begin(), mSections.end());
        std::sort(sorted.begin(), sorted.end(),
            [](const auto& a, const auto& b){
                return a.second.totalMicroseconds > b.second.totalMicroseconds;
            });

        int64_t wallTotal = 0;
        for (auto& s : sorted) wallTotal += s.second.totalMicroseconds;

        std::cout << "\n";
        std::cout << "╔══════════════════════════════════════════════════════════════════════════════╗\n";
        std::cout << "║                         SIMULATION PROFILER SUMMARY                         ║\n";
        std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";

        if (mTimestepCount > 0)
        {
            double avgStepMs = (double)mTimestepTotalUs / mTimestepCount / 1000.0;
            double maxStepMs = (double)mTimestepMaxUs / 1000.0;
            double totalStepS = (double)mTimestepTotalUs / 1e6;
            std::cout << "║  Timesteps: " << std::setw(8) << mTimestepCount
                      << "   Avg: " << std::fixed << std::setprecision(2) << std::setw(8) << avgStepMs << " ms"
                      << "   Max: " << std::setw(8) << maxStepMs << " ms"
                      << "   Total: " << std::setw(7) << std::setprecision(1) << totalStepS << " s"
                      << "  ║\n";
            std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
        }

        // Header
        std::cout << "║  " << std::left << std::setw(28) << "Section"
                  << std::right << std::setw(9) << "Total(ms)"
                  << std::setw(9) << "Calls"
                  << std::setw(10) << "Avg(µs)"
                  << std::setw(10) << "Max(µs)"
                  << std::setw(8) << "Pct%"
                  << "  ║\n";
        std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";

        for (auto& [name, data] : sorted)
        {
            double totalMs = data.totalMicroseconds / 1000.0;
            double avgUs = (data.callCount > 0) ? (double)data.totalMicroseconds / data.callCount : 0.0;
            double pct = (wallTotal > 0) ? 100.0 * data.totalMicroseconds / wallTotal : 0.0;

            std::cout << "║  " << std::left << std::setw(28) << name
                      << std::right << std::fixed
                      << std::setprecision(1) << std::setw(9) << totalMs
                      << std::setw(9) << data.callCount
                      << std::setprecision(0) << std::setw(10) << avgUs
                      << std::setw(10) << data.maxMicroseconds
                      << std::setprecision(1) << std::setw(7) << pct << "%"
                      << "  ║\n";
        }

        std::cout << "╠══════════════════════════════════════════════════════════════════════════════╣\n";
        double wallTotalMs = wallTotal / 1000.0;
        std::cout << "║  " << std::left << std::setw(28) << "TOTAL (instrumented)"
                  << std::right << std::fixed << std::setprecision(1)
                  << std::setw(9) << wallTotalMs
                  << std::setw(9) << ""
                  << std::setw(10) << ""
                  << std::setw(10) << ""
                  << std::setw(7) << "100.0" << "%"
                  << "  ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";
    }

    /** Write summary to a file (CSV). */
    void WriteSummaryCSV(const std::string& filepath) const
    {
        std::ofstream file(filepath);
        if (!file.is_open()) return;

        file << "section,total_us,calls,avg_us,max_us,pct\n";

        int64_t wallTotal = 0;
        for (auto& [name, data] : mSections)
            wallTotal += data.totalMicroseconds;

        for (auto& [name, data] : mSections)
        {
            double avgUs = (data.callCount > 0) ? (double)data.totalMicroseconds / data.callCount : 0.0;
            double pct = (wallTotal > 0) ? 100.0 * data.totalMicroseconds / wallTotal : 0.0;
            file << name << "," << data.totalMicroseconds << "," << data.callCount
                 << "," << (int64_t)avgUs << "," << data.maxMicroseconds
                 << "," << std::fixed << std::setprecision(2) << pct << "\n";
        }
        file.close();
    }

    /** Reset all counters. */
    void Reset()
    {
        mSections.clear();
        mRunning.clear();
        mTimestepCount = 0;
        mTimestepTotalUs = 0;
        mTimestepMaxUs = 0;
    }

private:
    using Clock = std::chrono::high_resolution_clock;
    using TimePoint = Clock::time_point;
    using Microseconds = std::chrono::microseconds;

    struct SectionData
    {
        int64_t totalMicroseconds = 0;
        int64_t callCount = 0;
        int64_t maxMicroseconds = 0;
    };

    std::map<std::string, SectionData> mSections;
    std::map<std::string, TimePoint> mRunning;

    int64_t mTimestepCount = 0;
    int64_t mTimestepTotalUs = 0;
    int64_t mTimestepMaxUs = 0;
    TimePoint mTimestepStart;

    SimProfiler() = default;
};

/**
 * RAII scoped timer — automatically starts on construction and stops on destruction.
 *
 * Usage:
 *   {
 *       ScopedTimer timer("MySection");
 *       ... work ...
 *   }  // automatically stops and records
 */
class ScopedTimer
{
public:
    explicit ScopedTimer(const std::string& name) : mName(name)
    {
        SimProfiler::Instance().Start(mName);
    }
    ~ScopedTimer()
    {
        SimProfiler::Instance().Stop(mName);
    }
private:
    std::string mName;
};

#endif /* SIMPROFILER_HPP_ */
