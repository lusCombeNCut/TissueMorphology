/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef LUMENPRESSUREFORCE_HPP_
#define LUMENPRESSUREFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimulationTime.hpp"
#include "SimProfiler.hpp"
#include "RingTopologyTracker.hpp"

/**
 * Hydrostatic lumen pressure force using pressure–volume work.
 *
 * Models the lumen as a fluid-filled cavity exerting uniform hydrostatic
 * pressure P on the apical (inner) surface of the epithelium. 
 *
 * Two modes are supported:
 *
 * 1. CONSTANT PRESSURE MODE (default, mUseTargetVolume=false):
 *    F_i = P × A_i × n̂_i
 *    where P is constant.
 *
 * 2. TARGET VOLUME MODE (mUseTargetVolume=true):
 *    Models an incompressible fluid with slowly growing volume.
 *    Target volume grows as: V₀(t) = V_initial × (1 + growth_rate × t)
 *    Pressure is: P = bulk_modulus × (V₀ - V_current) / V₀
 *    This creates feedback: if lumen is compressed below target, pressure increases.
 *
 * Surface area computation:
 *   2D (ring): A_i = ½(|x_L − x_i| + |x_R − x_i|) where L,R are ring
 *              neighbors from RingTopologyTracker (arc-length element).
 *   3D (sphere, no ring): A_i = 4πr²/N (equal Voronoi partition estimate).
 */
template<unsigned DIM>
class LumenPressureForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Hydrostatic lumen pressure (constant mode) or bulk modulus (target volume mode) */
    double mPressure;

    /** Center of the lumen (auto-tracked from centroid if mTrackCenter=true) */
    c_vector<double, DIM> mLumenCenter;

    /** Whether to auto-track lumen center from population centroid */
    bool mTrackCenter;

    /** Pointer to ring topology tracker (2D only, for computing arc lengths) */
    const RingTopologyTracker<DIM>* mpRingTracker;

    // ---- Target volume mode parameters ----
    
    /** Whether to use target volume mode (incompressible fluid) */
    bool mUseTargetVolume;
    
    /** Initial target lumen volume/area */
    double mInitialTargetVolume;
    
    /** Rate of volume growth (fraction per hour, e.g. 0.1 = 10% per hour) */
    double mVolumeGrowthRate;
    
    /** Bulk modulus - how stiff the fluid is (pressure per fractional compression) */
    double mBulkModulus;
    
    /** Current simulation time (updated each timestep) */
    double mCurrentTime;
    
    /** Whether initial volume has been set */
    bool mInitialVolumeSet;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mPressure;
        archive & mLumenCenter;
        archive & mTrackCenter;
        archive & mUseTargetVolume;
        archive & mInitialTargetVolume;
        archive & mVolumeGrowthRate;
        archive & mBulkModulus;
        archive & mCurrentTime;
        archive & mInitialVolumeSet;
    }

public:

    LumenPressureForce()
        : AbstractForce<DIM>(),
          mPressure(1.0),
          mTrackCenter(true),
          mpRingTracker(nullptr),
          mUseTargetVolume(false),
          mInitialTargetVolume(0.0),
          mVolumeGrowthRate(0.05),  // 5% per hour default
          mBulkModulus(10.0),
          mCurrentTime(0.0),
          mInitialVolumeSet(false)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mLumenCenter[i] = 0.0;
        }
    }

    virtual ~LumenPressureForce() {}

    /**
     * Compute the current lumen volume/area from cell positions.
     * In 2D, uses the shoelace formula for the polygon area.
     * In 3D, estimates volume from mean radius.
     */
    double ComputeCurrentLumenVolume(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        if constexpr (DIM == 2)
        {
            // Use shoelace formula for polygon area
            // First, collect all cell positions and sort by angle
            std::vector<std::pair<double, c_vector<double, 2>>> anglePos;
            
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                c_vector<double, 2> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double angle = std::atan2(pos[1] - mLumenCenter[1], pos[0] - mLumenCenter[0]);
                anglePos.push_back(std::make_pair(angle, pos));
            }
            
            // Sort by angle
            std::sort(anglePos.begin(), anglePos.end(),
                [](const std::pair<double, c_vector<double, 2>>& a, 
                   const std::pair<double, c_vector<double, 2>>& b) {
                    return a.first < b.first;
                });
            
            // Shoelace formula
            double area = 0.0;
            size_t n = anglePos.size();
            for (size_t i = 0; i < n; i++)
            {
                size_t j = (i + 1) % n;
                area += anglePos[i].second[0] * anglePos[j].second[1];
                area -= anglePos[j].second[0] * anglePos[i].second[1];
            }
            return std::abs(area) / 2.0;
        }
        else // DIM == 3
        {
            // Estimate volume from mean radius: V = (4/3)πr³
            double sum_r = 0.0;
            unsigned count = 0;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                sum_r += norm_2(pos - mLumenCenter);
                count++;
            }
            double mean_r = sum_r / count;
            return (4.0 / 3.0) * M_PI * mean_r * mean_r * mean_r;
        }
    }

    /**
     * Hydrostatic pressure force: F_i = P × A_i × n̂_i
     *
     * n̂_i points radially outward from the lumen centroid.
     * A_i depends on dimension and topology (see class docstring).
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        ScopedTimer _prof("LumenPressure");
        if (mTrackCenter)
        {
            UpdateLumenCenter(rCellPopulation);
        }

        // Compute effective pressure
        double effectivePressure = mPressure;
        
        if (mUseTargetVolume)
        {
            // Compute current lumen volume
            double currentVolume = ComputeCurrentLumenVolume(rCellPopulation);
            
            // Initialize target volume on first call
            if (!mInitialVolumeSet)
            {
                mInitialTargetVolume = currentVolume;
                mInitialVolumeSet = true;
            }
            
            // Target volume grows over time: V₀(t) = V_initial × (1 + growth_rate × t)
            double currentTime = SimulationTime::Instance()->GetTime();
            double targetVolume = mInitialTargetVolume * (1.0 + mVolumeGrowthRate * currentTime);
            
            // Pressure = bulk_modulus × (V_target - V_current) / V_target
            // Positive when compressed (V_current < V_target), pushing outward
            double volumeDeficit = (targetVolume - currentVolume) / targetVolume;
            effectivePressure = mBulkModulus * volumeDeficit;
            
            // Clamp to prevent negative pressure (lumen pulling inward)
            if (effectivePressure < 0.0)
            {
                effectivePressure = 0.0;
            }
        }

        // Detect vertex-based population (for force distribution)
        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        // Build index→position map for neighbor lookups in 2D ring mode
        std::map<unsigned, c_vector<double, DIM>> indexToPos;
        if (DIM == 2 && mpRingTracker != nullptr)
        {
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                indexToPos[idx] = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            }
        }

        // Count cells for 3D equal-partition area estimate
        unsigned numCells = rCellPopulation.GetNumRealCells();

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, DIM> disp = pos - mLumenCenter;
            double r = norm_2(disp);

            // Skip cells at the exact center (would give zero normal)
            if (r < 1e-6)
            {
                continue;
            }

            c_vector<double, DIM> n_hat = disp / r;  // outward radial unit normal

            // Compute effective apical surface area element A_i
            double area_element = 0.0;

            if (DIM == 2 && mpRingTracker != nullptr)
            {
                // 2D ring: arc length = ½(|x_L - x_i| + |x_R - x_i|)
                auto [leftIdx, rightIdx] = mpRingTracker->GetNeighbors(loc_index);
                if (indexToPos.count(leftIdx) && indexToPos.count(rightIdx))
                {
                    double dL = norm_2(indexToPos[leftIdx] - pos);
                    double dR = norm_2(indexToPos[rightIdx] - pos);
                    area_element = 0.5 * (dL + dR);
                }
                else
                {
                    // Fallback: equal partition of circumference
                    area_element = 2.0 * M_PI * r / numCells;
                }
            }
            else if (DIM == 2)
            {
                // 2D without ring tracker: equal partition of circumference
                area_element = 2.0 * M_PI * r / numCells;
            }
            else // DIM == 3
            {
                // 3D: equal Voronoi partition of sphere surface
                area_element = 4.0 * M_PI * r * r / numCells;
            }

            // F_i = P × A_i × n̂_i
            c_vector<double, DIM> force = effectivePressure * area_element * n_hat;

            // Apply force
            if (p_vertex_pop)
            {
                VertexElement<DIM, DIM>* p_element = p_vertex_pop->rGetMesh().GetElement(loc_index);
                unsigned n_nodes = p_element->GetNumNodes();
                c_vector<double, DIM> force_per_node = force / static_cast<double>(n_nodes);
                for (unsigned i = 0; i < n_nodes; i++)
                {
                    p_element->GetNode(i)->AddAppliedForceContribution(force_per_node);
                }
            }
            else
            {
                rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(force);
            }
        }
    }

    void UpdateLumenCenter(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> center = zero_vector<double>(DIM);
        unsigned count = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            center += rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            count++;
        }

        if (count > 0)
        {
            mLumenCenter = center / static_cast<double>(count);
        }
    }

    // ---- Setters/Getters ----

    void SetPressure(double pressure)
    {
        mPressure = pressure;
    }

    double GetPressure() const { return mPressure; }

    void SetLumenCenter(c_vector<double, DIM> center)
    {
        mLumenCenter = center;
        mTrackCenter = false;
    }

    c_vector<double, DIM> GetLumenCenter() const { return mLumenCenter; }

    void SetTrackCenter(bool track) { mTrackCenter = track; }
    bool GetTrackCenter() const { return mTrackCenter; }

    void SetRingTopologyTracker(const RingTopologyTracker<DIM>* pTracker)
    {
        mpRingTracker = pTracker;
    }

    // ---- Target Volume Mode (incompressible fluid) ----
    
    /** Enable target volume mode with growing lumen */
    void SetUseTargetVolume(bool use) { mUseTargetVolume = use; }
    bool GetUseTargetVolume() const { return mUseTargetVolume; }
    
    /** Set initial target volume (if not set, auto-detected from first timestep) */
    void SetInitialTargetVolume(double v) 
    { 
        mInitialTargetVolume = v; 
        mInitialVolumeSet = true;
    }
    double GetInitialTargetVolume() const { return mInitialTargetVolume; }
    
    /** Volume growth rate (fraction per hour, e.g. 0.1 = 10% per hour) */
    void SetVolumeGrowthRate(double rate) { mVolumeGrowthRate = rate; }
    double GetVolumeGrowthRate() const { return mVolumeGrowthRate; }
    
    /** Bulk modulus - stiffness of the incompressible fluid */
    void SetBulkModulus(double k) { mBulkModulus = k; }
    double GetBulkModulus() const { return mBulkModulus; }
    
    /** Update current simulation time (call from simulation modifier) */
    void SetCurrentTime(double t) { mCurrentTime = t; }
    double GetCurrentTime() const { return mCurrentTime; }

    // Keep old setters as no-ops for backward compatibility during transition
    void SetPressureStrength(double s) { mPressure = s; }
    void SetLumenEquilibriumRadius(double /*r*/) { /* no-op: no longer used */ }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Pressure>" << mPressure << "</Pressure>\n";
        *rParamsFile << "\t\t\t<TrackCenter>" << mTrackCenter << "</TrackCenter>\n";
        *rParamsFile << "\t\t\t<UseTargetVolume>" << mUseTargetVolume << "</UseTargetVolume>\n";
        *rParamsFile << "\t\t\t<VolumeGrowthRate>" << mVolumeGrowthRate << "</VolumeGrowthRate>\n";
        *rParamsFile << "\t\t\t<BulkModulus>" << mBulkModulus << "</BulkModulus>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenPressureForce)

#endif /*LUMENPRESSUREFORCE_HPP_*/
