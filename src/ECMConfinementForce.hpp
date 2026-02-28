/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

*/

#ifndef ECMCONFINEMENTFORCE_HPP_
#define ECMCONFINEMENTFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "DynamicECMField.hpp"
#include "SimulationTime.hpp"
#include "SimProfiler.hpp"
#include <boost/shared_ptr.hpp>

/**
 * ECM Confinement Force
 *
 * Replaces the simple radial BasementMembraneForce with a fiber-density-based
 * confinement that naturally supports budding through local ECM degradation.
 *
 * Physics:
 *   Each cell samples the ECM density at its position. The confining force
 *   pushes the cell inward (toward the centroid) proportional to the local
 *   ECM density:
 *
 *     F_i = -k × ρ(x_i) × r̂_i
 *
 *   where k is the confinement stiffness, ρ ∈ [0,1] is the local ECM density,
 *   and r̂_i points radially outward from the tissue centroid.
 *
 *   Cells also degrade the ECM at their basal surface (MMP secretion):
 *     ρ_{n+1} = max(0, ρ_n - k_deg × dt)
 *
 *   This creates a natural feedback loop:
 *     proliferation → crowding → cells push outward → degrade ECM → 
 *     local density drops → confinement weakens → budding
 *
 * The force also supports cell traction remodeling of ECM fiber orientation
 * when enabled.
 *
 * Parameters:
 *   mConfinementStiffness  — force magnitude per unit density
 *   mDegradationEnabled    — whether cells actively degrade ECM
 *   mRemodelingEnabled     — whether cell traction aligns ECM fibers
 */
template<unsigned DIM>
class ECMConfinementForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Shared pointer to the ECM field */
    boost::shared_ptr<DynamicECMField> mpECMField;

    /** Confinement stiffness (force per unit density) */
    double mConfinementStiffness;

    /** Whether cells degrade ECM at their position */
    bool mDegradationEnabled;

    /** Whether cell traction remodels ECM fibers */
    bool mRemodelingEnabled;

    /** Whether to auto-track centroid */
    bool mTrackCenter;

    /** The tissue centroid (auto-updated or manually set) */
    c_vector<double, DIM> mCenter;

    /** Last update time (for computing dt) */
    double mLastUpdateTime;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mConfinementStiffness;
        archive & mDegradationEnabled;
        archive & mRemodelingEnabled;
        archive & mTrackCenter;
        archive & mCenter;
    }

public:

    ECMConfinementForce()
        : AbstractForce<DIM>(),
          mConfinementStiffness(5.0),
          mDegradationEnabled(true),
          mRemodelingEnabled(true),
          mTrackCenter(true),
          mLastUpdateTime(0.0)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mCenter[i] = 0.0;
        }
    }

    virtual ~ECMConfinementForce() {}

    /**
     * Apply ECM confinement and degradation.
     *
     * For each cell:
     *   1. Sample ECM density at cell position → ρ
     *   2. Apply inward force F = -k × ρ × r̂
     *   3. Degrade ECM at cell position (MMP secretion)
     *   4. Optionally remodel ECM fibers via traction
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        ScopedTimer _prof("ECMConfinement");
        assert(DIM == 2);  // 2D only for now (uses DynamicECMField which is 2D)

        if (!mpECMField)
        {
            EXCEPTION("ECM field not set! Call SetECMField() before simulation.");
        }

        // Update centroid
        if (mTrackCenter)
        {
            c_vector<double, DIM> center = zero_vector<double>(DIM);
            unsigned count = 0;
            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                center += rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                count++;
            }
            if (count > 0) mCenter = center / static_cast<double>(count);
        }

        // Compute dt for degradation/remodeling
        double current_time = SimulationTime::Instance()->GetTime();
        double dt = current_time - mLastUpdateTime;
        if (dt < 1e-10) dt = 0.005;  // fallback
        mLastUpdateTime = current_time;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, DIM> disp = pos - mCenter;
            double r = norm_2(disp);

            if (r < 1e-6) continue;

            c_vector<double, DIM> r_hat = disp / r;  // outward unit radial

            // Sample ECM density at cell position
            double density = mpECMField->GetDensityAt(pos);

            // Store for visualization
            cell_iter->GetCellData()->SetItem("ecm_density", density);

            // Confinement force: pushes cell inward proportional to local ECM density
            // F = -k × ρ × r̂
            c_vector<double, DIM> force = -mConfinementStiffness * density * r_hat;
            rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(force);

            // Cell degrades ECM at its position (MMP secretion)
            if (mDegradationEnabled)
            {
                mpECMField->DegradeECM(pos, dt);
            }

            // Cell traction remodels ECM fiber orientation
            if (mRemodelingEnabled)
            {
                // Cell's net traction = the applied force (from all other forces)
                // Use radial direction as proxy for traction direction
                c_vector<double, 2> traction;
                traction[0] = r_hat[0];
                traction[1] = r_hat[1];
                mpECMField->ApplyCellTraction(pos, traction, dt);
            }
        }

        // Diffuse ECM field (cross-linking / smoothing)
        {
            ScopedTimer _diffProf("ECM::DiffuseECM");
            mpECMField->DiffuseECM(dt);
        }
    }

    // ---- Setters ----

    void SetECMField(boost::shared_ptr<DynamicECMField> pField)
    {
        mpECMField = pField;
    }

    boost::shared_ptr<DynamicECMField> GetECMField() const { return mpECMField; }

    void SetConfinementStiffness(double k)
    {
        mConfinementStiffness = k;
    }

    double GetConfinementStiffness() const { return mConfinementStiffness; }

    void SetDegradationEnabled(bool enabled) { mDegradationEnabled = enabled; }
    void SetRemodelingEnabled(bool enabled) { mRemodelingEnabled = enabled; }
    void SetTrackCenter(bool track) { mTrackCenter = track; }

    void SetCenter(c_vector<double, DIM> center)
    {
        mCenter = center;
        mTrackCenter = false;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ConfinementStiffness>" << mConfinementStiffness << "</ConfinementStiffness>\n";
        *rParamsFile << "\t\t\t<DegradationEnabled>" << mDegradationEnabled << "</DegradationEnabled>\n";
        *rParamsFile << "\t\t\t<RemodelingEnabled>" << mRemodelingEnabled << "</RemodelingEnabled>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ECMConfinementForce)

#endif /* ECMCONFINEMENTFORCE_HPP_ */
