/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

*/

#ifndef ECMCONFINEMENTFORCE3D_HPP_
#define ECMCONFINEMENTFORCE3D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "DynamicECMField3d.hpp"
#include "SimulationTime.hpp"
#include "SimProfiler.hpp"
#include <boost/shared_ptr.hpp>

/**
 * 3D ECM Confinement Force
 *
 * Three-dimensional analogue of ECMConfinementForce that uses DynamicECMField3d.
 * Provides fiber-density-based confinement that naturally supports budding
 * through local ECM degradation.
 *
 * Physics:
 *   F_i = -k × ρ(x_i) × r̂_i
 *
 *   where k is the confinement stiffness, ρ ∈ [0,1] is the local ECM density,
 *   and r̂_i points radially outward from the tissue centroid.
 *
 *   Cells degrade the ECM at their position (MMP secretion) and optionally
 *   remodel fiber orientation through traction forces.
 */
class ECMConfinementForce3d : public AbstractForce<3>
{
    friend class boost::serialization::access;

private:

    /** Shared pointer to the 3D ECM field */
    boost::shared_ptr<DynamicECMField3d> mpECMField;

    /** Confinement stiffness (force per unit density) */
    double mConfinementStiffness;

    /** Whether cells degrade ECM at their position */
    bool mDegradationEnabled;

    /** Whether cell traction remodels ECM fibers */
    bool mRemodelingEnabled;

    /** Whether to auto-track centroid */
    bool mTrackCenter;

    /** The tissue centroid (auto-updated or manually set) */
    c_vector<double, 3> mCenter;

    /** Last update time (for computing dt) */
    double mLastUpdateTime;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<3>>(*this);
        archive & mConfinementStiffness;
        archive & mDegradationEnabled;
        archive & mRemodelingEnabled;
        archive & mTrackCenter;
        archive & mCenter;
    }

public:

    ECMConfinementForce3d()
        : AbstractForce<3>(),
          mConfinementStiffness(5.0),
          mDegradationEnabled(true),
          mRemodelingEnabled(true),
          mTrackCenter(true),
          mLastUpdateTime(0.0)
    {
        mCenter[0] = 0.0;
        mCenter[1] = 0.0;
        mCenter[2] = 0.0;
    }

    virtual ~ECMConfinementForce3d() {}

    /**
     * Apply 3D ECM confinement and degradation.
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
    {
        ScopedTimer _prof("ECMConfinement3d");

        if (!mpECMField)
        {
            EXCEPTION("ECM field not set! Call SetECMField() before simulation.");
        }

        // Update centroid
        if (mTrackCenter)
        {
            c_vector<double, 3> center = zero_vector<double>(3);
            unsigned count = 0;
            for (typename AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
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
        if (dt < 1e-10) dt = 0.005;
        mLastUpdateTime = current_time;

        for (typename AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, 3> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, 3> disp = pos - mCenter;
            double r = norm_2(disp);

            if (r < 1e-6) continue;

            c_vector<double, 3> r_hat = disp / r;

            // Sample ECM density at cell position
            double density = mpECMField->GetDensityAt(pos);

            // Store for visualization
            cell_iter->GetCellData()->SetItem("ecm_density", density);

            // Confinement force: pushes cell inward proportional to local ECM density
            c_vector<double, 3> force = -mConfinementStiffness * density * r_hat;
            rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(force);

            // Cell degrades ECM at its position (MMP secretion)
            if (mDegradationEnabled)
            {
                mpECMField->DegradeECM(pos, dt);
            }

            // Cell traction remodels ECM fiber orientation
            if (mRemodelingEnabled)
            {
                mpECMField->ApplyCellTraction(pos, r_hat, dt);
            }
        }

        // Diffuse ECM field (cross-linking / smoothing)
        {
            ScopedTimer _diffProf("ECM3d::DiffuseECM");
            mpECMField->DiffuseECM(dt);
        }
    }

    // ---- Setters ----

    void SetECMField(boost::shared_ptr<DynamicECMField3d> pField)
    {
        mpECMField = pField;
    }

    boost::shared_ptr<DynamicECMField3d> GetECMField() const { return mpECMField; }

    void SetConfinementStiffness(double k)
    {
        mConfinementStiffness = k;
    }

    double GetConfinementStiffness() const { return mConfinementStiffness; }

    void SetDegradationEnabled(bool enabled) { mDegradationEnabled = enabled; }
    void SetRemodelingEnabled(bool enabled) { mRemodelingEnabled = enabled; }
    void SetTrackCenter(bool track) { mTrackCenter = track; }

    void SetCenter(c_vector<double, 3> center)
    {
        mCenter = center;
        mTrackCenter = false;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ConfinementStiffness>" << mConfinementStiffness << "</ConfinementStiffness>\n";
        *rParamsFile << "\t\t\t<DegradationEnabled>" << mDegradationEnabled << "</DegradationEnabled>\n";
        *rParamsFile << "\t\t\t<RemodelingEnabled>" << mRemodelingEnabled << "</RemodelingEnabled>\n";
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};

#endif /* ECMCONFINEMENTFORCE3D_HPP_ */
