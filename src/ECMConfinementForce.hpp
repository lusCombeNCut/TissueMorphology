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
 * ECM Confinement Force — agent-based viscoelastic model
 *
 * Models cell–ECM interaction by treating each ECM grid element as an
 * individual agent connected to nearby cells via linear springs. This
 * replaces the simple radial confinement force and provides:
 *
 *   - Tangential adhesion (prevents whole-organoid rotation/translation)
 *   - Density-dependent force (degraded ECM exerts weaker force)
 *   - Purely linear force law (cells can push through ECM without the
 *     exponential blow-up used in the standard cell–cell springs)
 *
 * Physics:
 *   For each cell i interacting with ECM element j within the cutoff:
 *
 *     F_ij = k × ρ_j × (d_ij − s₀) × d̂_ij
 *
 *   where:
 *     k      = spring stiffness (mConfinementStiffness)
 *     ρ_j    = ECM density at element j ∈ [0,1] — acts as a force multiplier
 *     d_ij   = distance between cell and ECM element
 *     s₀     = spring rest length (mEcmSpringRestLength)
 *     d̂_ij   = unit vector from cell to ECM element
 *
 *   When d > s₀: cell is pulled toward ECM element (adhesion)
 *   When d < s₀: cell is pushed away from ECM element (repulsion)
 *   Force is purely linear — no log/exp corrections.
 *
 *   Cells also degrade ECM at their position (MMP secretion) and
 *   optionally remodel ECM fiber orientation via traction.
 */
template<unsigned DIM>
class ECMConfinementForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Shared pointer to the ECM field */
    boost::shared_ptr<DynamicECMField> mpECMField;

    /** Spring stiffness for cell–ECM interaction (force per unit density per unit overlap) */
    double mConfinementStiffness;

    /** Rest length of the cell–ECM spring */
    double mEcmSpringRestLength;

    /** Cutoff distance for cell–ECM spring interactions */
    double mEcmInteractionCutoff;

    /** Whether cells degrade ECM at their position */
    bool mDegradationEnabled;

    /** Whether cell traction remodels ECM fibers */
    bool mRemodelingEnabled;

    /** Whether to auto-track centroid (only used for degradation/remodeling direction) */
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
        archive & mEcmSpringRestLength;
        archive & mEcmInteractionCutoff;
        archive & mDegradationEnabled;
        archive & mRemodelingEnabled;
        archive & mTrackCenter;
        archive & mCenter;
    }

public:

    ECMConfinementForce()
        : AbstractForce<DIM>(),
          mConfinementStiffness(5.0),
          mEcmSpringRestLength(0.0),
          mEcmInteractionCutoff(1.5),
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
     * Apply cell–ECM spring forces, degradation, and remodeling.
     *
     * For each cell:
     *   1. Find nearby ECM elements within the interaction cutoff
     *   2. Apply a linear spring force for each element:
     *        F = k × ρ × (d − s₀) × d̂   (toward ECM element)
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

        // Update centroid (used for degradation/remodeling direction)
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

        // Temporary vectors for nearby element queries (reused per cell)
        std::vector<c_vector<double, 2>> nearby_positions;
        std::vector<double> nearby_densities;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // ---- Agent-based cell–ECM spring forces ----
            // Find nearby ECM elements and apply linear springs
            mpECMField->GetNearbyElements(pos, mEcmInteractionCutoff,
                                          nearby_positions, nearby_densities);

            c_vector<double, DIM> total_ecm_force = zero_vector<double>(DIM);
            double max_density = 0.0;

            for (unsigned e = 0; e < nearby_positions.size(); e++)
            {
                c_vector<double, DIM> disp = nearby_positions[e] - pos;  // cell → ECM element
                double distance = norm_2(disp);

                if (distance < 1e-10) continue;  // skip ECM elements exactly at cell position

                c_vector<double, DIM> unit_disp = disp / distance;
                double density = nearby_densities[e];

                // Purely linear spring: F = k × ρ × (d − s₀) × d̂
                // d > s₀ → force toward ECM element (adhesion/attraction)
                // d < s₀ → force away from ECM element (repulsion)
                double overlap = distance - mEcmSpringRestLength;
                c_vector<double, DIM> spring_force = mConfinementStiffness * density * overlap * unit_disp;

                total_ecm_force += spring_force;

                if (density > max_density) max_density = density;
            }

            rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(total_ecm_force);

            // Store max local ECM density for visualization
            cell_iter->GetCellData()->SetItem("ecm_density", max_density);

            // ---- ECM degradation (MMP secretion at cell position) ----
            if (mDegradationEnabled)
            {
                mpECMField->DegradeECM(pos, dt);
            }

            // ---- ECM fiber remodeling via cell traction ----
            if (mRemodelingEnabled)
            {
                c_vector<double, DIM> disp_from_center = pos - mCenter;
                double r = norm_2(disp_from_center);
                if (r > 1e-6)
                {
                    c_vector<double, 2> traction;
                    traction[0] = disp_from_center[0] / r;
                    traction[1] = disp_from_center[1] / r;
                    mpECMField->ApplyCellTraction(pos, traction, dt);
                }
            }
        }

        // Diffuse ECM field (cross-linking / smoothing) — skip if diffusion coefficient is 0
        if (mpECMField->GetDiffusionCoeff() > 0.0)
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

    void SetEcmSpringRestLength(double restLength)
    {
        mEcmSpringRestLength = restLength;
    }

    double GetEcmSpringRestLength() const { return mEcmSpringRestLength; }

    void SetEcmInteractionCutoff(double cutoff)
    {
        mEcmInteractionCutoff = cutoff;
    }

    double GetEcmInteractionCutoff() const { return mEcmInteractionCutoff; }

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
        *rParamsFile << "\t\t\t<EcmSpringRestLength>" << mEcmSpringRestLength << "</EcmSpringRestLength>\n";
        *rParamsFile << "\t\t\t<EcmInteractionCutoff>" << mEcmInteractionCutoff << "</EcmInteractionCutoff>\n";
        *rParamsFile << "\t\t\t<DegradationEnabled>" << mDegradationEnabled << "</DegradationEnabled>\n";
        *rParamsFile << "\t\t\t<RemodelingEnabled>" << mRemodelingEnabled << "</RemodelingEnabled>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ECMConfinementForce)

#endif /* ECMCONFINEMENTFORCE_HPP_ */
