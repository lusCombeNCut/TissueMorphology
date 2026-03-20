/*
 * ViscoelasticGhostNodeEcmForce.hpp
 *
 * Force class for off-lattice ghost node ECM with true viscoelastic
 * constitutive law (generalized Maxwell / Standard Linear Solid).
 *
 * Handles:
 *   - Cell–ghost node spring forces (with density-weighted attraction)
 *   - Ghost–ghost viscoelastic spring forces (evolving rest lengths)
 *   - ECM degradation (MMP secretion)
 *   - Fibre remodeling (cell traction)
 *   - Periodic ghost node removal
 *
 * Constitutive model (Fertala et al. 2025):
 *   E(t) = E_0 + E_1 * exp(-t/tau)
 *
 * Discrete analogue via per-pair rest length evolution:
 *   F_ij   = (E_0 + E_1) * (d_ij - s_ij)
 *   ds/dt  = (d_ij - s_ij) / tau
 */
#ifndef VISCOELASTICGHOSTNODEECMFORCE_HPP_
#define VISCOELASTICGHOSTNODEECMFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "SimulationTime.hpp"
#include "SimProfiler.hpp"
#include <boost/shared_ptr.hpp>
#include <climits>

template<unsigned DIM>
class ViscoelasticGhostNodeEcmForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Shared pointer to the viscoelastic ghost node ECM field */
    boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> mpGhostField;

    /** Spring stiffness for cell–ghost node interaction */
    double mCellGhostStiffness;

    /** Rest length for cell–ghost node springs */
    double mCellGhostRestLength;

    /** Cutoff distance for cell–ghost node interactions */
    double mCellGhostCutoff;

    /** Whether cells degrade nearby ghost nodes (MMP secretion) */
    bool mDegradationEnabled;

    /** Whether cell traction remodels ghost node fibres */
    bool mRemodelingEnabled;

    /** Whether to auto-track tissue centroid */
    bool mTrackCenter;

    /** Tissue centroid (for radial traction direction) */
    c_vector<double, DIM> mCenter;

    /** Counter for periodic ghost node removal checks */
    unsigned mRemovalCheckInterval;
    unsigned mStepCounter;

    /** Last processed timestep */
    unsigned mLastStepProcessed;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mCellGhostStiffness;
        archive & mCellGhostRestLength;
        archive & mCellGhostCutoff;
        archive & mDegradationEnabled;
        archive & mRemodelingEnabled;
        archive & mTrackCenter;
        archive & mCenter;
    }

public:

    ViscoelasticGhostNodeEcmForce()
        : AbstractForce<DIM>(),
          mCellGhostStiffness(10.0),
          mCellGhostRestLength(0.0),
          mCellGhostCutoff(1.5),
          mDegradationEnabled(true),
          mRemodelingEnabled(true),
          mTrackCenter(true),
          mRemovalCheckInterval(100),
          mStepCounter(0),
          mLastStepProcessed(UINT_MAX)
    {
        for (unsigned i = 0; i < DIM; i++)
            mCenter[i] = 0.0;
    }

    virtual ~ViscoelasticGhostNodeEcmForce() {}

    /**
     * Apply cell–ghost node spring forces and manage viscoelastic ECM dynamics.
     *
     * This method:
     *   1. Computes cell–ghost node spring forces (density-weighted)
     *   2. Applies forces to cell nodes (distributed for vertex populations)
     *   3. Degrades ghost nodes near cells (MMP secretion)
     *   4. Remodels fibre orientation via cell traction
     *   5. Updates ghost–ghost viscoelastic spring forces + evolves rest lengths
     *   6. Periodically removes depleted ghost nodes
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        ScopedTimer _prof("ViscoelasticGhostNodeECM");

        if (!mpGhostField)
        {
            EXCEPTION("Viscoelastic ghost node ECM field not set! Call SetGhostField() before simulation.");
        }

        double dt = SimulationTime::Instance()->GetTimeStep();
        unsigned current_step = SimulationTime::Instance()->GetTimeStepsElapsed();
        bool new_step = (current_step != mLastStepProcessed);
        mLastStepProcessed = current_step;

        // Update tissue centroid
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

        // Clear ghost forces at the start of each timestep
        mpGhostField->ClearGhostForces();

        // Detect vertex-based population types
        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        MonolayerVertexBasedCellPopulation<DIM>* p_monolayer_pop =
            dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        std::vector<unsigned> nearby_ghosts;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> cell_pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // ── Cell–ghost node spring forces ──────────────────
            mpGhostField->GetNearbyGhostNodes(cell_pos, mCellGhostCutoff, nearby_ghosts);

            c_vector<double, DIM> total_force = zero_vector<double>(DIM);
            double max_density = 0.0;

            for (unsigned gn_idx : nearby_ghosts)
            {
                const ViscoelasticGhostNode<DIM>& gn = mpGhostField->GetNode(gn_idx);

                c_vector<double, DIM> disp = gn.position - cell_pos;
                double dist = norm_2(disp);
                if (dist < 1e-10) continue;

                c_vector<double, DIM> unit_disp = disp / dist;

                // Anisotropic modulation
                double aniso_factor = 1.0;
                if (gn.anisotropy > 0.0)
                {
                    double cos_angle = std::abs(inner_prod(unit_disp, gn.fibre_direction));
                    double a = 0.5;
                    aniso_factor = 1.0 - a * gn.anisotropy * (1.0 - cos_angle);
                }

                // Cell–ghost spring force (same regime as original)
                double overlap = dist - mCellGhostRestLength;
                double force_mag;

                if (mCellGhostRestLength > 1e-10)
                {
                    if (overlap <= 0)  // compression: log repulsion
                    {
                        double log_term = mCellGhostRestLength
                                          * std::log(1.0 + overlap / mCellGhostRestLength);
                        force_mag = mCellGhostStiffness * log_term * aniso_factor;
                    }
                    else  // extension: density-weighted attraction
                    {
                        double alpha = 5.0;
                        force_mag = mCellGhostStiffness * gn.density * overlap
                                    * std::exp(-alpha * overlap / mCellGhostRestLength) * aniso_factor;
                    }
                }
                else
                {
                    double alpha = 5.0;
                    force_mag = mCellGhostStiffness * gn.density * overlap
                                * std::exp(-alpha * overlap / mCellGhostCutoff) * aniso_factor;
                }

                c_vector<double, DIM> spring_force = force_mag * unit_disp;
                total_force += spring_force;

                // Newton's 3rd law reaction on ghost node
                mpGhostField->AddForceToGhostNode(gn_idx, -spring_force);

                if (gn.density > max_density) max_density = gn.density;
            }

            // Apply force to cell nodes
            if (p_monolayer_pop)
            {
                auto* p_element = p_monolayer_pop->GetElement(loc_index);
                unsigned n_nodes = p_element->GetNumNodes();
                c_vector<double, DIM> force_per_node = total_force / static_cast<double>(n_nodes);
                for (unsigned i = 0; i < n_nodes; i++)
                    p_element->GetNode(i)->AddAppliedForceContribution(force_per_node);
            }
            else if (p_vertex_pop)
            {
                VertexElement<DIM, DIM>* p_element = p_vertex_pop->rGetMesh().GetElement(loc_index);
                unsigned n_nodes = p_element->GetNumNodes();
                c_vector<double, DIM> force_per_node = total_force / static_cast<double>(n_nodes);
                for (unsigned i = 0; i < n_nodes; i++)
                    p_element->GetNode(i)->AddAppliedForceContribution(force_per_node);
            }
            else
            {
                rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(total_force);
            }
            cell_iter->GetCellData()->SetItem("ecm_density", max_density);

            // ── ECM degradation ────────────────────────────────
            if (new_step && mDegradationEnabled)
            {
                mpGhostField->DegradeNearby(cell_pos, nearby_ghosts, mCellGhostCutoff, dt);
            }

            // ── Fibre remodeling via cell traction ─────────────
            if (new_step && mRemodelingEnabled)
            {
                c_vector<double, DIM> disp_from_center = cell_pos - mCenter;
                double r = norm_2(disp_from_center);
                if (r > 1e-6)
                {
                    c_vector<double, DIM> traction = disp_from_center / r;
                    mpGhostField->RemodelNearby(cell_pos, nearby_ghosts, traction, mCellGhostCutoff, dt);
                }
            }
        }

        // ── Ghost–ghost viscoelastic dynamics ──────────────────
        if (new_step)
        {
            mpGhostField->UpdateGhostNodePositions(dt);
        }

        // ── Periodic removal of depleted ghost nodes ───────────
        if (new_step)
        {
            mStepCounter++;
            if (mStepCounter >= mRemovalCheckInterval)
            {
                mStepCounter = 0;
                unsigned removed = mpGhostField->RemoveDepletedNodes();
                if (removed > 0)
                {
                    std::cout << "  ViscoelasticGhostNodeECM: removed " << removed
                              << " depleted nodes (" << mpGhostField->GetNumActive()
                              << " remaining)" << std::endl;
                }
            }
        }
    }

    // ── Setters ────────────────────────────────────────────────

    void SetGhostField(boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> pField)
    {
        mpGhostField = pField;
    }

    boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> GetGhostField() const { return mpGhostField; }

    void SetCellGhostStiffness(double k) { mCellGhostStiffness = k; }
    double GetCellGhostStiffness() const { return mCellGhostStiffness; }

    void SetCellGhostRestLength(double s0) { mCellGhostRestLength = s0; }
    double GetCellGhostRestLength() const { return mCellGhostRestLength; }

    void SetCellGhostCutoff(double cutoff) { mCellGhostCutoff = cutoff; }
    double GetCellGhostCutoff() const { return mCellGhostCutoff; }

    void SetDegradationEnabled(bool enabled) { mDegradationEnabled = enabled; }
    void SetRemodelingEnabled(bool enabled) { mRemodelingEnabled = enabled; }
    void SetTrackCenter(bool track) { mTrackCenter = track; }
    void SetRemovalCheckInterval(unsigned interval) { mRemovalCheckInterval = interval; }

    void SetCenter(c_vector<double, DIM> center)
    {
        mCenter = center;
        mTrackCenter = false;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<CellGhostStiffness>" << mCellGhostStiffness << "</CellGhostStiffness>\n";
        *rParamsFile << "\t\t\t<CellGhostRestLength>" << mCellGhostRestLength << "</CellGhostRestLength>\n";
        *rParamsFile << "\t\t\t<CellGhostCutoff>" << mCellGhostCutoff << "</CellGhostCutoff>\n";
        *rParamsFile << "\t\t\t<DegradationEnabled>" << mDegradationEnabled << "</DegradationEnabled>\n";
        *rParamsFile << "\t\t\t<RemodelingEnabled>" << mRemodelingEnabled << "</RemodelingEnabled>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ViscoelasticGhostNodeEcmForce)

#endif // VISCOELASTICGHOSTNODEECMFORCE_HPP_
