/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

*/

#ifndef CELLPOLARITYFORCE_HPP_
#define CELLPOLARITYFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SimulationTime.hpp"

/**
 * A force class implementing cell polarity for maintaining epithelial monolayers.
 *
 * Based on the ya||a framework, polarity is represented as a unit vector p for each
 * cell, specified in spherical coordinates (theta, phi). The force drives cells
 * to form a single-layered structure where:
 *
 * 1. **Epithelial bending force** (U_Epi): Minimized when cell polarity is orthogonal
 *    to cell-cell connections, enforcing monolayer geometry.
 *
 * 2. **Tissue polarity force** (U_Pol): Minimized when neighboring polarities are
 *    parallel, providing tissue-level coherence.
 *
 * The polarity angles (theta, phi) are stored in cell data and evolved according
 * to forces derived from these potentials. The spatial force contribution ensures
 * cells physically arrange to satisfy polarity constraints.
 *
 * References:
 * - Marin-Riera et al. (2016) "Computational modeling of development..."
 * - ya||a framework: https://github.com/germannp/yalla
 */
template<unsigned DIM>
class CellPolarityForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Strength of epithelial bending force (monolayer enforcement) */
    double mBendingStrength;

    /** Strength of tissue polarity alignment force */
    double mPolarityAlignmentStrength;

    /** Maximum interaction distance for polarity forces */
    double mInteractionCutoff;

    /** Whether to initialize polarities pointing radially outward */
    bool mInitializeRadially;

    /** Small value to avoid division by zero */
    static constexpr double EPSILON = 1e-10;

    /**
     * Boost Serialization method for archiving/checkpointing.
     *
     * @param archive  the archive
     * @param version  the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mBendingStrength;
        archive & mPolarityAlignmentStrength;
        archive & mInteractionCutoff;
        archive & mInitializeRadially;
    }

    /**
     * Convert polarity angles (theta, phi) to a Cartesian unit vector.
     */
    c_vector<double, DIM> PolarityToVector(double theta, [[maybe_unused]] double phi) const
    {
        c_vector<double, DIM> p;
        if constexpr (DIM == 3)
        {
            p[0] = sin(theta) * cos(phi);
            p[1] = sin(theta) * sin(phi);
            p[2] = cos(theta);
        }
        else // DIM == 2
        {
            // In 2D, theta represents angle from x-axis (phi unused)
            p[0] = cos(theta);
            p[1] = sin(theta);
        }
        return p;
    }

    /**
     * Compute dot product between polarity vector p and unit direction r_hat.
     */
    double PolarityDotProduct(double theta_i, double phi_i,
                              double theta_r, double phi_r) const
    {
        if constexpr (DIM == 3)
        {
            return sin(theta_i) * sin(theta_r) * cos(phi_i - phi_r) +
                   cos(theta_i) * cos(theta_r);
        }
        else
        {
            return cos(theta_i - theta_r);
        }
    }

    /**
     * Compute the unidirectional polarization gradient (aligning p_i toward direction p).
     * Returns (dF_theta, dF_phi) for polarity evolution.
     */
    std::pair<double, double> UnidirectionalPolarizationGradient(
        double theta_i, double phi_i, double theta_r, double phi_r) const
    {
        double dF_theta = 0.0;
        double dF_phi = 0.0;

        if constexpr (DIM == 3)
        {
            dF_theta = cos(theta_i) * sin(theta_r) * cos(phi_i - phi_r) -
                       sin(theta_i) * cos(theta_r);

            double sin_theta_i = sin(theta_i);
            if (fabs(sin_theta_i) > EPSILON)
            {
                dF_phi = -sin(theta_r) * sin(phi_i - phi_r) / sin_theta_i;
            }
        }
        else
        {
            dF_theta = -sin(theta_i - theta_r);
        }

        return {dF_theta, dF_phi};
    }

public:

    /**
     * Constructor.
     */
    CellPolarityForce()
        : AbstractForce<DIM>(),
          mBendingStrength(0.3),
          mPolarityAlignmentStrength(0.1),
          mInteractionCutoff(3.0),
          mInitializeRadially(true)
    {
    }

    /**
     * Destructor.
     */
    virtual ~CellPolarityForce()
    {
    }

    /**
     * Initialize cell polarities if not already set.
     * If mInitializeRadially is true, polarities point radially outward from centroid.
     */
    void InitializePolarities(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            // Check if polarity already initialized
            bool has_theta = false;
            try
            {
                cell_iter->GetCellData()->GetItem("polarity_theta");
                has_theta = true;
            }
            catch (Exception&) {}

            if (!has_theta)
            {
                c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                c_vector<double, DIM> r = cell_location - centroid;
                double dist = norm_2(r);

                double theta = 0.0;
                double phi = 0.0;

                if (dist > EPSILON)
                {
                    if constexpr (DIM == 3)
                    {
                        // Spherical coords: theta = acos(z/r), phi = atan2(y, x)
                        theta = acos(r[2] / dist);
                        phi = atan2(r[1], r[0]);
                    }
                    else
                    {
                        theta = atan2(r[1], r[0]);
                    }
                }

                cell_iter->GetCellData()->SetItem("polarity_theta", theta);
                cell_iter->GetCellData()->SetItem("polarity_phi", phi);
            }
        }
    }

    /**
     * Add force contribution from cell polarity.
     *
     * Implements:
     * 1. Bending force from U_Epi = Σ(p_i · r_ij/r)^2/2
     *    - Keeps polarity orthogonal to cell-cell vectors (monolayer)
     * 2. Polarity alignment from U_Pol = -Σ(p_i · p_j)^2/2
     *    - Aligns neighboring polarities
     *
     * Also evolves polarity angles (theta, phi) stored in cell data.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        // Ensure polarities are initialized
        InitializePolarities(rCellPopulation);

        // Only works with node-based populations
        NodeBasedCellPopulation<DIM>* p_node_pop =
            dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation);

        if (!p_node_pop)
        {
            EXCEPTION("CellPolarityForce only works with NodeBasedCellPopulation");
        }

        // Get time step for polarity evolution
        double dt = SimulationTime::Instance()->GetTimeStep();

        // Store polarity updates to apply after iteration
        std::map<unsigned, std::pair<double, double>> polarity_updates;

        // Iterate over all cells
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index_i = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> location_i = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // Get polarity of cell i
            double theta_i = cell_iter->GetCellData()->GetItem("polarity_theta");
            double phi_i = cell_iter->GetCellData()->GetItem("polarity_phi");
            c_vector<double, DIM> p_i = PolarityToVector(theta_i, phi_i);

            // Accumulate forces and polarity derivatives
            c_vector<double, DIM> total_force = zero_vector<double>(DIM);
            double dTheta_i = 0.0;
            double dPhi_i = 0.0;

            // Get neighbors within cutoff
            std::set<unsigned> neighbours = p_node_pop->GetNodesWithinNeighbourhoodRadius(node_index_i, mInteractionCutoff);

            for (unsigned node_index_j : neighbours)
            {
                if (node_index_j == node_index_i) continue;

                CellPtr p_cell_j = rCellPopulation.GetCellUsingLocationIndex(node_index_j);
                c_vector<double, DIM> location_j = rCellPopulation.GetLocationOfCellCentre(p_cell_j);
                c_vector<double, DIM> r_ij = location_i - location_j;
                double dist = norm_2(r_ij);

                if (dist < EPSILON || dist > mInteractionCutoff) continue;

                // Get polarity of cell j
                double theta_j = p_cell_j->GetCellData()->GetItem("polarity_theta");
                double phi_j = p_cell_j->GetCellData()->GetItem("polarity_phi");
                c_vector<double, DIM> p_j = PolarityToVector(theta_j, phi_j);

                // Unit vector from j to i
                c_vector<double, DIM> r_hat = r_ij / dist;
                double theta_r, phi_r;
                if constexpr (DIM == 3)
                {
                    theta_r = acos(r_hat[2]);
                    phi_r = atan2(r_hat[1], r_hat[0]);
                }
                else
                {
                    theta_r = atan2(r_hat[1], r_hat[0]);
                    phi_r = 0.0;
                }

                // ================================================================
                // 1. EPITHELIAL BENDING FORCE (monolayer enforcement)
                //    U_Epi = Σ(p_i · r_ij/r)^2/2
                //    Force: ∆F = -(p_i · r_hat) * [p_i/r - (p_i · r_hat) * r_hat/r]
                // ================================================================
                double prod_i = inner_prod(p_i, r_hat);
                c_vector<double, DIM> bending_force_i = 
                    -prod_i / dist * p_i + prod_i * prod_i / (dist * dist) * r_ij;
                
                // Contribution from (p_j · r_ji)^2/2
                double prod_j = inner_prod(p_j, -r_hat);
                c_vector<double, DIM> bending_force_j = 
                    -prod_j / dist * p_j + prod_j * prod_j / (dist * dist) * r_ij;

                total_force += mBendingStrength * (bending_force_i + bending_force_j);

                // Polarity evolution from bending potential
                auto [dTheta_bend, dPhi_bend] = UnidirectionalPolarizationGradient(
                    theta_i, phi_i, theta_r, phi_r);
                dTheta_i -= mBendingStrength * prod_i * dTheta_bend;
                dPhi_i -= mBendingStrength * prod_i * dPhi_bend;

                // ================================================================
                // 2. TISSUE POLARITY ALIGNMENT
                //    U_Pol = -Σ(p_i · p_j)^2/2  (minimize when parallel)
                //    Bidirectional alignment
                // ================================================================
                double prod_pol = PolarityDotProduct(theta_i, phi_i, theta_j, phi_j);
                auto [dTheta_pol, dPhi_pol] = UnidirectionalPolarizationGradient(
                    theta_i, phi_i, theta_j, phi_j);
                
                // Gradient of (p_i · p_j)^2/2 gives prod_pol * gradient
                dTheta_i += mPolarityAlignmentStrength * prod_pol * dTheta_pol;
                dPhi_i += mPolarityAlignmentStrength * prod_pol * dPhi_pol;
            }

            // Apply spatial force
            rCellPopulation.GetNode(node_index_i)->AddAppliedForceContribution(total_force);

            // Store polarity update
            polarity_updates[node_index_i] = {theta_i + dTheta_i * dt, phi_i + dPhi_i * dt};
        }

        // Apply polarity updates
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            if (polarity_updates.count(node_index))
            {
                auto [new_theta, new_phi] = polarity_updates[node_index];
                
                // Clamp theta to [0, pi] for 3D, or normalize for 2D
                if constexpr (DIM == 3)
                {
                    new_theta = std::max(EPSILON, std::min(M_PI - EPSILON, new_theta));
                    // Wrap phi to [-pi, pi]
                    while (new_phi > M_PI) new_phi -= 2.0 * M_PI;
                    while (new_phi < -M_PI) new_phi += 2.0 * M_PI;
                }
                else
                {
                    // Wrap theta to [-pi, pi] in 2D
                    while (new_theta > M_PI) new_theta -= 2.0 * M_PI;
                    while (new_theta < -M_PI) new_theta += 2.0 * M_PI;
                }

                cell_iter->GetCellData()->SetItem("polarity_theta", new_theta);
                cell_iter->GetCellData()->SetItem("polarity_phi", new_phi);
            }
        }
    }

    /**
     * Overridden OutputForceParameters() method.
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<BendingStrength>" << mBendingStrength << "</BendingStrength>\n";
        *rParamsFile << "\t\t\t<PolarityAlignmentStrength>" << mPolarityAlignmentStrength << "</PolarityAlignmentStrength>\n";
        *rParamsFile << "\t\t\t<InteractionCutoff>" << mInteractionCutoff << "</InteractionCutoff>\n";
        *rParamsFile << "\t\t\t<InitializeRadially>" << mInitializeRadially << "</InitializeRadially>\n";

        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }

    // ================== SETTERS ==================

    void SetBendingStrength(double strength)
    {
        mBendingStrength = strength;
    }

    void SetPolarityAlignmentStrength(double strength)
    {
        mPolarityAlignmentStrength = strength;
    }

    void SetInteractionCutoff(double cutoff)
    {
        mInteractionCutoff = cutoff;
    }

    void SetInitializeRadially(bool initialize)
    {
        mInitializeRadially = initialize;
    }

    // ================== GETTERS ==================

    double GetBendingStrength() const
    {
        return mBendingStrength;
    }

    double GetPolarityAlignmentStrength() const
    {
        return mPolarityAlignmentStrength;
    }

    double GetInteractionCutoff() const
    {
        return mInteractionCutoff;
    }

    bool GetInitializeRadially() const
    {
        return mInitializeRadially;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellPolarityForce)

#endif // CELLPOLARITYFORCE_HPP_
