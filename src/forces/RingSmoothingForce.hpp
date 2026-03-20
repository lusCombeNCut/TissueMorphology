/*
 * RingSmoothingForce.hpp
 *
 * A discrete Laplacian smoothing force for 2D node-based ring simulations.
 *
 * For each cell, computes the midpoint of its two ring-adjacent neighbors
 * (from RingTopologyTracker) and applies a restoring force toward that
 * midpoint. This is discrete mean curvature flow: it smooths the ring
 * without pushing cells radially inward/outward.
 *
 * The force on cell i is:
 *   F_i = κ * (midpoint(left, right) - position_i)
 *
 * where κ is the smoothing stiffness. This:
 *   - Pulls kinked cells back in line with neighbors
 *   - Does NOT change the overall ring radius
 *   - Works correctly during budding (local smoothing, not global)
 *
 * Requires a RingTopologyTracker to define neighbors.
 */
#ifndef RINGSMOOTHINGFORCE_HPP_
#define RINGSMOOTHINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "RingTopologyTracker.hpp"
#include "StemCellProliferativeType.hpp"
#include "SimProfiler.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
class RingSmoothingForce : public AbstractForce<DIM>
{
private:

    /** Smoothing stiffness — how strongly cells are pulled toward neighbor midpoint */
    double mSmoothingStiffness;

    /** Per-cell-type stiffness multipliers (stem cells bend more, Paneth less) */
    double mStemScale;
    double mTransitScale;
    double mDiffScale;

    /** Pointer to the ring topology tracker (not owned, not serialized) */
    RingTopologyTracker<DIM>* mpRingTopology;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM>>(*this);
        archive & mSmoothingStiffness;
        archive & mStemScale;
        archive & mTransitScale;
        archive & mDiffScale;
        // mpRingTopology not serialized — must be re-attached after checkpoint
    }

public:

    RingSmoothingForce()
        : AbstractForce<DIM>(),
          mSmoothingStiffness(5.0),
          mStemScale(0.7),
          mTransitScale(1.0),
          mDiffScale(1.3),
          mpRingTopology(nullptr)
    {
    }

    virtual ~RingSmoothingForce() {}

    void SetSmoothingStiffness(double stiffness)
    {
        mSmoothingStiffness = stiffness;
    }

    double GetSmoothingStiffness() const
    {
        return mSmoothingStiffness;
    }

    /**
     * Set per-cell-type scaling factors for bending stiffness.
     * @param stemScale   multiplier for stem cells (< 1 = more compliant)
     * @param transitScale multiplier for transit-amplifying cells (1.0 = baseline)
     * @param diffScale   multiplier for differentiated/Paneth cells (> 1 = stiffer)
     */
    void SetCellTypeScaling(double stemScale, double transitScale, double diffScale)
    {
        mStemScale = stemScale;
        mTransitScale = transitScale;
        mDiffScale = diffScale;
    }

    void SetRingTopologyTracker(RingTopologyTracker<DIM>* pTracker)
    {
        mpRingTopology = pTracker;
    }

    /**
     * For each cell, apply a force toward the midpoint of its two ring neighbors.
     *
     * F_i = κ * (midpoint - cell_i)
     *
     * This is the discrete Laplacian of position on the ring, and drives
     * the curve toward uniform spacing (smoothness).
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override
    {
        ScopedTimer _prof("RingSmoothing");
        if (mpRingTopology == nullptr)
        {
            return;  // No topology — skip
        }

        // Build position map
        std::map<unsigned, c_vector<double, DIM>> positions;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            positions[idx] = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }

        // Apply smoothing force to each cell
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            if (!mpRingTopology->HasNeighbors(idx))
            {
                continue;
            }

            auto [left_idx, right_idx] = mpRingTopology->GetNeighbors(idx);

            if (left_idx == UINT_MAX || right_idx == UINT_MAX)
            {
                continue;
            }

            if (positions.find(left_idx) == positions.end() ||
                positions.find(right_idx) == positions.end())
            {
                continue;
            }

            c_vector<double, DIM> cell_pos = positions[idx];
            c_vector<double, DIM> left_pos = positions[left_idx];
            c_vector<double, DIM> right_pos = positions[right_idx];

            // Midpoint of the two neighbors
            c_vector<double, DIM> midpoint = 0.5 * (left_pos + right_pos);

            // Discrete Laplacian: displacement from cell to neighbor midpoint
            c_vector<double, DIM> displacement = midpoint - cell_pos;

            // Scale stiffness by cell type
            double scale = mTransitScale;  // default
            if (cell_iter->template HasCellProperty<StemCellProliferativeType>())
            {
                scale = mStemScale;
            }
            else if (cell_iter->template HasCellProperty<DifferentiatedCellProliferativeType>())
            {
                scale = mDiffScale;
            }

            // Force = stiffness * scale * displacement
            c_vector<double, DIM> force = mSmoothingStiffness * scale * displacement;

            rCellPopulation.GetNode(idx)->AddAppliedForceContribution(force);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<SmoothingStiffness>" << mSmoothingStiffness << "</SmoothingStiffness>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RingSmoothingForce)

#endif // RINGSMOOTHINGFORCE_HPP_
