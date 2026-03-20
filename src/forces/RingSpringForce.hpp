/*
 * RingSpringForce.hpp
 *
 * A spring force for 2D node-based monolayer simulations that ONLY applies
 * forces between topologically adjacent cells on the ring, as tracked by
 * a RingTopologyTracker.
 *
 * The standard GeneralisedLinearSpringForce applies springs between ALL
 * node pairs within the spatial interaction cutoff. On a ring this means
 * cells 2–3 positions apart can attract each other, crumpling the monolayer.
 *
 * This force instead iterates over each cell's two ring neighbors (left
 * and right) and applies the same spring law as GeneralisedLinearSpringForce.
 * Non-adjacent interactions are completely eliminated.
 *
 * Usage:
 *   - Requires a RingTopologyTracker to be added as a SimulationModifier
 *   - Set the tracker pointer via SetRingTopologyTracker() before simulation
 *   - Replaces GeneralisedLinearSpringForce for node2d runs
 */
#ifndef RINGSPRINGFORCE_HPP_
#define RINGSPRINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <cmath>

#include "GeneralisedLinearSpringForce.hpp"
#include "RingTopologyTracker.hpp"
#include "SimProfiler.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned DIM>
class RingSpringForce : public GeneralisedLinearSpringForce<DIM, DIM>
{
private:

    /** Pointer to the ring topology tracker (not owned, not serialized) */
    RingTopologyTracker<DIM>* mpRingTopology;

    /** Stiffness scale factor for Transit (TA) cells (multiplier relative to stem baseline) */
    double mTAStiffnessScale;

    /** Stiffness scale factor for Differentiated cells (multiplier relative to stem baseline) */
    double mDiffStiffnessScale;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<DIM, DIM>>(*this);
        // mpRingTopology is not serialized — must be re-attached after checkpoint
    }

public:

    RingSpringForce()
        : GeneralisedLinearSpringForce<DIM, DIM>(),
          mpRingTopology(nullptr),
          mTAStiffnessScale(1.0),
          mDiffStiffnessScale(1.0)
    {
    }

    virtual ~RingSpringForce() {}

    /**
     * Set the ring topology tracker. Must be called before simulation starts.
     * @param pTracker pointer to the RingTopologyTracker (not owned)
     */
    void SetRingTopologyTracker(RingTopologyTracker<DIM>* pTracker)
    {
        mpRingTopology = pTracker;
    }

    /**
     * Set the stiffness scale factor for Transit (TA) cells.
     * Stem cells are the baseline (scale = 1.0).
     * @param scale multiplicative factor applied to TA–TA or TA–other interactions
     */
    void SetTAStiffnessScale(double scale)
    {
        mTAStiffnessScale = scale;
    }

    /**
     * Set the stiffness scale factor for Differentiated cells.
     * Stem cells are the baseline (scale = 1.0).
     * @param scale multiplicative factor applied to Diff–Diff or Diff–other interactions
     */
    void SetDiffStiffnessScale(double scale)
    {
        mDiffStiffnessScale = scale;
    }

    /**
     * Override to return a per-cell-type stiffness multiplier.
     *
     * For each cell in the pair we look up its proliferative type:
     *   Stem → 1.0,  Transit → mTAStiffnessScale,  Differentiated → mDiffStiffnessScale.
     *
     * The pairwise multiplier is the geometric mean: sqrt(scaleA * scaleB).
     * This gives a symmetric, monotone blend:
     *   stem–stem  → 1,
     *   stem–TA    → sqrt(TAscale),
     *   TA–TA      → TAscale,
     *   TA–diff    → sqrt(TAscale * Diffscale),
     *   diff–diff  → Diffscale.
     */
    double VariableSpringConstantMultiplicationFactor(
        unsigned nodeAGlobalIndex,
        unsigned nodeBGlobalIndex,
        AbstractCellPopulation<DIM, DIM>& rCellPopulation,
        bool isCloserThanRestLength) override
    {
        // If both scales are 1.0, skip the lookup entirely
        if (mTAStiffnessScale == 1.0 && mDiffStiffnessScale == 1.0)
        {
            return 1.0;
        }

        auto GetScaleForNode = [&](unsigned nodeIndex) -> double
        {
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(nodeIndex);
            if (p_cell->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
            {
                return mTAStiffnessScale;
            }
            else if (p_cell->GetCellProliferativeType()->template IsType<DifferentiatedCellProliferativeType>())
            {
                return mDiffStiffnessScale;
            }
            return 1.0; // Stem or any other type
        };

        double scale_a = GetScaleForNode(nodeAGlobalIndex);
        double scale_b = GetScaleForNode(nodeBGlobalIndex);

        return std::sqrt(scale_a * scale_b);
    }

    /**
     * Override AddForceContribution to iterate ONLY over ring-adjacent pairs
     * instead of all spatial pairs.
     *
     * For each cell, we get its two neighbors from the RingTopologyTracker
     * and apply the GeneralisedLinearSpringForce spring law (via
     * CalculateForceBetweenNodes) to each pair. We process each pair once
     * by only computing when nodeA < nodeB.
     */
    void AddForceContribution(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer _prof("RingSpring");
        if (mpRingTopology == nullptr)
        {
            // No topology set — fall back to standard behaviour
            GeneralisedLinearSpringForce<DIM, DIM>::AddForceContribution(rCellPopulation);
            return;
        }

        // Iterate over all cells and apply forces only to ring neighbors
        for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_a = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            if (!mpRingTopology->HasNeighbors(node_a))
            {
                continue;
            }

            auto [left_idx, right_idx] = mpRingTopology->GetNeighbors(node_a);

            // Process each neighbor — only compute when node_a < neighbor to avoid doubles
            for (unsigned node_b : {left_idx, right_idx})
            {
                if (node_b == UINT_MAX)
                {
                    continue;
                }

                // Only compute each pair once (lower index drives the calculation)
                if (node_a >= node_b)
                {
                    continue;
                }

                // Use the parent class spring law to compute the force
                c_vector<double, DIM> force = this->CalculateForceBetweenNodes(
                    node_a, node_b, rCellPopulation);

                rCellPopulation.GetNode(node_a)->AddAppliedForceContribution(force);
                rCellPopulation.GetNode(node_b)->AddAppliedForceContribution(-force);
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<UsingRingTopology>"
                     << (mpRingTopology != nullptr ? "true" : "false")
                     << "</UsingRingTopology>\n";
        *rParamsFile << "\t\t\t<TAStiffnessScale>"
                     << mTAStiffnessScale
                     << "</TAStiffnessScale>\n";
        *rParamsFile << "\t\t\t<DiffStiffnessScale>"
                     << mDiffStiffnessScale
                     << "</DiffStiffnessScale>\n";
        GeneralisedLinearSpringForce<DIM, DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RingSpringForce)

#endif // RINGSPRINGFORCE_HPP_
