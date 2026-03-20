/*
 * SurfaceSpringForce.hpp
 *
 * A spring force for 3D node-based monolayer simulations that ONLY applies
 * forces between topologically adjacent cells, as tracked by a
 * SurfaceTopologyTracker.
 *
 * The standard GeneralisedLinearSpringForce applies springs between ALL
 * node pairs within the spatial interaction cutoff.  On a thin shell this
 * causes cells on opposite sides of a fold to attract each other, breaking
 * the monolayer.
 *
 * This force instead iterates over the topological adjacency graph and
 * applies the same spring law.  Non-adjacent interactions are eliminated.
 *
 * Usage:
 *   - Requires a SurfaceTopologyTracker to be added as a SimulationModifier
 *   - Set the tracker pointer via SetSurfaceTopologyTracker() before
 *     the simulation starts
 *   - Replaces GeneralisedLinearSpringForce for node3d runs
 */
#ifndef SURFACESPRINGFORCE_HPP_
#define SURFACESPRINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "GeneralisedLinearSpringForce.hpp"
#include "SurfaceTopologyTracker.hpp"
#include "SimProfiler.hpp"

template<unsigned DIM>
class SurfaceSpringForce : public GeneralisedLinearSpringForce<DIM, DIM>
{
private:

    /** Pointer to the surface topology tracker (not owned) */
    SurfaceTopologyTracker<DIM>* mpSurfaceTopology;

    /** Maximum force magnitude per pair (prevents blowup on division) */
    double mMaxForceMagnitude;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<
            GeneralisedLinearSpringForce<DIM, DIM>>(*this);
        // mpSurfaceTopology not serialized — re-attach after checkpoint
    }

public:

    SurfaceSpringForce()
        : GeneralisedLinearSpringForce<DIM, DIM>(),
          mpSurfaceTopology(nullptr),
          mMaxForceMagnitude(100.0)
    {
    }

    virtual ~SurfaceSpringForce() {}

    /**
     * Set the surface topology tracker.  Must be called before Solve().
     * @param pTracker  raw pointer (not owned by this class)
     */
    void SetSurfaceTopologyTracker(SurfaceTopologyTracker<DIM>* pTracker)
    {
        mpSurfaceTopology = pTracker;
    }

    /**
     * Override AddForceContribution to iterate ONLY over topologically
     * adjacent pairs instead of all spatial pairs.
     */
    void AddForceContribution(
        AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer timer("SurfaceSpring");
        if (mpSurfaceTopology == nullptr)
        {
            // No topology set — fall back to standard spatial behavior
            GeneralisedLinearSpringForce<DIM, DIM>::AddForceContribution(
                rCellPopulation);
            return;
        }

        // Iterate over every cell
        for (auto cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned node_a =
                rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            if (!mpSurfaceTopology->HasCell(node_a))
                continue;

            const std::set<unsigned>& neighbors =
                mpSurfaceTopology->GetNeighbors(node_a);

            for (unsigned node_b : neighbors)
            {
                // Only compute each pair once (lower index drives)
                if (node_a >= node_b)
                    continue;

                // Use the parent-class spring law
                c_vector<double, DIM> force =
                    this->CalculateForceBetweenNodes(
                        node_a, node_b, rCellPopulation);

                // Cap per-pair force to prevent blowups on cell division
                double fmag = norm_2(force);
                if (fmag > mMaxForceMagnitude)
                {
                    force *= mMaxForceMagnitude / fmag;
                }

                rCellPopulation.GetNode(node_a)
                    ->AddAppliedForceContribution(force);
                rCellPopulation.GetNode(node_b)
                    ->AddAppliedForceContribution(-force);
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<UsingSurfaceTopology>"
                     << (mpSurfaceTopology != nullptr ? "true" : "false")
                     << "</UsingSurfaceTopology>\n";
        GeneralisedLinearSpringForce<DIM, DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceSpringForce)

#endif // SURFACESPRINGFORCE_HPP_
