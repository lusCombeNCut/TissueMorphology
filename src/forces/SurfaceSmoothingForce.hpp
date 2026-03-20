/*
 * SurfaceSmoothingForce.hpp
 *
 * Discrete Laplacian smoothing force for 3D node-based monolayer
 * simulations.
 *
 * For each cell, computes the centroid of all its topological neighbors
 * (from SurfaceTopologyTracker) and applies a restoring force toward
 * that centroid.  This is the 3D generalization of RingSmoothingForce:
 *
 *   F_i = κ * (centroid(neighbors) − position_i)
 *
 * Effects:
 *   - Regularizes the mesh on the surface (prevents clustering / gaps)
 *   - Does NOT impose a preferred curvature
 *   - Works correctly during budding (purely local smoothing)
 *
 * Requires a SurfaceTopologyTracker to define the adjacency graph.
 */
#ifndef SURFACESMOOTHINGFORCE_HPP_
#define SURFACESMOOTHINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "SurfaceTopologyTracker.hpp"
#include "SimProfiler.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

template<unsigned DIM>
class SurfaceSmoothingForce : public AbstractForce<DIM>
{
private:

    double mSmoothingStiffness;

    /** Per-cell-type stiffness multipliers */
    double mStemScale;
    double mTransitScale;
    double mDiffScale;

    /** Maximum force magnitude per cell (prevents blowup on division) */
    double mMaxForceMagnitude;

    /** Pointer to the surface topology tracker (not owned) */
    SurfaceTopologyTracker<DIM>* mpSurfaceTopology;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM>>(*this);
        archive & mSmoothingStiffness;
        archive & mStemScale;
        archive & mTransitScale;
        archive & mDiffScale;
        archive & mMaxForceMagnitude;
    }

public:

    SurfaceSmoothingForce()
        : AbstractForce<DIM>(),
          mSmoothingStiffness(5.0),
          mStemScale(0.7),
          mTransitScale(1.0),
          mDiffScale(1.3),
          mMaxForceMagnitude(100.0),
          mpSurfaceTopology(nullptr)
    {
    }

    virtual ~SurfaceSmoothingForce() {}

    void SetSmoothingStiffness(double s)     { mSmoothingStiffness = s; }
    double GetSmoothingStiffness() const     { return mSmoothingStiffness; }

    void SetCellTypeScaling(double stem, double transit, double diff)
    {
        mStemScale    = stem;
        mTransitScale = transit;
        mDiffScale    = diff;
    }

    void SetSurfaceTopologyTracker(SurfaceTopologyTracker<DIM>* p)
    {
        mpSurfaceTopology = p;
    }

    void AddForceContribution(
        AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer timer("SurfaceSmoothing");
        if (mpSurfaceTopology == nullptr) return;

        for (auto cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx =
                rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            if (!mpSurfaceTopology->HasCell(idx)) continue;

            const std::set<unsigned>& neighbors =
                mpSurfaceTopology->GetNeighbors(idx);

            if (neighbors.empty()) continue;

            // Compute centroid of neighbors
            c_vector<double, DIM> centroid = zero_vector<double>(DIM);
            unsigned count = 0;

            for (unsigned nbr : neighbors)
            {
                try
                {
                    c_vector<double, DIM> nbrPos =
                        rCellPopulation.GetLocationOfCellCentre(
                            rCellPopulation.GetCellUsingLocationIndex(nbr));
                    centroid += nbrPos;
                    count++;
                }
                catch (Exception&)
                {
                    // Neighbor may have just been killed — skip
                }
            }

            if (count == 0) continue;
            centroid /= (double)count;

            c_vector<double, DIM> myPos =
                rCellPopulation.GetLocationOfCellCentre(*cell_iter);

            // Per-cell-type stiffness scaling
            double scale = mTransitScale;
            if (cell_iter->GetCellProliferativeType()
                    ->template IsType<StemCellProliferativeType>())
            {
                scale = mStemScale;
            }
            else if (cell_iter->GetCellProliferativeType()
                         ->template IsType<DifferentiatedCellProliferativeType>())
            {
                scale = mDiffScale;
            }

            c_vector<double, DIM> force =
                mSmoothingStiffness * scale * (centroid - myPos);

            // Cap force magnitude to prevent blowups during division events
            double forceMag = norm_2(force);
            if (forceMag > mMaxForceMagnitude)
            {
                force *= mMaxForceMagnitude / forceMag;
            }

            rCellPopulation.GetNode(idx)->AddAppliedForceContribution(force);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<SmoothingStiffness>" << mSmoothingStiffness
                     << "</SmoothingStiffness>\n";
        *rParamsFile << "\t\t\t<StemScale>" << mStemScale
                     << "</StemScale>\n";
        *rParamsFile << "\t\t\t<TransitScale>" << mTransitScale
                     << "</TransitScale>\n";
        *rParamsFile << "\t\t\t<DiffScale>" << mDiffScale
                     << "</DiffScale>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceSmoothingForce)

#endif // SURFACESMOOTHINGFORCE_HPP_
