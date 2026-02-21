/*

Copyright (c) 2005-2026, University of Oxford.
All rights reserved.

Inspired by: Drasdo, D. "Buckling Instabilities of One-Layered Growing Tissues"
             Phys. Rev. Lett. 84, 4244 (2000)

This force implements curvature-dependent bending energy to maintain monolayer
integrity in node-based cell simulations. It prevents cells from filling the
lumen by penalizing deviations from the expected surface curvature.

*/

#ifndef CURVATUREBENDINGFORCE_HPP_
#define CURVATUREBENDINGFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

// Forward declaration - include when using
template<unsigned DIM>
class RingTopologyTracker;

#include <vector>
#include <algorithm>
#include <cmath>

/**
 * A force class implementing Drasdo-style curvature-dependent bending energy
 * to enforce monolayer structure in node-based organoid models.
 *
 * The force has two components:
 *
 * 1. BENDING ENERGY (Drasdo 2000):
 *    For cells on the surface, computes discrete curvature from neighbors and
 *    applies restoring force when curvature deviates from target (1/R).
 *    E_bend = (κ_bend/2) * (κ_local - κ_target)²
 *
 * 2. LUMEN EXCLUSION:
 *    Strong radial repulsion for any cell that drifts inside the expected
 *    lumen radius. This catches cells that escape the bending penalty.
 *    F_excl = k_excl * (r_min - r) * r_hat  for r < r_min
 *
 * Together these enforce that cells stay in a single-layer shell around the lumen.
 */
template<unsigned DIM>
class CurvatureBendingForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Bending rigidity - stiffness of curvature penalty */
    double mBendingStiffness;

    /** Target curvature (= 1/target_radius for circular/spherical shell) */
    double mTargetCurvature;

    /** Target shell radius (derived from target curvature) */
    double mTargetRadius;

    /** Lumen exclusion strength - repels cells from entering lumen */
    double mLumenExclusionStrength;

    /** Minimum allowed distance from center (inner edge of shell) */
    double mMinRadiusFraction;  // as fraction of target radius

    /** Maximum allowed distance from center (outer edge of shell) */
    double mMaxRadiusFraction;  // as fraction of target radius

    /** Neighbor search cutoff for curvature calculation */
    double mNeighborCutoff;

    /** Center of the organoid (auto-tracked) */
    c_vector<double, DIM> mOrganoidCenter;

    /** Whether to auto-track center from population centroid */
    bool mTrackCenter;

    /** Optional ring topology tracker for stable neighbors (not serialized - must be re-attached) */
    RingTopologyTracker<DIM>* mpRingTopology;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mBendingStiffness;
        archive & mTargetCurvature;
        archive & mTargetRadius;
        archive & mLumenExclusionStrength;
        archive & mMinRadiusFraction;
        archive & mMaxRadiusFraction;
        archive & mNeighborCutoff;
        archive & mOrganoidCenter;
        archive & mTrackCenter;
    }

public:

    /**
     * Constructor with sensible defaults for organoid simulations.
     */
    CurvatureBendingForce()
        : AbstractForce<DIM>(),
          mBendingStiffness(5.0),
          mTargetCurvature(0.125),   // 1/8.0 for radius 8.0
          mTargetRadius(8.0),
          mLumenExclusionStrength(50.0),  // Strong repulsion from lumen
          mMinRadiusFraction(0.7),   // Cells should stay outside 70% of target radius
          mMaxRadiusFraction(1.5),   // Allow some outward budding
          mNeighborCutoff(2.5),
          mTrackCenter(true),
          mpRingTopology(nullptr)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mOrganoidCenter[i] = 0.0;
        }
    }

    virtual ~CurvatureBendingForce() {}

    // ========== Setters ==========

    void SetBendingStiffness(double stiffness)
    {
        mBendingStiffness = stiffness;
    }

    void SetTargetRadius(double radius)
    {
        mTargetRadius = radius;
        mTargetCurvature = 1.0 / radius;
    }

    void SetLumenExclusionStrength(double strength)
    {
        mLumenExclusionStrength = strength;
    }

    void SetMinRadiusFraction(double fraction)
    {
        mMinRadiusFraction = fraction;
    }

    void SetMaxRadiusFraction(double fraction)
    {
        mMaxRadiusFraction = fraction;
    }

    void SetNeighborCutoff(double cutoff)
    {
        mNeighborCutoff = cutoff;
    }

    void SetOrganoidCenter(const c_vector<double, DIM>& center)
    {
        mOrganoidCenter = center;
    }

    void SetTrackCenter(bool track)
    {
        mTrackCenter = track;
    }

    /**
     * Set the ring topology tracker for stable neighbor definitions.
     * When set, neighbors are taken from the tracker instead of the
     * cell population's GetNeighbouringNodeIndices().
     */
    void SetRingTopologyTracker(RingTopologyTracker<DIM>* pTracker)
    {
        mpRingTopology = pTracker;
    }

    // ========== Getters ==========

    double GetBendingStiffness() const { return mBendingStiffness; }
    double GetTargetRadius() const { return mTargetRadius; }
    double GetTargetCurvature() const { return mTargetCurvature; }
    double GetLumenExclusionStrength() const { return mLumenExclusionStrength; }

    /**
     * Update organoid center from cell population centroid.
     */
    void UpdateOrganoidCenter(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> centroid = zero_vector<double>(DIM);
        unsigned count = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            c_vector<double, DIM> loc = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            centroid += loc;
            count++;
        }

        if (count > 0)
        {
            mOrganoidCenter = centroid / (double)count;
        }
    }

    /**
     * Compute discrete curvature at a cell from its angular neighbors.
     *
     * For 2D: Uses the Menger curvature formula for three consecutive points:
     *   κ = 4 * Area(triangle) / (|P1-P2| * |P2-P3| * |P3-P1|)
     *
     * This identifies cells that have moved inward (higher curvature than expected)
     * or outward (lower curvature).
     */
    double ComputeLocalCurvature2D(
        const c_vector<double, 2>& cellPos,
        const c_vector<double, 2>& leftNeighbor,
        const c_vector<double, 2>& rightNeighbor) const
    {
        // Menger curvature: κ = 4A / (abc) where a,b,c are side lengths
        // and A is the area of the triangle

        double a = norm_2(leftNeighbor - cellPos);
        double b = norm_2(rightNeighbor - cellPos);
        double c = norm_2(rightNeighbor - leftNeighbor);

        // Avoid degenerate triangles
        if (a < 1e-6 || b < 1e-6 || c < 1e-6)
        {
            return mTargetCurvature;  // Return target if degenerate
        }

        // Signed area using cross product (2D)
        double crossZ = (leftNeighbor[0] - cellPos[0]) * (rightNeighbor[1] - cellPos[1])
                      - (leftNeighbor[1] - cellPos[1]) * (rightNeighbor[0] - cellPos[0]);
        double area = 0.5 * std::abs(crossZ);

        double curvature = 4.0 * area / (a * b * c);

        // Sign: positive if curving outward (convex), negative if inward (concave)
        // For a circular shell centered at origin, we want positive curvature
        // Check if cell is on the convex side
        c_vector<double, 2> midpoint = 0.5 * (leftNeighbor + rightNeighbor);
        c_vector<double, 2> toCell = cellPos - midpoint;
        c_vector<double, 2> toCenter = mOrganoidCenter - midpoint;

        // If cell is on opposite side of midpoint from center, curvature is "correct"
        double dotSign = toCell[0] * toCenter[0] + toCell[1] * toCenter[1];
        if (dotSign > 0)
        {
            // Cell is between neighbors and center - wrong side, high curvature penalty
            curvature = -curvature;
        }

        return curvature;
    }

    /**
     * Find the two angular neighbors of a cell (for 2D shell).
     * These are the neighbors that are roughly on either side when
     * moving around the shell circumference.
     */
    std::pair<unsigned, unsigned> FindAngularNeighbors2D(
        unsigned cellIndex,
        const c_vector<double, 2>& cellPos,
        const std::vector<std::pair<unsigned, c_vector<double, 2>>>& neighbors) const
    {
        if (neighbors.size() < 2)
        {
            return std::make_pair(UINT_MAX, UINT_MAX);
        }

        // Angular position of this cell relative to center
        c_vector<double, 2> relPos = cellPos - mOrganoidCenter;
        double cellAngle = std::atan2(relPos[1], relPos[0]);

        // Sort neighbors by their angular position
        std::vector<std::pair<double, unsigned>> angularNeighbors;
        for (const auto& nbr : neighbors)
        {
            c_vector<double, 2> nbrRel = nbr.second - mOrganoidCenter;
            double nbrAngle = std::atan2(nbrRel[1], nbrRel[0]);
            double angleDiff = nbrAngle - cellAngle;

            // Normalize to [-π, π]
            while (angleDiff > M_PI) angleDiff -= 2.0 * M_PI;
            while (angleDiff < -M_PI) angleDiff += 2.0 * M_PI;

            angularNeighbors.push_back(std::make_pair(angleDiff, nbr.first));
        }

        std::sort(angularNeighbors.begin(), angularNeighbors.end());

        // Find closest neighbor on each side (positive and negative angle)
        unsigned leftIdx = UINT_MAX;
        unsigned rightIdx = UINT_MAX;
        double minPosAngle = 2.0 * M_PI;
        double maxNegAngle = -2.0 * M_PI;

        for (const auto& an : angularNeighbors)
        {
            if (an.first > 0.01 && an.first < minPosAngle)
            {
                minPosAngle = an.first;
                rightIdx = an.second;
            }
            if (an.first < -0.01 && an.first > maxNegAngle)
            {
                maxNegAngle = an.first;
                leftIdx = an.second;
            }
        }

        return std::make_pair(leftIdx, rightIdx);
    }

    /**
     * Main force contribution method.
     *
     * Works with any centre-based population (NodeBased or MeshBased).
     * Uses the population's GetNeighbouringNodeIndices() which returns:
     *   - For MeshBasedCellPopulation: Delaunay mesh neighbors (topological)
     *   - For NodeBasedCellPopulation: distance-based neighbors (touching cells)
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override
    {
        // Cast to centre-based population (works for Node, Mesh, etc.)
        AbstractCentreBasedCellPopulation<DIM>* p_centre_pop =
            dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation);

        if (!p_centre_pop)
        {
            // Not a centre-based population - skip
            return;
        }

        // Update center if tracking enabled
        if (mTrackCenter)
        {
            UpdateOrganoidCenter(rCellPopulation);
        }

        double minRadius = mMinRadiusFraction * mTargetRadius;

        // Build position map for neighbor lookup
        std::map<unsigned, c_vector<double, DIM>> nodePositions;
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned nodeIdx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            nodePositions[nodeIdx] = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
        }

        // Apply forces to each cell
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned nodeIdx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> cellPos = nodePositions[nodeIdx];
            c_vector<double, DIM> displacement = cellPos - mOrganoidCenter;
            double radius = norm_2(displacement);

            // Unit radial vector (outward from center)
            c_vector<double, DIM> radialUnit = zero_vector<double>(DIM);
            if (radius > 1e-6)
            {
                radialUnit = displacement / radius;
            }
            else
            {
                // Cell at center - give it a random push outward
                radialUnit[0] = 1.0;
            }

            c_vector<double, DIM> force = zero_vector<double>(DIM);

            // ===== COMPONENT 1: Lumen Exclusion Force =====
            // Strong outward force if cell is inside the minimum radius
            if (radius < minRadius)
            {
                double penetration = minRadius - radius;
                // Quadratic penalty for deeper penetration
                double forceMag = mLumenExclusionStrength * penetration * (1.0 + penetration / minRadius);
                force += forceMag * radialUnit;
            }

            // ===== COMPONENT 2: Bending/Curvature Force (2D only) =====
            if constexpr (DIM == 2)
            {
                if (mBendingStiffness > 0)
                {
                    c_vector<double, 2> cellPos2D;
                    cellPos2D[0] = cellPos[0];
                    cellPos2D[1] = cellPos[1];

                    c_vector<double, 2> leftPos, rightPos;
                    bool hasNeighbors = false;

                    // Try ring topology tracker first (stable, position-independent neighbors)
                    if (mpRingTopology != nullptr && mpRingTopology->HasNeighbors(nodeIdx))
                    {
                        auto ringNeighbors = mpRingTopology->GetNeighbors(nodeIdx);
                        unsigned leftIdx = ringNeighbors.first;
                        unsigned rightIdx = ringNeighbors.second;

                        if (nodePositions.count(leftIdx) && nodePositions.count(rightIdx))
                        {
                            leftPos[0] = nodePositions[leftIdx][0];
                            leftPos[1] = nodePositions[leftIdx][1];
                            rightPos[0] = nodePositions[rightIdx][0];
                            rightPos[1] = nodePositions[rightIdx][1];
                            hasNeighbors = true;
                        }
                    }
                    
                    // Fallback to population-based neighbor lookup
                    if (!hasNeighbors)
                    {
                        std::set<unsigned> neighborIndices = rCellPopulation.GetNeighbouringNodeIndices(nodeIdx);

                        std::vector<std::pair<unsigned, c_vector<double, 2>>> neighbors;
                        for (unsigned nbrIdx : neighborIndices)
                        {
                            if (nodePositions.count(nbrIdx))
                            {
                                c_vector<double, 2> nbrPos;
                                nbrPos[0] = nodePositions[nbrIdx][0];
                                nbrPos[1] = nodePositions[nbrIdx][1];
                                neighbors.push_back(std::make_pair(nbrIdx, nbrPos));
                            }
                        }

                        if (neighbors.size() >= 2)
                        {
                            auto angularNbrs = FindAngularNeighbors2D(nodeIdx, cellPos2D, neighbors);

                            if (angularNbrs.first != UINT_MAX && angularNbrs.second != UINT_MAX)
                            {
                                for (const auto& nbr : neighbors)
                                {
                                    if (nbr.first == angularNbrs.first) leftPos = nbr.second;
                                    if (nbr.first == angularNbrs.second) rightPos = nbr.second;
                                }
                                hasNeighbors = true;
                            }
                        }
                    }

                    // Compute bending force if we have valid neighbors
                    if (hasNeighbors)
                    {
                        double localCurvature = ComputeLocalCurvature2D(cellPos2D, leftPos, rightPos);
                        double curvatureError = localCurvature - mTargetCurvature;

                        // Bending force: push radially to correct curvature
                        // Positive error (too curved inward) -> push outward
                        // Negative error (too flat/curved outward) -> push inward
                        double bendingForceMag = mBendingStiffness * curvatureError;

                        force += bendingForceMag * radialUnit;
                    }
                }
            }

            // Apply force to node
            rCellPopulation.GetNode(nodeIdx)->AddAppliedForceContribution(force);
        }
    }

    /**
     * Output parameters to file.
     */
    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<BendingStiffness>" << mBendingStiffness << "</BendingStiffness>\n";
        *rParamsFile << "\t\t\t<TargetRadius>" << mTargetRadius << "</TargetRadius>\n";
        *rParamsFile << "\t\t\t<TargetCurvature>" << mTargetCurvature << "</TargetCurvature>\n";
        *rParamsFile << "\t\t\t<LumenExclusionStrength>" << mLumenExclusionStrength << "</LumenExclusionStrength>\n";
        *rParamsFile << "\t\t\t<MinRadiusFraction>" << mMinRadiusFraction << "</MinRadiusFraction>\n";
        *rParamsFile << "\t\t\t<NeighborCutoff>" << mNeighborCutoff << "</NeighborCutoff>\n";

        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

// Explicit instantiation for 2D and 3D
template class CurvatureBendingForce<2>;
template class CurvatureBendingForce<3>;

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CurvatureBendingForce)

#endif // CURVATUREBENDINGFORCE_HPP_
