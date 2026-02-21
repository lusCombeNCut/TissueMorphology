/*
 * RingTopologyTracker.hpp
 *
 * Tracks the 1D ring topology of a monolayer organoid simulation.
 * 
 * For a ring of cells around a lumen, each cell has exactly 2 neighbors
 * (left and right). This tracker maintains a circular doubly-linked list
 * of cell IDs representing this topology.
 *
 * Key features:
 *   - Initializes based on angular position (only at t=0)
 *   - Updates topology when cells divide (inserts daughter next to parent)
 *   - Provides O(1) neighbor lookup
 *   - Does NOT depend on cell positions after initialization
 *
 * Usage:
 *   RingTopologyTracker tracker;
 *   tracker.InitializeFromAngularPositions(population, center);
 *   
 *   // Get neighbors for curvature calculation
 *   auto [left, right] = tracker.GetNeighbors(cellId);
 */
#ifndef RINGTOPOLOGYTRACKER_HPP_
#define RINGTOPOLOGYTRACKER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellPopulation.hpp"

#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

template<unsigned DIM>
class RingTopologyTracker : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:
    /** Maps cell location index -> {left neighbor index, right neighbor index} */
    std::map<unsigned, std::pair<unsigned, unsigned>> mNeighborMap;

    /** Center of the ring (for initial angular sorting) */
    c_vector<double, DIM> mCenter;

    /** Whether topology has been initialized */
    bool mInitialized;

    /** Track which cells existed at last update (for detecting new cells) */
    std::set<unsigned> mKnownCells;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM, DIM>>(*this);
        archive & mNeighborMap;
        archive & mCenter;
        archive & mInitialized;
        archive & mKnownCells;
    }

public:
    RingTopologyTracker()
        : AbstractCellBasedSimulationModifier<DIM, DIM>(),
          mInitialized(false)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mCenter[i] = 0.0;
        }
    }

    virtual ~RingTopologyTracker() {}

    void SetCenter(const c_vector<double, DIM>& center)
    {
        mCenter = center;
    }

    /**
     * Initialize the ring topology from angular positions.
     * Call this ONCE at the start of simulation.
     */
    void InitializeFromAngularPositions(AbstractCellPopulation<DIM, DIM>& rCellPopulation)
    {
        if (mInitialized)
        {
            return;  // Already initialized
        }

        // Collect all cell indices and their angular positions
        std::vector<std::pair<double, unsigned>> angularCells;

        for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            
            // Compute angle from center (only uses x,y for 2D)
            double angle = std::atan2(pos[1] - mCenter[1], pos[0] - mCenter[0]);
            angularCells.push_back(std::make_pair(angle, idx));
            mKnownCells.insert(idx);
        }

        // Sort by angle
        std::sort(angularCells.begin(), angularCells.end());

        // Build circular neighbor relationships
        unsigned n = angularCells.size();
        for (unsigned i = 0; i < n; i++)
        {
            unsigned thisIdx = angularCells[i].second;
            unsigned leftIdx = angularCells[(i + n - 1) % n].second;  // Previous in ring
            unsigned rightIdx = angularCells[(i + 1) % n].second;     // Next in ring

            mNeighborMap[thisIdx] = std::make_pair(leftIdx, rightIdx);
        }

        mInitialized = true;
        std::cout << "RingTopologyTracker: Initialized with " << n << " cells" << std::endl;
    }

    /**
     * Get the two neighbors of a cell.
     * Returns {leftNeighborIdx, rightNeighborIdx}
     */
    std::pair<unsigned, unsigned> GetNeighbors(unsigned cellIdx) const
    {
        auto it = mNeighborMap.find(cellIdx);
        if (it != mNeighborMap.end())
        {
            return it->second;
        }
        return std::make_pair(UINT_MAX, UINT_MAX);  // Not found
    }

    /**
     * Check if a cell has valid neighbors in the topology.
     */
    bool HasNeighbors(unsigned cellIdx) const
    {
        return mNeighborMap.find(cellIdx) != mNeighborMap.end();
    }

    /**
     * Insert a new cell (daughter) into the ring after a division.
     * The daughter is inserted between the parent and one of its neighbors.
     *
     * @param parentIdx  The dividing cell's index
     * @param daughterIdx  The new cell's index
     * @param insertRight  If true, insert between parent and right neighbor;
     *                     if false, insert between parent and left neighbor
     */
    void InsertDaughter(unsigned parentIdx, unsigned daughterIdx, bool insertRight = true)
    {
        auto it = mNeighborMap.find(parentIdx);
        if (it == mNeighborMap.end())
        {
            std::cerr << "Warning: Parent cell " << parentIdx << " not in topology" << std::endl;
            return;
        }

        unsigned leftOfParent = it->second.first;
        unsigned rightOfParent = it->second.second;

        if (insertRight)
        {
            // Insert daughter between parent and right neighbor
            // Before: ... <-> parent <-> rightNeighbor <-> ...
            // After:  ... <-> parent <-> daughter <-> rightNeighbor <-> ...

            // Daughter's neighbors
            mNeighborMap[daughterIdx] = std::make_pair(parentIdx, rightOfParent);

            // Update parent's right neighbor to daughter
            mNeighborMap[parentIdx].second = daughterIdx;

            // Update rightNeighbor's left neighbor to daughter
            if (mNeighborMap.find(rightOfParent) != mNeighborMap.end())
            {
                mNeighborMap[rightOfParent].first = daughterIdx;
            }
        }
        else
        {
            // Insert daughter between left neighbor and parent
            // Before: ... <-> leftNeighbor <-> parent <-> ...
            // After:  ... <-> leftNeighbor <-> daughter <-> parent <-> ...

            mNeighborMap[daughterIdx] = std::make_pair(leftOfParent, parentIdx);
            mNeighborMap[parentIdx].first = daughterIdx;

            if (mNeighborMap.find(leftOfParent) != mNeighborMap.end())
            {
                mNeighborMap[leftOfParent].second = daughterIdx;
            }
        }

        mKnownCells.insert(daughterIdx);
    }

    /**
     * Remove a cell from the topology (e.g., when it dies/sloughs).
     * Reconnects its left and right neighbors directly.
     */
    void RemoveCell(unsigned cellIdx)
    {
        auto it = mNeighborMap.find(cellIdx);
        if (it == mNeighborMap.end())
        {
            return;  // Not in topology
        }

        unsigned leftNeighbor = it->second.first;
        unsigned rightNeighbor = it->second.second;

        // Reconnect: leftNeighbor <-> rightNeighbor
        if (mNeighborMap.find(leftNeighbor) != mNeighborMap.end())
        {
            mNeighborMap[leftNeighbor].second = rightNeighbor;
        }
        if (mNeighborMap.find(rightNeighbor) != mNeighborMap.end())
        {
            mNeighborMap[rightNeighbor].first = leftNeighbor;
        }

        mNeighborMap.erase(cellIdx);
        mKnownCells.erase(cellIdx);
    }

    /**
     * Called at each timestep. Detects new cells (from division) and
     * removes dead cells from topology.
     */
    void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        // Initialize on first call
        if (!mInitialized)
        {
            InitializeFromAngularPositions(rCellPopulation);
            return;
        }

        // Find current cells
        std::set<unsigned> currentCells;
        std::map<unsigned, CellPtr> cellMap;

        for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            currentCells.insert(idx);
            cellMap[idx] = *cell_iter;
        }

        // Find new cells (born since last update)
        std::vector<unsigned> newCells;
        for (unsigned idx : currentCells)
        {
            if (mKnownCells.find(idx) == mKnownCells.end())
            {
                newCells.push_back(idx);
            }
        }

        // Find removed cells (died since last update)
        std::vector<unsigned> deadCells;
        for (unsigned idx : mKnownCells)
        {
            if (currentCells.find(idx) == currentCells.end())
            {
                deadCells.push_back(idx);
            }
        }

        // Handle new cells (from division)
        for (unsigned daughterIdx : newCells)
        {
            // Find the parent by checking which existing cell is closest
            c_vector<double, DIM> daughterPos = rCellPopulation.GetLocationOfCellCentre(cellMap[daughterIdx]);
            
            double minDist = DBL_MAX;
            unsigned parentIdx = UINT_MAX;

            for (unsigned existingIdx : mKnownCells)
            {
                if (cellMap.find(existingIdx) == cellMap.end()) continue;
                
                c_vector<double, DIM> existingPos = rCellPopulation.GetLocationOfCellCentre(cellMap[existingIdx]);
                double dist = norm_2(existingPos - daughterPos);
                
                if (dist < minDist)
                {
                    minDist = dist;
                    parentIdx = existingIdx;
                }
            }

            if (parentIdx != UINT_MAX)
            {
                // Insert daughter next to parent
                // Determine which side based on angular position
                c_vector<double, DIM> parentPos = rCellPopulation.GetLocationOfCellCentre(cellMap[parentIdx]);
                
                double parentAngle = std::atan2(parentPos[1] - mCenter[1], parentPos[0] - mCenter[0]);
                double daughterAngle = std::atan2(daughterPos[1] - mCenter[1], daughterPos[0] - mCenter[0]);
                
                // Normalize angle difference to [-π, π]
                double angleDiff = daughterAngle - parentAngle;
                while (angleDiff > M_PI) angleDiff -= 2.0 * M_PI;
                while (angleDiff < -M_PI) angleDiff += 2.0 * M_PI;

                // Insert on the appropriate side
                InsertDaughter(parentIdx, daughterIdx, angleDiff > 0);
            }
        }

        // Handle dead cells
        for (unsigned deadIdx : deadCells)
        {
            RemoveCell(deadIdx);
        }
    }

    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation, std::string outputDirectory) override
    {
        // Initialize topology at start of simulation
        InitializeFromAngularPositions(rCellPopulation);
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<RingTopologyTracker/>\n";
        AbstractCellBasedSimulationModifier<DIM, DIM>::OutputSimulationModifierParameters(rParamsFile);
    }
};

// Explicit instantiation
template class RingTopologyTracker<2>;
template class RingTopologyTracker<3>;

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RingTopologyTracker)

#endif // RINGTOPOLOGYTRACKER_HPP_
