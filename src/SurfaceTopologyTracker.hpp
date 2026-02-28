/*
 * SurfaceTopologyTracker.hpp
 *
 * Tracks the **graph topology** of a 3D node-based monolayer organoid.
 *
 * Unlike the 2D RingTopologyTracker (each cell has exactly 2 neighbors),
 * cells on a 3D surface have a variable number of neighbors (~6 on average
 * for a Delaunay triangulation on a sphere).
 *
 * The tracker maintains an undirected adjacency graph of cell-cell
 * connections.  This graph is:
 *   - Initialized once at t = 0 from spatial proximity (Delaunay-like)
 *   - Updated on cell **division**: daughter inherits a subset of the
 *     mother's neighbors (those closest to the daughter), and the
 *     daughter is also connected to the mother.
 *   - Updated on cell **death**: the dead cell is removed and its
 *     former neighbors are reconnected to each other if they fall
 *     within one graph hop of each other (depth-1 reconnection).
 *
 * Forces (springs, smoothing, etc.) query this tracker instead of
 * using the raw spatial neighbor list, ensuring that only topologically
 * valid connections carry forces.
 *
 * Key API:
 *   GetNeighbors(cellIdx)       → set<unsigned>   all neighbors
 *   HasNeighbor(a, b)           → bool            edge query
 *   GetNumNeighbors(cellIdx)    → unsigned         degree
 *   InsertDaughter(...)         → void            handle division
 *   RemoveCell(...)             → void            handle death
 */
#ifndef SURFACETOPOLOGYTRACKER_HPP_
#define SURFACETOPOLOGYTRACKER_HPP_

#include "AbstractCellBasedSimulationModifier.hpp"
#include "AbstractCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"

#include <map>
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <cfloat>

template<unsigned DIM>
class SurfaceTopologyTracker : public AbstractCellBasedSimulationModifier<DIM, DIM>
{
private:
    /* ------------------------------------------------------------------ */
    /*  Data                                                               */
    /* ------------------------------------------------------------------ */

    /** Adjacency list: cell index → set of neighbor indices */
    std::map<unsigned, std::set<unsigned>> mAdjacency;

    /** Center of mass used for initialization (radial direction) */
    c_vector<double, DIM> mCenter;

    /** The spatial cutoff used when building the initial graph */
    double mInitialCutoff;

    /** Whether the topology has been initialized */
    bool mInitialized;

    /** Snapshot of which cells were alive at the last update */
    std::set<unsigned> mKnownCells;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<
            AbstractCellBasedSimulationModifier<DIM, DIM>>(*this);
        archive & mAdjacency;
        archive & mCenter;
        archive & mInitialCutoff;
        archive & mInitialized;
        archive & mKnownCells;
    }

    /* ------------------------------------------------------------------ */
    /*  Helpers                                                            */
    /* ------------------------------------------------------------------ */

    /** Add an undirected edge a↔b */
    void AddEdge(unsigned a, unsigned b)
    {
        mAdjacency[a].insert(b);
        mAdjacency[b].insert(a);
    }

    /** Remove an undirected edge a↔b (if it exists) */
    void RemoveEdge(unsigned a, unsigned b)
    {
        auto it_a = mAdjacency.find(a);
        if (it_a != mAdjacency.end()) it_a->second.erase(b);

        auto it_b = mAdjacency.find(b);
        if (it_b != mAdjacency.end()) it_b->second.erase(a);
    }

public:

    /* ------------------------------------------------------------------ */
    /*  Construction                                                       */
    /* ------------------------------------------------------------------ */

    SurfaceTopologyTracker()
        : AbstractCellBasedSimulationModifier<DIM, DIM>(),
          mInitialCutoff(3.0),
          mInitialized(false)
    {
        for (unsigned d = 0; d < DIM; d++)
            mCenter[d] = 0.0;
    }

    virtual ~SurfaceTopologyTracker() {}

    void SetCenter(const c_vector<double, DIM>& center) { mCenter = center; }
    void SetInitialCutoff(double cutoff) { mInitialCutoff = cutoff; }

    /* ------------------------------------------------------------------ */
    /*  Initialization                                                     */
    /* ------------------------------------------------------------------ */

    /**
     * Build the initial neighbor graph from spatial proximity.
     *
     * Two cells are connected if their distance is less than
     * mInitialCutoff.  This is equivalent to the spatial neighbor
     * list that GeneralisedLinearSpringForce would use, but it is
     * frozen into a topological graph that persists across time.
     */
    void InitializeFromSpatialProximity(
        AbstractCellPopulation<DIM, DIM>& rCellPopulation)
    {
        if (mInitialized) return;

        // Collect all cell indices and positions
        std::vector<unsigned> indices;
        std::map<unsigned, c_vector<double, DIM>> positions;

        for (auto cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            c_vector<double, DIM> pos =
                rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            indices.push_back(idx);
            positions[idx] = pos;
            mKnownCells.insert(idx);
            mAdjacency[idx];  // ensure entry exists (may have 0 neighbors)
        }

        // Build edges for all pairs within cutoff
        for (unsigned i = 0; i < indices.size(); i++)
        {
            for (unsigned j = i + 1; j < indices.size(); j++)
            {
                unsigned a = indices[i];
                unsigned b = indices[j];
                double dist = norm_2(positions[a] - positions[b]);
                if (dist < mInitialCutoff)
                {
                    AddEdge(a, b);
                }
            }
        }

        mInitialized = true;
        std::cout << "SurfaceTopologyTracker: Initialized "
                  << indices.size() << " cells, "
                  << CountEdges() << " edges (cutoff="
                  << mInitialCutoff << ")" << std::endl;
    }

    /* ------------------------------------------------------------------ */
    /*  Queries                                                            */
    /* ------------------------------------------------------------------ */

    /** Get all neighbors of a cell. */
    const std::set<unsigned>& GetNeighbors(unsigned cellIdx) const
    {
        static const std::set<unsigned> empty;
        auto it = mAdjacency.find(cellIdx);
        return (it != mAdjacency.end()) ? it->second : empty;
    }

    /** Check whether two cells are topologically connected. */
    bool HasNeighbor(unsigned a, unsigned b) const
    {
        auto it = mAdjacency.find(a);
        if (it == mAdjacency.end()) return false;
        return it->second.count(b) > 0;
    }

    /** Check whether a cell is tracked. */
    bool HasCell(unsigned cellIdx) const
    {
        return mAdjacency.find(cellIdx) != mAdjacency.end();
    }

    /** Get the number of neighbors (graph degree) of a cell. */
    unsigned GetNumNeighbors(unsigned cellIdx) const
    {
        auto it = mAdjacency.find(cellIdx);
        return (it != mAdjacency.end()) ? it->second.size() : 0;
    }

    /** Count total number of undirected edges. */
    unsigned CountEdges() const
    {
        unsigned total = 0;
        for (const auto& entry : mAdjacency)
            total += entry.second.size();
        return total / 2;
    }

    /* ------------------------------------------------------------------ */
    /*  Division handling                                                   */
    /* ------------------------------------------------------------------ */

    /**
     * Handle cell division: insert a daughter cell into the graph.
     *
     * The daughter inherits connections to the subset of the mother's
     * neighbors that are spatially closer to the daughter than to the
     * mother.  An edge mother↔daughter is always created.
     *
     * @param motherIdx     Location index of the mother cell.
     * @param daughterIdx   Location index of the new daughter cell.
     * @param rCellPopulation  Reference to the cell population (for positions).
     */
    void InsertDaughter(unsigned motherIdx,
                        unsigned daughterIdx,
                        AbstractCellPopulation<DIM, DIM>& rCellPopulation)
    {
        c_vector<double, DIM> motherPos =
            rCellPopulation.GetLocationOfCellCentre(
                rCellPopulation.GetCellUsingLocationIndex(motherIdx));
        c_vector<double, DIM> daughterPos =
            rCellPopulation.GetLocationOfCellCentre(
                rCellPopulation.GetCellUsingLocationIndex(daughterIdx));

        // Get mother's current neighbors
        std::set<unsigned> motherNeighbors = GetNeighbors(motherIdx);  // copy

        // Partition neighbors: those closer to daughter get transferred
        for (unsigned nbr : motherNeighbors)
        {
            if (mAdjacency.find(nbr) == mAdjacency.end()) continue;

            c_vector<double, DIM> nbrPos =
                rCellPopulation.GetLocationOfCellCentre(
                    rCellPopulation.GetCellUsingLocationIndex(nbr));

            double distToMother   = norm_2(nbrPos - motherPos);
            double distToDaughter = norm_2(nbrPos - daughterPos);

            if (distToDaughter < distToMother)
            {
                // Transfer this edge: remove mother↔nbr, add daughter↔nbr
                RemoveEdge(motherIdx, nbr);
                AddEdge(daughterIdx, nbr);
            }
            else
            {
                // Neighbor stays with mother, but also connect to daughter
                // (both mother and daughter are close — shared neighbor)
                AddEdge(daughterIdx, nbr);
            }
        }

        // Always connect mother ↔ daughter
        AddEdge(motherIdx, daughterIdx);

        mKnownCells.insert(daughterIdx);
    }

    /* ------------------------------------------------------------------ */
    /*  Death handling                                                      */
    /* ------------------------------------------------------------------ */

    /**
     * Remove a dead cell from the graph and reconnect its former
     * neighbors to each other (depth-1 reconnection).
     *
     * Before removal the dead cell's neighbor set N is collected.
     * After erasing the dead node, every pair (a, b) in N that was
     * at graph distance 2 through the dead cell is now connected
     * directly, maintaining surface integrity.
     *
     * @param deadIdx  Location index of the dead/sloughed cell.
     */
    void RemoveCell(unsigned deadIdx)
    {
        auto it = mAdjacency.find(deadIdx);
        if (it == mAdjacency.end()) return;

        // Copy the neighbor set before erasing
        std::set<unsigned> deadNeighbors = it->second;

        // Remove all edges involving the dead cell
        for (unsigned nbr : deadNeighbors)
        {
            auto nbr_it = mAdjacency.find(nbr);
            if (nbr_it != mAdjacency.end())
                nbr_it->second.erase(deadIdx);
        }
        mAdjacency.erase(deadIdx);
        mKnownCells.erase(deadIdx);

        // Depth-1 reconnection: connect all former-neighbor pairs
        // that are not already connected
        std::vector<unsigned> nbrVec(deadNeighbors.begin(), deadNeighbors.end());
        for (unsigned i = 0; i < nbrVec.size(); i++)
        {
            for (unsigned j = i + 1; j < nbrVec.size(); j++)
            {
                if (!HasNeighbor(nbrVec[i], nbrVec[j]))
                {
                    AddEdge(nbrVec[i], nbrVec[j]);
                }
            }
        }
    }

    /* ------------------------------------------------------------------ */
    /*  SimulationModifier interface                                        */
    /* ------------------------------------------------------------------ */

    void SetupSolve(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                    std::string outputDirectory) override
    {
        InitializeFromSpatialProximity(rCellPopulation);
    }

    /**
     * Called every timestep.  Detects new cells (from division) and
     * removed cells (from death/sloughing) and updates the graph.
     */
    void UpdateAtEndOfTimeStep(
        AbstractCellPopulation<DIM, DIM>& rCellPopulation) override
    {
        ScopedTimer timer("SurfaceTopology");
        if (!mInitialized)
        {
            InitializeFromSpatialProximity(rCellPopulation);
            return;
        }

        // ---- Gather current live cells ---------------------------------
        std::set<unsigned> currentCells;
        std::map<unsigned, CellPtr> cellMap;

        for (auto cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            currentCells.insert(idx);
            cellMap[idx] = *cell_iter;
        }

        // ---- Detect new cells (born since last update) -----------------
        std::vector<unsigned> newCells;
        for (unsigned idx : currentCells)
        {
            if (mKnownCells.find(idx) == mKnownCells.end())
                newCells.push_back(idx);
        }

        // ---- Detect dead cells (gone since last update) ----------------
        std::vector<unsigned> deadCells;
        for (unsigned idx : mKnownCells)
        {
            if (currentCells.find(idx) == currentCells.end())
                deadCells.push_back(idx);
        }

        // ---- Handle births ---------------------------------------------
        for (unsigned daughterIdx : newCells)
        {
            // Find the most likely mother: the closest existing cell
            c_vector<double, DIM> daughterPos =
                rCellPopulation.GetLocationOfCellCentre(cellMap[daughterIdx]);

            double minDist = DBL_MAX;
            unsigned motherIdx = UINT_MAX;

            for (unsigned existingIdx : mKnownCells)
            {
                if (cellMap.find(existingIdx) == cellMap.end()) continue;

                c_vector<double, DIM> existingPos =
                    rCellPopulation.GetLocationOfCellCentre(cellMap[existingIdx]);
                double dist = norm_2(existingPos - daughterPos);

                if (dist < minDist)
                {
                    minDist = dist;
                    motherIdx = existingIdx;
                }
            }

            if (motherIdx != UINT_MAX)
            {
                InsertDaughter(motherIdx, daughterIdx, rCellPopulation);
            }
        }

        // ---- Handle deaths ---------------------------------------------
        for (unsigned deadIdx : deadCells)
        {
            RemoveCell(deadIdx);
        }
    }

    void OutputSimulationModifierParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<InitialCutoff>" << mInitialCutoff
                     << "</InitialCutoff>\n";
        AbstractCellBasedSimulationModifier<DIM, DIM>::
            OutputSimulationModifierParameters(rParamsFile);
    }
};

// Explicit instantiation
template class SurfaceTopologyTracker<2>;
template class SurfaceTopologyTracker<3>;

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SurfaceTopologyTracker)

#endif // SURFACETOPOLOGYTRACKER_HPP_
