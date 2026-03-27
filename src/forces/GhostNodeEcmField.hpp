/*
 * GhostNodeEcmField.hpp
 *
 * Off-lattice agent-based ECM field using "ghost nodes".
 * Templated on spatial dimension (DIM = 2 or 3).
 *
 * Each ghost node is a movable ECM particle carrying density, fibre
 * orientation, and spring connectivity to its neighbours.
 *
 * Replaces the fixed-grid DynamicECMField with a particle-based approach
 * inspired by Carrasco-Mantis et al. (2024) doi:10.1016/j.compbiomed.2024.109559
 *
 *   - Ghost nodes can move under viscoelastic spring forces (ECM-ECM)
 *   - Cells interact with nearby ghost nodes via springs (cell-ECM)
 *   - Fibre direction modulates anisotropic resistance to movement
 *   - Ghost nodes are removed when density drops below a threshold
 */
#ifndef GHOSTNODEECMFIELD_HPP_
#define GHOSTNODEECMFIELD_HPP_

#include <vector>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <set>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include "UblasVectorInclude.hpp"
#include "RandomNumberGenerator.hpp"
#include "OutputFileHandler.hpp"

/**
 * A single ghost node representing a discrete ECM element.
 */
template<unsigned DIM>
struct GhostNode
{
    unsigned id;                              ///< Unique identifier
    c_vector<double, DIM> position;           ///< Current position (movable)
    c_vector<double, DIM> velocity;           ///< Current velocity
    c_vector<double, DIM> force;              ///< Accumulated force for this timestep
    double density;                           ///< ECM density in [0,1]
    c_vector<double, DIM> fibre_direction;    ///< Fibre orientation (unit vector)
    double anisotropy;                        ///< Degree of fibre alignment in [0,1]
    bool is_active;                           ///< False if removed (density below threshold)
    std::vector<unsigned> neighbours;         ///< Indices of connected ghost node neighbours
};

/**
 * Off-lattice ECM field built from movable ghost nodes.
 *
 * Initialization places ghost nodes on a grid (square/hex for 2D, cubic for 3D)
 * then connects each to its neighbours. Subsequent dynamics are fully
 * off-lattice: nodes move, degrade, remodel, and can be removed.
 */
template<unsigned DIM>
class GhostNodeEcmField
{
private:

    /** All ghost nodes (indexed by id) */
    std::vector<GhostNode<DIM>> mNodes;

    /** Number of active (non-removed) nodes */
    unsigned mNumActive;

    /** Ghost-ghost spring stiffness */
    double mGhostGhostStiffness;

    /** Ghost-ghost spring rest length (initial spacing) */
    double mGhostRestLength;

    /** Ghost-ghost damping coefficient */
    double mGhostDamping;

    /** ECM degradation rate per cell per unit time */
    double mDegradationRate;

    /** Density threshold below which ghost nodes are removed */
    double mRemovalThreshold;

    /** Fibre remodeling rate (how quickly fibres align with traction) */
    double mFibreRemodelingRate;

    /** Anisotropy strength: controls fibre-direction bias on forces */
    double mAnisotropyStrength;

    /** Grid spacing used during initialization */
    double mInitialSpacing;

    /** Domain bounds (tracked for boundary expansion) */
    double mDomainMin[3];
    double mDomainMax[3];

    /** Grid type and fibre pattern (for generating new nodes) */
    std::string mGridType;
    std::string mFibrePattern;

    /** Spatial hash grid for efficient neighbor queries */
    double mSpatialHashCellSize;
    int mHashN[3];        ///< Hash grid dimensions [Nx, Ny, Nz]
    double mHashMin[3];   ///< Hash grid origin [xMin, yMin, zMin]
    std::vector<std::vector<unsigned>> mSpatialHash;

    /** Compute flat bucket index from cell coordinates */
    inline int GetHashIndex(int cx, int cy, int cz = 0) const
    {
        if (DIM == 2)
            return cx * mHashN[1] + cy;
        else
            return (cx * mHashN[1] + cy) * mHashN[2] + cz;
    }

    /** Total number of spatial hash buckets */
    inline int GetTotalBuckets() const
    {
        if (DIM == 2) return mHashN[0] * mHashN[1];
        else return mHashN[0] * mHashN[1] * mHashN[2];
    }

public:

    /**
     * Construct a 2D ghost node ECM field.
     *
     * @param fibrePattern  Initial fibre orientation: "radial", "random", "parallel"
     * @param spacing       Spacing between ghost nodes
     * @param xMin,xMax     Domain extent in x
     * @param yMin,yMax     Domain extent in y
     * @param gridType      "square" or "hex"
     */
    GhostNodeEcmField(const std::string& fibrePattern,
                      double spacing,
                      double xMin, double xMax,
                      double yMin, double yMax,
                      const std::string& gridType = "square")
        : mNumActive(0),
          mGhostGhostStiffness(5.0),
          mGhostRestLength(spacing),
          mGhostDamping(1.0),
          mDegradationRate(0.02),
          mRemovalThreshold(0.01),
          mFibreRemodelingRate(0.1),
          mAnisotropyStrength(0.5),
          mInitialSpacing(spacing),
          mGridType(gridType),
          mFibrePattern(fibrePattern),
          mSpatialHashCellSize(spacing * 2.0)
    {
        mDomainMin[0] = xMin; mDomainMax[0] = xMax;
        mDomainMin[1] = yMin; mDomainMax[1] = yMax;
        mDomainMin[2] = 0.0;  mDomainMax[2] = 0.0;

        mHashN[0] = static_cast<int>(std::ceil((xMax - xMin) / mSpatialHashCellSize)) + 1;
        mHashN[1] = static_cast<int>(std::ceil((yMax - yMin) / mSpatialHashCellSize)) + 1;
        mHashN[2] = 1;
        mHashMin[0] = xMin;
        mHashMin[1] = yMin;
        mHashMin[2] = 0.0;
        mSpatialHash.resize(GetTotalBuckets());

        if (gridType == "hex")
        {
            GenerateHexGrid(spacing, xMin, xMax, yMin, yMax);
        }
        else
        {
            GenerateSquareGrid(spacing, xMin, xMax, yMin, yMax);
        }

        InitialiseFibres(fibrePattern);
        RebuildSpatialHash();
        double neighbour_cutoff = spacing * 1.6;
        BuildNeighbourConnectivity(neighbour_cutoff);

        std::cout << "GhostNodeEcmField<" << DIM << ">: created " << mNodes.size()
                  << " ghost nodes (" << gridType << " grid, spacing="
                  << spacing << ")" << std::endl;
    }

    /**
     * Construct a 3D ghost node ECM field.
     *
     * @param fibrePattern  Initial fibre orientation: "radial", "random", "parallel"
     * @param spacing       Spacing between ghost nodes
     * @param xMin,xMax     Domain extent in x
     * @param yMin,yMax     Domain extent in y
     * @param zMin,zMax     Domain extent in z
     * @param gridType      "cubic"
     */
    GhostNodeEcmField(const std::string& fibrePattern,
                      double spacing,
                      double xMin, double xMax,
                      double yMin, double yMax,
                      double zMin, double zMax,
                      const std::string& gridType = "cubic")
        : mNumActive(0),
          mGhostGhostStiffness(5.0),
          mGhostRestLength(spacing),
          mGhostDamping(1.0),
          mDegradationRate(0.02),
          mRemovalThreshold(0.01),
          mFibreRemodelingRate(0.1),
          mAnisotropyStrength(0.5),
          mInitialSpacing(spacing),
          mGridType(gridType),
          mFibrePattern(fibrePattern),
          mSpatialHashCellSize(spacing * 2.0)
    {
        mDomainMin[0] = xMin; mDomainMax[0] = xMax;
        mDomainMin[1] = yMin; mDomainMax[1] = yMax;
        mDomainMin[2] = zMin; mDomainMax[2] = zMax;

        mHashN[0] = static_cast<int>(std::ceil((xMax - xMin) / mSpatialHashCellSize)) + 1;
        mHashN[1] = static_cast<int>(std::ceil((yMax - yMin) / mSpatialHashCellSize)) + 1;
        mHashN[2] = static_cast<int>(std::ceil((zMax - zMin) / mSpatialHashCellSize)) + 1;
        mHashMin[0] = xMin;
        mHashMin[1] = yMin;
        mHashMin[2] = zMin;
        mSpatialHash.resize(GetTotalBuckets());

        GenerateCubicGrid(spacing, xMin, xMax, yMin, yMax, zMin, zMax);

        InitialiseFibres(fibrePattern);
        RebuildSpatialHash();
        double neighbour_cutoff = spacing * 1.6;
        BuildNeighbourConnectivity(neighbour_cutoff);

        std::cout << "GhostNodeEcmField<" << DIM << ">: created " << mNodes.size()
                  << " ghost nodes (" << gridType << " grid, spacing="
                  << spacing << ")" << std::endl;
    }

    // -- Accessors ------------------------------------------------

    unsigned GetNumNodes() const { return mNodes.size(); }
    unsigned GetNumActive() const { return mNumActive; }

    const GhostNode<DIM>& GetNode(unsigned i) const { return mNodes[i]; }
    GhostNode<DIM>& rGetNode(unsigned i) { return mNodes[i]; }

    const std::vector<GhostNode<DIM>>& GetNodes() const { return mNodes; }

    // -- Force accumulation ----------------------------------------

    /**
     * Clear all ghost node forces. Call at the start of each timestep
     * before accumulating external + internal forces.
     */
    void ClearGhostForces()
    {
        for (auto& node : mNodes)
        {
            if (node.is_active)
                node.force = zero_vector<double>(DIM);
        }
    }

    /**
     * Add an external force to a ghost node (e.g. reaction from cell-ghost spring).
     */
    void AddForceToGhostNode(unsigned index, const c_vector<double, DIM>& force)
    {
        mNodes[index].force += force;
    }

    // -- Parameter setters -----------------------------------------

    void SetGhostGhostStiffness(double k) { mGhostGhostStiffness = k; }
    double GetGhostGhostStiffness() const { return mGhostGhostStiffness; }

    void SetGhostRestLength(double s0) { mGhostRestLength = s0; }
    double GetGhostRestLength() const { return mGhostRestLength; }

    void SetGhostDamping(double eta) { mGhostDamping = eta; }
    double GetGhostDamping() const { return mGhostDamping; }

    void SetDegradationRate(double rate) { mDegradationRate = rate; }
    double GetDegradationRate() const { return mDegradationRate; }

    void SetRemovalThreshold(double thresh) { mRemovalThreshold = thresh; }
    double GetRemovalThreshold() const { return mRemovalThreshold; }

    void SetFibreRemodelingRate(double rate) { mFibreRemodelingRate = rate; }
    void SetAnisotropyStrength(double a) { mAnisotropyStrength = a; }

    double GetDomainMin(unsigned dim) const { return mDomainMin[dim]; }
    double GetDomainMax(unsigned dim) const { return mDomainMax[dim]; }
    double GetInitialSpacing() const { return mInitialSpacing; }

    // -- Core methods ---------------------------------------------

    /**
     * Clear ECM density inside a sphere (3D) or circle (2D).
     */
    void ClearDensityInsideRadius(c_vector<double, DIM> center, double radius)
    {
        for (auto& node : mNodes)
        {
            double dist = norm_2(node.position - center);
            if (dist < radius)
            {
                node.density = 0.0;
                node.is_active = false;
            }
        }
        UpdateActiveCount();
        RebuildSpatialHash();
    }

    /**
     * Remove ghost nodes that are closer than clearance to any of the given points.
     * This produces a tight-fitting hole that matches the actual cell geometry
     * rather than a circular approximation.
     */
    void ClearDensityNearPoints(const std::vector<c_vector<double, DIM>>& points,
                                double clearance)
    {
        double clearance_sq = clearance * clearance;
        for (auto& node : mNodes)
        {
            if (!node.is_active) continue;
            for (const auto& pt : points)
            {
                double dist_sq = inner_prod(node.position - pt, node.position - pt);
                if (dist_sq < clearance_sq)
                {
                    node.density = 0.0;
                    node.is_active = false;
                    break;
                }
            }
        }
        UpdateActiveCount();
        RebuildSpatialHash();
    }

    /**
     * Find all active ghost nodes within cutoff of a position.
     */
    void GetNearbyGhostNodes(const c_vector<double, DIM>& pos,
                             double cutoff,
                             std::vector<unsigned>& rIndices) const
    {
        rIndices.clear();
        double cutoff_sq = cutoff * cutoff;

        int cx = static_cast<int>((pos[0] - mHashMin[0]) / mSpatialHashCellSize);
        int cy = static_cast<int>((pos[1] - mHashMin[1]) / mSpatialHashCellSize);
        int range = static_cast<int>(std::ceil(cutoff / mSpatialHashCellSize));

        if (DIM == 2)
        {
            for (int ix = cx - range; ix <= cx + range; ix++)
            {
                for (int iy = cy - range; iy <= cy + range; iy++)
                {
                    if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1]) continue;
                    const auto& bucket = mSpatialHash[GetHashIndex(ix, iy)];
                    for (unsigned idx : bucket)
                    {
                        c_vector<double, DIM> disp = mNodes[idx].position - pos;
                        double dist_sq = inner_prod(disp, disp);
                        if (dist_sq < cutoff_sq && dist_sq > 1e-20)
                        {
                            rIndices.push_back(idx);
                        }
                    }
                }
            }
        }
        else // DIM == 3
        {
            int cz = static_cast<int>((pos[2] - mHashMin[2]) / mSpatialHashCellSize);
            for (int ix = cx - range; ix <= cx + range; ix++)
            {
                for (int iy = cy - range; iy <= cy + range; iy++)
                {
                    for (int iz = cz - range; iz <= cz + range; iz++)
                    {
                        if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1] ||
                            iz < 0 || iz >= mHashN[2]) continue;
                        const auto& bucket = mSpatialHash[GetHashIndex(ix, iy, iz)];
                        for (unsigned idx : bucket)
                        {
                            c_vector<double, DIM> disp = mNodes[idx].position - pos;
                            double dist_sq = inner_prod(disp, disp);
                            if (dist_sq < cutoff_sq && dist_sq > 1e-20)
                            {
                                rIndices.push_back(idx);
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Compute ghost-ghost spring forces and update ghost node positions.
     *
     * Force on ghost node i from neighbour j (Chaste spring regime):
     *   Compression (d < s0): F = k * s0 * log(1 + (d-s0)/s0)
     *   Extension   (d > s0): F = k * (d-s0) * exp(-5(d-s0)/s0)
     *
     * Fibre-direction anisotropy modulates the force:
     *   F_ij *= (1 - a + a * |cos(angle between displacement and fibre)|)
     *
     * Position update (overdamped): x_new = x_old + (F / eta) * dt
     */
    void UpdateGhostNodePositions(double dt)
    {
        // NOTE: forces should already contain external forces (cell->ghost
        // reactions) accumulated via AddForceToGhostNode().
        // ClearGhostForces() must be called at the start of each timestep.

        // Accumulate ghost-ghost spring forces using Newton's 3rd law:
        // each pair (i,j) is processed once, force applied to both.
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            GhostNode<DIM>& node_i = mNodes[i];
            if (!node_i.is_active) continue;

            for (unsigned n_idx : node_i.neighbours)
            {
                // Only process each pair once (i < j)
                if (n_idx <= i) continue;

                GhostNode<DIM>& node_j = mNodes[n_idx];
                if (!node_j.is_active) continue;

                c_vector<double, DIM> disp = node_j.position - node_i.position;
                double dist = norm_2(disp);
                if (dist < 1e-10) continue;

                c_vector<double, DIM> unit_disp = disp / dist;

                // Chaste spring regime
                double overlap = dist - mGhostRestLength;
                double spring_mag;
                if (overlap <= 0)  // compression: log repulsion
                {
                    spring_mag = mGhostGhostStiffness * mGhostRestLength
                                 * std::log(1.0 + overlap / mGhostRestLength);
                }
                else  // extension: exponential decay
                {
                    double alpha = 5.0;
                    spring_mag = mGhostGhostStiffness * overlap
                                 * std::exp(-alpha * overlap / mGhostRestLength);
                }

                // Anisotropic modulation: average the anisotropy of both nodes
                double aniso_factor = 1.0;
                if (mAnisotropyStrength > 0.0)
                {
                    double avg_aniso = 0.5 * (node_i.anisotropy + node_j.anisotropy);
                    if (avg_aniso > 0.0)
                    {
                        // Use average fibre direction for the pair
                        c_vector<double, DIM> avg_fibre = node_i.fibre_direction + node_j.fibre_direction;
                        double fn = norm_2(avg_fibre);
                        if (fn > 1e-10)
                        {
                            avg_fibre /= fn;
                            double cos_angle_to_fibre = std::abs(inner_prod(unit_disp, avg_fibre));
                            aniso_factor = 1.0 - mAnisotropyStrength * avg_aniso
                                           * (1.0 - cos_angle_to_fibre);
                        }
                    }
                }

                c_vector<double, DIM> force_ij = spring_mag * aniso_factor * unit_disp;
                node_i.force += force_ij;
                node_j.force -= force_ij;
            }
        }

        // Update positions (overdamped dynamics)
        double inv_damping = 1.0 / mGhostDamping;
        double scale = inv_damping * dt;
        for (auto& node : mNodes)
        {
            if (!node.is_active) continue;
            node.position += node.force * scale;
        }

        RebuildSpatialHash();
    }

    /**
     * Degrade ECM at a cell position.
     */
    void DegradeAtPosition(const c_vector<double, DIM>& cellPos,
                           double cutoff, double dt)
    {
        std::vector<unsigned> nearby;
        GetNearbyGhostNodes(cellPos, cutoff, nearby);
        DegradeNearby(cellPos, nearby, cutoff, dt);
    }

    /**
     * Degrade ECM using pre-computed nearby ghost node indices.
     * Avoids redundant spatial hash query when caller already has the list.
     */
    void DegradeNearby(const c_vector<double, DIM>& cellPos,
                       const std::vector<unsigned>& nearbyIndices,
                       double cutoff, double dt)
    {
        double inv_cutoff = 1.0 / cutoff;
        for (unsigned idx : nearbyIndices)
        {
            GhostNode<DIM>& gn = mNodes[idx];
            double dist = norm_2(gn.position - cellPos);
            double weight = std::max(0.0, 1.0 - dist * inv_cutoff);
            gn.density = std::max(0.0, gn.density - mDegradationRate * weight * dt);
        }
    }

    /**
     * Remove ghost nodes whose density has fallen below the removal threshold.
     */
    unsigned RemoveDepletedNodes()
    {
        unsigned removed_count = 0;

        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            GhostNode<DIM>& node = mNodes[i];
            if (!node.is_active) continue;

            if (node.density < mRemovalThreshold)
            {
                node.is_active = false;
                node.density = 0.0;
                removed_count++;

                for (unsigned n_idx : node.neighbours)
                {
                    auto& n_list = mNodes[n_idx].neighbours;
                    n_list.erase(
                        std::remove(n_list.begin(), n_list.end(), i),
                        n_list.end());
                }
                node.neighbours.clear();
            }
        }

        if (removed_count > 0)
        {
            UpdateActiveCount();
            RebuildSpatialHash();
        }

        return removed_count;
    }

    /**
     * Expand the ghost node domain if any cell position is within a threshold
     * of the current domain boundary. New nodes are generated in a strip
     * extending the domain by one spacing width in each direction that needs
     * expansion. Neighbours are connected via spatial hash lookup.
     *
     * @param cellPositions  All current cell centroid positions
     * @param threshold      Distance from boundary that triggers expansion
     * @return Number of new ghost nodes added
     */
    unsigned ExpandBoundary(const std::vector<c_vector<double, DIM>>& cellPositions,
                            double threshold)
    {
        // Determine which faces need expansion
        bool expandLo[3] = {false, false, false};
        bool expandHi[3] = {false, false, false};

        for (const auto& pos : cellPositions)
        {
            for (unsigned d = 0; d < DIM; d++)
            {
                if (pos[d] - mDomainMin[d] < threshold) expandLo[d] = true;
                if (mDomainMax[d] - pos[d] < threshold) expandHi[d] = true;
            }
        }

        bool any_expansion = false;
        for (unsigned d = 0; d < DIM; d++)
        {
            if (expandLo[d] || expandHi[d]) { any_expansion = true; break; }
        }
        if (!any_expansion) return 0;

        // Compute new domain bounds (extend by one spacing in each flagged direction)
        double newMin[3], newMax[3];
        for (unsigned d = 0; d < DIM; d++)
        {
            newMin[d] = expandLo[d] ? mDomainMin[d] - mInitialSpacing : mDomainMin[d];
            newMax[d] = expandHi[d] ? mDomainMax[d] + mInitialSpacing : mDomainMax[d];
        }

        // Generate candidate node positions on the grid within the NEW domain
        // but OUTSIDE the old domain (i.e. only in the expansion strips)
        std::vector<c_vector<double, DIM>> new_positions;
        double spacing = mInitialSpacing;

        if (DIM == 2)
        {
            if (mGridType == "hex")
            {
                double row_spacing = spacing * std::sqrt(3.0) / 2.0;
                int row = 0;
                for (double y = newMin[1]; y <= newMax[1]; y += row_spacing)
                {
                    double x_offset = (row % 2 == 1) ? spacing * 0.5 : 0.0;
                    for (double x = newMin[0] + x_offset; x <= newMax[0]; x += spacing)
                    {
                        // Only add if outside old domain
                        if (x >= mDomainMin[0] && x <= mDomainMax[0] &&
                            y >= mDomainMin[1] && y <= mDomainMax[1])
                            continue;
                        c_vector<double, DIM> pos = zero_vector<double>(DIM);
                        pos[0] = x; pos[1] = y;
                        new_positions.push_back(pos);
                    }
                    row++;
                }
            }
            else // square
            {
                for (double y = newMin[1]; y <= newMax[1]; y += spacing)
                {
                    for (double x = newMin[0]; x <= newMax[0]; x += spacing)
                    {
                        if (x >= mDomainMin[0] && x <= mDomainMax[0] &&
                            y >= mDomainMin[1] && y <= mDomainMax[1])
                            continue;
                        c_vector<double, DIM> pos = zero_vector<double>(DIM);
                        pos[0] = x; pos[1] = y;
                        new_positions.push_back(pos);
                    }
                }
            }
        }
        else // DIM == 3
        {
            for (double z = newMin[2]; z <= newMax[2]; z += spacing)
            {
                for (double y = newMin[1]; y <= newMax[1]; y += spacing)
                {
                    for (double x = newMin[0]; x <= newMax[0]; x += spacing)
                    {
                        if (x >= mDomainMin[0] && x <= mDomainMax[0] &&
                            y >= mDomainMin[1] && y <= mDomainMax[1] &&
                            z >= mDomainMin[2] && z <= mDomainMax[2])
                            continue;
                        c_vector<double, DIM> pos = zero_vector<double>(DIM);
                        pos[0] = x; pos[1] = y; pos[2] = z;
                        new_positions.push_back(pos);
                    }
                }
            }
        }

        if (new_positions.empty()) return 0;

        // Create new ghost nodes
        unsigned first_new_id = mNodes.size();
        for (unsigned i = 0; i < new_positions.size(); i++)
        {
            GhostNode<DIM> gn;
            gn.id = first_new_id + i;
            gn.position = new_positions[i];
            gn.velocity = zero_vector<double>(DIM);
            gn.force = zero_vector<double>(DIM);
            gn.density = 1.0;
            gn.is_active = true;
            gn.anisotropy = 0.0;
            gn.fibre_direction = zero_vector<double>(DIM);

            // Radial fibre pattern relative to origin
            double r = norm_2(gn.position);
            if (r > 1e-10)
            {
                gn.fibre_direction = gn.position / r;
                gn.anisotropy = 0.5;
            }
            else
            {
                gn.fibre_direction[0] = 1.0;
            }

            mNodes.push_back(gn);
        }

        unsigned num_added = new_positions.size();

        // Update domain bounds
        for (unsigned d = 0; d < DIM; d++)
        {
            mDomainMin[d] = newMin[d];
            mDomainMax[d] = newMax[d];
        }

        // Expand spatial hash to cover new domain
        mHashMin[0] = newMin[0];
        mHashMin[1] = newMin[1];
        mHashN[0] = static_cast<int>(std::ceil((newMax[0] - newMin[0]) / mSpatialHashCellSize)) + 1;
        mHashN[1] = static_cast<int>(std::ceil((newMax[1] - newMin[1]) / mSpatialHashCellSize)) + 1;
        if (DIM == 3)
        {
            mHashMin[2] = newMin[2];
            mHashN[2] = static_cast<int>(std::ceil((newMax[2] - newMin[2]) / mSpatialHashCellSize)) + 1;
        }
        mSpatialHash.clear();
        mSpatialHash.resize(GetTotalBuckets());
        RebuildSpatialHash();

        // Build neighbour connectivity for new nodes (connect to both new and existing nodes)
        double neighbour_cutoff = mInitialSpacing * 1.6;
        double cutoff_sq = neighbour_cutoff * neighbour_cutoff;
        int range = static_cast<int>(std::ceil(neighbour_cutoff / mSpatialHashCellSize));

        for (unsigned ni = first_new_id; ni < mNodes.size(); ni++)
        {
            const auto& pos_i = mNodes[ni].position;
            int cx = static_cast<int>((pos_i[0] - mHashMin[0]) / mSpatialHashCellSize);
            int cy = static_cast<int>((pos_i[1] - mHashMin[1]) / mSpatialHashCellSize);

            if (DIM == 2)
            {
                for (int ix = cx - range; ix <= cx + range; ix++)
                {
                    for (int iy = cy - range; iy <= cy + range; iy++)
                    {
                        if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1]) continue;
                        const auto& bucket = mSpatialHash[GetHashIndex(ix, iy)];
                        for (unsigned j : bucket)
                        {
                            if (j == ni) continue;
                            if (!mNodes[j].is_active) continue;
                            // Avoid duplicate links: only connect if j < ni (new-new)
                            // or j is an existing node (j < first_new_id)
                            // For new-new pairs, process once via j < ni
                            if (j >= first_new_id && j >= ni) continue;
                            c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                            double dist_sq = inner_prod(disp, disp);
                            if (dist_sq < cutoff_sq)
                            {
                                mNodes[ni].neighbours.push_back(j);
                                mNodes[j].neighbours.push_back(ni);
                            }
                        }
                    }
                }
            }
            else // DIM == 3
            {
                int cz = static_cast<int>((pos_i[2] - mHashMin[2]) / mSpatialHashCellSize);
                for (int ix = cx - range; ix <= cx + range; ix++)
                {
                    for (int iy = cy - range; iy <= cy + range; iy++)
                    {
                        for (int iz = cz - range; iz <= cz + range; iz++)
                        {
                            if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1]
                                || iz < 0 || iz >= mHashN[2]) continue;
                            const auto& bucket = mSpatialHash[GetHashIndex(ix, iy, iz)];
                            for (unsigned j : bucket)
                            {
                                if (j == ni) continue;
                                if (!mNodes[j].is_active) continue;
                                if (j >= first_new_id && j >= ni) continue;
                                c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                                double dist_sq = inner_prod(disp, disp);
                                if (dist_sq < cutoff_sq)
                                {
                                    mNodes[ni].neighbours.push_back(j);
                                    mNodes[j].neighbours.push_back(ni);
                                }
                            }
                        }
                    }
                }
            }
        }

        UpdateActiveCount();

        std::cout << "  GhostNodeECM: expanded boundary, added " << num_added
                  << " new nodes (" << mNumActive << " total active)" << std::endl;

        return num_added;
    }

    /**
     * Apply cell traction to remodel fibre orientation.
     * Fibres near the cell rotate toward the traction direction.
     */
    void ApplyCellTraction(const c_vector<double, DIM>& cellPos,
                           const c_vector<double, DIM>& traction,
                           double cutoff, double dt)
    {
        std::vector<unsigned> nearby;
        GetNearbyGhostNodes(cellPos, cutoff, nearby);
        RemodelNearby(cellPos, nearby, traction, cutoff, dt);
    }

    /**
     * Remodel fibre orientation using pre-computed nearby ghost node indices.
     * Avoids redundant spatial hash query when caller already has the list.
     */
    void RemodelNearby(const c_vector<double, DIM>& cellPos,
                       const std::vector<unsigned>& nearbyIndices,
                       const c_vector<double, DIM>& traction,
                       double cutoff, double dt)
    {
        double inv_cutoff = 1.0 / cutoff;
        for (unsigned idx : nearbyIndices)
        {
            GhostNode<DIM>& gn = mNodes[idx];
            double dist = norm_2(gn.position - cellPos);
            double weight = std::max(0.0, 1.0 - dist * inv_cutoff);

            // Fibres are unoriented: pick the closest hemisphere
            c_vector<double, DIM> target = traction;
            if (inner_prod(gn.fibre_direction, target) < 0)
                target = -target;

            // Rotate fibre_direction toward traction (linear blend + renormalize)
            double blend = mFibreRemodelingRate * weight * dt;
            gn.fibre_direction += blend * (target - gn.fibre_direction);
            double fn = norm_2(gn.fibre_direction);
            if (fn > 1e-10)
                gn.fibre_direction /= fn;

            // Increase anisotropy toward 1.0
            gn.anisotropy = std::min(1.0, gn.anisotropy + blend);
        }
    }

    /**
     * Write ghost node field to a VTP (VTK PolyData) file for ParaView.
     */
    void WriteOutput(const std::string& filePath, double time) const
    {
        std::ofstream file(filePath);
        if (!file.is_open()) return;

        unsigned num_active = 0;
        for (const auto& n : mNodes) if (n.is_active) num_active++;

        unsigned num_edges = 0;
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (!mNodes[i].is_active) continue;
            for (unsigned n_idx : mNodes[i].neighbours)
            {
                if (n_idx > i && mNodes[n_idx].is_active) num_edges++;
            }
        }

        std::map<unsigned, unsigned> index_map;
        unsigned contiguous = 0;
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (mNodes[i].is_active)
                index_map[i] = contiguous++;
        }

        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"PolyData\" version=\"0.1\">\n";
        file << "<PolyData>\n";
        file << "<Piece NumberOfPoints=\"" << num_active
             << "\" NumberOfVerts=\"0\" NumberOfLines=\"" << num_edges
             << "\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";

        // Points (VTK always uses 3 components)
        file << "<Points>\n";
        file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& n : mNodes)
        {
            if (!n.is_active) continue;
            for (unsigned d = 0; d < DIM; d++)
                file << n.position[d] << " ";
            for (unsigned d = DIM; d < 3; d++)
                file << "0 ";
            file << "\n";
        }
        file << "</DataArray>\n</Points>\n";

        // Lines (ECM-ECM bonds)
        file << "<Lines>\n";
        file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (!mNodes[i].is_active) continue;
            for (unsigned n_idx : mNodes[i].neighbours)
            {
                if (n_idx > i && mNodes[n_idx].is_active)
                {
                    file << index_map[i] << " " << index_map[n_idx] << "\n";
                }
            }
        }
        file << "</DataArray>\n";
        file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        unsigned offset = 0;
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (!mNodes[i].is_active) continue;
            for (unsigned n_idx : mNodes[i].neighbours)
            {
                if (n_idx > i && mNodes[n_idx].is_active)
                {
                    offset += 2;
                    file << offset << "\n";
                }
            }
        }
        file << "</DataArray>\n</Lines>\n";

        // Point data
        file << "<PointData>\n";

        // Density
        file << "<DataArray type=\"Float64\" Name=\"density\" format=\"ascii\">\n";
        for (const auto& n : mNodes) if (n.is_active) file << n.density << "\n";
        file << "</DataArray>\n";

        // Fibre direction as 3D vector
        file << "<DataArray type=\"Float64\" Name=\"fibre_direction\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (const auto& n : mNodes)
        {
            if (!n.is_active) continue;
            for (unsigned d = 0; d < DIM; d++)
                file << n.fibre_direction[d] << " ";
            for (unsigned d = DIM; d < 3; d++)
                file << "0 ";
            file << "\n";
        }
        file << "</DataArray>\n";

        // Anisotropy
        file << "<DataArray type=\"Float64\" Name=\"anisotropy\" format=\"ascii\">\n";
        for (const auto& n : mNodes) if (n.is_active) file << n.anisotropy << "\n";
        file << "</DataArray>\n";

        // Node ID
        file << "<DataArray type=\"Int32\" Name=\"id\" format=\"ascii\">\n";
        for (const auto& n : mNodes) if (n.is_active) file << n.id << "\n";
        file << "</DataArray>\n";

        file << "</PointData>\n";

        file << "</Piece>\n</PolyData>\n</VTKFile>\n";
        file.close();
    }

private:

    void GenerateSquareGrid(double spacing, double xMin, double xMax,
                            double yMin, double yMax)
    {
        unsigned id = 0;
        for (double y = yMin; y <= yMax; y += spacing)
        {
            for (double x = xMin; x <= xMax; x += spacing)
            {
                GhostNode<DIM> gn;
                gn.id = id++;
                gn.position = zero_vector<double>(DIM);
                gn.position[0] = x;
                gn.position[1] = y;
                gn.velocity = zero_vector<double>(DIM);
                gn.force = zero_vector<double>(DIM);
                // Clipped Gaussian: mean=0.95, std=0.025, range [0.9, 1.0]
                double raw_density = RandomNumberGenerator::Instance()->NormalRandomDeviate(0.95, 0.025);
                gn.density = std::max(0.9, std::min(1.0, raw_density));
                gn.fibre_direction = zero_vector<double>(DIM);
                gn.fibre_direction[0] = 1.0;
                gn.anisotropy = 0.0;
                gn.is_active = true;
                mNodes.push_back(gn);
            }
        }
        mNumActive = mNodes.size();
    }

    void GenerateHexGrid(double spacing, double xMin, double xMax,
                         double yMin, double yMax)
    {
        double row_spacing = spacing * std::sqrt(3.0) / 2.0;
        unsigned id = 0;
        int row = 0;
        for (double y = yMin; y <= yMax; y += row_spacing)
        {
            double x_offset = (row % 2 == 1) ? spacing * 0.5 : 0.0;
            for (double x = xMin + x_offset; x <= xMax; x += spacing)
            {
                GhostNode<DIM> gn;
                gn.id = id++;
                gn.position = zero_vector<double>(DIM);
                gn.position[0] = x;
                gn.position[1] = y;
                gn.velocity = zero_vector<double>(DIM);
                gn.force = zero_vector<double>(DIM);
                // Clipped Gaussian: mean=0.95, std=0.025, range [0.9, 1.0]
                double raw_density = RandomNumberGenerator::Instance()->NormalRandomDeviate(0.95, 0.025);
                gn.density = std::max(0.9, std::min(1.0, raw_density));
                gn.fibre_direction = zero_vector<double>(DIM);
                gn.fibre_direction[0] = 1.0;
                gn.anisotropy = 0.0;
                gn.is_active = true;
                mNodes.push_back(gn);
            }
            row++;
        }
        mNumActive = mNodes.size();
    }

    void GenerateCubicGrid(double spacing,
                           double xMin, double xMax,
                           double yMin, double yMax,
                           double zMin, double zMax)
    {
        unsigned id = 0;
        for (double z = zMin; z <= zMax; z += spacing)
        {
            for (double y = yMin; y <= yMax; y += spacing)
            {
                for (double x = xMin; x <= xMax; x += spacing)
                {
                    GhostNode<DIM> gn;
                    gn.id = id++;
                    gn.position = zero_vector<double>(DIM);
                    gn.position[0] = x;
                    gn.position[1] = y;
                    gn.position[2] = z;
                    gn.velocity = zero_vector<double>(DIM);
                    gn.force = zero_vector<double>(DIM);
                    // Clipped Gaussian: mean=0.95, std=0.025, range [0.9, 1.0]
                    double raw_density = RandomNumberGenerator::Instance()->NormalRandomDeviate(0.95, 0.025);
                    gn.density = std::max(0.9, std::min(1.0, raw_density));
                    gn.fibre_direction = zero_vector<double>(DIM);
                    gn.fibre_direction[0] = 1.0;
                    gn.anisotropy = 0.0;
                    gn.is_active = true;
                    mNodes.push_back(gn);
                }
            }
        }
        mNumActive = mNodes.size();
    }

    void InitialiseFibres(const std::string& pattern)
    {
        for (auto& node : mNodes)
        {
            if (pattern == "radial")
            {
                double r = norm_2(node.position);
                if (r > 1e-10)
                {
                    node.fibre_direction = node.position / r;
                }
                else
                {
                    node.fibre_direction = zero_vector<double>(DIM);
                    node.fibre_direction[0] = 1.0;
                }
                node.anisotropy = 0.5;
            }
            else if (pattern == "parallel")
            {
                node.fibre_direction = zero_vector<double>(DIM);
                node.fibre_direction[0] = 1.0;
                node.anisotropy = 1.0;
            }
            else  // "random"
            {
                if (DIM == 2)
                {
                    double angle = RandomNumberGenerator::Instance()->ranf() * 2.0 * M_PI;
                    node.fibre_direction[0] = std::cos(angle);
                    node.fibre_direction[1] = std::sin(angle);
                }
                else  // DIM == 3
                {
                    double theta = std::acos(1.0 - 2.0 * RandomNumberGenerator::Instance()->ranf());
                    double phi = RandomNumberGenerator::Instance()->ranf() * 2.0 * M_PI;
                    node.fibre_direction[0] = std::sin(theta) * std::cos(phi);
                    node.fibre_direction[1] = std::sin(theta) * std::sin(phi);
                    node.fibre_direction[2] = std::cos(theta);
                }
                node.anisotropy = 0.0;
            }
        }
    }

    void BuildNeighbourConnectivity(double cutoff)
    {
        // Use spatial hash for O(N) neighbour search instead of O(N²) brute-force.
        // RebuildSpatialHash() must be called before this method.
        double cutoff_sq = cutoff * cutoff;
        int range = static_cast<int>(std::ceil(cutoff / mSpatialHashCellSize));

        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            const auto& pos_i = mNodes[i].position;
            int cx = static_cast<int>((pos_i[0] - mHashMin[0]) / mSpatialHashCellSize);
            int cy = static_cast<int>((pos_i[1] - mHashMin[1]) / mSpatialHashCellSize);

            if (DIM == 2)
            {
                for (int ix = cx - range; ix <= cx + range; ix++)
                {
                    for (int iy = cy - range; iy <= cy + range; iy++)
                    {
                        if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1]) continue;
                        const auto& bucket = mSpatialHash[GetHashIndex(ix, iy)];
                        for (unsigned j : bucket)
                        {
                            if (j <= i) continue;  // avoid duplicates
                            c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                            double dist_sq = inner_prod(disp, disp);
                            if (dist_sq < cutoff_sq)
                            {
                                mNodes[i].neighbours.push_back(j);
                                mNodes[j].neighbours.push_back(i);
                            }
                        }
                    }
                }
            }
            else // DIM == 3
            {
                int cz = static_cast<int>((pos_i[2] - mHashMin[2]) / mSpatialHashCellSize);
                for (int ix = cx - range; ix <= cx + range; ix++)
                {
                    for (int iy = cy - range; iy <= cy + range; iy++)
                    {
                        for (int iz = cz - range; iz <= cz + range; iz++)
                        {
                            if (ix < 0 || ix >= mHashN[0] || iy < 0 || iy >= mHashN[1]
                                || iz < 0 || iz >= mHashN[2]) continue;
                            const auto& bucket = mSpatialHash[GetHashIndex(ix, iy, iz)];
                            for (unsigned j : bucket)
                            {
                                if (j <= i) continue;
                                c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                                double dist_sq = inner_prod(disp, disp);
                                if (dist_sq < cutoff_sq)
                                {
                                    mNodes[i].neighbours.push_back(j);
                                    mNodes[j].neighbours.push_back(i);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void UpdateActiveCount()
    {
        mNumActive = 0;
        for (const auto& n : mNodes)
        {
            if (n.is_active) mNumActive++;
        }
    }

    void RebuildSpatialHash()
    {
        for (auto& bucket : mSpatialHash)
            bucket.clear();

        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (!mNodes[i].is_active) continue;
            int cx = static_cast<int>((mNodes[i].position[0] - mHashMin[0]) / mSpatialHashCellSize);
            int cy = static_cast<int>((mNodes[i].position[1] - mHashMin[1]) / mSpatialHashCellSize);

            if (DIM == 2)
            {
                if (cx >= 0 && cx < mHashN[0] && cy >= 0 && cy < mHashN[1])
                {
                    mSpatialHash[GetHashIndex(cx, cy)].push_back(i);
                }
            }
            else // DIM == 3
            {
                int cz = static_cast<int>((mNodes[i].position[2] - mHashMin[2]) / mSpatialHashCellSize);
                if (cx >= 0 && cx < mHashN[0] && cy >= 0 && cy < mHashN[1]
                    && cz >= 0 && cz < mHashN[2])
                {
                    mSpatialHash[GetHashIndex(cx, cy, cz)].push_back(i);
                }
            }
        }
    }
};

#endif // GHOSTNODEECMFIELD_HPP_
