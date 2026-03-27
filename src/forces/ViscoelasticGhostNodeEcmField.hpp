/*
 * ViscoelasticGhostNodeEcmField.hpp
 *
 * Off-lattice agent-based ECM field with Standard Linear Solid (SLS/Zener)
 * constitutive law. Templated on spatial dimension (DIM = 2 or 3).
 *
 * Each ghost node carries density, fibre orientation, and spring connectivity.
 * Ghost–ghost springs implement a true SLS with two per-pair rest lengths:
 *
 *   E(t) = E_0 + E_1 * exp(-t/tau)
 *
 * Discrete SLS split-force formulation:
 *   Force: F_ij = E_0*(d_ij - s0_ij) + E_1*(d_ij - s_ij(t))
 *   where s0_ij = permanent rest length (initial distance, never changes)
 *         s_ij  = evolving dashpot length (relaxes toward current distance)
 *
 * Rest length evolution uses the exact exponential integrator:
 *   s_ij^{n+1} = d_ij + (s_ij^n - d_ij) * exp(-dt/tau)
 *
 * This gives: instantaneous stiffness = E_0+E_1, long-time stiffness = E_0 != 0
 * (solid-like; force never decays to zero).
 *
 * Based on the constitutive framework of:
 *   Fertala et al. (2025) doi:10.1101/2025.07.02.662292
 *   "Scale-Specific Viscoelastic Characterization of Hydrogels"
 *
 * Parameter values calibrated for cell-diameter length scale (~10-20 um)
 * using microgel AFM stress relaxation data from the above reference.
 */
#ifndef VISCOELASTICGHOSTNODEECMFIELD_HPP_
#define VISCOELASTICGHOSTNODEECMFIELD_HPP_

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
 * A single ghost node representing a discrete ECM element with
 * viscoelastic per-pair rest length storage.
 */
template<unsigned DIM>
struct ViscoelasticGhostNode
{
    unsigned id;                              ///< Unique identifier
    c_vector<double, DIM> position;           ///< Current position (movable)
    c_vector<double, DIM> force;              ///< Accumulated force for this timestep
    double density;                           ///< ECM density in [0,1]
    c_vector<double, DIM> fibre_direction;    ///< Fibre orientation (unit vector)
    double anisotropy;                        ///< Degree of fibre alignment in [0,1]
    bool is_active;                           ///< False if removed (density below threshold)
    std::vector<unsigned> neighbours;         ///< Indices of connected ghost node neighbours
    std::vector<double> rest_lengths;         ///< Per-pair evolving dashpot lengths s_ij(t)
    std::vector<double> initial_rest_lengths; ///< Per-pair permanent rest lengths s0_ij (never change)
};

/**
 * Off-lattice ECM field with Standard Linear Solid (SLS/Zener) constitutive law.
 *
 * Each ghost–ghost spring pair (i,j) maintains two rest lengths:
 *   s0_ij   — permanent rest length (initial distance; represents covalent crosslinks)
 *   s_ij(t) — evolving dashpot length (relaxes toward current distance on timescale tau;
 *              represents reversible physical network rearrangements)
 *
 * SLS split-force:
 *   F = E_0*(d_ij - s0_ij) + E_1*(d_ij - s_ij(t))
 *
 *   - Instantaneous stiffness: E_0 + E_1
 *   - Long-time stiffness:     E_0  (non-zero; solid-like)
 *   - Relaxation time:         tau = eta_VE / E_1
 */
template<unsigned DIM>
class ViscoelasticGhostNodeEcmField
{
private:

    /** All ghost nodes (indexed by id) */
    std::vector<ViscoelasticGhostNode<DIM>> mNodes;

    /** Number of active (non-removed) nodes */
    unsigned mNumActive;

    // ── Viscoelastic constitutive parameters ──────────────────

    /** Relaxed (long-time, equilibrium) modulus E_0 [Pa in Chaste force units] */
    double mRelaxedStiffness;    // E_0

    /** Relaxation modulus E_1: viscous contribution that decays [Pa] */
    double mRelaxationModulus;   // E_1

    /** Viscoelastic relaxation time tau = eta_VE / E_1 [time units] */
    double mRelaxationTime;      // tau

    // ── Equation-of-motion drag (NOT a constitutive parameter) ──

    /** Drag coefficient for overdamped position update */
    double mGhostDamping;

    // ── ECM biology parameters ────────────────────────────────

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

    // ── Spatial hash for neighbour queries ─────────────────────

    double mSpatialHashCellSize;
    int mHashN[3];
    double mHashMin[3];
    std::vector<std::vector<unsigned>> mSpatialHash;

    inline int GetHashIndex(int cx, int cy, int cz = 0) const
    {
        if (DIM == 2)
            return cx * mHashN[1] + cy;
        else
            return (cx * mHashN[1] + cy) * mHashN[2] + cz;
    }

    inline int GetTotalBuckets() const
    {
        if (DIM == 2) return mHashN[0] * mHashN[1];
        else return mHashN[0] * mHashN[1] * mHashN[2];
    }

public:

    /**
     * Construct a 2D viscoelastic ghost node ECM field.
     */
    ViscoelasticGhostNodeEcmField(const std::string& fibrePattern,
                                  double spacing,
                                  double xMin, double xMax,
                                  double yMin, double yMax,
                                  const std::string& gridType = "square")
        : mNumActive(0),
          mRelaxedStiffness(5.0),
          mRelaxationModulus(1.0),
          mRelaxationTime(1.0),
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
            GenerateHexGrid(spacing, xMin, xMax, yMin, yMax);
        else
            GenerateSquareGrid(spacing, xMin, xMax, yMin, yMax);

        InitialiseFibres(fibrePattern);
        RebuildSpatialHash();
        double neighbour_cutoff = spacing * 1.6;
        BuildNeighbourConnectivity(neighbour_cutoff);

        std::cout << "ViscoelasticGhostNodeEcmField<" << DIM << ">: created "
                  << mNodes.size() << " ghost nodes (" << gridType
                  << " grid, spacing=" << spacing << ")" << std::endl;
    }

    /**
     * Construct a 3D viscoelastic ghost node ECM field.
     */
    ViscoelasticGhostNodeEcmField(const std::string& fibrePattern,
                                  double spacing,
                                  double xMin, double xMax,
                                  double yMin, double yMax,
                                  double zMin, double zMax,
                                  const std::string& gridType = "cubic")
        : mNumActive(0),
          mRelaxedStiffness(5.0),
          mRelaxationModulus(1.0),
          mRelaxationTime(1.0),
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

        std::cout << "ViscoelasticGhostNodeEcmField<" << DIM << ">: created "
                  << mNodes.size() << " ghost nodes (" << gridType
                  << " grid, spacing=" << spacing << ")" << std::endl;
    }

    // ── Accessors ────────────────────────────────────────────

    unsigned GetNumNodes() const { return mNodes.size(); }
    unsigned GetNumActive() const { return mNumActive; }

    const ViscoelasticGhostNode<DIM>& GetNode(unsigned i) const { return mNodes[i]; }
    ViscoelasticGhostNode<DIM>& rGetNode(unsigned i) { return mNodes[i]; }

    const std::vector<ViscoelasticGhostNode<DIM>>& GetNodes() const { return mNodes; }

    // ── Force accumulation ────────────────────────────────────

    void ClearGhostForces()
    {
        for (auto& node : mNodes)
        {
            if (node.is_active)
                node.force = zero_vector<double>(DIM);
        }
    }

    void AddForceToGhostNode(unsigned index, const c_vector<double, DIM>& force)
    {
        mNodes[index].force += force;
    }

    // ── Viscoelastic parameter setters ────────────────────────

    /** Set the relaxed (equilibrium) stiffness E_0 */
    void SetRelaxedStiffness(double E0) { mRelaxedStiffness = E0; }
    double GetRelaxedStiffness() const { return mRelaxedStiffness; }

    /** Set the relaxation modulus E_1 (viscous contribution) */
    void SetRelaxationModulus(double E1) { mRelaxationModulus = E1; }
    double GetRelaxationModulus() const { return mRelaxationModulus; }

    /** Set the viscoelastic relaxation time tau */
    void SetRelaxationTime(double tau) { mRelaxationTime = tau; }
    double GetRelaxationTime() const { return mRelaxationTime; }

    /** Get the instantaneous (unrelaxed) stiffness E_u = E_0 + E_1 */
    double GetInstantaneousStiffness() const { return mRelaxedStiffness + mRelaxationModulus; }

    // ── EoM drag setter (distinct from constitutive parameters) ──

    void SetGhostDamping(double eta) { mGhostDamping = eta; }
    double GetGhostDamping() const { return mGhostDamping; }

    // ── Biology parameter setters ─────────────────────────────

    void SetDegradationRate(double rate) { mDegradationRate = rate; }
    double GetDegradationRate() const { return mDegradationRate; }

    void SetRemovalThreshold(double thresh) { mRemovalThreshold = thresh; }
    double GetRemovalThreshold() const { return mRemovalThreshold; }

    void SetFibreRemodelingRate(double rate) { mFibreRemodelingRate = rate; }
    void SetAnisotropyStrength(double a) { mAnisotropyStrength = a; }

    double GetDomainMin(unsigned dim) const { return mDomainMin[dim]; }
    double GetDomainMax(unsigned dim) const { return mDomainMax[dim]; }
    double GetInitialSpacing() const { return mInitialSpacing; }

    // ── Core methods ─────────────────────────────────────────

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
     * Compute ghost–ghost SLS spring forces, evolve dashpot rest lengths,
     * and update positions.
     *
     * For each connected pair (i,j):
     *   1. Compute current distance d_ij
     *   2. SLS split-force:
     *        F = [E_0*(d_ij - s0_ij) + E_1*(d_ij - s_ij)] * A_ij * d_hat_ij
     *      where s0_ij is permanent and s_ij(t) is the evolving dashpot arm
     *   3. Evolve dashpot rest length (exact exponential integrator):
     *        s_ij^{n+1} = d_ij + (s_ij^n - d_ij) * exp(-dt/tau)
     *      Unconditionally stable for all dt > 0, tau > 0.
     *   4. Update position (overdamped EoM): x += F_total/eta_drag * dt
     */
    void UpdateGhostNodePositions(double dt)
    {
        // SLS: E_0 arm and E_1 (Maxwell) arm applied separately via split-force

        // Accumulate ghost-ghost viscoelastic spring forces (Newton's 3rd law)
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            ViscoelasticGhostNode<DIM>& node_i = mNodes[i];
            if (!node_i.is_active) continue;

            for (unsigned nb = 0; nb < node_i.neighbours.size(); nb++)
            {
                unsigned j = node_i.neighbours[nb];
                if (j <= i) continue;  // each pair processed once

                ViscoelasticGhostNode<DIM>& node_j = mNodes[j];
                if (!node_j.is_active) continue;

                c_vector<double, DIM> disp = node_j.position - node_i.position;
                double dist = norm_2(disp);
                if (dist < 1e-10) continue;

                c_vector<double, DIM> unit_disp = disp / dist;

                // Per-pair rest lengths: permanent (s0) and evolving dashpot (s)
                double s_ij  = node_i.rest_lengths[nb];
                double s0_ij = node_i.initial_rest_lengths[nb];

                // Find the reciprocal index (j's entry pointing to i)
                // to keep rest lengths symmetric
                unsigned nb_recip = FindReciprocalIndex(j, i);

                // SLS split-force: F = E_0*(d - s0) + E_1*(d - s)
                double force_mag = mRelaxedStiffness   * (dist - s0_ij)
                                 + mRelaxationModulus  * (dist - s_ij);

                // Anisotropic modulation
                double aniso_factor = 1.0;
                if (mAnisotropyStrength > 0.0)
                {
                    double avg_aniso = 0.5 * (node_i.anisotropy + node_j.anisotropy);
                    if (avg_aniso > 0.0)
                    {
                        c_vector<double, DIM> avg_fibre = node_i.fibre_direction + node_j.fibre_direction;
                        double fn = norm_2(avg_fibre);
                        if (fn > 1e-10)
                        {
                            avg_fibre /= fn;
                            double cos_angle = std::abs(inner_prod(unit_disp, avg_fibre));
                            aniso_factor = 1.0 - mAnisotropyStrength * avg_aniso
                                           * (1.0 - cos_angle);
                        }
                    }
                }

                force_mag *= aniso_factor;

                c_vector<double, DIM> force_ij = force_mag * unit_disp;
                node_i.force += force_ij;
                node_j.force -= force_ij;

                // ── Evolve dashpot rest length (exact exponential integrator) ──
                // Solves ds/dt = (d - s)/tau exactly over the interval dt.
                // s^{n+1} = d + (s^n - d) * exp(-dt/tau)  [unconditionally stable]
                double new_s = dist + (s_ij - dist) * std::exp(-dt / mRelaxationTime);

                node_i.rest_lengths[nb] = new_s;
                node_j.rest_lengths[nb_recip] = new_s;
            }
        }
        
        // Update positions (overdamped EoM — separate from constitutive law)
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
     * Degrade ECM using pre-computed nearby ghost node indices.
     */
    void DegradeNearby(const c_vector<double, DIM>& cellPos,
                       const std::vector<unsigned>& nearbyIndices,
                       double cutoff, double dt)
    {
        double inv_cutoff = 1.0 / cutoff;
        for (unsigned idx : nearbyIndices)
        {
            ViscoelasticGhostNode<DIM>& gn = mNodes[idx];
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
            ViscoelasticGhostNode<DIM>& node = mNodes[i];
            if (!node.is_active) continue;

            if (node.density < mRemovalThreshold)
            {
                node.is_active = false;
                node.density = 0.0;
                removed_count++;

                // Remove i from all its neighbours' lists (and rest length arrays)
                for (unsigned nb = 0; nb < node.neighbours.size(); nb++)
                {
                    unsigned n_idx = node.neighbours[nb];
                    auto& n_list  = mNodes[n_idx].neighbours;
                    auto& r_list  = mNodes[n_idx].rest_lengths;
                    auto& r0_list = mNodes[n_idx].initial_rest_lengths;
                    for (unsigned k = 0; k < n_list.size(); k++)
                    {
                        if (n_list[k] == i)
                        {
                            n_list.erase(n_list.begin() + k);
                            r_list.erase(r_list.begin() + k);
                            r0_list.erase(r0_list.begin() + k);
                            break;
                        }
                    }
                }
                node.neighbours.clear();
                node.rest_lengths.clear();
                node.initial_rest_lengths.clear();
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
     * expansion. Neighbours are connected via spatial hash lookup, and per-pair
     * rest lengths are initialised to the actual distance.
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

        // Compute new domain bounds
        double newMin[3], newMax[3];
        for (unsigned d = 0; d < DIM; d++)
        {
            newMin[d] = expandLo[d] ? mDomainMin[d] - mInitialSpacing : mDomainMin[d];
            newMax[d] = expandHi[d] ? mDomainMax[d] + mInitialSpacing : mDomainMax[d];
        }

        // Generate candidate node positions in expansion strips only
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
            ViscoelasticGhostNode<DIM> gn;
            gn.id = first_new_id + i;
            gn.position = new_positions[i];
            gn.force = zero_vector<double>(DIM);
            gn.density = 1.0;
            gn.is_active = true;
            gn.anisotropy = 0.0;
            gn.fibre_direction = zero_vector<double>(DIM);

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

        // Expand spatial hash
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

        // Build neighbour connectivity for new nodes with per-pair rest length init
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
                            if (j >= first_new_id && j >= ni) continue;
                            c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                            double dist_sq = inner_prod(disp, disp);
                            if (dist_sq < cutoff_sq)
                            {
                                double initial_rest = std::sqrt(dist_sq);
                                mNodes[ni].neighbours.push_back(j);
                                mNodes[ni].rest_lengths.push_back(initial_rest);
                                mNodes[ni].initial_rest_lengths.push_back(initial_rest);
                                mNodes[j].neighbours.push_back(ni);
                                mNodes[j].rest_lengths.push_back(initial_rest);
                                mNodes[j].initial_rest_lengths.push_back(initial_rest);
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
                                    double initial_rest = std::sqrt(dist_sq);
                                    mNodes[ni].neighbours.push_back(j);
                                    mNodes[ni].rest_lengths.push_back(initial_rest);
                                    mNodes[ni].initial_rest_lengths.push_back(initial_rest);
                                    mNodes[j].neighbours.push_back(ni);
                                    mNodes[j].rest_lengths.push_back(initial_rest);
                                    mNodes[j].initial_rest_lengths.push_back(initial_rest);
                                }
                            }
                        }
                    }
                }
            }
        }

        UpdateActiveCount();

        std::cout << "  ViscoelasticGhostNodeECM: expanded boundary, added " << num_added
                  << " new nodes (" << mNumActive << " total active)" << std::endl;

        return num_added;
    }

    /**
     * Remodel fibre orientation using pre-computed nearby ghost node indices.
     */
    void RemodelNearby(const c_vector<double, DIM>& cellPos,
                       const std::vector<unsigned>& nearbyIndices,
                       const c_vector<double, DIM>& traction,
                       double cutoff, double dt)
    {
        double inv_cutoff = 1.0 / cutoff;
        for (unsigned idx : nearbyIndices)
        {
            ViscoelasticGhostNode<DIM>& gn = mNodes[idx];
            double dist = norm_2(gn.position - cellPos);
            double weight = std::max(0.0, 1.0 - dist * inv_cutoff);

            c_vector<double, DIM> target = traction;
            if (inner_prod(gn.fibre_direction, target) < 0)
                target = -target;

            double blend = mFibreRemodelingRate * weight * dt;
            gn.fibre_direction += blend * (target - gn.fibre_direction);
            double fn = norm_2(gn.fibre_direction);
            if (fn > 1e-10)
                gn.fibre_direction /= fn;

            gn.anisotropy = std::min(1.0, gn.anisotropy + blend);
        }
    }

    /**
     * Write ghost node field to VTP file for ParaView,
     * including per-node average rest length strain for visualisation.
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
            for (unsigned nb = 0; nb < mNodes[i].neighbours.size(); nb++)
            {
                unsigned j = mNodes[i].neighbours[nb];
                if (j > i && mNodes[j].is_active) num_edges++;
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

        // Points
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
            for (unsigned nb = 0; nb < mNodes[i].neighbours.size(); nb++)
            {
                unsigned j = mNodes[i].neighbours[nb];
                if (j > i && mNodes[j].is_active)
                    file << index_map[i] << " " << index_map[j] << "\n";
            }
        }
        file << "</DataArray>\n";
        file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        unsigned offset = 0;
        for (unsigned i = 0; i < mNodes.size(); i++)
        {
            if (!mNodes[i].is_active) continue;
            for (unsigned nb = 0; nb < mNodes[i].neighbours.size(); nb++)
            {
                unsigned j = mNodes[i].neighbours[nb];
                if (j > i && mNodes[j].is_active)
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

        // Fibre direction
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

        // Average rest length strain: mean of (s_ij - s_initial) / s_initial
        // Positive = rest length has grown (material has "relaxed" / yielded)
        // Negative = rest length has shrunk
        file << "<DataArray type=\"Float64\" Name=\"rest_length_strain\" format=\"ascii\">\n";
        for (const auto& n : mNodes)
        {
            if (!n.is_active) continue;
            double avg_strain = 0.0;
            if (!n.rest_lengths.empty())
            {
                for (unsigned nb = 0; nb < n.rest_lengths.size(); nb++)
                {
                    double s0 = n.initial_rest_lengths[nb];
                    if (s0 > 1e-10)
                        avg_strain += (n.rest_lengths[nb] - s0) / s0;
                }
                avg_strain /= static_cast<double>(n.rest_lengths.size());
            }
            file << avg_strain << "\n";
        }
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
                ViscoelasticGhostNode<DIM> gn;
                gn.id = id++;
                gn.position = zero_vector<double>(DIM);
                gn.position[0] = x;
                gn.position[1] = y;
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
                ViscoelasticGhostNode<DIM> gn;
                gn.id = id++;
                gn.position = zero_vector<double>(DIM);
                gn.position[0] = x;
                gn.position[1] = y;
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
                    ViscoelasticGhostNode<DIM> gn;
                    gn.id = id++;
                    gn.position = zero_vector<double>(DIM);
                    gn.position[0] = x;
                    gn.position[1] = y;
                    gn.position[2] = z;
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
                    node.fibre_direction = node.position / r;
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
                else
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

    /**
     * Build neighbour connectivity and initialise per-pair rest lengths
     * to the initial spacing.
     */
    void BuildNeighbourConnectivity(double cutoff)
    {
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
                            if (j <= i) continue;
                            c_vector<double, DIM> disp = mNodes[j].position - pos_i;
                            double dist_sq = inner_prod(disp, disp);
                            if (dist_sq < cutoff_sq)
                            {
                                double initial_rest = std::sqrt(dist_sq);
                                mNodes[i].neighbours.push_back(j);
                                mNodes[i].rest_lengths.push_back(initial_rest);
                                mNodes[i].initial_rest_lengths.push_back(initial_rest);
                                mNodes[j].neighbours.push_back(i);
                                mNodes[j].rest_lengths.push_back(initial_rest);
                                mNodes[j].initial_rest_lengths.push_back(initial_rest);
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
                                    double initial_rest = std::sqrt(dist_sq);
                                    mNodes[i].neighbours.push_back(j);
                                    mNodes[i].rest_lengths.push_back(initial_rest);
                                    mNodes[i].initial_rest_lengths.push_back(initial_rest);
                                    mNodes[j].neighbours.push_back(i);
                                    mNodes[j].rest_lengths.push_back(initial_rest);
                                    mNodes[j].initial_rest_lengths.push_back(initial_rest);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * Find the index in node j's neighbour list that points back to node i.
     * Required to keep per-pair rest lengths symmetric.
     */
    unsigned FindReciprocalIndex(unsigned j, unsigned i) const
    {
        const auto& n_list = mNodes[j].neighbours;
        for (unsigned k = 0; k < n_list.size(); k++)
        {
            if (n_list[k] == i) return k;
        }
        // Should never happen for a well-formed connectivity graph
        assert(false);
        return 0;
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
                    mSpatialHash[GetHashIndex(cx, cy)].push_back(i);
            }
            else // DIM == 3
            {
                int cz = static_cast<int>((mNodes[i].position[2] - mHashMin[2]) / mSpatialHashCellSize);
                if (cx >= 0 && cx < mHashN[0] && cy >= 0 && cy < mHashN[1]
                    && cz >= 0 && cz < mHashN[2])
                    mSpatialHash[GetHashIndex(cx, cy, cz)].push_back(i);
            }
        }
    }
};

#endif // VISCOELASTICGHOSTNODEECMFIELD_HPP_
