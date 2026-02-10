/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef DYNAMICECMFIELD3D_HPP_
#define DYNAMICECMFIELD3D_HPP_

#include <map>
#include <tuple>
#include <cmath>
#include <fstream>
#include "UblasVectorInclude.hpp"
#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * 3D Dynamic ECM Field
 * 
 * Extension of DynamicECMField to three dimensions. In 3D, fiber orientation
 * is represented as a unit vector (not just an angle), enabling full 3D
 * contact guidance, degradation, remodeling and deposition.
 * 
 * The ECM is stored on a structured 3D Cartesian grid. Each voxel stores:
 *   - fiber_direction: unit vector giving local collagen fiber alignment
 *   - density: scalar in [0,1] representing ECM protein concentration
 *   - anisotropy: scalar in [0,1] representing degree of fiber alignment
 * 
 * Key biological processes modeled:
 *   1. Contact guidance: cells preferentially migrate along fiber direction
 *   2. Degradation: cells secrete MMPs that break down ECM (density decreases)
 *   3. Remodeling: cell traction forces rotate fibers toward pulling direction
 *   4. Deposition: cells deposit new ECM aligned with their motion
 *   5. Diffusion: ECM properties smooth over space (collagen cross-linking)
 * 
 * References:
 *   Painter 2009; Metzcar et al. 2025
 */
class DynamicECMField3d
{
private:
    /** Grid spacing in micrometers */
    double mGridSpacing;
    
    /** Domain bounds */
    double mXMin, mXMax, mYMin, mYMax, mZMin, mZMax;
    
    /** ECM remodeling rate (how fast cells align ECM) */
    double mRemodelingRate;
    
    /** ECM degradation rate per unit time */
    double mDegradationRate;
    
    /** ECM deposition rate */
    double mDepositionRate;
    
    /** Initial ECM orientation type */
    std::string mInitialECMType;
    
    /** Grid cell structure storing ECM properties per voxel */
    struct ECMGridCell
    {
        c_vector<double, 3> fiber_direction;  ///< Unit vector for fiber alignment
        double density;                        ///< ECM density [0, 1]
        double anisotropy;                     ///< Degree of alignment [0, 1]
        
        ECMGridCell()
          : density(1.0),
            anisotropy(1.0)
        {
            fiber_direction[0] = 0.0;
            fiber_direction[1] = 0.0;
            fiber_direction[2] = 1.0;  // Default: radial (z-aligned)
        }
    };
    
    /** Grid storage: map from (i,j,k) to ECM properties */
    typedef std::tuple<int, int, int> GridIndex;
    std::map<GridIndex, ECMGridCell> mGrid;
    
    /** Grid dimensions */
    int mNx, mNy, mNz;
    
    /** Convert position to grid indices */
    GridIndex GetGridIndex(const c_vector<double, 3>& position) const
    {
        int i = (int)floor((position[0] - mXMin) / mGridSpacing);
        int j = (int)floor((position[1] - mYMin) / mGridSpacing);
        int k = (int)floor((position[2] - mZMin) / mGridSpacing);
        // Clamp to valid range
        i = std::max(0, std::min(i, mNx - 1));
        j = std::max(0, std::min(j, mNy - 1));
        k = std::max(0, std::min(k, mNz - 1));
        return std::make_tuple(i, j, k);
    }
    
    /** Generate a random unit vector on the sphere */
    c_vector<double, 3> RandomUnitVector() const
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        // Marsaglia method for uniform distribution on sphere
        double v1, v2, s;
        do {
            v1 = 2.0 * p_gen->ranf() - 1.0;
            v2 = 2.0 * p_gen->ranf() - 1.0;
            s = v1 * v1 + v2 * v2;
        } while (s >= 1.0 || s < 1e-10);
        
        c_vector<double, 3> result;
        result[0] = 2.0 * v1 * sqrt(1.0 - s);
        result[1] = 2.0 * v2 * sqrt(1.0 - s);
        result[2] = 1.0 - 2.0 * s;
        return result;
    }
    
    /** Initialize grid with chosen pattern */
    void InitializeGrid()
    {
        mNx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
        mNy = (int)ceil((mYMax - mYMin) / mGridSpacing) + 1;
        mNz = (int)ceil((mZMax - mZMin) / mGridSpacing) + 1;
        
        for (int i = 0; i < mNx; i++)
        {
            for (int j = 0; j < mNy; j++)
            {
                for (int k = 0; k < mNz; k++)
                {
                    double x = mXMin + i * mGridSpacing;
                    double y = mYMin + j * mGridSpacing;
                    double z = mZMin + k * mGridSpacing;
                    
                    ECMGridCell cell;
                    cell.density = 1.0;
                    
                    if (mInitialECMType == "random")
                    {
                        cell.fiber_direction = RandomUnitVector();
                        cell.anisotropy = 0.3;  // Low alignment for random ECM
                    }
                    else if (mInitialECMType == "radial")
                    {
                        // Fibers point radially outward from domain center
                        c_vector<double, 3> center;
                        center[0] = 0.5 * (mXMin + mXMax);
                        center[1] = 0.5 * (mYMin + mYMax);
                        center[2] = 0.5 * (mZMin + mZMax);
                        
                        c_vector<double, 3> radial;
                        radial[0] = x - center[0];
                        radial[1] = y - center[1];
                        radial[2] = z - center[2];
                        
                        double mag = norm_2(radial);
                        if (mag > 1e-10)
                        {
                            cell.fiber_direction = radial / mag;
                        }
                        else
                        {
                            cell.fiber_direction[0] = 0.0;
                            cell.fiber_direction[1] = 0.0;
                            cell.fiber_direction[2] = 1.0;
                        }
                        cell.anisotropy = 0.8;  // High alignment for radial
                    }
                    else if (mInitialECMType == "circumferential")
                    {
                        // Fibers tangential to spheres centered at domain center
                        // (perpendicular to radial)
                        c_vector<double, 3> center;
                        center[0] = 0.5 * (mXMin + mXMax);
                        center[1] = 0.5 * (mYMin + mYMax);
                        center[2] = 0.5 * (mZMin + mZMax);
                        
                        c_vector<double, 3> radial;
                        radial[0] = x - center[0];
                        radial[1] = y - center[1];
                        radial[2] = z - center[2];
                        
                        double mag = norm_2(radial);
                        if (mag > 1e-10)
                        {
                            radial /= mag;
                            // Cross product with z-axis to get tangent
                            c_vector<double, 3> z_axis;
                            z_axis[0] = 0.0; z_axis[1] = 0.0; z_axis[2] = 1.0;
                            
                            c_vector<double, 3> tangent;
                            tangent[0] = radial[1] * z_axis[2] - radial[2] * z_axis[1];
                            tangent[1] = radial[2] * z_axis[0] - radial[0] * z_axis[2];
                            tangent[2] = radial[0] * z_axis[1] - radial[1] * z_axis[0];
                            
                            double tmag = norm_2(tangent);
                            if (tmag > 1e-10)
                            {
                                cell.fiber_direction = tangent / tmag;
                            }
                        }
                        cell.anisotropy = 0.8;
                    }
                    else  // Default: aligned along z
                    {
                        cell.fiber_direction[0] = 0.0;
                        cell.fiber_direction[1] = 0.0;
                        cell.fiber_direction[2] = 1.0;
                        cell.anisotropy = 1.0;
                    }
                    
                    mGrid[std::make_tuple(i, j, k)] = cell;
                }
            }
        }
    }

public:
    /**
     * Constructor
     * 
     * @param ecmType Initial fiber pattern: "random", "radial", "circumferential", "aligned"
     * @param gridSpacing Voxel size in micrometers
     * @param xMin, xMax, yMin, yMax, zMin, zMax Domain bounds
     */
    DynamicECMField3d(std::string ecmType = "random",
                      double gridSpacing = 25.0,
                      double xMin = -100.0, double xMax = 100.0,
                      double yMin = -100.0, double yMax = 100.0,
                      double zMin = -100.0, double zMax = 100.0)
        : mGridSpacing(gridSpacing),
          mXMin(xMin), mXMax(xMax),
          mYMin(yMin), mYMax(yMax),
          mZMin(zMin), mZMax(zMax),
          mRemodelingRate(0.05),
          mDegradationRate(0.001),
          mDepositionRate(0.0005),
          mInitialECMType(ecmType),
          mNx(0), mNy(0), mNz(0)
    {
        InitializeGrid();
    }
    
    // ========================================================================
    // Query methods
    // ========================================================================
    
    /**
     * Get ECM fiber direction at a 3D position (unit vector).
     */
    c_vector<double, 3> GetFiberDirectionAt(const c_vector<double, 3>& position) const
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            return it->second.fiber_direction;
        }
        c_vector<double, 3> default_dir;
        default_dir[0] = 0.0; default_dir[1] = 0.0; default_dir[2] = 1.0;
        return default_dir;
    }
    
    /**
     * Get ECM density at a 3D position.
     */
    double GetDensityAt(const c_vector<double, 3>& position) const
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            return it->second.density;
        }
        return 1.0;
    }
    
    /**
     * Get ECM anisotropy at a 3D position.
     */
    double GetAnisotropyAt(const c_vector<double, 3>& position) const
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            return it->second.anisotropy;
        }
        return 1.0;
    }
    
    // ========================================================================
    // Cell-ECM interaction: degradation, remodeling, deposition
    // ========================================================================
    
    /**
     * Degrade ECM at the cell's location (MMP secretion).
     * 
     * rho_{n+1} = max(0, rho_n - k_deg * dt)
     * When density drops below 0.5, anisotropy decays (fiber breakdown).
     */
    void DegradeECM(const c_vector<double, 3>& position, double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it == mGrid.end()) return;
        
        it->second.density -= mDegradationRate * dt;
        it->second.density = std::max(0.0, it->second.density);
        
        // Degradation breaks fiber alignment
        if (it->second.density < 0.5)
        {
            it->second.anisotropy *= 0.99;
        }
    }
    
    /**
     * Remodel ECM due to cell traction forces.
     * 
     * Cells pull ECM fibers toward their traction direction:
     *   f_new = f_old + k_remodel * min(1, |F|/F_ref) * (f_cell - (f_cell.f_old)*f_old) * dt
     * then re-normalize.
     */
    void ApplyCellTraction(const c_vector<double, 3>& position,
                           const c_vector<double, 3>& tractionForce,
                           double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it == mGrid.end()) return;
        
        double force_magnitude = norm_2(tractionForce);
        if (force_magnitude < 1e-10) return;
        
        // Cell's traction direction (unit vector)
        c_vector<double, 3> cell_dir = tractionForce / force_magnitude;
        
        // Current fiber direction
        c_vector<double, 3>& fiber = it->second.fiber_direction;
        
        // Rotation: move fiber toward cell_dir using Rodrigues-like update
        // Project cell_dir onto plane perpendicular to current fiber
        double dot = inner_prod(cell_dir, fiber);
        c_vector<double, 3> perp = cell_dir - dot * fiber;
        double perp_mag = norm_2(perp);
        
        if (perp_mag > 1e-10)
        {
            double effective_remodeling = mRemodelingRate * std::min(1.0, force_magnitude / 10.0);
            fiber += effective_remodeling * (perp / perp_mag) * dt;
            
            // Re-normalize to unit vector
            double fmag = norm_2(fiber);
            if (fmag > 1e-10)
            {
                fiber /= fmag;
            }
        }
        
        // Traction increases anisotropy
        it->second.anisotropy += 0.01 * dt;
        it->second.anisotropy = std::min(1.0, it->second.anisotropy);
    }
    
    /**
     * Deposit new ECM at the cell's location.
     * 
     * Increases density and blends deposited fiber direction with existing.
     */
    void DepositECM(const c_vector<double, 3>& position,
                    const c_vector<double, 3>& depositionDirection,
                    double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it == mGrid.end()) return;
        
        double old_density = it->second.density;
        it->second.density += mDepositionRate * dt;
        it->second.density = std::min(1.0, it->second.density);
        
        // Blend new fiber direction (weighted by relative deposit amount)
        if (old_density < 0.9)
        {
            double weight = mDepositionRate * dt / (old_density + mDepositionRate * dt);
            c_vector<double, 3>& fiber = it->second.fiber_direction;
            fiber = (1.0 - weight) * fiber + weight * depositionDirection;
            
            double fmag = norm_2(fiber);
            if (fmag > 1e-10)
            {
                fiber /= fmag;
            }
        }
    }
    
    /**
     * Diffuse ECM properties spatially (6-neighbor stencil in 3D).
     */
    void DiffuseECM(double dt)
    {
        std::map<GridIndex, ECMGridCell> new_grid = mGrid;
        double diffusion_coeff = 0.01;
        
        for (auto& entry : mGrid)
        {
            int i = std::get<0>(entry.first);
            int j = std::get<1>(entry.first);
            int k = std::get<2>(entry.first);
            
            // Average with 6 neighbors
            c_vector<double, 3> avg_dir = entry.second.fiber_direction;
            double avg_density = entry.second.density;
            int count = 1;
            
            std::vector<GridIndex> neighbors = {
                std::make_tuple(i-1,j,k), std::make_tuple(i+1,j,k),
                std::make_tuple(i,j-1,k), std::make_tuple(i,j+1,k),
                std::make_tuple(i,j,k-1), std::make_tuple(i,j,k+1)
            };
            
            for (auto& n_idx : neighbors)
            {
                auto nit = mGrid.find(n_idx);
                if (nit != mGrid.end())
                {
                    avg_dir += nit->second.fiber_direction;
                    avg_density += nit->second.density;
                    count++;
                }
            }
            
            avg_dir /= (double)count;
            avg_density /= (double)count;
            
            // Diffuse density
            new_grid[entry.first].density = entry.second.density +
                diffusion_coeff * dt * (avg_density - entry.second.density);
            
            // Diffuse fiber direction (and re-normalize)
            c_vector<double, 3> new_dir = entry.second.fiber_direction +
                diffusion_coeff * dt * (avg_dir - entry.second.fiber_direction);
            double nmag = norm_2(new_dir);
            if (nmag > 1e-10)
            {
                new_grid[entry.first].fiber_direction = new_dir / nmag;
            }
        }
        
        mGrid = new_grid;
    }
    
    // ========================================================================
    // Setters for initialization
    // ========================================================================
    
    /** Set fiber direction at a specific position */
    void SetFiberDirection(const c_vector<double, 3>& position,
                           const c_vector<double, 3>& direction)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            double mag = norm_2(direction);
            if (mag > 1e-10)
            {
                it->second.fiber_direction = direction / mag;
            }
        }
    }
    
    /** Set density at a specific position */
    void SetDensity(const c_vector<double, 3>& position, double density)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            it->second.density = std::max(0.0, std::min(1.0, density));
        }
    }
    
    /** Parameter setters */
    void SetRemodelingRate(double rate) { mRemodelingRate = rate; }
    void SetDegradationRate(double rate) { mDegradationRate = rate; }
    void SetDepositionRate(double rate) { mDepositionRate = rate; }
    
    /** Getters for grid metadata */
    double GetGridSpacing() const { return mGridSpacing; }
    int GetNx() const { return mNx; }
    int GetNy() const { return mNy; }
    int GetNz() const { return mNz; }
    double GetXMin() const { return mXMin; }
    double GetYMin() const { return mYMin; }
    double GetZMin() const { return mZMin; }
    double GetXMax() const { return mXMax; }
    double GetYMax() const { return mYMax; }
    double GetZMax() const { return mZMax; }
    
    // ========================================================================
    // VTK output for ParaView visualization
    // ========================================================================
    
    /**
     * Write ECM grid to VTI (VTK ImageData XML) file.
     */
    void WriteToVTI(const std::string& filename, double time) const
    {
        std::ofstream vtifile(filename);
        if (!vtifile.is_open())
        {
            EXCEPTION("Could not open file: " + filename);
        }
        
        vtifile << "<?xml version=\"1.0\"?>\n";
        vtifile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        vtifile << "  <ImageData WholeExtent=\"0 " << (mNx-1) << " 0 " << (mNy-1)
                << " 0 " << (mNz-1) << "\" ";
        vtifile << "Origin=\"" << mXMin << " " << mYMin << " " << mZMin << "\" ";
        vtifile << "Spacing=\"" << mGridSpacing << " " << mGridSpacing << " " << mGridSpacing << "\">\n";
        vtifile << "    <Piece Extent=\"0 " << (mNx-1) << " 0 " << (mNy-1)
                << " 0 " << (mNz-1) << "\">\n";
        vtifile << "      <PointData Vectors=\"ecm_fiber_direction\">\n";
        
        // ECM density
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_density\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int k = 0; k < mNz; k++)
            for (int j = 0; j < mNy; j++)
                for (int i = 0; i < mNx; i++)
                {
                    auto it = mGrid.find(std::make_tuple(i, j, k));
                    vtifile << (it != mGrid.end() ? it->second.density : 0.0) << " ";
                }
        vtifile << "\n        </DataArray>\n";
        
        // ECM anisotropy
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_anisotropy\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int k = 0; k < mNz; k++)
            for (int j = 0; j < mNy; j++)
                for (int i = 0; i < mNx; i++)
                {
                    auto it = mGrid.find(std::make_tuple(i, j, k));
                    vtifile << (it != mGrid.end() ? it->second.anisotropy : 0.0) << " ";
                }
        vtifile << "\n        </DataArray>\n";
        
        // ECM fiber direction (3-component vector)
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_fiber_direction\" "
                << "NumberOfComponents=\"3\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int k = 0; k < mNz; k++)
            for (int j = 0; j < mNy; j++)
                for (int i = 0; i < mNx; i++)
                {
                    auto it = mGrid.find(std::make_tuple(i, j, k));
                    if (it != mGrid.end())
                    {
                        const auto& fd = it->second.fiber_direction;
                        vtifile << fd[0] << " " << fd[1] << " " << fd[2] << " ";
                    }
                    else
                    {
                        vtifile << "0 0 0 ";
                    }
                }
        vtifile << "\n        </DataArray>\n";
        
        vtifile << "      </PointData>\n";
        vtifile << "    </Piece>\n";
        vtifile << "  </ImageData>\n";
        vtifile << "</VTKFile>\n";
        
        vtifile.close();
    }
};

#endif /* DYNAMICECMFIELD3D_HPP_ */
