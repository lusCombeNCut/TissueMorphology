/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

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

#ifndef DYNAMICECMFIELD_HPP_
#define DYNAMICECMFIELD_HPP_

#include <map>
#include <cmath>
#include <fstream>
#include "UblasVectorInclude.hpp"
#include "Exception.hpp"

/**
 * Dynamic ECM Field
 * 
 * Stores ECM properties on a spatial grid and allows cell-mediated remodeling.
 * 
 * Key features:
 * - Grid-based storage of ECM fiber orientation and density
 * - Cells mechanically align ECM through traction forces
 * - Cells degrade ECM as they migrate
 * - Cells deposit new ECM
 * - Interpolation for smooth field queries
 */
class DynamicECMField
{
private:
    /** Grid spacing in micrometers */
    double mGridSpacing;
    
    /** Domain bounds */
    double mXMin, mXMax, mYMin, mYMax;
    
    /** ECM remodeling rate (how fast cells align ECM) */
    double mRemodelingRate;
    
    /** ECM degradation rate per unit time */
    double mDegradationRate;
    
    /** ECM deposition rate */
    double mDepositionRate;
    
    /** Initial ECM orientation type */
    std::string mInitialECMType;
    
    /** Grid cell structure */
    struct ECMGridCell
    {
        double fiber_angle;      // Orientation in radians [0, 2π)
        double density;          // ECM density [0, 1], 1 = maximum
        double anisotropy;       // Degree of alignment [0, 1], 1 = perfectly aligned
        
        ECMGridCell() : fiber_angle(0.0), density(1.0), anisotropy(1.0) {}
    };
    
    /** Grid storage: map from (i,j) to ECM properties */
    std::map<std::pair<int, int>, ECMGridCell> mGrid;
    
    /** Convert position to grid indices */
    std::pair<int, int> GetGridIndex(c_vector<double, 2> position) const
    {
        int i = (int)floor((position[0] - mXMin) / mGridSpacing);
        int j = (int)floor((position[1] - mYMin) / mGridSpacing);
        return std::make_pair(i, j);
    }
    
    /** Initialize grid with pattern */
    void InitializeGrid()
    {
        int nx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
        int ny = (int)ceil((mYMax - mYMin) / mGridSpacing) + 1;
        
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                double x = mXMin + i * mGridSpacing;
                // double y = mYMin + j * mGridSpacing;  // Not used yet
                
                ECMGridCell cell;
                cell.density = 1.0;
                cell.anisotropy = 1.0;
                
                // Set initial orientation based on type
                if (mInitialECMType == "random")
                {
                    // Random will be set per-cell later
                    cell.fiber_angle = 0.0;
                    cell.anisotropy = 0.3;  // Less aligned for random
                }
                else if (mInitialECMType == "parallel")
                {
                    cell.fiber_angle = 0.0;  // Horizontal
                }
                else if (mInitialECMType == "perpendicular")
                {
                    cell.fiber_angle = M_PI / 2.0;  // Vertical
                }
                else if (mInitialECMType == "mixed")
                {
                    double center_start = (mXMax - mXMin) / 3.0 + mXMin;
                    double center_end = 2.0 * (mXMax - mXMin) / 3.0 + mXMin;
                    
                    if (x > center_start && x < center_end)
                    {
                        cell.fiber_angle = M_PI / 2.0;  // Perpendicular in center
                    }
                    else
                    {
                        cell.fiber_angle = 0.0;  // Parallel on sides
                    }
                }
                
                mGrid[std::make_pair(i, j)] = cell;
            }
        }
    }
    
    /** Calculate angle difference handling wraparound */
    double AngleDifference(double angle1, double angle2) const
    {
        double diff = angle1 - angle2;
        
        // Normalize to [-π, π]
        while (diff > M_PI) diff -= 2.0 * M_PI;
        while (diff < -M_PI) diff += 2.0 * M_PI;
        
        return diff;
    }

public:
    /**
     * Constructor
     */
    DynamicECMField(std::string ecmType = "random",
                   double gridSpacing = 25.0,
                   double xMin = 0.0, double xMax = 600.0,
                   double yMin = 0.0, double yMax = 1000.0)
        : mGridSpacing(gridSpacing),
          mXMin(xMin), mXMax(xMax),
          mYMin(yMin), mYMax(yMax),
          mRemodelingRate(0.05),      // Moderate remodeling
          mDegradationRate(0.001),    // Slow degradation
          mDepositionRate(0.0005),    // Slower deposition
          mInitialECMType(ecmType)
    {
        InitializeGrid();
    }
    
    /**
     * Get ECM fiber orientation at a position
     * Uses bilinear interpolation for smooth field
     */
    double GetFiberAngleAt(c_vector<double, 2> position) const
    {
        auto idx = GetGridIndex(position);
        
        // Simple nearest-neighbor for now (could upgrade to bilinear)
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            return it->second.fiber_angle;
        }
        
        return 0.0;  // Default
    }
    
    /**
     * Get ECM density at a position
     */
    double GetDensityAt(c_vector<double, 2> position) const
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
     * Get ECM anisotropy at a position
     */
    double GetAnisotropyAt(c_vector<double, 2> position) const
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        if (it != mGrid.end())
        {
            return it->second.anisotropy;
        }
        return 1.0;
    }
    
    /**
     * Update ECM based on cell traction
     * 
     * Cells mechanically remodel ECM by pulling fibers in their direction of motion
     * 
     * @param position Cell position
     * @param tractionForce Cell's traction force vector
     * @param dt Time step
     */
    void ApplyCellTraction(c_vector<double, 2> position, 
                          c_vector<double, 2> tractionForce,
                          double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        
        if (it == mGrid.end())
            return;
        
        // Calculate cell's pulling direction
        double force_magnitude = norm_2(tractionForce);
        if (force_magnitude < 1e-10)
            return;
        
        double cell_angle = atan2(tractionForce[1], tractionForce[0]);
        
        // Gradually align ECM toward cell's pulling direction
        double current_angle = it->second.fiber_angle;
        double angle_diff = AngleDifference(cell_angle, current_angle);
        
        // Stronger forces cause faster remodeling
        double effective_remodeling = mRemodelingRate * std::min(1.0, force_magnitude / 10.0);
        
        it->second.fiber_angle += effective_remodeling * angle_diff * dt;
        
        // Normalize angle to [0, 2π)
        while (it->second.fiber_angle < 0.0) it->second.fiber_angle += 2.0 * M_PI;
        while (it->second.fiber_angle >= 2.0 * M_PI) it->second.fiber_angle -= 2.0 * M_PI;
        
        // Cell traction increases anisotropy (creates aligned fibers)
        it->second.anisotropy += 0.01 * dt;
        it->second.anisotropy = std::min(1.0, it->second.anisotropy);
    }
    
    /**
     * Apply ECM degradation
     * 
     * Cells secrete MMPs (matrix metalloproteinases) that degrade ECM
     * 
     * @param position Cell position
     * @param dt Time step
     */
    void DegradeECM(c_vector<double, 2> position, double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        
        if (it == mGrid.end())
            return;
        
        // Reduce density
        it->second.density -= mDegradationRate * dt;
        it->second.density = std::max(0.0, it->second.density);
        
        // Degradation reduces anisotropy (breaks alignment)
        if (it->second.density < 0.5)
        {
            it->second.anisotropy *= 0.99;  // Gradual loss of alignment
        }
    }
    
    /**
     * Apply ECM deposition
     * 
     * Cells secrete ECM components (collagen, fibronectin, etc.)
     * 
     * @param position Cell position
     * @param depositionAngle Angle at which to deposit ECM
     * @param dt Time step
     */
    void DepositECM(c_vector<double, 2> position, double depositionAngle, double dt)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        
        if (it == mGrid.end())
            return;
        
        // Increase density
        double old_density = it->second.density;
        it->second.density += mDepositionRate * dt;
        it->second.density = std::min(1.0, it->second.density);
        
        // Newly deposited ECM aligned with cell's orientation
        if (old_density < 0.9)  // Only if not already saturated
        {
            double current_angle = it->second.fiber_angle;
            double angle_diff = AngleDifference(depositionAngle, current_angle);
            
            // Weighted average based on relative amounts
            double weight = mDepositionRate * dt / (old_density + mDepositionRate * dt);
            it->second.fiber_angle += weight * angle_diff;
        }
    }
    
    /**
     * Set ECM fiber angle directly (for initialization)
     */
    void SetFiberAngle(c_vector<double, 2> position, double angle)
    {
        auto idx = GetGridIndex(position);
        auto it = mGrid.find(idx);
        
        if (it != mGrid.end())
        {
            it->second.fiber_angle = angle;
        }
    }
    
    /**
     * Diffuse ECM properties to neighbors (smoothing)
     * 
     * @param dt Time step
     */
    void DiffuseECM(double dt)
    {
        // Create temporary copy for simultaneous update
        std::map<std::pair<int, int>, ECMGridCell> new_grid = mGrid;
        
        double diffusion_coeff = 0.01;  // Small diffusion
        
        for (auto& entry : mGrid)
        {
            auto idx = entry.first;
            int i = idx.first;
            int j = idx.second;
            
            // Average with neighbors (simple 4-neighbor stencil)
            double avg_cos = cos(entry.second.fiber_angle);
            double avg_sin = sin(entry.second.fiber_angle);
            double avg_density = entry.second.density;
            int count = 1;
            
            std::vector<std::pair<int,int>> neighbors = {
                {i-1, j}, {i+1, j}, {i, j-1}, {i, j+1}
            };
            
            for (auto& n_idx : neighbors)
            {
                auto it = mGrid.find(n_idx);
                if (it != mGrid.end())
                {
                    avg_cos += cos(it->second.fiber_angle);
                    avg_sin += sin(it->second.fiber_angle);
                    avg_density += it->second.density;
                    count++;
                }
            }
            
            avg_cos /= count;
            avg_sin /= count;
            avg_density /= count;
            
            // Update with diffused values
            double avg_angle = atan2(avg_sin, avg_cos);
            
            new_grid[idx].fiber_angle = entry.second.fiber_angle + 
                                       diffusion_coeff * dt * AngleDifference(avg_angle, entry.second.fiber_angle);
            new_grid[idx].density = entry.second.density + 
                                   diffusion_coeff * dt * (avg_density - entry.second.density);
        }
        
        mGrid = new_grid;
    }
    
    /**
     * Write ECM grid to VTI (XML ImageData) file for visualization with PVD
     */
    void WriteToVTI(const std::string& filename, double time)
    {
        std::ofstream vtifile;
        vtifile.open(filename);
        
        if (!vtifile.is_open())
        {
            EXCEPTION("Could not open file: " + filename);
        }
        
        // Grid dimensions
        int nx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
        int ny = (int)ceil((mYMax - mYMin) / mGridSpacing) + 1;
        
        // VTI header (XML format)
        vtifile << "<?xml version=\"1.0\"?>\n";
        vtifile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        vtifile << "  <ImageData WholeExtent=\"0 " << (nx-1) << " 0 " << (ny-1) << " 0 0\" ";
        vtifile << "Origin=\"" << mXMin << " " << mYMin << " 0\" ";
        vtifile << "Spacing=\"" << mGridSpacing << " " << mGridSpacing << " 1\">\n";
        vtifile << "    <Piece Extent=\"0 " << (nx-1) << " 0 " << (ny-1) << " 0 0\">\n";
        vtifile << "      <PointData Vectors=\"ecm_orientation\">\n";
        
        // ECM density (scalar)
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_density\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    vtifile << it->second.density << " ";
                }
                else
                {
                    vtifile << "0.0 ";
                }
            }
        }
        vtifile << "\n        </DataArray>\n";
        
        // ECM anisotropy (scalar)
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_anisotropy\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    vtifile << it->second.anisotropy << " ";
                }
                else
                {
                    vtifile << "0.0 ";
                }
            }
        }
        vtifile << "\n        </DataArray>\n";
        
        // ECM orientation (vector - 3 components)
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_orientation\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        vtifile << "          ";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    double angle = it->second.fiber_angle;
                    vtifile << cos(angle) << " " << sin(angle) << " 0.0 ";
                }
                else
                {
                    vtifile << "0.0 0.0 0.0 ";
                }
            }
        }
        vtifile << "\n        </DataArray>\n";
        
        // Close tags
        vtifile << "      </PointData>\n";
        vtifile << "    </Piece>\n";
        vtifile << "  </ImageData>\n";
        vtifile << "</VTKFile>\n";
        
        vtifile.close();
    }
    
    /**
     * Write ECM grid to VTK file for visualization (legacy format)
     */
    void WriteToVTK(const std::string& filename, double time)
    {
        std::ofstream vtkfile;
        vtkfile.open(filename);
        
        if (!vtkfile.is_open())
        {
            EXCEPTION("Could not open file: " + filename);
        }
        
        // VTK header
        vtkfile << "# vtk DataFile Version 3.0\n";
        vtkfile << "ECM Grid at time " << time << "\n";
        vtkfile << "ASCII\n";
        vtkfile << "DATASET STRUCTURED_POINTS\n";
        
        // Grid dimensions
        int nx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
        int ny = (int)ceil((mYMax - mYMin) / mGridSpacing) + 1;
        vtkfile << "DIMENSIONS " << nx << " " << ny << " 1\n";
        vtkfile << "ORIGIN " << mXMin << " " << mYMin << " 0.0\n";
        vtkfile << "SPACING " << mGridSpacing << " " << mGridSpacing << " 1.0\n";
        
        // Point data
        vtkfile << "POINT_DATA " << (nx * ny) << "\n";
        
        // ECM density
        vtkfile << "SCALARS ecm_density float 1\n";
        vtkfile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    vtkfile << it->second.density << "\n";
                }
                else
                {
                    vtkfile << "0.0\n";
                }
            }
        }
        
        // ECM anisotropy
        vtkfile << "SCALARS ecm_anisotropy float 1\n";
        vtkfile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    vtkfile << it->second.anisotropy << "\n";
                }
                else
                {
                    vtkfile << "0.0\n";
                }
            }
        }
        
        // ECM fiber orientation (as vectors)
        vtkfile << "VECTORS ecm_orientation float\n";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                auto it = mGrid.find(std::make_pair(i, j));
                if (it != mGrid.end())
                {
                    double angle = it->second.fiber_angle;
                    vtkfile << cos(angle) << " " << sin(angle) << " 0.0\n";
                }
                else
                {
                    vtkfile << "0.0 0.0 0.0\n";
                }
            }
        }
        
        vtkfile.close();
    }
    
    /** Setters for parameters */
    void SetRemodelingRate(double rate) { mRemodelingRate = rate; }
    void SetDegradationRate(double rate) { mDegradationRate = rate; }
    void SetDepositionRate(double rate) { mDepositionRate = rate; }
};

#endif /* DYNAMICECMFIELD_HPP_ */
