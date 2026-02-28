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

#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <sstream>
#include "UblasVectorInclude.hpp"
#include "Exception.hpp"

/**
 * Dynamic ECM Field
 * 
 * Stores ECM properties on a spatial grid and allows cell-mediated remodeling.
 * 
 * Performance note: the grid is stored as a flat std::vector with row-major
 * indexing (index = j * mNx + i) for O(1) access and cache-friendly iteration.
 * A second buffer vector is used for double-buffered diffusion (swap instead
 * of copy).
 * 
 * Key features:
 * - Grid-based storage of ECM fiber orientation and density
 * - Cells mechanically align ECM through traction forces
 * - Cells degrade ECM as they migrate
 * - Cells deposit new ECM
 * - Supports square and hexagonal grid layouts
 */
class DynamicECMField
{
private:
    /** Grid spacing in micrometers (distance between hex/square cell centres in a row) */
    double mGridSpacing;
    
    /** Row spacing for hex grids (= mGridSpacing * sqrt(3)/2); equals mGridSpacing for square */
    double mRowSpacing;
    
    /** Domain bounds */
    double mXMin, mXMax, mYMin, mYMax;
    
    /** Grid dimensions */
    int mNx, mNy;
    
    /** Grid type: "square" or "hex" */
    std::string mGridType;
    
    /** Fast bool flag for grid type (avoids string compare in hot path) */
    bool mIsHex;
    
    /** ECM remodeling rate (how fast cells align ECM) */
    double mRemodelingRate;
    
    /** ECM degradation rate per unit time */
    double mDegradationRate;
    
    /** ECM deposition rate */
    double mDepositionRate;
    
    /** Diffusion coefficient for density smoothing */
    double mDiffusionCoeff;
    
    /** Initial ECM orientation type */
    std::string mInitialECMType;
    
    /** Grid cell structure */
    struct ECMGridCell
    {
        double fiber_angle;      // Orientation in radians [0, 2π)
        double density;          // ECM density [0, 1], 1 = maximum
        double anisotropy;       // Degree of alignment [0, 1], 1 = perfectly aligned
        double cos_angle;        // Cached cos(fiber_angle) for diffusion & output
        double sin_angle;        // Cached sin(fiber_angle) for diffusion & output
        
        ECMGridCell() : fiber_angle(0.0), density(1.0), anisotropy(1.0),
                        cos_angle(1.0), sin_angle(0.0) {}
        
        /** Update cached trig values. Call after changing fiber_angle. */
        void UpdateTrigCache()
        {
            cos_angle = std::cos(fiber_angle);
            sin_angle = std::sin(fiber_angle);
        }
    };
    
    /** Grid storage: flat row-major vector, index = j * mNx + i */
    std::vector<ECMGridCell> mGrid;
    
    /** Second buffer for double-buffered diffusion */
    std::vector<ECMGridCell> mDiffuseBuffer;
    
    // ---- inline helpers ----
    
    /** Flat row-major index from (i,j). No bounds check. */
    inline int FlatIndex(int i, int j) const
    {
        return j * mNx + i;
    }
    
    /** Is (i,j) inside the grid? */
    inline bool InBounds(int i, int j) const
    {
        return i >= 0 && i < mNx && j >= 0 && j < mNy;
    }
    
    /** Convert position to grid (i,j), clamped to valid range. Returns flat index. */
    int GetGridIndexFlat(c_vector<double, 2> position) const
    {
        int i, j;
        if (mIsHex)
        {
            j = (int)floor((position[1] - mYMin) / mRowSpacing);
            j = std::max(0, std::min(j, mNy - 1));
            double xOffset = (j % 2 == 1) ? mGridSpacing * 0.5 : 0.0;
            i = (int)floor((position[0] - mXMin - xOffset) / mGridSpacing);
            i = std::max(0, std::min(i, mNx - 1));
        }
        else
        {
            i = (int)floor((position[0] - mXMin) / mGridSpacing);
            j = (int)floor((position[1] - mYMin) / mGridSpacing);
            i = std::max(0, std::min(i, mNx - 1));
            j = std::max(0, std::min(j, mNy - 1));
        }
        return FlatIndex(i, j);
    }
    
    /**
     * Get the world-space position of a grid cell centre.
     */
    c_vector<double, 2> GetCellPosition(int i, int j) const
    {
        c_vector<double, 2> pos;
        if (mIsHex)
        {
            double xOffset = (j % 2 == 1) ? mGridSpacing * 0.5 : 0.0;
            pos[0] = mXMin + i * mGridSpacing + xOffset;
            pos[1] = mYMin + j * mRowSpacing;
        }
        else
        {
            pos[0] = mXMin + i * mGridSpacing;
            pos[1] = mYMin + j * mGridSpacing;
        }
        return pos;
    }
    
    /** Initialize grid with pattern */
    void InitializeGrid()
    {
        mIsHex = (mGridType == "hex");
        
        if (mIsHex)
        {
            mRowSpacing = mGridSpacing * std::sqrt(3.0) / 2.0;
            mNx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
            mNy = (int)ceil((mYMax - mYMin) / mRowSpacing) + 1;
        }
        else
        {
            mRowSpacing = mGridSpacing;
            mNx = (int)ceil((mXMax - mXMin) / mGridSpacing) + 1;
            mNy = (int)ceil((mYMax - mYMin) / mGridSpacing) + 1;
        }
        
        // Allocate flat grids
        mGrid.resize(mNx * mNy);
        mDiffuseBuffer.resize(mNx * mNy);
        
        for (int j = 0; j < mNy; j++)
        {
            for (int i = 0; i < mNx; i++)
            {
                c_vector<double, 2> cellPos = GetCellPosition(i, j);
                double x = cellPos[0];
                
                ECMGridCell cell;
                cell.density = 1.0;
                cell.anisotropy = 1.0;
                
                if (mInitialECMType == "random")
                {
                    cell.fiber_angle = 0.0;
                    cell.anisotropy = 0.3;
                }
                else if (mInitialECMType == "parallel")
                {
                    cell.fiber_angle = 0.0;
                }
                else if (mInitialECMType == "perpendicular")
                {
                    cell.fiber_angle = M_PI / 2.0;
                }
                else if (mInitialECMType == "mixed")
                {
                    double center_start = (mXMax - mXMin) / 3.0 + mXMin;
                    double center_end = 2.0 * (mXMax - mXMin) / 3.0 + mXMin;
                    
                    if (x > center_start && x < center_end)
                    {
                        cell.fiber_angle = M_PI / 2.0;
                    }
                    else
                    {
                        cell.fiber_angle = 0.0;
                    }
                }
                
                cell.UpdateTrigCache();
                mGrid[FlatIndex(i, j)] = cell;
            }
        }
    }
    
    /** Calculate angle difference handling wraparound */
    double AngleDifference(double angle1, double angle2) const
    {
        double diff = angle1 - angle2;
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
                   double yMin = 0.0, double yMax = 1000.0,
                   std::string gridType = "square")
        : mGridSpacing(gridSpacing),
          mRowSpacing(gridSpacing),
          mXMin(xMin), mXMax(xMax),
          mYMin(yMin), mYMax(yMax),
          mNx(0), mNy(0),
          mGridType(gridType),
          mIsHex(gridType == "hex"),
          mRemodelingRate(0.05),
          mDegradationRate(0.001),
          mDepositionRate(0.0005),
          mDiffusionCoeff(0.01),
          mInitialECMType(ecmType)
    {
        InitializeGrid();
    }
    
    /**
     * Get ECM fiber orientation at a position
     */
    double GetFiberAngleAt(c_vector<double, 2> position) const
    {
        return mGrid[GetGridIndexFlat(position)].fiber_angle;
    }
    
    /**
     * Get ECM density at a position
     */
    double GetDensityAt(c_vector<double, 2> position) const
    {
        return mGrid[GetGridIndexFlat(position)].density;
    }
    
    /**
     * Get ECM anisotropy at a position
     */
    double GetAnisotropyAt(c_vector<double, 2> position) const
    {
        return mGrid[GetGridIndexFlat(position)].anisotropy;
    }
    
    /**
     * Update ECM based on cell traction
     */
    void ApplyCellTraction(c_vector<double, 2> position, 
                          c_vector<double, 2> tractionForce,
                          double dt)
    {
        double force_magnitude = norm_2(tractionForce);
        if (force_magnitude < 1e-10)
            return;
        
        ECMGridCell& cell = mGrid[GetGridIndexFlat(position)];
        
        double cell_angle = atan2(tractionForce[1], tractionForce[0]);
        double angle_diff = AngleDifference(cell_angle, cell.fiber_angle);
        double effective_remodeling = mRemodelingRate * std::min(1.0, force_magnitude / 10.0);
        
        cell.fiber_angle += effective_remodeling * angle_diff * dt;
        
        // Normalize angle to [0, 2π)
        while (cell.fiber_angle < 0.0) cell.fiber_angle += 2.0 * M_PI;
        while (cell.fiber_angle >= 2.0 * M_PI) cell.fiber_angle -= 2.0 * M_PI;
        
        cell.UpdateTrigCache();
        cell.anisotropy = std::min(1.0, cell.anisotropy + 0.01 * dt);
    }
    
    /**
     * Apply ECM degradation
     */
    void DegradeECM(c_vector<double, 2> position, double dt)
    {
        ECMGridCell& cell = mGrid[GetGridIndexFlat(position)];
        
        cell.density -= mDegradationRate * dt;
        cell.density = std::max(0.0, cell.density);
        
        if (cell.density < 0.5)
        {
            cell.anisotropy *= 0.99;
        }
    }
    
    /**
     * Apply ECM deposition
     */
    void DepositECM(c_vector<double, 2> position, double depositionAngle, double dt)
    {
        ECMGridCell& cell = mGrid[GetGridIndexFlat(position)];
        
        double old_density = cell.density;
        cell.density = std::min(1.0, cell.density + mDepositionRate * dt);
        
        if (old_density < 0.9)
        {
            double angle_diff = AngleDifference(depositionAngle, cell.fiber_angle);
            double weight = mDepositionRate * dt / (old_density + mDepositionRate * dt);
            cell.fiber_angle += weight * angle_diff;
            cell.UpdateTrigCache();
        }
    }
    
    /**
     * Set ECM fiber angle directly (for initialization)
     */
    void SetFiberAngle(c_vector<double, 2> position, double angle)
    {
        ECMGridCell& cell = mGrid[GetGridIndexFlat(position)];
        cell.fiber_angle = angle;
        cell.UpdateTrigCache();
    }
    
    /**
     * Diffuse ECM properties to neighbours (smoothing).
     * Uses double-buffered swap instead of copying the entire grid.
     * Neighbour stencil: 4 for square, 6 for hex (odd-r offset).
     * 
     * @param dt Time step
     */
    void DiffuseECM(double dt)
    {
        for (int j = 0; j < mNy; j++)
        {
            for (int i = 0; i < mNx; i++)
            {
                int idx = FlatIndex(i, j);
                const ECMGridCell& self = mGrid[idx];
                
                // Use cached cos/sin — no trig calls per cell
                double avg_cos = self.cos_angle;
                double avg_sin = self.sin_angle;
                double avg_density = self.density;
                int count = 1;
                
                // --- accumulate neighbours inline ---
                // Left / right (always present for both grid types)
                if (i > 0)
                {
                    const ECMGridCell& n = mGrid[idx - 1];
                    avg_cos += n.cos_angle;
                    avg_sin += n.sin_angle;
                    avg_density += n.density;
                    count++;
                }
                if (i < mNx - 1)
                {
                    const ECMGridCell& n = mGrid[idx + 1];
                    avg_cos += n.cos_angle;
                    avg_sin += n.sin_angle;
                    avg_density += n.density;
                    count++;
                }
                
                if (mIsHex)
                {
                    // Hex 6-neighbour (odd-row offset)
                    if (j % 2 == 0)
                    {
                        if (j > 0)
                        {
                            if (i > 0) { const ECMGridCell& n = mGrid[FlatIndex(i-1,j-1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                            { const ECMGridCell& n = mGrid[FlatIndex(i,j-1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                        }
                        if (j < mNy - 1)
                        {
                            if (i > 0) { const ECMGridCell& n = mGrid[FlatIndex(i-1,j+1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                            { const ECMGridCell& n = mGrid[FlatIndex(i,j+1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                        }
                    }
                    else
                    {
                        if (j > 0)
                        {
                            { const ECMGridCell& n = mGrid[FlatIndex(i,j-1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                            if (i < mNx-1) { const ECMGridCell& n = mGrid[FlatIndex(i+1,j-1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                        }
                        if (j < mNy - 1)
                        {
                            { const ECMGridCell& n = mGrid[FlatIndex(i,j+1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                            if (i < mNx-1) { const ECMGridCell& n = mGrid[FlatIndex(i+1,j+1)]; avg_cos+=n.cos_angle; avg_sin+=n.sin_angle; avg_density+=n.density; count++; }
                        }
                    }
                }
                else
                {
                    // Square: up / down
                    if (j > 0)
                    {
                        const ECMGridCell& n = mGrid[idx - mNx];
                        avg_cos += n.cos_angle;
                        avg_sin += n.sin_angle;
                        avg_density += n.density;
                        count++;
                    }
                    if (j < mNy - 1)
                    {
                        const ECMGridCell& n = mGrid[idx + mNx];
                        avg_cos += n.cos_angle;
                        avg_sin += n.sin_angle;
                        avg_density += n.density;
                        count++;
                    }
                }
                
                double inv_count = 1.0 / (double)count;
                avg_cos *= inv_count;
                avg_sin *= inv_count;
                avg_density *= inv_count;
                
                double avg_angle = atan2(avg_sin, avg_cos);
                
                ECMGridCell& dest = mDiffuseBuffer[idx];
                dest.fiber_angle = self.fiber_angle +
                    mDiffusionCoeff * dt * AngleDifference(avg_angle, self.fiber_angle);
                dest.density = self.density +
                    mDiffusionCoeff * dt * (avg_density - self.density);
                dest.anisotropy = self.anisotropy;  // anisotropy not diffused
                dest.UpdateTrigCache();  // refresh cache for next step
            }
        }
        
        // Swap buffers (O(1) pointer swap, no copy)
        mGrid.swap(mDiffuseBuffer);
    }
    
    // ========================================================================
    // VTK output
    // ========================================================================
    
    /**
     * Write ECM grid to VTI (XML ImageData) file for visualization with PVD.
     * Only valid for square grids.
     */
    void WriteToVTI(const std::string& filename, double time)
    {
        std::ofstream vtifile;
        vtifile.open(filename);
        
        if (!vtifile.is_open())
        {
            EXCEPTION("Could not open file: " + filename);
        }
        
        vtifile << "<?xml version=\"1.0\"?>\n";
        vtifile << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
        vtifile << "  <ImageData WholeExtent=\"0 " << (mNx-1) << " 0 " << (mNy-1) << " 0 0\" ";
        vtifile << "Origin=\"" << mXMin << " " << mYMin << " 0\" ";
        vtifile << "Spacing=\"" << mGridSpacing << " " << mGridSpacing << " 1\">\n";
        vtifile << "    <Piece Extent=\"0 " << (mNx-1) << " 0 " << (mNy-1) << " 0 0\">\n";
        vtifile << "      <PointData Vectors=\"ecm_orientation\">\n";
        
        // ECM density
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_density\" format=\"ascii\">\n";
        for (int j = 0; j < mNy; j++)
        {
            vtifile << "          ";
            for (int i = 0; i < mNx; i++)
            {
                vtifile << mGrid[FlatIndex(i, j)].density << " ";
            }
            vtifile << "\n";
        }
        vtifile << "        </DataArray>\n";
        
        // ECM anisotropy
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_anisotropy\" format=\"ascii\">\n";
        for (int j = 0; j < mNy; j++)
        {
            vtifile << "          ";
            for (int i = 0; i < mNx; i++)
            {
                vtifile << mGrid[FlatIndex(i, j)].anisotropy << " ";
            }
            vtifile << "\n";
        }
        vtifile << "        </DataArray>\n";
        
        // ECM orientation (3-component vector) — use cached trig
        vtifile << "        <DataArray type=\"Float32\" Name=\"ecm_orientation\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (int j = 0; j < mNy; j++)
        {
            vtifile << "          ";
            for (int i = 0; i < mNx; i++)
            {
                const ECMGridCell& cell = mGrid[FlatIndex(i, j)];
                vtifile << cell.cos_angle << " " << cell.sin_angle << " 0.0 ";
            }
            vtifile << "\n";
        }
        vtifile << "        </DataArray>\n";
        
        vtifile << "      </PointData>\n";
        vtifile << "    </Piece>\n";
        vtifile << "  </ImageData>\n";
        vtifile << "</VTKFile>\n";
        
        vtifile.close();
    }
    
    /**
     * Write ECM grid to VTU (XML UnstructuredGrid) file.
     * Used for hex grids. Each cell is a VTK_POLYGON (hexagon) with 6 vertices.
     */
    void WriteToVTU(const std::string& filename, double time)
    {
        std::ofstream vtufile(filename, std::ios::binary);
        
        if (!vtufile.is_open())
        {
            EXCEPTION("Could not open file: " + filename);
        }
        
        const int nCells  = mNx * mNy;
        const int nPoints = 6 * nCells;
        
        double r = mGridSpacing / std::sqrt(3.0);
        double hex_dx[6], hex_dy[6];
        for (int v = 0; v < 6; v++)
        {
            double angle = M_PI / 6.0 + v * M_PI / 3.0;
            hex_dx[v] = r * std::cos(angle);
            hex_dy[v] = r * std::sin(angle);
        }
        
        // Block sizes (bytes)
        const uint32_t pointsBytes  = nPoints * 3 * sizeof(float);
        const uint32_t connBytes    = nCells  * 6 * sizeof(int32_t);
        const uint32_t offsetsBytes = nCells      * sizeof(int32_t);
        const uint32_t typesBytes   = nCells      * sizeof(uint8_t);
        const uint32_t densityBytes = nCells      * sizeof(float);
        const uint32_t anisoBytes   = nCells      * sizeof(float);
        const uint32_t orientBytes  = nCells  * 3 * sizeof(float);
        const uint32_t hdr = sizeof(uint32_t);   // 4-byte block header
        
        // Appended-data offsets (each block = header + data)
        const uint32_t off0 = 0;
        const uint32_t off1 = off0 + hdr + pointsBytes;
        const uint32_t off2 = off1 + hdr + connBytes;
        const uint32_t off3 = off2 + hdr + offsetsBytes;
        const uint32_t off4 = off3 + hdr + typesBytes;
        const uint32_t off5 = off4 + hdr + densityBytes;
        const uint32_t off6 = off5 + hdr + anisoBytes;
        
        // --- XML header (written as text into a binary stream) ---
        std::ostringstream xml;
        xml << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\""
            << " byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
            << "  <UnstructuredGrid>\n"
            << "    <Piece NumberOfPoints=\"" << nPoints
            << "\" NumberOfCells=\"" << nCells << "\">\n"
            << "      <Points>\n"
            << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\""
            << " format=\"appended\" offset=\"" << off0 << "\"/>\n"
            << "      </Points>\n"
            << "      <Cells>\n"
            << "        <DataArray type=\"Int32\" Name=\"connectivity\""
            << " format=\"appended\" offset=\"" << off1 << "\"/>\n"
            << "        <DataArray type=\"Int32\" Name=\"offsets\""
            << " format=\"appended\" offset=\"" << off2 << "\"/>\n"
            << "        <DataArray type=\"UInt8\" Name=\"types\""
            << " format=\"appended\" offset=\"" << off3 << "\"/>\n"
            << "      </Cells>\n"
            << "      <CellData Scalars=\"ecm_density\" Vectors=\"ecm_orientation\">\n"
            << "        <DataArray type=\"Float32\" Name=\"ecm_density\""
            << " format=\"appended\" offset=\"" << off4 << "\"/>\n"
            << "        <DataArray type=\"Float32\" Name=\"ecm_anisotropy\""
            << " format=\"appended\" offset=\"" << off5 << "\"/>\n"
            << "        <DataArray type=\"Float32\" Name=\"ecm_orientation\""
            << " NumberOfComponents=\"3\" format=\"appended\" offset=\"" << off6 << "\"/>\n"
            << "      </CellData>\n"
            << "    </Piece>\n"
            << "  </UnstructuredGrid>\n"
            << "  <AppendedData encoding=\"raw\">\n_";
        std::string xmlStr = xml.str();
        vtufile.write(xmlStr.data(), xmlStr.size());
        
        // --- Binary data blocks ---
        
        // Block 0: Points (6 vertices per hex cell)
        {
            std::vector<float> pts(nPoints * 3);
            int pidx = 0;
            for (int j = 0; j < mNy; j++)
            {
                for (int i = 0; i < mNx; i++)
                {
                    c_vector<double, 2> ctr = GetCellPosition(i, j);
                    for (int v = 0; v < 6; v++)
                    {
                        pts[pidx++] = static_cast<float>(ctr[0] + hex_dx[v]);
                        pts[pidx++] = static_cast<float>(ctr[1] + hex_dy[v]);
                        pts[pidx++] = 0.0f;
                    }
                }
            }
            vtufile.write(reinterpret_cast<const char*>(&pointsBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(pts.data()), pointsBytes);
        }
        
        // Block 1: Connectivity  (0,1,2,3,4,5  6,7,8,9,10,11 ...)
        {
            std::vector<int32_t> conn(nCells * 6);
            for (int c = 0; c < nCells; c++)
            {
                int base = 6 * c;
                conn[6*c]   = base;     conn[6*c+1] = base+1;
                conn[6*c+2] = base+2;   conn[6*c+3] = base+3;
                conn[6*c+4] = base+4;   conn[6*c+5] = base+5;
            }
            vtufile.write(reinterpret_cast<const char*>(&connBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(conn.data()), connBytes);
        }
        
        // Block 2: Offsets  (6, 12, 18, ...)
        {
            std::vector<int32_t> offs(nCells);
            for (int c = 0; c < nCells; c++) offs[c] = 6 * (c + 1);
            vtufile.write(reinterpret_cast<const char*>(&offsetsBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(offs.data()), offsetsBytes);
        }
        
        // Block 3: Cell types  (all 7 = VTK_POLYGON)
        {
            std::vector<uint8_t> types(nCells, 7);
            vtufile.write(reinterpret_cast<const char*>(&typesBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(types.data()), typesBytes);
        }
        
        // Block 4: ecm_density
        {
            std::vector<float> buf(nCells);
            for (int k = 0; k < nCells; k++) buf[k] = static_cast<float>(mGrid[k].density);
            vtufile.write(reinterpret_cast<const char*>(&densityBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(buf.data()), densityBytes);
        }
        
        // Block 5: ecm_anisotropy
        {
            std::vector<float> buf(nCells);
            for (int k = 0; k < nCells; k++) buf[k] = static_cast<float>(mGrid[k].anisotropy);
            vtufile.write(reinterpret_cast<const char*>(&anisoBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(buf.data()), anisoBytes);
        }
        
        // Block 6: ecm_orientation (use cached cos/sin – no trig calls)
        {
            std::vector<float> buf(nCells * 3);
            for (int k = 0; k < nCells; k++)
            {
                buf[3*k]   = static_cast<float>(mGrid[k].cos_angle);
                buf[3*k+1] = static_cast<float>(mGrid[k].sin_angle);
                buf[3*k+2] = 0.0f;
            }
            vtufile.write(reinterpret_cast<const char*>(&orientBytes), hdr);
            vtufile.write(reinterpret_cast<const char*>(buf.data()), orientBytes);
        }
        
        // Close XML
        const char* closing = "\n  </AppendedData>\n</VTKFile>\n";
        vtufile.write(closing, std::strlen(closing));
        vtufile.close();
    }
    
    /**
     * Unified output method — picks VTI for square grids, VTU for hex grids.
     */
    void WriteOutput(const std::string& filename, double time)
    {
        if (mGridType == "hex")
        {
            WriteToVTU(filename, time);
        }
        else
        {
            WriteToVTI(filename, time);
        }
    }
    
    /**
     * Get the VTK file extension appropriate for this grid type.
     */
    std::string GetOutputExtension() const
    {
        return (mGridType == "hex") ? ".vtu" : ".vti";
    }
    
    /**
     * Get the grid type.
     */
    std::string GetGridType() const
    {
        return mGridType;
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
        
        vtkfile << "# vtk DataFile Version 3.0\n";
        vtkfile << "ECM Grid at time " << time << "\n";
        vtkfile << "ASCII\n";
        vtkfile << "DATASET STRUCTURED_POINTS\n";
        vtkfile << "DIMENSIONS " << mNx << " " << mNy << " 1\n";
        vtkfile << "ORIGIN " << mXMin << " " << mYMin << " 0.0\n";
        vtkfile << "SPACING " << mGridSpacing << " " << mGridSpacing << " 1.0\n";
        vtkfile << "POINT_DATA " << (mNx * mNy) << "\n";
        
        vtkfile << "SCALARS ecm_density float 1\n";
        vtkfile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < mNy; j++)
            for (int i = 0; i < mNx; i++)
                vtkfile << mGrid[FlatIndex(i, j)].density << "\n";
        
        vtkfile << "SCALARS ecm_anisotropy float 1\n";
        vtkfile << "LOOKUP_TABLE default\n";
        for (int j = 0; j < mNy; j++)
            for (int i = 0; i < mNx; i++)
                vtkfile << mGrid[FlatIndex(i, j)].anisotropy << "\n";
        
        vtkfile << "VECTORS ecm_orientation float\n";
        for (int j = 0; j < mNy; j++)
            for (int i = 0; i < mNx; i++)
            {
                double angle = mGrid[FlatIndex(i, j)].fiber_angle;
                vtkfile << cos(angle) << " " << sin(angle) << " 0.0\n";
            }
        
        vtkfile.close();
    }
    
    /** Setters for parameters */
    void SetRemodelingRate(double rate) { mRemodelingRate = rate; }
    void SetDegradationRate(double rate) { mDegradationRate = rate; }
    void SetDepositionRate(double rate) { mDepositionRate = rate; }
    void SetDiffusionCoeff(double coeff) { mDiffusionCoeff = coeff; }

    /**
     * Clear ECM density inside a given radius from a center point.
     * Sets density=0 for all voxels whose center lies within the radius.
     */
    void ClearDensityInsideRadius(c_vector<double, 2> center, double radius)
    {
        double r2 = radius * radius;
        for (int j = 0; j < mNy; j++)
        {
            for (int i = 0; i < mNx; i++)
            {
                c_vector<double, 2> pos = GetCellPosition(i, j);
                double dx = pos[0] - center[0];
                double dy = pos[1] - center[1];
                if (dx * dx + dy * dy < r2)
                {
                    mGrid[FlatIndex(i, j)].density = 0.0;
                }
            }
        }
    }
};

#endif /* DYNAMICECMFIELD_HPP_ */
