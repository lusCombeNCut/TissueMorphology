/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef LUMENPRESSUREFORCE_HPP_
#define LUMENPRESSUREFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <algorithm>
#include <vector>
#include <cmath>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimProfiler.hpp"

/**
 * Hydrostatic lumen pressure force: F_i = P ∇_i V_lumen.
 *
 * Models the lumen as a fluid-filled cavity exerting uniform hydrostatic
 * pressure P on the apical (inner) surface of the epithelium.
 *
 * The force on each cell/node is the gradient of the enclosed
 * area (2D) or volume (3D) with respect to that cell's position,
 * scaled by the constant pressure P.
 *
 * 2D: Cells are sorted by angle around the centroid to form a polygon.
 *     The gradient is the exact shoelace-formula derivative:
 *       ∂A/∂x_i = ½(y_{i+1} − y_{i−1})
 *       ∂A/∂y_i = ½(x_{i−1} − x_{i+1})
 *
 * 3D (node-based): Approximate volume gradient assuming a roughly
 *     spherical arrangement:  ∇_i V ≈ (4πr²/N) n̂_i
 *
 * For vertex-based populations the per-cell force is distributed
 * equally to all nodes of the element.
 */
template<unsigned DIM>
class LumenPressureForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Constant hydrostatic lumen pressure */
    double mPressure;

    /** Center of the lumen (auto-tracked from centroid if mTrackCenter=true) */
    c_vector<double, DIM> mLumenCenter;

    /** Whether to auto-track lumen center from population centroid */
    bool mTrackCenter;

    /** Whether to store force vectors in CellData for visualization */
    bool mWriteForce;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mPressure;
        archive & mLumenCenter;
        archive & mTrackCenter;
        archive & mWriteForce;
    }

public:

    LumenPressureForce()
        : AbstractForce<DIM>(),
          mPressure(1.0),
          mTrackCenter(true),
          mWriteForce(false)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mLumenCenter[i] = 0.0;
        }
    }

    virtual ~LumenPressureForce() {}

    /**
     * Hydrostatic pressure force: F_i = P ∇_i V_lumen.
     *
     * In 2D the lumen is the polygon formed by cell centres (sorted by
     * angle around the centroid).  The gradient is the exact derivative
     * of the shoelace area formula.
     *
     * In 3D the lumen volume gradient is approximated as
     * (4πr²/N) n̂_i for a roughly spherical cell arrangement.
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        ScopedTimer _prof("LumenPressure");
        if (mTrackCenter)
        {
            UpdateLumenCenter(rCellPopulation);
        }

        if constexpr (DIM == 2)
        {
            // Detect vertex-based population (for force distribution)
            VertexBasedCellPopulation<DIM>* p_vertex_pop =
                dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

            // ── 2D: exact shoelace area gradient ──────────────────────
            // Collect cell centres and sort by angle to obtain the
            // polygon vertex ordering (CCW).
            struct CellInfo
            {
                double angle;
                unsigned locIndex;
                c_vector<double, 2> position;
            };
            std::vector<CellInfo> cells;
            cells.reserve(rCellPopulation.GetNumRealCells());

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                unsigned idx = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, 2> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                double angle = std::atan2(pos[1] - mLumenCenter[1],
                                          pos[0] - mLumenCenter[0]);
                cells.push_back({angle, idx, pos});
            }

            std::sort(cells.begin(), cells.end(),
                [](const CellInfo& a, const CellInfo& b) { return a.angle < b.angle; });

            const unsigned N = cells.size();
            for (unsigned i = 0; i < N; i++)
            {
                const unsigned iPrev = (i + N - 1) % N;
                const unsigned iNext = (i + 1) % N;

                // Shoelace area gradient w.r.t. vertex i:
                //   ∂A/∂x_i = ½(y_{i+1} − y_{i−1})
                //   ∂A/∂y_i = ½(x_{i−1} − x_{i+1})
                c_vector<double, 2> grad;
                grad[0] = 0.5 * (cells[iNext].position[1] - cells[iPrev].position[1]);
                grad[1] = 0.5 * (cells[iPrev].position[0] - cells[iNext].position[0]);

                c_vector<double, 2> force = mPressure * grad;

                unsigned locIndex = cells[i].locIndex;

                // Store force in CellData for visualization
                if (mWriteForce)
                {
                    CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(locIndex);
                    pCell->GetCellData()->SetItem("lumen_force_x", force[0]);
                    pCell->GetCellData()->SetItem("lumen_force_y", force[1]);
                    if constexpr (DIM == 3)
                    {
                        pCell->GetCellData()->SetItem("lumen_force_z", 0.0);
                    }
                }

                if (p_vertex_pop)
                {
                    VertexElement<DIM, DIM>* p_element =
                        p_vertex_pop->rGetMesh().GetElement(locIndex);
                    unsigned n_nodes = p_element->GetNumNodes();
                    c_vector<double, 2> force_per_node =
                        force / static_cast<double>(n_nodes);
                    for (unsigned j = 0; j < n_nodes; j++)
                    {
                        p_element->GetNode(j)->AddAppliedForceContribution(force_per_node);
                    }
                }
                else
                {
                    rCellPopulation.GetNode(locIndex)->AddAppliedForceContribution(force);
                }
            }
        }
        else // DIM == 3
        {
            // ── 3D: approximate volume gradient ──────────────────────
            // For a roughly spherical arrangement of N cells:
            //   ∇_i V ≈ (4πr_i²/N) n̂_i
            const unsigned numCells = rCellPopulation.GetNumRealCells();

            for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End(); ++cell_iter)
            {
                unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                c_vector<double, DIM> pos = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                c_vector<double, DIM> disp = pos - mLumenCenter;
                double r = norm_2(disp);

                if (r < 1e-6)
                {
                    continue;
                }

                c_vector<double, DIM> n_hat = disp / r;
                double area_element = 4.0 * M_PI * r * r / numCells;
                c_vector<double, DIM> force = mPressure * area_element * n_hat;

                // Store force in CellData for visualization
                if (mWriteForce)
                {
                    CellPtr pCell = rCellPopulation.GetCellUsingLocationIndex(loc_index);
                    pCell->GetCellData()->SetItem("lumen_force_x", force[0]);
                    pCell->GetCellData()->SetItem("lumen_force_y", force[1]);
                    pCell->GetCellData()->SetItem("lumen_force_z", force[2]);
                }

                rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(force);
            }
        }
    }

    void UpdateLumenCenter(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> center = zero_vector<double>(DIM);
        unsigned count = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End(); ++cell_iter)
        {
            center += rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            count++;
        }

        if (count > 0)
        {
            mLumenCenter = center / static_cast<double>(count);
        }
    }

    // ---- Setters/Getters ----

    void SetPressure(double pressure)
    {
        mPressure = pressure;
    }

    double GetPressure() const { return mPressure; }

    void SetLumenCenter(c_vector<double, DIM> center)
    {
        mLumenCenter = center;
        mTrackCenter = false;
    }

    c_vector<double, DIM> GetLumenCenter() const { return mLumenCenter; }

    void SetTrackCenter(bool track) { mTrackCenter = track; }
    bool GetTrackCenter() const { return mTrackCenter; }

    void SetWriteForce(bool write) { mWriteForce = write; }
    bool GetWriteForce() const { return mWriteForce; }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Pressure>" << mPressure << "</Pressure>\n";
        *rParamsFile << "\t\t\t<TrackCenter>" << mTrackCenter << "</TrackCenter>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenPressureForce)

#endif /*LUMENPRESSUREFORCE_HPP_*/
