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

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * A force class modelling hydrostatic/osmotic pressure from the organoid lumen.
 *
 * In a real organoid the lumen is maintained by apical fluid secretion and
 * tight-junction-mediated osmotic pressure. This force provides a simple
 * outward radial pressure on all cells, pushing them away from the organoid
 * centre. Combined with the inward BasementMembraneForce it creates the
 * apical (inner) / basal (outer) tension that maintains the epithelial shell.
 *
 * The force acts only on cells within the organoid's equilibrium radius. Cells
 * that have budded outward beyond that radius feel no additional lumen pressure
 * (the bud may form its own micro-lumen but that's not modelled here).
 *
 * Biology:
 *   - Apical surface faces the lumen (inner)
 *   - Basal surface faces the ECM (outer)
 *   - Lumen pressure ↔ cell proliferative pressure → budding instability
 *
 * Parameters:
 *   mPressureStrength    — magnitude of outward radial force per cell
 *   mLumenEquilibriumRadius — the "relaxed" lumen surface radius
 *   mLumenCenter         — center of the lumen (auto-tracked from centroid)
 */
template<unsigned DIM>
class LumenPressureForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Pressure strength (outward force magnitude per unit displacement) */
    double mPressureStrength;

    /** Equilibrium radius of the lumen surface */
    double mLumenEquilibriumRadius;

    /** Center of the lumen */
    c_vector<double, DIM> mLumenCenter;

    /** Whether to auto-track lumen center from population centroid */
    bool mTrackCenter;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mPressureStrength;
        archive & mLumenEquilibriumRadius;
        archive & mLumenCenter;
        archive & mTrackCenter;
    }

public:

    /**
     * Constructor.
     */
    LumenPressureForce()
        : AbstractForce<DIM>(),
          mPressureStrength(1.0),
          mLumenEquilibriumRadius(5.0),
          mTrackCenter(true)
    {
        for (unsigned i = 0; i < DIM; i++)
        {
            mLumenCenter[i] = 0.0;
        }
    }

    virtual ~LumenPressureForce() {}

    /**
     * Add lumen pressure force contribution.
     *
     * For each cell, compute the outward radial force from the lumen center.
     * The force pushes cells outward toward the equilibrium shell radius:
     *   F = pressure * (R_eq - r) * r_hat   for r < R_eq
     *   F = 0                                for r >= R_eq
     *
     * where r is the cell's distance from lumen center, R_eq is the
     * equilibrium radius, and r_hat is the outward unit radial vector.
     *
     * This creates a "pressure vessel" effect: cells inside the shell
     * are pushed outward, maintaining lumen openness.
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        // Optionally update lumen center from cell population centroid
        if (mTrackCenter)
        {
            UpdateLumenCenter(rCellPopulation);
        }

        // Detect vertex-based population (for force distribution)
        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            // Use cell centre (works for both node-based and vertex-based)
            c_vector<double, DIM> location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            c_vector<double, DIM> displacement = location - mLumenCenter;
            double distance = norm_2(displacement);

            // Avoid division by zero for cells at the exact center
            if (distance < 1e-6)
            {
                // Apply a small random outward nudge to break symmetry
                c_vector<double, DIM> nudge;
                for (unsigned d = 0; d < DIM; d++)
                {
                    nudge[d] = (static_cast<double>(rand()) / RAND_MAX - 0.5) * 0.01;
                }
                unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                if (p_vertex_pop)
                {
                    VertexElement<DIM, DIM>* p_element = p_vertex_pop->rGetMesh().GetElement(loc_index);
                    unsigned n_nodes = p_element->GetNumNodes();
                    c_vector<double, DIM> nudge_per_node = nudge / static_cast<double>(n_nodes);
                    for (unsigned i = 0; i < n_nodes; i++)
                    {
                        p_element->GetNode(i)->AddAppliedForceContribution(nudge_per_node);
                    }
                }
                else
                {
                    rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(nudge);
                }
                continue;
            }

            c_vector<double, DIM> unit_radial = displacement / distance;

            // Apply outward pressure only to cells inside the equilibrium shell
            if (distance < mLumenEquilibriumRadius)
            {
                double deficit = mLumenEquilibriumRadius - distance;
                c_vector<double, DIM> force = mPressureStrength * deficit * unit_radial;

                // Apply force: distribute across element vertices for vertex pops
                unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                if (p_vertex_pop)
                {
                    VertexElement<DIM, DIM>* p_element = p_vertex_pop->rGetMesh().GetElement(loc_index);
                    unsigned n_nodes = p_element->GetNumNodes();
                    c_vector<double, DIM> force_per_node = force / static_cast<double>(n_nodes);
                    for (unsigned i = 0; i < n_nodes; i++)
                    {
                        p_element->GetNode(i)->AddAppliedForceContribution(force_per_node);
                    }
                }
                else
                {
                    rCellPopulation.GetNode(loc_index)->AddAppliedForceContribution(force);
                }
            }
        }
    }

    /**
     * Update lumen center from population centroid.
     */
    void UpdateLumenCenter(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> center = zero_vector<double>(DIM);
        unsigned count = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
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

    void SetPressureStrength(double strength)
    {
        assert(strength >= 0.0);
        mPressureStrength = strength;
    }

    double GetPressureStrength() const { return mPressureStrength; }

    void SetLumenEquilibriumRadius(double radius)
    {
        assert(radius > 0.0);
        mLumenEquilibriumRadius = radius;
    }

    double GetLumenEquilibriumRadius() const { return mLumenEquilibriumRadius; }

    void SetLumenCenter(c_vector<double, DIM> center)
    {
        mLumenCenter = center;
        mTrackCenter = false;  // Manual center disables auto-tracking
    }

    c_vector<double, DIM> GetLumenCenter() const { return mLumenCenter; }

    void SetTrackCenter(bool track) { mTrackCenter = track; }

    bool GetTrackCenter() const { return mTrackCenter; }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<PressureStrength>" << mPressureStrength << "</PressureStrength>\n";
        *rParamsFile << "\t\t\t<LumenEquilibriumRadius>" << mLumenEquilibriumRadius << "</LumenEquilibriumRadius>\n";
        *rParamsFile << "\t\t\t<TrackCenter>" << mTrackCenter << "</TrackCenter>\n";
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenPressureForce)

#endif /*LUMENPRESSUREFORCE_HPP_*/
