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

#ifndef APICALCONSTRICTIONFORCE_HPP_
#define APICALCONSTRICTIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "SimProfiler.hpp"

/**
 * A force class to model apical constriction in epithelial tissues.
 *
 * This force causes cells to reduce their apical surface area, creating
 * a wedge shape that drives tissue invagination and folding.
 * 
 * The force acts radially inward on cells marked as "apical", simulating
 * actomyosin contraction at the apical surface.
 */
template<unsigned DIM>
class ApicalConstrictionForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** Strength of apical constriction */
    double mConstrictionStrength;

    /** Target reduction factor for apical area (0-1, where 0.5 = 50% reduction) */
    double mTargetReduction;

    /** Reference area for each cell (set at initialization) */
    std::map<unsigned, double> mReferenceAreas;

    /**
     * Boost Serialization method for archiving/checkpointing.
     *
     * @param archive  the archive
     * @param version  the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mConstrictionStrength;
        archive & mTargetReduction;
        archive & mReferenceAreas;
    }

public:

    /**
     * Constructor.
     */
    ApicalConstrictionForce()
        : AbstractForce<DIM>(),
          mConstrictionStrength(5.0),
          mTargetReduction(0.5)
    {
    }

    /**
     * Destructor.
     */
    virtual ~ApicalConstrictionForce()
    {
    }

    /**
     * Add force contribution from apical constriction.
     *
     * For each cell marked as "apical", applies an inward radial force
     * proportional to the difference between current and target area.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        ScopedTimer timer("ApicalConstriction");
        // Initialize reference areas if not done
        if (mReferenceAreas.empty())
        {
            InitializeReferenceAreas(rCellPopulation);
        }

        // Calculate population centroid for radial direction
        c_vector<double, DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

        // Detect vertex-based population (for force distribution)
        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        // Apply constriction force to apical cells
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            // Check if cell is apical
            bool is_apical = false;
            try
            {
                double apical_marker = cell_iter->GetCellData()->GetItem("is_apical");
                is_apical = (apical_marker > 0.5);
            }
            catch (Exception&)
            {
                // If not marked, assume not apical
                continue;
            }

            if (!is_apical)
            {
                continue;
            }

            unsigned loc_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

            // Get cell's current area (approximation based on local density)
            double current_area = EstimateCellArea(rCellPopulation, loc_index);
            
            // Get reference area
            double reference_area = mReferenceAreas[loc_index];
            
            // Calculate target area
            double target_area = reference_area * (1.0 - mTargetReduction);

            // Force magnitude proportional to area difference
            double area_difference = current_area - target_area;
            
            if (area_difference > 0.0)  // Only apply if cell needs to constrict
            {
                // Direction: toward centroid (apical constriction)
                c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
                c_vector<double, DIM> displacement = cell_location - centroid;
                double distance = norm_2(displacement);
                
                if (distance > 1e-6)  // Avoid division by zero
                {
                    c_vector<double, DIM> unit_radial = displacement / distance;
                    
                    // Force: inward (negative radial), proportional to excess area
                    c_vector<double, DIM> force = -mConstrictionStrength * area_difference * unit_radial;

                    // Apply force: distribute across element vertices for vertex pops
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
    }

    /**
     * Initialize reference areas for all cells.
     * Called once at the start of simulation.
     *
     * @param rCellPopulation reference to the cell population
     */
    void InitializeReferenceAreas(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            double area = EstimateCellArea(rCellPopulation, node_index);
            mReferenceAreas[node_index] = area;
        }
    }

    /**
     * Estimate cell area based on Voronoi tessellation or neighbor distances.
     * 
     * For node-based populations, uses average squared distance to neighbors.
     * For mesh-based populations, uses element area.
     *
     * @param rCellPopulation reference to the cell population
     * @param nodeIndex index of the cell/node
     * @return estimated area
     */
    double EstimateCellArea(AbstractCellPopulation<DIM>& rCellPopulation, unsigned locIndex)
    {
        // For vertex-based: use actual element area
        VertexBasedCellPopulation<DIM>* p_vertex_pop =
            dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        if (p_vertex_pop)
        {
            return p_vertex_pop->rGetMesh().GetVolumeOfElement(locIndex);
        }

        // For node-based: area ~ π * r_neighbor^2
        Node<DIM>* p_node = rCellPopulation.GetNode(locIndex);
        c_vector<double, DIM> node_location = p_node->rGetLocation();
        
        // Find neighboring nodes
        std::set<unsigned> neighbor_indices = rCellPopulation.GetNeighbouringNodeIndices(locIndex);
        
        if (neighbor_indices.empty())
        {
            return 1.0;  // Default area
        }
        
        double total_distance = 0.0;
        unsigned count = 0;
        
        for (std::set<unsigned>::iterator iter = neighbor_indices.begin();
             iter != neighbor_indices.end();
             ++iter)
        {
            Node<DIM>* p_neighbor = rCellPopulation.GetNode(*iter);
            c_vector<double, DIM> neighbor_location = p_neighbor->rGetLocation();
            double distance = norm_2(neighbor_location - node_location);
            total_distance += distance;
            count++;
        }
        
        double average_distance = total_distance / static_cast<double>(count);
        
        // Estimate area as circle with radius = average neighbor distance
        if (DIM == 2)
        {
            return M_PI * average_distance * average_distance;
        }
        else if (DIM == 3)
        {
            // For 3D, use surface area approximation
            return 4.0 * M_PI * average_distance * average_distance;
        }
        
        return 1.0;
    }

    /**
     * Set the constriction strength.
     *
     * @param strength force per unit area difference
     */
    void SetConstrictionStrength(double strength)
    {
        assert(strength >= 0.0);
        mConstrictionStrength = strength;
    }

    /**
     * Get the constriction strength.
     *
     * @return the constriction strength
     */
    double GetConstrictionStrength() const
    {
        return mConstrictionStrength;
    }

    /**
     * Set the target reduction factor.
     *
     * @param reduction factor between 0 and 1 (e.g., 0.5 for 50% reduction)
     */
    void SetTargetReduction(double reduction)
    {
        assert(reduction >= 0.0 && reduction <= 1.0);
        mTargetReduction = reduction;
    }

    /**
     * Get the target reduction factor.
     *
     * @return the target reduction
     */
    double GetTargetReduction() const
    {
        return mTargetReduction;
    }

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ConstrictionStrength>" << mConstrictionStrength << "</ConstrictionStrength>\n";
        *rParamsFile << "\t\t\t<TargetReduction>" << mTargetReduction << "</TargetReduction>\n";

        // Call method on direct parent class
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ApicalConstrictionForce)

#endif /*APICALCONSTRICTIONFORCE_HPP_*/
