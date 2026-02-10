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

#ifndef BASEMENTMEMBRANEFORCE_HPP_
#define BASEMENTMEMBRANEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "SimulationTime.hpp"

/**
 * A force class to model basement membrane interactions in organoid formation.
 *
 * This force applies a restoring force toward a central region, simulating
 * the constraining effect of a basement membrane with variable stiffness.
 * The force magnitude depends on the distance from the organoid center and
 * the basement membrane stiffness parameter stored in cell data.
 */
template<unsigned DIM>
class BasementMembraneForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** The center of the organoid (basement membrane center) */
    c_vector<double, DIM> mOrganoidCenter;

    /** The radius beyond which basement membrane forces apply */
    double mBasementMembraneRadius;

    /** Initial basement membrane radius (for degradation calculations) */
    double mInitialRadius;

    /** Default force strength */
    double mForceStrength;

    /** ECM degradation rate (radius increase per unit time) */
    double mEcmDegradationRate;

    /** Maximum basement membrane radius (biological limit) */
    double mMaxRadius;

    /** Whether ECM degradation is enabled */
    bool mEcmDegradationEnabled;

    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  the archive
     * @param version  the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mOrganoidCenter;
        archive & mBasementMembraneRadius;
        archive & mInitialRadius;
        archive & mForceStrength;
        archive & mEcmDegradationRate;
        archive & mMaxRadius;
        archive & mEcmDegradationEnabled;
    }

public:

    /**
     * Constructor.
     */
    BasementMembraneForce()
        : AbstractForce<DIM>(),
          mBasementMembraneRadius(2.0),
          mInitialRadius(2.0),
          mForceStrength(1.0),
          mEcmDegradationRate(0.1),  // Default: 0.1 radius units per time unit
          mMaxRadius(10.0),          // Default: 5x initial radius max
          mEcmDegradationEnabled(false)
    {
        // Initialize center to origin
        for (unsigned i = 0; i < DIM; i++)
        {
            mOrganoidCenter[i] = 0.0;
        }
    }

    /**
     * Destructor.
     */
    virtual ~BasementMembraneForce()
    {
    }

    /**
     * Add force to each node.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        // Update organoid center based on current cell positions
        UpdateOrganoidCenter(rCellPopulation);
        
        // Update basement membrane radius due to ECM degradation
        UpdateBasementMembraneRadius();

        // Apply basement membrane forces to cells outside the radius
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

            c_vector<double, DIM> node_location = p_node->rGetLocation();
            c_vector<double, DIM> displacement = node_location - mOrganoidCenter;
            double distance_from_center = norm_2(displacement);

            // Apply force if cell is beyond basement membrane radius
            if (distance_from_center > mBasementMembraneRadius)
            {
                // Get basement membrane stiffness from cell data
                double stiffness = mForceStrength;
                try
                {
                    stiffness = cell_iter->GetCellData()->GetItem("basement_membrane_stiffness");
                }
                catch (Exception& e)
                {
                    // Use default if not found
                }

                // Calculate restoring force toward the boundary
                double excess_distance = distance_from_center - mBasementMembraneRadius;
                c_vector<double, DIM> unit_displacement = displacement / distance_from_center;
                c_vector<double, DIM> force = -stiffness * excess_distance * unit_displacement;

                p_node->AddAppliedForceContribution(force);
            }
        }
    }

    /**
     * Update the organoid center based on current cell positions.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateOrganoidCenter(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        c_vector<double, DIM> center = zero_vector<double>(DIM);
        unsigned num_cells = 0;

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

            center += p_node->rGetLocation();
            num_cells++;
        }

        if (num_cells > 0)
        {
            mOrganoidCenter = center / static_cast<double>(num_cells);
        }
    }

    /**
     * @return the organoid center
     */
    c_vector<double, DIM> GetOrganoidCenter()
    {
        return mOrganoidCenter;
    }

    /**
     * Set the organoid center.
     *
     * @param center the new center location
     */
    void SetOrganoidCenter(c_vector<double, DIM> center)
    {
        mOrganoidCenter = center;
    }

    /**
     * @return the basement membrane radius
     */
    double GetBasementMembraneRadius()
    {
        return mBasementMembraneRadius;
    }

    /**
     * Set the basement membrane radius.
     *
     * @param radius the new radius
     */
    void SetBasementMembraneRadius(double radius)
    {
        assert(radius > 0.0);
        mBasementMembraneRadius = radius;
    }

    /**
     * @return the force strength
     */
    double GetForceStrength()
    {
        return mForceStrength;
    }

    /**
     * Set the force strength.
     *
     * @param strength the new force strength
     */
    void SetForceStrength(double strength)
    {
        mForceStrength = strength;
    }

    /**
     * Set target radius (convenience method for 3D organoid studies).
     *
     * @param targetRadius the target organoid radius
     */
    void SetTargetRadius(double targetRadius)
    {
        SetBasementMembraneRadius(targetRadius);
    }
    
    /**
     * Set basement membrane parameter (alias for force strength).
     *
     * @param parameter the basement membrane stiffness parameter
     */
    void SetBasementMembraneParameter(double parameter)
    {
        mForceStrength = parameter;
    }

    /**
     * Enable ECM degradation with specified rate.
     *
     * @param degradationRate rate of radius increase per unit time
     * @param maxRadius maximum allowed radius (prevents infinite expansion)
     */
    void EnableEcmDegradation(double degradationRate = 0.1, double maxRadius = 10.0)
    {
        mEcmDegradationEnabled = true;
        mEcmDegradationRate = degradationRate;
        mMaxRadius = maxRadius;
        mInitialRadius = mBasementMembraneRadius;  // Store initial value
    }

    /**
     * Disable ECM degradation (radius remains constant).
     */
    void DisableEcmDegradation()
    {
        mEcmDegradationEnabled = false;
    }

    /**
     * Update basement membrane radius due to ECM degradation.
     * Called each timestep if degradation is enabled.
     */
    void UpdateBasementMembraneRadius()
    {
        if (mEcmDegradationEnabled)
        {
            // Get current simulation time
            double current_time = SimulationTime::Instance()->GetTime();
            
            // Calculate new radius: R(t) = R_0 + rate * t
            double new_radius = mInitialRadius + mEcmDegradationRate * current_time;
            
            // Clamp to maximum radius
            mBasementMembraneRadius = std::min(new_radius, mMaxRadius);
        }
    }

    /**
     * Get ECM degradation rate.
     *
     * @return the current degradation rate
     */
    double GetEcmDegradationRate() const
    {
        return mEcmDegradationRate;
    }

    /**
     * Get maximum radius.
     *
     * @return the maximum allowed radius
     */
    double GetMaxRadius() const
    {
        return mMaxRadius;
    }

    /**
     * Check if ECM degradation is enabled.
     *
     * @return true if degradation is active
     */
    bool IsEcmDegradationEnabled() const
    {
        return mEcmDegradationEnabled;
    }

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<BasementMembraneRadius>" << mBasementMembraneRadius << "</BasementMembraneRadius>\n";
        *rParamsFile << "\t\t\t<InitialRadius>" << mInitialRadius << "</InitialRadius>\n";
        *rParamsFile << "\t\t\t<ForceStrength>" << mForceStrength << "</ForceStrength>\n";
        *rParamsFile << "\t\t\t<EcmDegradationRate>" << mEcmDegradationRate << "</EcmDegradationRate>\n";
        *rParamsFile << "\t\t\t<MaxRadius>" << mMaxRadius << "</MaxRadius>\n";
        *rParamsFile << "\t\t\t<EcmDegradationEnabled>" << mEcmDegradationEnabled << "</EcmDegradationEnabled>\n";

        // Call method on direct parent class
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BasementMembraneForce)

#endif /*BASEMENTMEMBRANEFORCE_HPP_*/