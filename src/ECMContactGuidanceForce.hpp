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

#ifndef ECMCONTACTGUIDANCEFORCE_HPP_
#define ECMCONTACTGUIDANCEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * ECM Contact Guidance Force
 * 
 * Part of the three-factor morphogenesis model from Metzcar et al. 2025.
 * This force represents structured ECM that guides cell migration through
 * contact guidance (cells align their motion with ECM fiber orientation).
 * 
 * Implements Section 3.1: "Invasive Cellular Front Pushing into ECM"
 * Based on Painter 2009's work on tumor invasion into structured matrices.
 * 
 * Key parameters:
 * - Base cell speed: v0 (default 1.25 µm/min)
 * - ECM sensitivity: s (how strongly cells follow ECM, default 1.0)
 * - Anisotropy: a (0=random walk, 1=fully aligned, default 1.0)
 * - ECM orientation pattern: random, parallel, perpendicular, or mixed
 * 
 * The force magnitude represents the cell's intrinsic motility modulated by
 * ECM structure. Cells preferentially migrate along ECM fibers.
 */
template<unsigned DIM>
class ECMContactGuidanceForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:

    /** ECM orientation pattern type */
    std::string mECMOrientationType;
    
    /** Base cell migration speed (µm/min) */
    double mBaseSpeed;
    
    /** ECM sensitivity parameter (0-1+) */
    double mECMSensitivity;
    
    /** Degree of anisotropy (0=isotropic/random, 1=fully aligned) */
    double mAnisotropy;
    
    /** Domain width for mixed ECM pattern (µm) */
    double mDomainWidth;
    
    /** Use persistent random walk (cells remember previous direction) */
    bool mUsePersistence;
    
    /** Persistence time (how long cells maintain direction) */
    double mPersistenceTime;
    
    /**
     * Boost Serialization method for archiving/checkpointing.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mECMOrientationType;
        archive & mBaseSpeed;
        archive & mECMSensitivity;
        archive & mAnisotropy;
        archive & mDomainWidth;
        archive & mUsePersistence;
        archive & mPersistenceTime;
    }
    
    /**
     * Get ECM fiber orientation angle at a given position
     * 
     * @param position The spatial position to query
     * @return Angle in radians [0, 2π]
     */
    double GetECMAngleAtPosition(c_vector<double, DIM> position)
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        if (mECMOrientationType == "random")
        {
            // Random orientation - isotropic ECM
            return p_gen->ranf() * 2.0 * M_PI;
        }
        else if (mECMOrientationType == "parallel")
        {
            // Parallel to domain edge (horizontal, θ = 0)
            // This restricts invasion perpendicular to the front
            return 0.0;
        }
        else if (mECMOrientationType == "perpendicular")
        {
            // Perpendicular to domain edge (vertical, θ = π/2)
            // This facilitates invasion into the domain
            return M_PI / 2.0;
        }
        else if (mECMOrientationType == "mixed")
        {
            // Mixed: perpendicular in center region, parallel on sides
            // Creates heterogeneous invasion pattern
            double x = position[0];
            double center_start = mDomainWidth / 3.0;
            double center_end = 2.0 * mDomainWidth / 3.0;
            
            if (x > center_start && x < center_end)
            {
                return M_PI / 2.0;  // Perpendicular (fast invasion)
            }
            else
            {
                return 0.0;  // Parallel (slow invasion)
            }
        }
        
        return 0.0;  // Default to horizontal
    }

public:

    /**
     * Constructor.
     */
    ECMContactGuidanceForce()
        : AbstractForce<DIM>(),
          mECMOrientationType("random"),
          mBaseSpeed(1.25),        // µm/min (from Metzcar 2025)
          mECMSensitivity(1.0),    // Full ECM sensitivity
          mAnisotropy(1.0),        // Highly anisotropic ECM
          mDomainWidth(600.0),     // Default domain width (µm)
          mUsePersistence(false),
          mPersistenceTime(10.0)   // 10 minutes
    {
    }
    
    /**
     * Destructor.
     */
    virtual ~ECMContactGuidanceForce()
    {
    }
    
    /**
     * Set ECM orientation pattern type
     * 
     * @param type "random", "parallel", "perpendicular", or "mixed"
     */
    void SetECMOrientationType(std::string type)
    {
        mECMOrientationType = type;
    }
    
    /**
     * Set base cell migration speed
     * 
     * @param speed Base speed in µm/min
     */
    void SetBaseSpeed(double speed)
    {
        mBaseSpeed = speed;
    }
    
    /**
     * Set ECM sensitivity parameter
     * 
     * @param sensitivity How strongly cells follow ECM (0-1+)
     */
    void SetECMSensitivity(double sensitivity)
    {
        mECMSensitivity = sensitivity;
    }
    
    /**
     * Set degree of ECM anisotropy
     * 
     * @param anisotropy 0=isotropic (random walk), 1=fully aligned
     */
    void SetAnisotropy(double anisotropy)
    {
        mAnisotropy = anisotropy;
    }
    
    /**
     * Set domain width for mixed ECM pattern
     * 
     * @param width Domain width in µm
     */
    void SetDomainWidth(double width)
    {
        mDomainWidth = width;
    }
    
    /**
     * Apply ECM contact guidance force to cells
     * 
     * Cells migrate with a preferred direction based on local ECM orientation.
     * Key improvements for realism:
     * 1. Bidirectional fiber movement (cells can move along fiber in either direction)
     * 2. Perpendicular diffusion (cells can move across fibers with reduced mobility)
     * 3. Random walk component (even in aligned ECM)
     * 
     * Also stores ECM orientation information in cell data for visualization.
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, DIM> position = p_node->rGetLocation();
            
            // Get ECM fiber orientation at this position
            double ecm_angle = GetECMAngleAtPosition(position);
            
            // ECM fiber direction vector (unit vector along fiber)
            c_vector<double, DIM> ecm_fiber_direction = zero_vector<double>(DIM);
            ecm_fiber_direction[0] = cos(ecm_angle);
            if (DIM > 1)
            {
                ecm_fiber_direction[1] = sin(ecm_angle);
            }
            
            // Store ECM orientation in cell data for visualization
            cell_iter->GetCellData()->SetItem("ecm_angle", ecm_angle);
            cell_iter->GetCellData()->SetItem("ecm_orientation_x", ecm_fiber_direction[0]);
            if (DIM > 1)
            {
                cell_iter->GetCellData()->SetItem("ecm_orientation_y", ecm_fiber_direction[1]);
            }
            if (DIM > 2)
            {
                cell_iter->GetCellData()->SetItem("ecm_orientation_z", 0.0);
            }
            
            // ========================================
            // REALISTIC CELL MIGRATION ON ECM FIBERS
            // ========================================
            
            // 1. Choose direction along fiber (bidirectional - cells can go either way)
            //    This is critical: fibers don't have arrows, just orientation!
            double fiber_direction_sign = (p_gen->ranf() < 0.5) ? -1.0 : 1.0;
            c_vector<double, DIM> fiber_motion = fiber_direction_sign * ecm_fiber_direction;
            
            // 2. Add perpendicular motion component (reduced by anisotropy)
            //    Cells can still move across fibers, just more slowly
            c_vector<double, DIM> perpendicular_direction = zero_vector<double>(DIM);
            if (DIM > 1)
            {
                // Perpendicular to fiber in 2D: rotate by 90 degrees
                perpendicular_direction[0] = -ecm_fiber_direction[1];
                perpendicular_direction[1] = ecm_fiber_direction[0];
            }
            
            // Random component for perpendicular motion
            double perp_magnitude = p_gen->ranf() * 2.0 - 1.0;  // [-1, 1]
            c_vector<double, DIM> perpendicular_motion = perp_magnitude * perpendicular_direction;
            
            // 3. Add isotropic random walk component
            //    Even in highly aligned ECM, cells have intrinsic random motility
            c_vector<double, DIM> random_direction = zero_vector<double>(DIM);
            double random_angle = p_gen->ranf() * 2.0 * M_PI;
            random_direction[0] = cos(random_angle);
            if (DIM > 1)
            {
                random_direction[1] = sin(random_angle);
            }
            
            // 4. Combine components with anisotropy weighting
            //    High anisotropy (a→1): strong fiber alignment, weak perpendicular/random
            //    Low anisotropy (a→0): isotropic random walk
            //
            // Mobility ratio from Metzcar: parallel/perpendicular ≈ 2-5 for aligned ECM
            double perpendicular_mobility = (1.0 - mAnisotropy) * 0.3;  // Reduced cross-fiber motion
            double random_mobility = (1.0 - mAnisotropy) * 0.5;        // Random walk component
            
            c_vector<double, DIM> migration_direction = 
                mAnisotropy * fiber_motion +                    // Along fiber (bidirectional)
                perpendicular_mobility * perpendicular_motion + // Across fiber (reduced)
                random_mobility * random_direction;             // Random walk
            
            // Normalize direction vector
            double magnitude = norm_2(migration_direction);
            if (magnitude > 1e-10)
            {
                migration_direction /= magnitude;
            }
            
            // Store actual migration direction in cell data
            cell_iter->GetCellData()->SetItem("migration_direction_x", migration_direction[0]);
            if (DIM > 1)
            {
                cell_iter->GetCellData()->SetItem("migration_direction_y", migration_direction[1]);
            }
            if (DIM > 2)
            {
                cell_iter->GetCellData()->SetItem("migration_direction_z", 0.0);
            }
            
            // Apply force with magnitude = base_speed * ECM_sensitivity
            // This represents the cell's motility biased by ECM structure
            c_vector<double, DIM> force = mBaseSpeed * mECMSensitivity * migration_direction;
            
            p_node->AddAppliedForceContribution(force);
        }
    }
    
    /**
     * Output force parameters to file
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ECMOrientationType>" << mECMOrientationType << "</ECMOrientationType>\n";
        *rParamsFile << "\t\t\t<BaseSpeed>" << mBaseSpeed << "</BaseSpeed>\n";
        *rParamsFile << "\t\t\t<ECMSensitivity>" << mECMSensitivity << "</ECMSensitivity>\n";
        *rParamsFile << "\t\t\t<Anisotropy>" << mAnisotropy << "</Anisotropy>\n";
        *rParamsFile << "\t\t\t<DomainWidth>" << mDomainWidth << "</DomainWidth>\n";
        
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#endif /* ECMCONTACTGUIDANCEFORCE_HPP_ */
