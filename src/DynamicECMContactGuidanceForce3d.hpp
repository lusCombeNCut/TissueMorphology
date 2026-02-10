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

#ifndef DYNAMICECMCONTACTGUIDANCEFORCE3D_HPP_
#define DYNAMICECMCONTACTGUIDANCEFORCE3D_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"
#include "DynamicECMField3d.hpp"

/**
 * 3D Dynamic ECM Contact Guidance Force
 * 
 * Cells migrate through a 3D ECM field where local fiber direction,
 * density, and anisotropy determine the bias and speed of cell migration.
 * 
 * Force on cell i at position r_i:
 * 
 *   F_i = v0 * s * sqrt(rho(r_i)) * n_migration
 * 
 * where n_migration is a weighted combination of:
 *   1. Fiber-aligned motion:    alpha_eff * (±f_hat)
 *   2. Perpendicular diffusion: (1 - alpha_eff) * 0.3 * xi * f_perp
 *   3. Random walk:             (1 - alpha_eff) * 0.5 * r_rand
 * 
 * The force also drives ECM remodeling:
 *   - Cells degrade ECM at their location (MMP activity)
 *   - Cell traction reorients ECM fibers toward pulling direction
 *   - Cells optionally deposit new ECM aligned with their motion
 * 
 * References: Painter 2009; Metzcar et al. 2025
 */
class DynamicECMContactGuidanceForce3d : public AbstractForce<3>
{
    friend class boost::serialization::access;

private:
    /** Pointer to the 3D ECM field */
    boost::shared_ptr<DynamicECMField3d> mpECMField;
    
    /** Base cell migration speed (µm/min) */
    double mBaseSpeed;
    
    /** ECM sensitivity multiplier */
    double mECMSensitivity;
    
    /** Enable/disable ECM remodeling by cell traction */
    bool mEnableRemodeling;
    
    /** Enable/disable ECM degradation by cells */
    bool mEnableDegradation;
    
    /** Enable/disable ECM deposition by cells */
    bool mEnableDeposition;
    
    /** Time of last force evaluation (for dt calculation) */
    double mLastUpdateTime;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<3> >(*this);
    }

public:
    /**
     * Constructor.
     */
    DynamicECMContactGuidanceForce3d()
        : AbstractForce<3>(),
          mBaseSpeed(1.25),
          mECMSensitivity(1.0),
          mEnableRemodeling(true),
          mEnableDegradation(true),
          mEnableDeposition(false),
          mLastUpdateTime(0.0)
    {
    }
    
    /**
     * Destructor.
     */
    ~DynamicECMContactGuidanceForce3d() {}
    
    // ========================================================================
    // Setters
    // ========================================================================
    
    void SetECMField(boost::shared_ptr<DynamicECMField3d> pECMField) { mpECMField = pECMField; }
    void SetBaseSpeed(double speed) { mBaseSpeed = speed; }
    void SetECMSensitivity(double sensitivity) { mECMSensitivity = sensitivity; }
    void SetEnableRemodeling(bool enable) { mEnableRemodeling = enable; }
    void SetEnableDegradation(bool enable) { mEnableDegradation = enable; }
    void SetEnableDeposition(bool enable) { mEnableDeposition = enable; }
    
    /**
     * Compute 3D ECM contact guidance forces on all cells.
     * 
     * For each cell:
     *   1. Query local ECM fiber direction, density, anisotropy
     *   2. Compute composite migration direction (fiber + perp + random)
     *   3. Scale by density-dependent effective speed
     *   4. Apply force and record traction for ECM remodeling
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation)
    {
        if (!mpECMField)
        {
            EXCEPTION("ECM field not set! Call SetECMField() before simulation.");
        }
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double current_time = SimulationTime::Instance()->GetTime();
        double dt = current_time - mLastUpdateTime;
        if (dt < 1e-10) dt = 1.0;
        
        // Collect tractions for second-pass remodeling
        std::vector<std::pair<c_vector<double, 3>, c_vector<double, 3>>> cell_tractions;
        
        for (typename AbstractCellPopulation<3>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<3>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, 3> position = p_node->rGetLocation();
            
            // ----- Query ECM at cell position -----
            c_vector<double, 3> ecm_fiber = mpECMField->GetFiberDirectionAt(position);
            double ecm_density = mpECMField->GetDensityAt(position);
            double ecm_anisotropy = mpECMField->GetAnisotropyAt(position);
            
            // Store for visualization
            cell_iter->GetCellData()->SetItem("ecm_density", ecm_density);
            cell_iter->GetCellData()->SetItem("ecm_fiber_x", ecm_fiber[0]);
            cell_iter->GetCellData()->SetItem("ecm_fiber_y", ecm_fiber[1]);
            cell_iter->GetCellData()->SetItem("ecm_fiber_z", ecm_fiber[2]);
            
            // ----- 1. Bidirectional fiber-aligned motion -----
            double fiber_sign = (p_gen->ranf() < 0.5) ? -1.0 : 1.0;
            c_vector<double, 3> fiber_motion = fiber_sign * ecm_fiber;
            
            // ----- 2. Perpendicular diffusion -----
            // Generate a random vector perpendicular to the fiber direction
            // by cross-product with a random vector, then normalize
            c_vector<double, 3> random_vec;
            random_vec[0] = p_gen->ranf() - 0.5;
            random_vec[1] = p_gen->ranf() - 0.5;
            random_vec[2] = p_gen->ranf() - 0.5;
            
            // Cross product: perp = fiber × random_vec
            c_vector<double, 3> perp_direction;
            perp_direction[0] = ecm_fiber[1] * random_vec[2] - ecm_fiber[2] * random_vec[1];
            perp_direction[1] = ecm_fiber[2] * random_vec[0] - ecm_fiber[0] * random_vec[2];
            perp_direction[2] = ecm_fiber[0] * random_vec[1] - ecm_fiber[1] * random_vec[0];
            
            double perp_mag = norm_2(perp_direction);
            if (perp_mag > 1e-10)
            {
                perp_direction /= perp_mag;
            }
            else
            {
                // Fiber was parallel to random_vec; try again with different axis
                perp_direction[0] = -ecm_fiber[1];
                perp_direction[1] = ecm_fiber[0];
                perp_direction[2] = 0.0;
                perp_mag = norm_2(perp_direction);
                if (perp_mag > 1e-10) perp_direction /= perp_mag;
            }
            
            double perp_magnitude = p_gen->ranf() * 2.0 - 1.0;
            c_vector<double, 3> perp_motion = perp_magnitude * perp_direction;
            
            // ----- 3. Random walk (uniform on sphere) -----
            c_vector<double, 3> random_direction;
            double v1, v2, s;
            do {
                v1 = 2.0 * p_gen->ranf() - 1.0;
                v2 = 2.0 * p_gen->ranf() - 1.0;
                s = v1 * v1 + v2 * v2;
            } while (s >= 1.0 || s < 1e-10);
            random_direction[0] = 2.0 * v1 * sqrt(1.0 - s);
            random_direction[1] = 2.0 * v2 * sqrt(1.0 - s);
            random_direction[2] = 1.0 - 2.0 * s;
            
            // ----- 4. Combine using effective anisotropy -----
            double effective_anisotropy = ecm_anisotropy * ecm_density;
            double perpendicular_mobility = (1.0 - effective_anisotropy) * 0.3;
            double random_mobility = (1.0 - effective_anisotropy) * 0.5;
            
            c_vector<double, 3> migration_direction =
                effective_anisotropy * fiber_motion +
                perpendicular_mobility * perp_motion +
                random_mobility * random_direction;
            
            double mag = norm_2(migration_direction);
            if (mag > 1e-10)
            {
                migration_direction /= mag;
            }
            
            // Store migration direction for visualization
            cell_iter->GetCellData()->SetItem("migration_dir_x", migration_direction[0]);
            cell_iter->GetCellData()->SetItem("migration_dir_y", migration_direction[1]);
            cell_iter->GetCellData()->SetItem("migration_dir_z", migration_direction[2]);
            
            // ----- Apply force -----
            double effective_speed = mBaseSpeed * mECMSensitivity * sqrt(ecm_density);
            c_vector<double, 3> force = effective_speed * migration_direction;
            
            p_node->AddAppliedForceContribution(force);
            
            // Store traction for ECM remodeling
            cell_tractions.push_back(std::make_pair(position, force));
            
            // ECM degradation (MMP activity at cell location)
            if (mEnableDegradation)
            {
                mpECMField->DegradeECM(position, dt);
            }
            
            // ECM deposition (optional)
            if (mEnableDeposition)
            {
                mpECMField->DepositECM(position, migration_direction, dt);
            }
        }
        
        // Second pass: ECM remodeling from cell traction
        if (mEnableRemodeling)
        {
            for (auto& traction : cell_tractions)
            {
                mpECMField->ApplyCellTraction(traction.first, traction.second, dt);
            }
            mpECMField->DiffuseECM(dt);
        }
        
        mLastUpdateTime = current_time;
    }
    
    /**
     * Output force parameters (required by AbstractForce interface).
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<BaseSpeed>" << mBaseSpeed << "</BaseSpeed>\n";
        *rParamsFile << "\t\t\t<ECMSensitivity>" << mECMSensitivity << "</ECMSensitivity>\n";
        *rParamsFile << "\t\t\t<EnableRemodeling>" << mEnableRemodeling << "</EnableRemodeling>\n";
        *rParamsFile << "\t\t\t<EnableDegradation>" << mEnableDegradation << "</EnableDegradation>\n";
        *rParamsFile << "\t\t\t<EnableDeposition>" << mEnableDeposition << "</EnableDeposition>\n";
        
        AbstractForce<3>::OutputForceParameters(rParamsFile);
    }
};

// Serialization
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(DynamicECMContactGuidanceForce3d)

#endif /* DYNAMICECMCONTACTGUIDANCEFORCE3D_HPP_ */
