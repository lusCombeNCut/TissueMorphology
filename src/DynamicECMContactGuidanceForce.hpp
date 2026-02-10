/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

*/

#ifndef DYNAMICECMCONTACTGUIDANCEFORCE_HPP_
#define DYNAMICECMCONTACTGUIDANCEFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "DynamicECMField.hpp"
#include "SimulationTime.hpp"

/**
 * Dynamic ECM Contact Guidance Force
 * 
 * Improved version with:
 * 1. Dynamic ECM remodeling (cells align ECM through traction)
 * 2. ECM degradation and deposition
 * 3. Bidirectional fiber movement
 * 4. Perpendicular diffusion
 * 5. Random walk component
 */
template<unsigned DIM>
class DynamicECMContactGuidanceForce : public AbstractForce<DIM>
{
    friend class boost::serialization::access;

private:
    /** Pointer to dynamic ECM field */
    boost::shared_ptr<DynamicECMField> mpECMField;
    
    /** Base cell migration speed (µm/min) */
    double mBaseSpeed;
    
    /** ECM sensitivity parameter (0-1+) */
    double mECMSensitivity;
    
    /** Whether to enable ECM remodeling */
    bool mEnableRemodeling;
    
    /** Whether to enable ECM degradation */
    bool mEnableDegradation;
    
    /** Whether to enable ECM deposition */
    bool mEnableDeposition;
    
    /** Last update time for ECM field */
    double mLastUpdateTime;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        // Note: mpECMField not serialized - will need to recreate
        archive & mBaseSpeed;
        archive & mECMSensitivity;
        archive & mEnableRemodeling;
        archive & mEnableDegradation;
        archive & mEnableDeposition;
    }

public:
    /**
     * Constructor
     */
    DynamicECMContactGuidanceForce()
        : AbstractForce<DIM>(),
          mBaseSpeed(1.25),
          mECMSensitivity(1.0),
          mEnableRemodeling(true),
          mEnableDegradation(true),
          mEnableDeposition(false),
          mLastUpdateTime(0.0)
    {
        // ECM field will be set externally
    }
    
    /**
     * Destructor
     */
    virtual ~DynamicECMContactGuidanceForce()
    {
    }
    
    /**
     * Set the ECM field
     */
    void SetECMField(boost::shared_ptr<DynamicECMField> pECMField)
    {
        mpECMField = pECMField;
    }
    
    void SetBaseSpeed(double speed) { mBaseSpeed = speed; }
    void SetECMSensitivity(double sensitivity) { mECMSensitivity = sensitivity; }
    void SetEnableRemodeling(bool enable) { mEnableRemodeling = enable; }
    void SetEnableDegradation(bool enable) { mEnableDegradation = enable; }
    void SetEnableDeposition(bool enable) { mEnableDeposition = enable; }
    
    /**
     * Apply force and update ECM field
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
    {
        if (!mpECMField)
        {
            EXCEPTION("ECM field not set! Call SetECMField() before simulation.");
        }
        
        // Only works in 2D for now
        assert(DIM == 2);
        
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double current_time = SimulationTime::Instance()->GetTime();
        double dt = current_time - mLastUpdateTime;
        
        if (dt < 1e-10)
            dt = 1.0;  // Default timestep
        
        // First pass: apply forces and collect traction
        std::vector<std::pair<c_vector<double, 2>, c_vector<double, 2>>> cell_tractions;
        
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);
            c_vector<double, DIM> position = p_node->rGetLocation();
            
            // Get ECM properties at cell position
            double ecm_angle = mpECMField->GetFiberAngleAt(position);
            double ecm_density = mpECMField->GetDensityAt(position);
            double ecm_anisotropy = mpECMField->GetAnisotropyAt(position);
            
            // ECM fiber direction (unit vector)
            c_vector<double, DIM> ecm_fiber = zero_vector<double>(DIM);
            ecm_fiber[0] = cos(ecm_angle);
            ecm_fiber[1] = sin(ecm_angle);
            
            // Store for visualization
            cell_iter->GetCellData()->SetItem("ecm_angle", ecm_angle);
            cell_iter->GetCellData()->SetItem("ecm_density", ecm_density);
            cell_iter->GetCellData()->SetItem("ecm_orientation_x", ecm_fiber[0]);
            cell_iter->GetCellData()->SetItem("ecm_orientation_y", ecm_fiber[1]);
            
            // === REALISTIC CELL MIGRATION ===
            
            // 1. Bidirectional fiber movement
            double fiber_sign = (p_gen->ranf() < 0.5) ? -1.0 : 1.0;
            c_vector<double, DIM> fiber_motion = fiber_sign * ecm_fiber;
            
            // 2. Perpendicular diffusion
            c_vector<double, DIM> perp_direction;
            perp_direction[0] = -ecm_fiber[1];
            perp_direction[1] = ecm_fiber[0];
            double perp_magnitude = p_gen->ranf() * 2.0 - 1.0;
            c_vector<double, DIM> perp_motion = perp_magnitude * perp_direction;
            
            // 3. Random walk
            double random_angle = p_gen->ranf() * 2.0 * M_PI;
            c_vector<double, DIM> random_direction;
            random_direction[0] = cos(random_angle);
            random_direction[1] = sin(random_angle);
            
            // 4. Combine with anisotropy from ECM density and alignment
            double effective_anisotropy = ecm_anisotropy * ecm_density;
            double perpendicular_mobility = (1.0 - effective_anisotropy) * 0.3;
            double random_mobility = (1.0 - effective_anisotropy) * 0.5;
            
            c_vector<double, DIM> migration_direction = 
                effective_anisotropy * fiber_motion +
                perpendicular_mobility * perp_motion +
                random_mobility * random_direction;
            
            // Normalize
            double mag = norm_2(migration_direction);
            if (mag > 1e-10)
            {
                migration_direction /= mag;
            }
            
            // Store migration direction
            cell_iter->GetCellData()->SetItem("migration_direction_x", migration_direction[0]);
            cell_iter->GetCellData()->SetItem("migration_direction_y", migration_direction[1]);
            
            // Apply force (reduced in low-density ECM)
            double effective_speed = mBaseSpeed * mECMSensitivity * sqrt(ecm_density);
            c_vector<double, DIM> force = effective_speed * migration_direction;
            
            p_node->AddAppliedForceContribution(force);
            
            // Store traction for ECM remodeling
            cell_tractions.push_back(std::make_pair(position, force));
            
            // ECM degradation (cells degrade as they migrate)
            if (mEnableDegradation)
            {
                mpECMField->DegradeECM(position, dt);
            }
            
            // ECM deposition (optional - cells lay down ECM)
            if (mEnableDeposition)
            {
                double deposition_angle = atan2(migration_direction[1], migration_direction[0]);
                mpECMField->DepositECM(position, deposition_angle, dt);
            }
        }
        
        // Second pass: ECM remodeling from cell traction
        if (mEnableRemodeling)
        {
            for (auto& traction : cell_tractions)
            {
                mpECMField->ApplyCellTraction(traction.first, traction.second, dt);
            }
            
            // Diffuse ECM to smooth field
            mpECMField->DiffuseECM(dt);  // DISABLED: Testing without diffusion
        }
        
        mLastUpdateTime = current_time;
    }
    
    /**
     * Output force parameters
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<BaseSpeed>" << mBaseSpeed << "</BaseSpeed>\n";
        *rParamsFile << "\t\t\t<ECMSensitivity>" << mECMSensitivity << "</ECMSensitivity>\n";
        *rParamsFile << "\t\t\t<EnableRemodeling>" << mEnableRemodeling << "</EnableRemodeling>\n";
        *rParamsFile << "\t\t\t<EnableDegradation>" << mEnableDegradation << "</EnableDegradation>\n";
        *rParamsFile << "\t\t\t<EnableDeposition>" << mEnableDeposition << "</EnableDeposition>\n";
        
        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#endif /* DYNAMICECMCONTACTGUIDANCEFORCE_HPP_ */
