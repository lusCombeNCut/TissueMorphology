/*
 * ECMWiringHelpers.hpp
 *
 * Shared helpers for configuring ECM fields and forces.
 * Each helper sets parameters on an already-constructed object,
 * keeping the dimension-specific constructor in the runner.
 */
#ifndef ECMWIRINGHELPERS_HPP_
#define ECMWIRINGHELPERS_HPP_

#include "GhostNodeEcmField.hpp"
#include "GhostNodeEcmForce.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "ViscoelasticGhostNodeEcmForce.hpp"
#include "ECMConfinementForce.hpp"
#include "DynamicECMField.hpp"
#include "CryptBuddingParams.hpp"

// ===================================================================
// Ghost node ECM (elastic) — field configuration
// ===================================================================

template<unsigned DIM>
void ConfigureGhostField(boost::shared_ptr<GhostNodeEcmField<DIM>> pField,
                          const CryptBuddingParams& p, double gnRest)
{
    pField->SetGhostGhostStiffness(p.cellGhostSpringStiffness);
    pField->SetGhostDamping(p.ghostDamping);
    pField->SetDegradationRate(p.ecmDegradationRate);
    pField->SetRemovalThreshold(p.ghostRemovalThreshold);
    pField->SetFibreRemodelingRate(p.ghostFibreRemodelingRate);
    pField->SetAnisotropyStrength(p.ghostAnisotropyStrength);
    pField->SetGhostRestLength(gnRest);
}

// ===================================================================
// Ghost node ECM (elastic) — force configuration
// ===================================================================

template<unsigned DIM>
void ConfigureGhostForce(boost::shared_ptr<GhostNodeEcmForce<DIM>> pForce,
                          boost::shared_ptr<GhostNodeEcmField<DIM>> pField,
                          const CryptBuddingParams& p)
{
    pForce->SetGhostField(pField);
    pForce->SetCellGhostStiffness(p.cellGhostSpringStiffness);
    pForce->SetCellGhostRestLength(p.ghostCellGhostRestLength);
    pForce->SetCellGhostCutoff(p.ghostCellGhostCutoff);
    pForce->SetDegradationEnabled(true);
    pForce->SetRemodelingEnabled(p.enableEcmGuidance);
    pForce->SetTrackCenter(true);
    pForce->SetBoundaryExpansionEnabled(p.enableGhostBoundaryExpansion);
    pForce->SetRemovalCheckInterval(p.ghostRemovalCheckInterval);
}

// ===================================================================
// Viscoelastic ghost node ECM — field configuration
// ===================================================================

template<unsigned DIM>
void ConfigureViscoelasticField(boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> pField,
                                 const CryptBuddingParams& p)
{
    pField->SetRelaxedStiffness(p.ghostE0);
    pField->SetRelaxationModulus(p.ghostE1);
    pField->SetRelaxationTime(p.ghostRelaxationTime);
    pField->SetGhostDamping(p.ghostDamping);
    pField->SetDegradationRate(p.ecmDegradationRate);
    pField->SetRemovalThreshold(p.ghostRemovalThreshold);
    pField->SetFibreRemodelingRate(p.ghostFibreRemodelingRate);
    pField->SetAnisotropyStrength(p.ghostAnisotropyStrength);
}

// ===================================================================
// Viscoelastic ghost node ECM — force configuration
// ===================================================================

template<unsigned DIM>
void ConfigureViscoelasticForce(boost::shared_ptr<ViscoelasticGhostNodeEcmForce<DIM>> pForce,
                                 boost::shared_ptr<ViscoelasticGhostNodeEcmField<DIM>> pField,
                                 const CryptBuddingParams& p)
{
    pForce->SetGhostField(pField);
    pForce->SetCellGhostStiffness(p.cellGhostSpringStiffness);
    pForce->SetCellGhostRestLength(p.ghostCellGhostRestLength);
    pForce->SetCellGhostCutoff(p.ghostCellGhostCutoff);
    pForce->SetDegradationEnabled(true);
    pForce->SetRemodelingEnabled(p.enableEcmGuidance);
    pForce->SetTrackCenter(true);
    pForce->SetBoundaryExpansionEnabled(p.enableGhostBoundaryExpansion);
    pForce->SetRemovalCheckInterval(p.ghostRemovalCheckInterval);
}

// ===================================================================
// Grid-based ECM confinement — force configuration (2D only)
// ===================================================================

inline void ConfigureGridECMForce(boost::shared_ptr<ECMConfinementForce<2>> pForce,
                                   boost::shared_ptr<DynamicECMField> pField,
                                   const CryptBuddingParams& p)
{
    pForce->SetECMField(pField);
    pForce->SetConfinementStiffness(p.cellGhostSpringStiffness);
    pForce->SetEcmSpringRestLength(p.ecmSpringRestLength);
    pForce->SetEcmInteractionCutoff(p.ecmInteractionCutoff);
    pForce->SetDegradationEnabled(true);
    pForce->SetRemodelingEnabled(p.enableEcmGuidance);
    pForce->SetTrackCenter(true);
}

#endif // ECMWIRINGHELPERS_HPP_
