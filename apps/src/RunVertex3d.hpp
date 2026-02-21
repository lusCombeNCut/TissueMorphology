/*
 * RunVertex3d.hpp — 3D Vertex-Based (OrganoidChaste) crypt budding model runner
 *
 * Uses the finite-thickness monolayer vertex mesh from OrganoidChaste.
 * Per-cell-type surface tensions are set via mutation states:
 *   WildType   → Stem cells    (softer: gamma × gammaStemScale)
 *   ApcOneHit  → Transit-amplifying (baseline: gamma × gammaTransitScale)
 *   ApcTwoHit  → Differentiated/Paneth (stiffer: gamma × gammaDiffScale)
 *
 * Two-phase simulation:
 *   Phase 1: Relaxation at small dt (all cells differentiated, no lumen pressure)
 *   Phase 2: Growth at dtGrow with lumen pressure enabled
 */
#ifndef RUNVERTEX3D_HPP_
#define RUNVERTEX3D_HPP_

#include <string>
#include <map>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"

#include "MutableMonolayerVertexMesh.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"
#include "FiniteThicknessRandomizedSphereMeshGenerator.hpp"
#include "FiniteThicknessSimulation3d.hpp"

#include "ContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "ApcOneHitCellMutationState.hpp"
#include "ApcTwoHitCellMutationState.hpp"

#include "SurfaceTensionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "LumenPressureSubForce.hpp"

#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"

#include "VolumeTrackingModifier.hpp"
#include "GeometricalTargetVolumeModifier.hpp"
#include "CellVolumesWriter.hpp"
#include "CellThicknessWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "FaceTypeWriter.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"

void RunVertex3d(const CryptBuddingParams& p, const std::string& outputDir)
{
    RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

    // ------------------------------------------------------------------
    // Mesh generation
    // ------------------------------------------------------------------
    double height = 2.0 / 3.0 / sqrt(3.0)
                    * cbrt((9.0 / 2.0) * (9.0 / 2.0))
                    * cbrt((p.gammaApical + p.gammaBasal) / p.gammaLateral
                           * (p.gammaApical + p.gammaBasal) / p.gammaLateral)
                    * 1.0;
    double t1_length = 0.66 / cbrt(3.0 * 3.0
                       * (1.0 + p.gammaLateral) / p.gammaLateral);

    std::cout << "  Cell height:    " << height << std::endl;
    std::cout << "  T1 threshold:   " << t1_length << std::endl;

    FiniteThicknessRandomizedSphereMeshGenerator generator(
        p.numCells3dVertex, t1_length, 0.001, height, p.sphereRadius3dVertex);
    MutableMonolayerVertexMesh<3, 3>* pMesh = generator.GetMesh();
    pMesh->SetProtorosetteFormationProbability(1.0);
    pMesh->SetProtorosetteResolutionProbabilityPerTimestep(0.1);

    if (pMesh->GetNumElements() != p.numCells3dVertex)
        std::cerr << "WARNING: mesh has " << pMesh->GetNumElements()
                  << " elements, expected " << p.numCells3dVertex << std::endl;

    double outerR = p.sphereRadius3dVertex + height;
    double avgVol = 4.0 / 3.0 * M_PI
                    * (outerR * outerR * outerR
                       - p.sphereRadius3dVertex * p.sphereRadius3dVertex * p.sphereRadius3dVertex)
                    / p.numCells3dVertex;

    // ------------------------------------------------------------------
    // Cell creation with per-cell-type mutation states
    // ------------------------------------------------------------------
    // Use different mutation states to encode cell type for SurfaceTensionForce:
    //   WildType   = stem cells     (softer, gamma × gammaStemScale)
    //   ApcOneHit  = transit-amplifying (baseline, gamma × gammaTransitScale)
    //   ApcTwoHit  = differentiated/Paneth (stiffer, gamma × gammaDiffScale)
    MAKE_PTR(WildTypeCellMutationState, p_stem_mut);
    MAKE_PTR(ApcOneHitCellMutationState, p_ta_mut);
    MAKE_PTR(ApcTwoHitCellMutationState, p_diff_mut);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

    std::vector<CellPtr> cells;

    for (unsigned i = 0; i < p.numCells3dVertex; i++)
    {
        ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
        p_cycle->SetDimension(3);
        p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        p_cycle->SetEquilibriumVolume(avgVol);

        // Uniform random cell type assignment across the organoid surface
        double u = RandomNumberGenerator::Instance()->ranf();
        boost::shared_ptr<AbstractCellMutationState> p_mut;
        boost::shared_ptr<AbstractCellProliferativeType> p_type;
        double type_id;

        if (u < p.stemFraction)
        {
            p_mut = p_stem_mut;
            p_type = p_stem;
            type_id = 0.0;
        }
        else if (u < p.stemFraction + p.transitFraction)
        {
            p_mut = p_ta_mut;
            p_type = p_ta;
            type_id = 1.0;
        }
        else
        {
            p_mut = p_diff_mut;
            p_type = p_diff;
            type_id = 2.0;
        }

        CellPtr p_cell(new Cell(p_mut, p_cycle));
        p_cell->SetCellProliferativeType(p_type);
        p_cell->GetCellData()->SetItem("cell_type_id", type_id);

        p_cell->SetBirthTime(-RandomNumberGenerator::Instance()->ranf() * 10.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->GetCellData()->SetItem("volume", avgVol);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);
        cells.push_back(p_cell);
    }

    // ------------------------------------------------------------------
    // Population
    // ------------------------------------------------------------------
    MonolayerVertexBasedCellPopulation<3> population(*pMesh, cells);
    population.SetOutputCellRearrangementLocations(false);
    population.SetRestrictVertexMovementBoolean(false);
    population.SetDoInitialVolumeRelaxation(true);

    population.AddCellWriter<CellVolumesWriter>();
    population.AddCellWriter<CellThicknessWriter>();
    population.AddCellWriter<CellProliferativeTypesWriter>();
    population.AddCellWriter<CellAgesWriter>();
    population.AddFaceWriter<FaceTypeWriter>();

    // ------------------------------------------------------------------
    // Simulator & forces
    // ------------------------------------------------------------------
    FiniteThicknessSimulation3d simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);

    // Surface tension force — per-cell-type (differential adhesion) or uniform
    MAKE_PTR(SurfaceTensionForce<3>, p_tension);
    if (p.enableDifferentialAdhesion)
    {
        // Per-cell-type tensions via mutation state
        std::map<boost::shared_ptr<AbstractCellMutationState>, std::array<double, 3>> mut_tension_map;
        mut_tension_map[p_stem_mut] = {{p.gammaApical * p.gammaStemScale,
                                        p.gammaBasal  * p.gammaStemScale,
                                        p.gammaLateral * p.gammaStemScale}};
        mut_tension_map[p_ta_mut]   = {{p.gammaApical * p.gammaTransitScale,
                                        p.gammaBasal  * p.gammaTransitScale,
                                        p.gammaLateral * p.gammaTransitScale}};
        mut_tension_map[p_diff_mut] = {{p.gammaApical * p.gammaDiffScale,
                                        p.gammaBasal  * p.gammaDiffScale,
                                        p.gammaLateral * p.gammaDiffScale}};
        p_tension->SetSurfaceTensionParametersByMutation(mut_tension_map);
        p_tension->UpdateSurfaceTensionsByMutation(&population);

        std::cout << "  Surface tensions (apical, basal, lateral) — DIFFERENTIAL:" << std::endl;
        std::cout << "    Stem:    " << p.gammaApical * p.gammaStemScale << ", "
                  << p.gammaBasal * p.gammaStemScale << ", "
                  << p.gammaLateral * p.gammaStemScale << std::endl;
        std::cout << "    Transit: " << p.gammaApical * p.gammaTransitScale << ", "
                  << p.gammaBasal * p.gammaTransitScale << ", "
                  << p.gammaLateral * p.gammaTransitScale << std::endl;
        std::cout << "    Paneth:  " << p.gammaApical * p.gammaDiffScale << ", "
                  << p.gammaBasal * p.gammaDiffScale << ", "
                  << p.gammaLateral * p.gammaDiffScale << std::endl;
    }
    else
    {
        // Uniform tensions — all cells get the same gamma values
        p_tension->CreateSurfaceTensionParametersForCells(p.gammaApical,
                                                          p.gammaBasal,
                                                          p.gammaLateral,
                                                          pMesh);

        std::cout << "  Surface tensions (apical, basal, lateral) — UNIFORM:" << std::endl;
        std::cout << "    All:     " << p.gammaApical << ", "
                  << p.gammaBasal << ", "
                  << p.gammaLateral << std::endl;
    }
    p_tension->SetSimulatedAnnealingParameters(0.0, 50.0, 0.0);
    p_tension->SetSimulationInstance(&simulator);
    simulator.AddForce(p_tension);

    // NOTE: LumenPressureSubForce is added in Phase 2 only.
    // Adding it during relaxation destabilises the mesh → DIVERGED_ITS.

    MAKE_PTR(BasementMembraneForce<3>, p_bm);
    p_bm->SetBasementMembraneParameter(p.bmStiffnessVertex);
    p_bm->SetTargetRadius(p.organoidRadius3d + height + p.bmOffset3dVertex);
    simulator.AddForce(p_bm);

    // ------------------------------------------------------------------
    // ECM guidance field (optional)
    // ------------------------------------------------------------------
    boost::shared_ptr<DynamicECMField3d> pEcmField;
    if (p.enableEcmGuidance)
    {
        pEcmField.reset(new DynamicECMField3d(
            "radial", p.ecmGridSpacing,
            -p.ecmDomainHalf, p.ecmDomainHalf,
            -p.ecmDomainHalf, p.ecmDomainHalf,
            -p.ecmDomainHalf, p.ecmDomainHalf));
        pEcmField->SetDegradationRate(0.002);
        pEcmField->SetRemodelingRate(0.05);
        pEcmField->SetDepositionRate(0.0003);
    }

    // ------------------------------------------------------------------
    // Modifiers
    // ------------------------------------------------------------------
    MAKE_PTR(VolumeTrackingModifier<3>, p_vol);
    simulator.AddSimulationModifier(p_vol);

    MAKE_PTR_ARGS(GeometricalTargetVolumeModifier<3>, p_gvol, (&population));
    p_gvol->SetGrowthDuration(0.0);
    p_gvol->SetT1AdaptationDuration(0.100);
    p_gvol->SetReferenceTargetVolume(avgVol);
    simulator.AddSimulationModifier(p_gvol);

    double totalSimTime = p.enableRelaxation ? (p.relaxationTime + p.endTime) : p.endTime;
    boost::shared_ptr<CryptBuddingSummaryModifier<3>> p_summary(
        new CryptBuddingSummaryModifier<3>(p.ecmStiffness, p.samplingMultiple,
                                           totalSimTime, avgVol));
    simulator.AddSimulationModifier(p_summary);

    // Note: no sloughing killers for vertex3d because
    // MutableMonolayerVertexMesh does not support DeleteElementPriorToReMesh in 3D.

    // ------------------------------------------------------------------
    // Two-phase solve
    // ------------------------------------------------------------------
    if (p.enableRelaxation)
    {
        // Phase 1: Relaxation — all cells differentiated, no lumen pressure
        std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
        for (AbstractCellPopulation<3>::Iterator it = population.Begin();
             it != population.End(); ++it)
        {
            origTypes[*it] = it->GetCellProliferativeType();
            it->SetCellProliferativeType(p_diff);
        }

        simulator.SetEndTime(p.relaxationTime);
        std::cout << "--- Phase 1: Relaxation (" << p.relaxationTime << " hours) ---" << std::endl;
        simulator.Solve();
        std::cout << "Relaxation complete. Cells: " << population.GetNumRealCells() << std::endl;

        // Restore original cell types
        for (AbstractCellPopulation<3>::Iterator it = population.Begin();
             it != population.End(); ++it)
        {
            if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
        }

        // Phase 2: Growth — use dtGrow (reduced for stability)
        double dt_grow = p.dtGrow;
        unsigned sampling_grow = static_cast<unsigned>(std::max(1.0, 1.0 / dt_grow / 10.0));

        simulator.SetEndTime(p.relaxationTime + p.endTime);
        simulator.SetSamplingTimestepMultiple(sampling_grow);
        simulator.SetDt(dt_grow);

        p_tension->SetSimulatedAnnealingParameters(0.003, 1900000.0, 1.0);
        p_tension->SetSimulationInstance(&simulator);
        p_tension->SetPerformActiveT1Swaps(true);
        p_tension->SetT1TransitionParameters(2.0, false);

        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);

        // Add lumen pressure in Phase 2 (not during relaxation)
        if (p.enableLumenPressure)
        {
            MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_sub, (p.lumenPressure));
            simulator.AddForce(p_lumen_sub);
        }

        // Add ECM guidance in Phase 2
        if (p.enableEcmGuidance && pEcmField)
        {
            MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm);
            p_ecm->SetECMField(pEcmField);
            p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
            p_ecm->SetECMSensitivity(1.0);
            p_ecm->SetEnableDegradation(true);
            p_ecm->SetEnableRemodeling(true);
            p_ecm->SetEnableDeposition(false);
            simulator.AddForce(p_ecm);
        }

        std::cout << "--- Phase 2: Growth (" << p.endTime << " hours, dt=" << dt_grow << ") ---" << std::endl;
        simulator.Solve();
    }
    else
    {
        // No relaxation — still need to configure growth-phase settings
        double dt_grow = p.dtGrow;
        unsigned sampling_grow = static_cast<unsigned>(std::max(1.0, 1.0 / dt_grow / 10.0));

        simulator.SetEndTime(p.endTime);
        simulator.SetSamplingTimestepMultiple(sampling_grow);
        simulator.SetDt(dt_grow);

        p_tension->SetSimulatedAnnealingParameters(0.003, 1900000.0, 1.0);
        p_tension->SetSimulationInstance(&simulator);
        p_tension->SetPerformActiveT1Swaps(true);
        p_tension->SetT1TransitionParameters(2.0, false);

        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);

        // Add lumen pressure
        if (p.enableLumenPressure)
        {
            MAKE_PTR_ARGS(LumenPressureSubForce<3>, p_lumen_sub, (p.lumenPressure));
            simulator.AddForce(p_lumen_sub);
        }

        // Add ECM guidance
        if (p.enableEcmGuidance && pEcmField)
        {
            MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm);
            p_ecm->SetECMField(pEcmField);
            p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
            p_ecm->SetECMSensitivity(1.0);
            p_ecm->SetEnableDegradation(true);
            p_ecm->SetEnableRemodeling(true);
            p_ecm->SetEnableDeposition(false);
            simulator.AddForce(p_ecm);
        }

        std::cout << "--- Single-phase Growth (" << p.endTime << " hours, dt=" << dt_grow << ") ---" << std::endl;
        simulator.Solve();
    }

    unsigned final_cells = population.GetNumRealCells();
    std::cout << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
    if (final_cells == 0)
        EXCEPTION("Simulation ended with zero cells");
}

#endif // RUNVERTEX3D_HPP_
