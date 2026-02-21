/*
 * RunNode3d.hpp — 3D Node-Based crypt budding model runner
 */
#ifndef RUNNODE3D_HPP_
#define RUNNODE3D_HPP_

#include <string>
#include <map>
#include <iostream>
#include <cmath>

#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"

#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"

#include "ContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "CellPolarityForce.hpp"

#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"

#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPolarityWriter.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"

void RunNode3d(const CryptBuddingParams& p, const std::string& outputDir)
{
    RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

    std::vector<Node<3>*> nodes;
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double golden = (1.0 + sqrt(5.0)) / 2.0;

    for (unsigned i = 0; i < p.numCells3dNode; i++)
    {
        double theta = 2.0 * M_PI * i / golden;
        double phi = acos(1.0 - 2.0 * (i + 0.5) / p.numCells3dNode);
        double r = p.organoidRadius3d + (p_gen->ranf() - 0.5) * p.shellThickness3d;
        double x = r * sin(phi) * cos(theta);
        double y = r * sin(phi) * sin(theta);
        double z = r * cos(phi);
        nodes.push_back(new Node<3>(i, false, x, y, z));
    }

    NodesOnlyMesh<3> mesh;
    mesh.ConstructNodesWithoutMesh(nodes, p.interactionCutoff3d);
    for (unsigned i = 0; i < nodes.size(); i++) delete nodes[i];

    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

    for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
    {
        ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
        p_cycle->SetDimension(3);
        p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        p_cycle->SetEquilibriumVolume(1.0);

        CellPtr p_cell(new Cell(p_state, p_cycle));
        AssignCellTypeByFraction(p_cell, 0.0, p_stem, p_ta, p_diff,
                                 p.stemFraction, p.transitFraction);

        p_cell->SetBirthTime(-p_gen->ranf() * 18.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->GetCellData()->SetItem("volume", 1.0);
        p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessNode);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);

        // Initialize cell polarity (spherical coords pointing radially outward)
        // Position on sphere gives natural polarity direction
        double theta_cell = acos(1.0 - 2.0 * (i + 0.5) / p.numCells3dNode);  // phi from Fibonacci
        double phi_cell = 2.0 * M_PI * i / golden;  // theta from Fibonacci (recalculated)
        // Polarity theta = angle from z-axis (same as position theta)
        // Polarity phi = azimuthal angle (same as position phi)
        p_cell->GetCellData()->SetItem("polarity_theta", theta_cell);
        p_cell->GetCellData()->SetItem("polarity_phi", phi_cell);

        cells.push_back(p_cell);
    }

    NodeBasedCellPopulation<3> population(mesh, cells);
    population.SetAbsoluteMovementThreshold(50.0);
    population.AddCellWriter<CellIdWriter>();
    population.AddCellWriter<CellAgesWriter>();
    population.AddCellWriter<CellPolarityWriter>();
    population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    OffLatticeSimulation<3> simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
    simulator.SetEndTime(p.endTime);

    if (p.enableDifferentialAdhesion)
    {
        MAKE_PTR(DifferentialAdhesionForce<3>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.interactionCutoff3d);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        p_spring->SetApicalApicalAdhesion(p.apicalApicalAdhesion);
        p_spring->SetBasalBasalAdhesion(p.basalBasalAdhesion);
        p_spring->SetApicalBasalAdhesion(p.apicalBasalAdhesion);
        simulator.AddForce(p_spring);
    }
    else
    {
        MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.interactionCutoff3d);
        simulator.AddForce(p_spring);
    }

    MAKE_PTR(BasementMembraneForce<3>, p_bm);
    p_bm->SetBasementMembraneParameter(p.bmStiffnessNode);
    p_bm->SetTargetRadius(p.bmRadius3d);
    p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);
    simulator.AddForce(p_bm);

    if (p.enableLumenPressure)
    {
        MAKE_PTR(LumenPressureForce<3>, p_lumen);
        p_lumen->SetPressureStrength(p.lumenPressure);
        p_lumen->SetLumenEquilibriumRadius(p.lumenEqRadius3d);
        p_lumen->SetTrackCenter(true);
        simulator.AddForce(p_lumen);
    }

    if (p.enableApicalConstriction)
    {
        MAKE_PTR(ApicalConstrictionForce<3>, p_ac);
        p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
        simulator.AddForce(p_ac);
    }

    // Cell polarity force for monolayer maintenance (ya||a-style)
    if (p.enableCellPolarity)
    {
        MAKE_PTR(CellPolarityForce<3>, p_polarity);
        p_polarity->SetBendingStrength(p.polarityBendingStrength);
        p_polarity->SetPolarityAlignmentStrength(p.polarityAlignmentStrength);
        p_polarity->SetInteractionCutoff(p.interactionCutoff3d);
        p_polarity->SetInitializeRadially(true);
        simulator.AddForce(p_polarity);
    }

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

        MAKE_PTR(DynamicECMContactGuidanceForce3d, p_ecm);
        p_ecm->SetECMField(pEcmField);
        p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
        p_ecm->SetECMSensitivity(1.0);
        p_ecm->SetEnableDegradation(true);
        p_ecm->SetEnableRemodeling(true);
        p_ecm->SetEnableDeposition(false);
        simulator.AddForce(p_ecm);
    }

    if (p.enableSloughing)
    {
        AddBoundingBoxKillers<3>(simulator, population,
                                p.organoidRadius3d * p.sloughRadiusFactor);
    }

    MAKE_PTR(VolumeTrackingModifier<3>, p_vol);
    simulator.AddSimulationModifier(p_vol);

    double totalSimTime = p.enableRelaxation ? (p.relaxationTime + p.endTime) : p.endTime;
    boost::shared_ptr<CryptBuddingSummaryModifier<3>> p_summary(
        new CryptBuddingSummaryModifier<3>(p.ecmStiffness, p.samplingMultiple,
                                           totalSimTime));
    simulator.AddSimulationModifier(p_summary);

    if (p.enableRelaxation)
    {
        std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
        for (AbstractCellPopulation<3>::Iterator it = population.Begin();
             it != population.End(); ++it)
        {
            origTypes[*it] = it->GetCellProliferativeType();
            it->SetCellProliferativeType(p_diff);
        }

        simulator.SetEndTime(p.relaxationTime);
        std::cout << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
        simulator.Solve();

        for (AbstractCellPopulation<3>::Iterator it = population.Begin();
             it != population.End(); ++it)
        {
            if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
        }

        simulator.SetEndTime(p.relaxationTime + p.endTime);
        std::cout << "--- Phase 2: Growth (" << p.endTime << "h) ---" << std::endl;
        simulator.Solve();
    }
    else
    {
        simulator.Solve();
    }

    unsigned final_cells = population.GetNumRealCells();
    std::cout << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells << std::endl;
    if (final_cells == 0)
        EXCEPTION("Simulation ended with zero cells");
}

#endif // RUNNODE3D_HPP_
