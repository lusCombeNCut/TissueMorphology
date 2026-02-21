/*
 * RunNode2d.hpp — 2D Node-Based crypt budding model runner
 */
#ifndef RUNNODE2D_HPP_
#define RUNNODE2D_HPP_

#include <string>
#include <map>
#include <iostream>

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
#include "CurvatureBendingForce.hpp"
#include "CellPolarityForce.hpp"
#include "RingTopologyTracker.hpp"

#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPolarityWriter.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"

void RunNode2d(const CryptBuddingParams& p, const std::string& outputDir)
{
    RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

    std::vector<Node<2>*> nodes;
    for (unsigned i = 0; i < p.numCells2dNode; i++)
    {
        double theta = 2.0 * M_PI * i / p.numCells2dNode;
        double r_noise = p.organoidRadius2d +
            (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.3;
        double x = r_noise * cos(theta);
        double y = r_noise * sin(theta);
        nodes.push_back(new Node<2>(i, false, x, y));
    }

    NodesOnlyMesh<2> mesh;
    mesh.ConstructNodesWithoutMesh(nodes, p.interactionCutoff2d);
    for (unsigned i = 0; i < nodes.size(); i++) delete nodes[i];

    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
    {
        ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
        p_cycle->SetDimension(2);
        p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        p_cycle->SetEquilibriumVolume(1.0);

        CellPtr p_cell(new Cell(p_state, p_cycle));
        AssignCellTypeByFraction(p_cell, 0.0, p_stem, p_ta, p_diff,
                                 p.stemFraction, p.transitFraction);

        p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->GetCellData()->SetItem("volume", 1.0);
        p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessNode);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);

        // Initialize cell polarity (2D: angle from x-axis pointing radially outward)
        double theta_cell = 2.0 * M_PI * i / p.numCells2dNode;
        p_cell->GetCellData()->SetItem("polarity_theta", theta_cell);
        p_cell->GetCellData()->SetItem("polarity_phi", 0.0);  // unused in 2D

        cells.push_back(p_cell);
    }

    NodeBasedCellPopulation<2> population(mesh, cells);
    population.SetAbsoluteMovementThreshold(50.0);
    population.AddCellWriter<CellIdWriter>();
    population.AddCellWriter<CellAgesWriter>();
    population.AddCellWriter<CellVolumesWriter>();
    population.AddCellWriter<CellPolarityWriter>();
    population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    OffLatticeSimulation<2> simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
    simulator.SetEndTime(p.endTime);

    if (p.enableDifferentialAdhesion)
    {
        MAKE_PTR(DifferentialAdhesionForce<2>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.springCutoff);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        p_spring->SetApicalApicalAdhesion(p.apicalApicalAdhesion);
        p_spring->SetBasalBasalAdhesion(p.basalBasalAdhesion);
        p_spring->SetApicalBasalAdhesion(p.apicalBasalAdhesion);
        simulator.AddForce(p_spring);
    }
    else
    {
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.springCutoff);
        simulator.AddForce(p_spring);
    }

    MAKE_PTR(BasementMembraneForce<2>, p_bm);
    p_bm->SetBasementMembraneParameter(p.bmStiffnessNode);
    p_bm->SetBasementMembraneRadius(p.bmRadius2d);
    c_vector<double, 2> center2d = zero_vector<double>(2);
    p_bm->SetOrganoidCenter(center2d);
    p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius2d);
    simulator.AddForce(p_bm);

    if (p.enableLumenPressure)
    {
        MAKE_PTR(LumenPressureForce<2>, p_lumen);
        p_lumen->SetPressureStrength(p.lumenPressure);
        p_lumen->SetLumenEquilibriumRadius(p.lumenEqRadius2d);
        p_lumen->SetTrackCenter(true);
        simulator.AddForce(p_lumen);
    }

    if (p.enableApicalConstriction)
    {
        MAKE_PTR(ApicalConstrictionForce<2>, p_ac);
        p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
        simulator.AddForce(p_ac);
    }

    // Ring topology tracker - maintains stable circular neighbors
    // Neighbors are determined by angular position at start, then updated only on division
    boost::shared_ptr<RingTopologyTracker<2>> p_ring_tracker(new RingTopologyTracker<2>());
    p_ring_tracker->SetCenter(center2d);
    simulator.AddSimulationModifier(p_ring_tracker);

    // Curvature bending force (Drasdo 2000) to enforce monolayer structure
    if (p.enableCurvatureBending)
    {
        MAKE_PTR(CurvatureBendingForce<2>, p_curv);
        p_curv->SetBendingStiffness(p.bendingStiffness);
        p_curv->SetTargetRadius(p.organoidRadius2d);
        p_curv->SetLumenExclusionStrength(p.lumenExclusionStrength);
        p_curv->SetMinRadiusFraction(p.minRadiusFraction);
        p_curv->SetNeighborCutoff(p.interactionCutoff2d);
        p_curv->SetTrackCenter(true);
        // Use ring topology for stable neighbors (not position-based)
        p_curv->SetRingTopologyTracker(p_ring_tracker.get());
        simulator.AddForce(p_curv);
    }

    // Cell polarity force for monolayer maintenance (ya||a-style)
    if (p.enableCellPolarity)
    {
        MAKE_PTR(CellPolarityForce<2>, p_polarity);
        p_polarity->SetBendingStrength(p.polarityBendingStrength);
        p_polarity->SetPolarityAlignmentStrength(p.polarityAlignmentStrength);
        p_polarity->SetInteractionCutoff(p.interactionCutoff2d);
        p_polarity->SetInitializeRadially(true);
        simulator.AddForce(p_polarity);
    }

    if (p.enableSloughing)
    {
        AddBoundingBoxKillers<2>(simulator, population,
                                p.organoidRadius2d * p.sloughRadiusFactor);
    }

    MAKE_PTR(VolumeTrackingModifier<2>, p_vol);
    simulator.AddSimulationModifier(p_vol);

    double totalSimTime = p.enableRelaxation ? (p.relaxationTime + p.endTime) : p.endTime;
    boost::shared_ptr<CryptBuddingSummaryModifier<2>> p_summary(
        new CryptBuddingSummaryModifier<2>(p.ecmStiffness, p.samplingMultiple,
                                           totalSimTime));
    simulator.AddSimulationModifier(p_summary);

    if (p.enableRelaxation)
    {
        std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
        for (AbstractCellPopulation<2>::Iterator it = population.Begin();
             it != population.End(); ++it)
        {
            origTypes[*it] = it->GetCellProliferativeType();
            it->SetCellProliferativeType(p_diff);
        }

        simulator.SetEndTime(p.relaxationTime);
        std::cout << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
        simulator.Solve();

        for (AbstractCellPopulation<2>::Iterator it = population.Begin();
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

#endif // RUNNODE2D_HPP_
