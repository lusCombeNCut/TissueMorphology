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

#include "UniformContactInhibitionCellCycleModel.hpp"
#include "UniformContactInhibitionGenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "RingSpringForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "ECMConfinementForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"
#include "ContinuousPvdModifier.hpp"
#include "RingOutlineWriter.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "CellPolarityForce.hpp"
#include "RingSmoothingForce.hpp"
#include "RingTopologyTracker.hpp"

#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPolarityWriter.hpp"
#include "TangentialCentreBasedDivisionRule.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"
#include "SimProfiler.hpp"

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
        // Create cell cycle model - generational or simple based on config
        AbstractCellCycleModel* p_cycle_base;
        if (p.enableGenerationalCascade)
        {
            auto* p_cycle = new UniformContactInhibitionGenerationalCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(1.0);
            p_cycle->SetTotalCycleMin(p.stemCycleMin);
            p_cycle->SetTotalCycleMax(p.stemCycleMax);
            p_cycle->SetTransitCycleRatio(p.taCycleRatio);
            p_cycle->SetMaxTransitGenerations(p.maxTransitGenerations);
            p_cycle_base = p_cycle;
        }
        else
        {
            auto* p_cycle = new UniformContactInhibitionCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(1.0);
            p_cycle->SetTotalCycleMin(p.stemCycleMin);
            p_cycle->SetTotalCycleMax(p.stemCycleMax);
            p_cycle->SetTransitCycleRatio(p.taCycleRatio);
            p_cycle_base = p_cycle;
        }

        CellPtr p_cell(new Cell(p_state, p_cycle_base));
        AssignCellTypeByFraction(p_cell, 0.0, p_stem, p_ta, p_diff,
                                 p.stemFraction, p.transitFraction);

        // For generational model, set initial generation based on cell type
        if (p.enableGenerationalCascade)
        {
            auto* p_gen_cycle = static_cast<UniformContactInhibitionGenerationalCellCycleModel*>(p_cycle_base);
            if (p_cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
                p_gen_cycle->SetGeneration(0);
            else if (p_cell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
                p_gen_cycle->SetGeneration(1);  // TA cells start at generation 1
            else
                p_gen_cycle->SetGeneration(p.maxTransitGenerations + 1);  // Already differentiated
        }

        p_cell->SetBirthTime(-p_gen->ranf() * p.stemCycleMax);
        p_cell->InitialiseCellCycleModel();
        p_cell->GetCellData()->SetItem("volume", 1.0);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);

        // Initialize cell polarity (2D: angle from x-axis pointing radially outward)
        double theta_cell = 2.0 * M_PI * i / p.numCells2dNode;
        p_cell->GetCellData()->SetItem("polarity_theta", theta_cell);
        p_cell->GetCellData()->SetItem("polarity_phi", 0.0);  // unused in 2D

        cells.push_back(p_cell);
    }

    NodeBasedCellPopulation<2> population(mesh, cells);
    population.SetAbsoluteMovementThreshold(50.0);

    // Tangential division rule: daughters placed side-by-side along the ring,
    // not radially (which would stack them and break the monolayer)
    boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2>> p_div_rule(
        new TangentialCentreBasedDivisionRule<2,2>());
    population.SetCentreBasedDivisionRule(p_div_rule);

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

    // Ring topology tracker — must be created before forces that use it
    // Tracks left/right neighbors on the ring; updates only on birth/death
    boost::shared_ptr<RingTopologyTracker<2>> p_ring_tracker(new RingTopologyTracker<2>());
    c_vector<double, 2> center2d = zero_vector<double>(2);
    p_ring_tracker->SetCenter(center2d);
    simulator.AddSimulationModifier(p_ring_tracker);

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
    else if (p.useTopologyBasedSprings)
    {
        // RingSpringForce: same spring law as GeneralisedLinearSpringForce,
        // but only applies forces between ring-adjacent cells (left & right neighbors).
        // Non-adjacent interactions are completely eliminated.
        MAKE_PTR(RingSpringForce<2>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.springCutoff);
        p_spring->SetRingTopologyTracker(p_ring_tracker.get());
        simulator.AddForce(p_spring);
    }
    else
    {
        // Distance-threshold springs: standard GeneralisedLinearSpringForce
        // applies springs between ALL node pairs within the cutoff distance.
        boost::shared_ptr<GeneralisedLinearSpringForce<2, 2>> p_spring(
            new GeneralisedLinearSpringForce<2, 2>());
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.springCutoff);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        simulator.AddForce(p_spring);
    }

    // ECM field: fiber-based confinement replacing simple radial BM force
    // Grid centered on organoid, spanning ecmMaxRadius in each direction
    double ecm_half = p.organoidRadius2d * p.ecmMaxRadiusFraction;
    boost::shared_ptr<DynamicECMField> p_ecm_field(
        new DynamicECMField("radial", p.ecmGridSpacing,
                           -ecm_half, ecm_half, -ecm_half, ecm_half,
                           p.ecmGridType));
    p_ecm_field->SetDegradationRate(p.ecmDegradationRate);
    p_ecm_field->SetDiffusionCoeff(p.ecmDiffusionCoeff);

    // Clear ECM density inside the organoid — no matrix in the lumen or
    // epithelial layer; ECM only exists outside the cell ring.
    // Use a slightly larger radius so cells sit inside the ECM boundary
    // rather than right on top of it.
    double ecm_clear_radius = p.organoidRadius2d + 0.5 * p.ecmGridSpacing;
    p_ecm_field->ClearDensityInsideRadius(center2d, ecm_clear_radius);

    if (p.enableEcmConfinement)
    {
        MAKE_PTR(ECMConfinementForce<2>, p_ecm);
        p_ecm->SetECMField(p_ecm_field);
        p_ecm->SetConfinementStiffness(p.ecmStiffness);
        p_ecm->SetDegradationEnabled(true);
        p_ecm->SetRemodelingEnabled(p.enableEcmGuidance);
        p_ecm->SetTrackCenter(true);
        simulator.AddForce(p_ecm);
    }

    // ECM field VTK writer for ParaView visualization
    boost::shared_ptr<ECMFieldWriter<2>> p_ecm_writer(
        new ECMFieldWriter<2>(p_ecm_field, p.samplingMultiple));
    simulator.AddSimulationModifier(p_ecm_writer);

    // Ring outline writer — outputs .vtp polydata with line-segment
    // connectivity so ParaView can render the organoid boundary
    boost::shared_ptr<RingOutlineWriter<2>> p_outline_writer(
        new RingOutlineWriter<2>(p_ring_tracker.get(), p.samplingMultiple));
    simulator.AddSimulationModifier(p_outline_writer);

    if (p.enableLumenPressure)
    {
        MAKE_PTR(LumenPressureForce<2>, p_lumen);
        p_lumen->SetPressure(p.lumenPressure);
        p_lumen->SetRingTopologyTracker(p_ring_tracker.get());
        p_lumen->SetTrackCenter(true);
        simulator.AddForce(p_lumen);
    }

    if (p.enableApicalConstriction)
    {
        MAKE_PTR(ApicalConstrictionForce<2>, p_ac);
        p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
        simulator.AddForce(p_ac);
    }

    // Ring smoothing force — discrete Laplacian that pulls each cell toward
    // the midpoint of its two ring neighbors. Smooths kinks without pushing
    // cells radially (which was breaking the monolayer).
    if (p.enableCurvatureBending)
    {
        MAKE_PTR(RingSmoothingForce<2>, p_smooth);
        p_smooth->SetSmoothingStiffness(p.bendingStiffness);
        p_smooth->SetCellTypeScaling(p.gammaStemScale, p.gammaTransitScale, p.gammaDiffScale);
        p_smooth->SetRingTopologyTracker(p_ring_tracker.get());
        simulator.AddForce(p_smooth);
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

    // Continuous PVD shadow copies (keeps .pvd files valid for ParaView during simulation)
    if (p.enableContinuousPvd)
    {
        boost::shared_ptr<ContinuousPvdModifier<2>> p_pvd(
            new ContinuousPvdModifier<2>(p.samplingMultiple));
        simulator.AddSimulationModifier(p_pvd);
    }

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

    // Print profiling summary
    SimProfiler::Instance().PrintSummary();

    // Also save CSV for further analysis
    OutputFileHandler prof_handler(outputDir, false);
    std::string profCsvPath = prof_handler.GetOutputDirectoryFullPath() + "profiler.csv";
    SimProfiler::Instance().WriteSummaryCSV(profCsvPath);

    if (final_cells == 0)
        EXCEPTION("Simulation ended with zero cells");
}

#endif // RUNNODE2D_HPP_
