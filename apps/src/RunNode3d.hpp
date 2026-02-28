/*
 * RunNode3d.hpp — 3D Node-Based crypt budding model runner
 */
#ifndef RUNNODE3D_HPP_
#define RUNNODE3D_HPP_

#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>

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
#include "DifferentialAdhesionForce.hpp"
#include "BasementMembraneForce.hpp"
#include "ContinuousPvdModifier.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "CellPolarityForce.hpp"

#include "SurfaceTopologyTracker.hpp"
#include "SurfaceSpringForce.hpp"
#include "SurfaceSmoothingForce.hpp"
#include "SurfaceMeshWriter.hpp"

#include "DynamicECMContactGuidanceForce3d.hpp"
#include "DynamicECMField3d.hpp"
#include "ECMConfinementForce3d.hpp"
#include "ECMFieldWriter3d.hpp"

#include "VolumeTrackingModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPolarityWriter.hpp"
#include "TangentialCentreBasedDivisionRule.hpp"

#include "TimedForce.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"
#include "SimProfiler.hpp"
#include "OutputFileHandler.hpp"

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
        // Create cell cycle model - generational or simple based on config
        AbstractCellCycleModel* p_cycle_base;
        if (p.enableGenerationalCascade)
        {
            auto* p_cycle = new UniformContactInhibitionGenerationalCellCycleModel();
            p_cycle->SetDimension(3);
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
            p_cycle->SetDimension(3);
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

    // Tangential division rule: daughters placed side-by-side on the sphere surface,
    // not radially (which would stack them and break the monolayer)
    boost::shared_ptr<AbstractCentreBasedDivisionRule<3,3>> p_div_rule(
        new TangentialCentreBasedDivisionRule<3,3>());
    population.SetCentreBasedDivisionRule(p_div_rule);

    population.AddCellWriter<CellIdWriter>();
    population.AddCellWriter<CellAgesWriter>();
    population.AddCellWriter<CellPolarityWriter>();
    population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    OffLatticeSimulation<3> simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
    simulator.SetEndTime(p.endTime);

    // Surface topology tracker — maintains the cell adjacency graph on the
    // 3D surface.  Must be created before any force that depends on it.
    // Initial connectivity is built from spatial proximity at t = 0;
    // thereafter it updates topologically on birth/death only.
    boost::shared_ptr<SurfaceTopologyTracker<3>> p_surface_tracker(
        new SurfaceTopologyTracker<3>());
    c_vector<double, 3> center3d = zero_vector<double>(3);
    p_surface_tracker->SetCenter(center3d);
    p_surface_tracker->SetInitialCutoff(p.interactionCutoff3d);
    simulator.AddSimulationModifier(p_surface_tracker);

    if (p.useTopologyBasedSprings)
    {
        // SurfaceSpringForce: same spring law as GeneralisedLinearSpringForce,
        // but only between topologically connected cells (from the surface
        // adjacency graph). Non-adjacent spatial neighbors are ignored.
        auto p_spring = boost::make_shared<SurfaceSpringForce<3>>();
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.interactionCutoff3d);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        p_spring->SetSurfaceTopologyTracker(p_surface_tracker.get());
        simulator.AddForce(p_spring);
    }
    else
    {
        // Distance-threshold springs: standard GeneralisedLinearSpringForce
        // applies springs between ALL node pairs within the cutoff distance.
        auto p_spring = boost::make_shared<GeneralisedLinearSpringForce<3, 3>>();
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.interactionCutoff3d);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        auto p_timed_spring = boost::make_shared<TimedForce<3>>(p_spring, "DistanceSpring");
        simulator.AddForce(p_timed_spring);
    }

    {
        auto p_bm = boost::make_shared<BasementMembraneForce<3>>();
        p_bm->SetBasementMembraneParameter(p.bmStiffnessNode);
        p_bm->SetTargetRadius(p.bmRadius3d);
        p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.ecmMaxRadius3d);
        simulator.AddForce(p_bm);
    }

    if (p.enableLumenPressure)
    {
        auto p_lumen = boost::make_shared<LumenPressureForce<3>>();
        p_lumen->SetPressure(p.lumenPressure);
        p_lumen->SetTrackCenter(true);
        simulator.AddForce(p_lumen);
    }

    if (p.enableApicalConstriction)
    {
        auto p_ac = boost::make_shared<ApicalConstrictionForce<3>>();
        p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
        simulator.AddForce(p_ac);
    }

    // Cell polarity force for monolayer maintenance (ya||a-style)
    if (p.enableCellPolarity)
    {
        auto p_polarity = boost::make_shared<CellPolarityForce<3>>();
        p_polarity->SetBendingStrength(p.polarityBendingStrength);
        p_polarity->SetPolarityAlignmentStrength(p.polarityAlignmentStrength);
        p_polarity->SetInteractionCutoff(p.interactionCutoff3d);
        p_polarity->SetInitializeRadially(true);
        simulator.AddForce(p_polarity);
    }

    // Surface smoothing force — discrete Laplacian that pulls each cell
    // toward the centroid of its topological neighbors, regularizing
    // the mesh on the sphere surface.
    if (p.enableCurvatureBending)
    {
        auto p_smooth = boost::make_shared<SurfaceSmoothingForce<3>>();
        p_smooth->SetSmoothingStiffness(p.bendingStiffness);
        p_smooth->SetCellTypeScaling(p.gammaStemScale, p.gammaTransitScale, p.gammaDiffScale);
        p_smooth->SetSurfaceTopologyTracker(p_surface_tracker.get());
        simulator.AddForce(p_smooth);
    }

    // ------------------------------------------------------------------
    // ECM field (shared between confinement and guidance forces)
    // ------------------------------------------------------------------
    boost::shared_ptr<DynamicECMField3d> pEcmField;
    if (p.enableEcmConfinement || p.enableEcmGuidance)
    {
        pEcmField.reset(new DynamicECMField3d(
            "radial", p.ecmGridSpacing,
            -p.ecmDomainHalf, p.ecmDomainHalf,
            -p.ecmDomainHalf, p.ecmDomainHalf,
            -p.ecmDomainHalf, p.ecmDomainHalf));
        pEcmField->SetDegradationRate(p.ecmDegradationRate);
        pEcmField->SetDiffusionCoeff(p.ecmDiffusionCoeff);
        pEcmField->SetRemodelingRate(0.05);
        pEcmField->SetDepositionRate(0.0003);

        // Clear ECM density inside the organoid — no matrix in the lumen;
        // ECM only exists outside the cell shell.
        double ecm_clear_radius = p.organoidRadius3d + 0.5 * p.ecmGridSpacing;
        pEcmField->ClearDensityInsideRadius(center3d, ecm_clear_radius);
    }

    if (p.enableEcmConfinement && pEcmField)
    {
        auto p_ecm_conf = boost::make_shared<ECMConfinementForce3d>();
        p_ecm_conf->SetECMField(pEcmField);
        p_ecm_conf->SetConfinementStiffness(p.ecmStiffness);
        p_ecm_conf->SetDegradationEnabled(true);
        p_ecm_conf->SetRemodelingEnabled(p.enableEcmGuidance);
        p_ecm_conf->SetTrackCenter(true);
        simulator.AddForce(p_ecm_conf);
    }

    if (p.enableEcmGuidance && pEcmField)
    {
        auto p_ecm = boost::make_shared<DynamicECMContactGuidanceForce3d>();
        p_ecm->SetECMField(pEcmField);
        p_ecm->SetBaseSpeed(p.ecmBaseSpeed);
        p_ecm->SetECMSensitivity(1.0);
        p_ecm->SetEnableDegradation(true);
        p_ecm->SetEnableRemodeling(true);
        p_ecm->SetEnableDeposition(false);
        simulator.AddForce(p_ecm);
    }

    // ECM field VTI writer for ParaView visualization
    if (pEcmField)
    {
        boost::shared_ptr<ECMFieldWriter3d> p_ecm_writer(
            new ECMFieldWriter3d(pEcmField, p.samplingMultiple));
        simulator.AddSimulationModifier(p_ecm_writer);
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

    // Surface mesh writer — outputs triangulated surface as VTP files
    // for ParaView visualization (similar to RingOutlineWriter for 2D)
    boost::shared_ptr<SurfaceMeshWriter<3>> p_mesh_writer(
        new SurfaceMeshWriter<3>(p_surface_tracker.get(), p.samplingMultiple));
    simulator.AddSimulationModifier(p_mesh_writer);

    // Continuous PVD shadow copies (keeps .pvd files valid for ParaView during simulation)
    if (p.enableContinuousPvd)
    {
        boost::shared_ptr<ContinuousPvdModifier<3>> p_pvd(
            new ContinuousPvdModifier<3>(p.samplingMultiple));
        simulator.AddSimulationModifier(p_pvd);
    }

    auto wall_t0 = std::chrono::steady_clock::now();
    SimProfiler::Instance().Reset();

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

    auto wall_t1 = std::chrono::steady_clock::now();
    double wall_secs = std::chrono::duration<double>(wall_t1 - wall_t0).count();

    unsigned final_cells = population.GetNumRealCells();
    std::cout << "\nSIMULATION COMPLETE  |  Final cells: " << final_cells
              << "  |  Wall time: " << std::fixed << std::setprecision(1)
              << wall_secs << "s" << std::endl;

    // Print unified profiling summary
    SimProfiler::Instance().PrintSummary();

    // Also save CSV for further analysis
    OutputFileHandler prof_handler(outputDir, false);
    std::string profCsvPath = prof_handler.GetOutputDirectoryFullPath() + "profiler.csv";
    SimProfiler::Instance().WriteSummaryCSV(profCsvPath);

    if (final_cells == 0)
        EXCEPTION("Simulation ended with zero cells");
}

#endif // RUNNODE3D_HPP_
