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

#include "WildTypeCellMutationState.hpp"
#include "TACellMutationState.hpp"
#include "PanethCellMutationState.hpp"
#include "EnterocyteCellMutationState.hpp"

#include "GeneralisedLinearSpringForce.hpp"
#include "RingSpringForce.hpp"
#include "DifferentialAdhesionForce.hpp"
#include "ECMConfinementForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"
#include "GhostNodeEcmField.hpp"
#include "GhostNodeEcmForce.hpp"
#include "GhostNodeEcmWriter.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "ViscoelasticGhostNodeEcmForce.hpp"
#include "ViscoelasticGhostNodeEcmWriter.hpp"
#include "RingOutlineWriter.hpp"
#include "RingSmoothingForce.hpp"
#include "RingTopologyTracker.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellPolarityWriter.hpp"
#include "CellContactInhibitionStatusWriter.hpp"
#include "CellTrueMutationStateWriter.hpp"
#include "TangentialCentreBasedDivisionRule.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"
#include "SimProfiler.hpp"

#include "CellSetupHelpers.hpp"
#include "ECMWiringHelpers.hpp"
#include "SimulationHelpers.hpp"

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
    MAKE_PTR(TACellMutationState, p_ta_mut);
    MAKE_PTR(PanethCellMutationState, p_paneth_mut);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    for (unsigned i = 0; i < mesh.GetNumNodes(); i++)
    {
        AbstractCellCycleModel* p_cycle_base = CreateCellCycleModel(p, 2, 1.0);

        CellPtr p_cell(new Cell(p_state, p_cycle_base));

        if (p.enableStochasticFourType)
        {
            AssignStochasticFourType(p_cell, p, p_state, p_ta_mut, p_paneth_mut, p_ta, p_diff);
        }
        else
        {
            AssignCellTypeByFraction(p_cell, 0.0, p_stem, p_ta, p_diff,
                                     p.stemFraction, p.transitFraction);
        }

        SetInitialGeneration(p_cell, p_cycle_base, p);

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
    population.AddCellWriter<CellContactInhibitionStatusWriter>();
    population.AddCellWriter<CellTrueMutationStateWriter>();
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

    // Wire the division rule to the tracker so it can record exact parent indices
    // at division time, preventing spurious long-range spring connections.
    static_cast<TangentialCentreBasedDivisionRule<2,2>*>(p_div_rule.get())
        ->SetRingTopologyTracker(p_ring_tracker.get());

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
        if (p.useTopologyBasedSprings)
        {
            p_spring->SetRingTopologyTracker(p_ring_tracker.get());
        }
        p_spring->SetTAStiffnessScale(p.springStiffnessTAScale);
        p_spring->SetDiffStiffnessScale(p.springStiffnessDiffScale);
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
        p_spring->SetTAStiffnessScale(p.springStiffnessTAScale);
        p_spring->SetDiffStiffnessScale(p.springStiffnessDiffScale);
        simulator.AddForce(p_spring);
    }
    else
    {
        // Distance-threshold springs: standard spring law applied
        // between ALL node pairs within the cutoff distance.
        // Uses RingSpringForce without a tracker (falls back to
        // GeneralisedLinearSpringForce behaviour) so that per-cell-type
        // stiffness scaling is still available.
        MAKE_PTR(RingSpringForce<2>, p_spring);
        p_spring->SetMeinekeSpringStiffness(p.springStiffness);
        p_spring->SetCutOffLength(p.springCutoff);
        p_spring->SetMeinekeDivisionRestingSpringLength(0.5);
        p_spring->SetMeinekeSpringGrowthDuration(1.0);
        p_spring->SetTAStiffnessScale(p.springStiffnessTAScale);
        p_spring->SetDiffStiffnessScale(p.springStiffnessDiffScale);
        simulator.AddForce(p_spring);
    }

    // ECM field: fiber-based confinement replacing simple radial BM force
    // Grid centered on organoid, spanning ecmMaxRadius in each direction
    double ecm_half = p.organoidRadius2d * p.ecmMaxRadiusFraction;

    if (p.enableGhostNodeECM)
    {
        // ── Ghost node ECM (off-lattice particle-based) ──────────
        double gn_spacing = p.ecmGridSpacing;
        double gn_rest = (p.ghostSpringRestLength > 0.0) ? p.ghostSpringRestLength : gn_spacing;

        if (p.enableViscoelasticECM)
        {
            // ── Viscoelastic constitutive model ──────────────────
            boost::shared_ptr<ViscoelasticGhostNodeEcmField<2>> p_ve_field(
                new ViscoelasticGhostNodeEcmField<2>("random", gn_spacing, -ecm_half, ecm_half, -ecm_half, ecm_half,
                                      p.ecmGridType));
            ConfigureViscoelasticField(p_ve_field, p);

            double ecm_clear_radius = p.organoidRadius2d + 0.5 * gn_spacing;
            p_ve_field->ClearDensityInsideRadius(center2d, ecm_clear_radius);

            if (p.enableEcmConfinement)
            {
                MAKE_PTR(ViscoelasticGhostNodeEcmForce<2>, p_ve_force);
                ConfigureViscoelasticForce(p_ve_force, p_ve_field, p);
                simulator.AddForce(p_ve_force);
            }

            boost::shared_ptr<ViscoelasticGhostNodeEcmWriter<2>> p_ve_writer(
                new ViscoelasticGhostNodeEcmWriter<2>(p_ve_field, p.samplingMultiple));
            simulator.AddSimulationModifier(p_ve_writer);
        }
        else
        {
            // ── Original elastic ghost node ECM ──────────────────
            boost::shared_ptr<GhostNodeEcmField<2>> p_ghost_field(
                new GhostNodeEcmField<2>("random", gn_spacing, -ecm_half, ecm_half, -ecm_half, ecm_half,
                                      p.ecmGridType));
            ConfigureGhostField(p_ghost_field, p, gn_rest);

            double ecm_clear_radius = p.organoidRadius2d + 0.5 * gn_spacing;
            p_ghost_field->ClearDensityInsideRadius(center2d, ecm_clear_radius);

            if (p.enableEcmConfinement)
            {
                MAKE_PTR(GhostNodeEcmForce<2>, p_gn_force);
                ConfigureGhostForce(p_gn_force, p_ghost_field, p);
                simulator.AddForce(p_gn_force);
            }

            boost::shared_ptr<GhostNodeEcmWriter<2>> p_gn_writer(
                new GhostNodeEcmWriter<2>(p_ghost_field, p.samplingMultiple));
            simulator.AddSimulationModifier(p_gn_writer);
        }
    }
    else
    {
        // ── Grid-based ECM (original) ────────────────────────────
        boost::shared_ptr<DynamicECMField> p_ecm_field(
            new DynamicECMField("random", p.ecmGridSpacing,
                               -ecm_half, ecm_half, -ecm_half, ecm_half,
                               p.ecmGridType));
        p_ecm_field->SetDegradationRate(p.ecmDegradationRate);
        p_ecm_field->SetDiffusionCoeff(p.ecmDiffusionCoeff);

        double ecm_clear_radius = p.organoidRadius2d + 0.5 * p.ecmGridSpacing;
        p_ecm_field->ClearDensityInsideRadius(center2d, ecm_clear_radius);

        if (p.enableEcmConfinement)
        {
            MAKE_PTR(ECMConfinementForce<2>, p_ecm);
            ConfigureGridECMForce(p_ecm, p_ecm_field, p);
            simulator.AddForce(p_ecm);
        }

        boost::shared_ptr<ECMFieldWriter<2>> p_ecm_writer(
            new ECMFieldWriter<2>(p_ecm_field, p.samplingMultiple));
        simulator.AddSimulationModifier(p_ecm_writer);
    }

    // Ring outline writer — outputs .vtp polydata with line-segment
    // connectivity so ParaView can render the organoid boundary
    boost::shared_ptr<RingOutlineWriter<2>> p_outline_writer(
        new RingOutlineWriter<2>(p_ring_tracker.get(), p.samplingMultiple));
    simulator.AddSimulationModifier(p_outline_writer);

    WireLumenPressure<2>(simulator, population, p);

    WireApicalConstriction<2>(simulator, p);

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

    WireCellPolarity<2>(simulator, p, p.interactionCutoff2d);

    if (p.enableSloughing)
    {
        AddBoundingBoxKillers<2>(simulator, population,
                                p.organoidRadius2d * p.sloughRadiusFactor);
    }

    double totalSimTime = p.enableRelaxation ? (p.relaxationTime + p.endTime) : p.endTime;
    WireCommonModifiers<2>(simulator, p, totalSimTime);

    if (p.enableRelaxation)
    {
        RunRelaxationPhase<2>(simulator, population, p);
    }
    else
    {
        simulator.Solve();
    }

    PrintSimulationSummary<2>(population, outputDir);
}

#endif // RUNNODE2D_HPP_
