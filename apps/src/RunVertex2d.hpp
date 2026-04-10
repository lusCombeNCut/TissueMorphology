/*
 * RunVertex2d.hpp — 2D Vertex-Based crypt budding model runner
 */
#ifndef RUNVERTEX2D_HPP_
#define RUNVERTEX2D_HPP_

#include <string>
#include <map>
#include <iostream>

#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"
#include "ChastePoint.hpp"

#include "MutableVertexMesh.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"

#include "WildTypeCellMutationState.hpp"
#include "TACellMutationState.hpp"
#include "PanethCellMutationState.hpp"
#include "EnterocyteCellMutationState.hpp"

#include "FastNagaiHondaForce.hpp"
#include "LocalTangentVertexBasedDivisionRule.hpp"

#include "ECMConfinementForce.hpp"
#include "NonNeighbourCellRepulsionForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"
#include "GhostNodeEcmField.hpp"
#include "GhostNodeEcmForce.hpp"
#include "GhostNodeEcmWriter.hpp"
#include "ViscoelasticGhostNodeEcmField.hpp"
#include "ViscoelasticGhostNodeEcmForce.hpp"
#include "ViscoelasticGhostNodeEcmWriter.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellContactInhibitionStatusWriter.hpp"
#include "CellTrueMutationStateWriter.hpp"
#include "CryptBuddingParams.hpp"
#include "CryptBuddingUtils.hpp"
#include "TimedForce.hpp"

#include "CellSetupHelpers.hpp"
#include "ECMWiringHelpers.hpp"
#include "SimulationHelpers.hpp"

void RunVertex2d(const CryptBuddingParams& p, const std::string& outputDir)
{
    RandomNumberGenerator::Instance()->Reseed(p.randomSeed);

    boost::shared_ptr<MutableVertexMesh<2,2>> pMesh =
        MakeAnnularVertexMesh(p.numCells2dVertex, p.innerRadius2d, p.outerRadius2d,
                              p.t1Threshold2d, p.t2Threshold2d);

    for (unsigned i = 0; i < pMesh->GetNumNodes(); i++)
    {
        c_vector<double, 2> pos = pMesh->GetNode(i)->rGetLocation();
        double r = norm_2(pos);
        if (r > 1e-6)
        {
            double noise = (RandomNumberGenerator::Instance()->ranf() - 0.5) * 0.2;
            pos += noise * (pos / r);
            ChastePoint<2> pt(pos[0], pos[1]);
            pMesh->GetNode(i)->SetPoint(pt);
        }
    }

    double dtheta = 2.0 * M_PI / static_cast<double>(p.numCells2dVertex);
    double target_area = 0.5 * std::sin(dtheta)
        * (p.outerRadius2d * p.outerRadius2d - p.innerRadius2d * p.innerRadius2d);

    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(TACellMutationState, p_ta_mut);
    MAKE_PTR(PanethCellMutationState, p_paneth_mut);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    for (unsigned i = 0; i < pMesh->GetNumElements(); i++)
    {
        AbstractCellCycleModel* p_cycle_base = CreateCellCycleModel(p, 2, target_area);

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
        p_cell->GetCellData()->SetItem("target area", target_area);
        p_cell->GetCellData()->SetItem("volume", target_area);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);
        cells.push_back(p_cell);
    }

    VertexBasedCellPopulation<2> population(*pMesh, cells);

    // Use local tangent division rule to ensure daughters stay side-by-side
    // within the monolayer. This computes local surface direction from neighbors
    // rather than global centroid, so it works correctly for budding shapes.
    boost::shared_ptr<AbstractVertexBasedDivisionRule<2>> p_div_rule(
        new LocalTangentVertexBasedDivisionRule<2>());
    population.SetVertexBasedDivisionRule(p_div_rule);

    population.AddCellWriter<CellIdWriter>();
    population.AddCellWriter<CellAgesWriter>();
    population.AddCellWriter<CellVolumesWriter>();
    population.AddCellWriter<CellContactInhibitionStatusWriter>();
    population.AddCellWriter<CellTrueMutationStateWriter>();
    population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    OffLatticeSimulation<2> simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
    simulator.SetEndTime(p.endTime);

    auto p_nh = boost::make_shared<FastNagaiHondaForce<2>>();
    // Note: Deformation energy should be ~100 (Nagai-Honda default) for area stability.
    // cellGhostSpringStiffness is used for ECM field stiffness, not the vertex deformation energy.
    p_nh->SetDeformationEnergyParameter(100.0);
    p_nh->SetMembraneSurfaceEnergyParameter(p.nhMembraneSurface);
    if (p.enableDifferentialAdhesion)
    {
        // Per-cell-type adhesion matrix
        p_nh->SetAdhesionMatrix(
            p.nhStemStemAdhesion,    p.nhStemTransitAdhesion, p.nhStemDiffAdhesion,
            p.nhTransitTransitAdhesion, p.nhTransitDiffAdhesion,
            p.nhDiffDiffAdhesion);
        p_nh->SetBoundaryAdhesion(0, p.nhStemBoundaryAdhesion);
        p_nh->SetBoundaryAdhesion(1, p.nhTransitBoundaryAdhesion);
        p_nh->SetBoundaryAdhesion(2, p.nhDiffBoundaryAdhesion);
    }
    else
    {
        p_nh->SetUniformCellCellAdhesion(p.nhCellCellAdhesion);
        p_nh->SetUniformBoundaryAdhesion(p.nhBoundaryAdhesion);
    }
    auto p_timed_nh = boost::make_shared<TimedForce<2>>(p_nh, "FastNagaiHonda");
    simulator.AddForce(p_timed_nh);

    c_vector<double, 2> center2d = zero_vector<double>(2);

    // ---- ECM field: fiber-based confinement (replaces BasementMembraneForce) ----
    double ecm_half = p.organoidRadius2d * p.ecmMaxRadiusFraction;

    if (p.enableGhostNodeECM)
    {
        // ── Ghost node ECM (off-lattice particle-based) ──────────
        double gn_spacing = p.ghostGridSpacing;
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

    if (p.enableNonNeighbourCellRepulsion)
    {
        MAKE_PTR(NonNeighbourCellRepulsionForce<2>, p_repulsion);
        p_repulsion->SetRepulsionStiffness(p.nonNeighbourRepulsionStiffness);
        p_repulsion->SetInteractionCutoff(p.nonNeighbourRepulsionCutoff);
        simulator.AddForce(p_repulsion);
    }

    WireLumenPressure<2>(simulator, population, p);

    WireApicalConstriction<2>(simulator, p);

    if (p.enableSloughing)
    {
        AddBoundingBoxKillers<2>(simulator, population,
                                p.outerRadius2d * p.sloughRadiusFactor);
    }

    MAKE_PTR(SimpleTargetAreaModifier<2>, p_area);
    p_area->SetReferenceTargetArea(target_area);
    simulator.AddSimulationModifier(p_area);

    WireCommonModifiers<2>(simulator, p, p.endTime);

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

#endif // RUNVERTEX2D_HPP_
