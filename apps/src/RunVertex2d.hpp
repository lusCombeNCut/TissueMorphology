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

#include "UniformContactInhibitionCellCycleModel.hpp"
#include "UniformContactInhibitionGenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "FastNagaiHondaForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "LocalTangentVertexBasedDivisionRule.hpp"
#include "ContinuousPvdModifier.hpp"

#include "ECMConfinementForce.hpp"
#include "DynamicECMField.hpp"
#include "ECMFieldWriter.hpp"

#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"
#include "SimProfiler.hpp"
#include "TimedForce.hpp"
#include "OutputFileHandler.hpp"

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
    double target_area = 0.5 * dtheta
        * (p.outerRadius2d * p.outerRadius2d - p.innerRadius2d * p.innerRadius2d);

    std::vector<CellPtr> cells;
    MAKE_PTR(WildTypeCellMutationState, p_state);
    MAKE_PTR(StemCellProliferativeType, p_stem);
    MAKE_PTR(TransitCellProliferativeType, p_ta);
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    for (unsigned i = 0; i < pMesh->GetNumElements(); i++)
    {
        // Create cell cycle model - generational or simple based on config
        AbstractCellCycleModel* p_cycle_base;
        if (p.enableGenerationalCascade)
        {
            auto* p_cycle = new UniformContactInhibitionGenerationalCellCycleModel();
            p_cycle->SetDimension(2);
            p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
            p_cycle->SetEquilibriumVolume(target_area);
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
            p_cycle->SetEquilibriumVolume(target_area);
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
    population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();

    OffLatticeSimulation<2> simulator(population);
    simulator.SetOutputDirectory(outputDir);
    simulator.SetDt(p.dt);
    simulator.SetSamplingTimestepMultiple(p.samplingMultiple);
    simulator.SetEndTime(p.endTime);

    auto p_nh = boost::make_shared<FastNagaiHondaForce<2>>();
    // Note: Deformation energy should be ~100 (Nagai-Honda default) for area stability.
    // ecmStiffness is used for ECM field stiffness, not the vertex deformation energy.
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
    boost::shared_ptr<DynamicECMField> p_ecm_field(
        new DynamicECMField("radial", p.ecmGridSpacing,
                           -ecm_half, ecm_half, -ecm_half, ecm_half,
                           p.ecmGridType));
    p_ecm_field->SetDegradationRate(p.ecmDegradationRate);
    p_ecm_field->SetDiffusionCoeff(p.ecmDiffusionCoeff);

    // Clear ECM density inside the organoid
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

    if (p.enableLumenPressure)
    {
        MAKE_PTR(LumenPressureForce<2>, p_lumen);
        p_lumen->SetPressure(p.lumenPressure);
        p_lumen->SetTrackCenter(true);
        
        // Configure target volume (incompressible fluid) mode if enabled
        if (p.lumenUseTargetVolume)
        {
            p_lumen->SetUseTargetVolume(true);
            p_lumen->SetVolumeGrowthRate(p.lumenVolumeGrowthRate);
            p_lumen->SetBulkModulus(p.lumenBulkModulus);
        }
        
        simulator.AddForce(p_lumen);
    }

    if (p.enableApicalConstriction)
    {
        MAKE_PTR(ApicalConstrictionForce<2>, p_ac);
        p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
        simulator.AddForce(p_ac);
    }

    if (p.enableSloughing)
    {
        AddBoundingBoxKillers<2>(simulator, population,
                                p.outerRadius2d * p.sloughRadiusFactor);
    }

    MAKE_PTR(SimpleTargetAreaModifier<2>, p_area);
    p_area->SetReferenceTargetArea(target_area);
    simulator.AddSimulationModifier(p_area);

    MAKE_PTR(VolumeTrackingModifier<2>, p_vol);
    simulator.AddSimulationModifier(p_vol);

    boost::shared_ptr<CryptBuddingSummaryModifier<2>> p_summary(
        new CryptBuddingSummaryModifier<2>(p.ecmStiffness, p.samplingMultiple, p.endTime));
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

#endif // RUNVERTEX2D_HPP_
