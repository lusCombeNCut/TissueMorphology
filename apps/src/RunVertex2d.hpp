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

#include "ContactInhibitionCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

#include "NagaiHondaForce.hpp"
#include "BasementMembraneForce.hpp"
#include "LumenPressureForce.hpp"
#include "ApicalConstrictionForce.hpp"
#include "RadialVertexBasedDivisionRule.hpp"

#include "VolumeTrackingModifier.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"

#include "CryptBuddingParams.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingUtils.hpp"

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
        ContactInhibitionCellCycleModel* p_cycle = new ContactInhibitionCellCycleModel();
        p_cycle->SetDimension(2);
        p_cycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        p_cycle->SetEquilibriumVolume(target_area);

        CellPtr p_cell(new Cell(p_state, p_cycle));
        AssignCellTypeByFraction(p_cell, 0.0, p_stem, p_ta, p_diff,
                                 p.stemFraction, p.transitFraction);

        p_cell->SetBirthTime(-p_gen->ranf() * 12.0);
        p_cell->InitialiseCellCycleModel();
        p_cell->GetCellData()->SetItem("target area", target_area);
        p_cell->GetCellData()->SetItem("volume", target_area);
        p_cell->GetCellData()->SetItem("basement_membrane_stiffness", p.bmStiffnessVertex);
        p_cell->GetCellData()->SetItem("is_apical", 1.0);
        cells.push_back(p_cell);
    }

    VertexBasedCellPopulation<2> population(*pMesh, cells);

    // Use radial division rule to ensure daughters stay side-by-side
    // within the monolayer (prevents radial stacking after T1 swaps)
    boost::shared_ptr<AbstractVertexBasedDivisionRule<2>> p_div_rule(
        new RadialVertexBasedDivisionRule<2>());
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

    MAKE_PTR(NagaiHondaForce<2>, p_nh);
    p_nh->SetNagaiHondaDeformationEnergyParameter(p.ecmStiffness);
    p_nh->SetNagaiHondaMembraneSurfaceEnergyParameter(p.nhMembraneSurface);
    p_nh->SetNagaiHondaCellCellAdhesionEnergyParameter(p.nhCellCellAdhesion);
    p_nh->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(p.nhBoundaryAdhesion);
    simulator.AddForce(p_nh);

    MAKE_PTR(BasementMembraneForce<2>, p_bm);
    p_bm->SetBasementMembraneParameter(p.bmStiffnessVertex);
    p_bm->SetBasementMembraneRadius(p.bmRadius2d);
    c_vector<double, 2> center2d = zero_vector<double>(2);
    p_bm->SetOrganoidCenter(center2d);
    p_bm->EnableEcmDegradation(p.ecmDegradationRate, p.outerRadius2d * 4.0);
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

#endif // RUNVERTEX2D_HPP_
