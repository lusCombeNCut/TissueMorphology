/*
 * SimulationHelpers.hpp
 *
 * Shared helpers for simulation phases (relaxation), common modifiers,
 * and post-simulation cleanup. Used by all Run*.hpp runners.
 */
#ifndef SIMULATIONHELPERS_HPP_
#define SIMULATIONHELPERS_HPP_

#include <map>
#include <iostream>
#include <string>

#include "AbstractCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "VolumeTrackingModifier.hpp"
#include "ContinuousPvdModifier.hpp"
#include "CryptBuddingSummaryModifier.hpp"
#include "CryptBuddingParams.hpp"
#include "SimProfiler.hpp"
#include "OutputFileHandler.hpp"

// ===================================================================
// Standard two-phase relaxation (Node2d, Vertex2d, Node3d)
//
// Phase 1: All cells set to differentiated (no proliferation), solve
// Phase 2: Restore original types, solve growth phase
//
// NOTE: Not used by Vertex3d which has a custom two-phase with
// delayed force addition.
// ===================================================================

template<unsigned DIM>
void RunRelaxationPhase(OffLatticeSimulation<DIM>& rSimulator,
                         AbstractCellPopulation<DIM>& rPopulation,
                         const CryptBuddingParams& p)
{
    MAKE_PTR(DifferentiatedCellProliferativeType, p_diff);

    std::map<CellPtr, boost::shared_ptr<AbstractCellProperty>> origTypes;
    for (typename AbstractCellPopulation<DIM>::Iterator it = rPopulation.Begin();
         it != rPopulation.End(); ++it)
    {
        origTypes[*it] = it->GetCellProliferativeType();
        it->SetCellProliferativeType(p_diff);
    }

    rSimulator.SetEndTime(p.relaxationTime);
    std::cout << "--- Phase 1: Relaxation (" << p.relaxationTime << "h) ---" << std::endl;
    rSimulator.Solve();

    for (typename AbstractCellPopulation<DIM>::Iterator it = rPopulation.Begin();
         it != rPopulation.End(); ++it)
    {
        if (origTypes.count(*it)) it->SetCellProliferativeType(origTypes[*it]);
    }

    rSimulator.SetEndTime(p.relaxationTime + p.endTime);
    std::cout << "--- Phase 2: Growth (" << p.endTime << "h) ---" << std::endl;
    rSimulator.Solve();
}

// ===================================================================
// Wire common modifiers (volume tracking, summary, continuous PVD)
// ===================================================================

template<unsigned DIM>
void WireCommonModifiers(OffLatticeSimulation<DIM>& rSimulator,
                          const CryptBuddingParams& p,
                          double totalSimTime)
{
    MAKE_PTR(VolumeTrackingModifier<DIM>, p_vol);
    rSimulator.AddSimulationModifier(p_vol);

    boost::shared_ptr<CryptBuddingSummaryModifier<DIM>> p_summary(
        new CryptBuddingSummaryModifier<DIM>(p.GetPrimaryEcmParam(), p.samplingMultiple,
                                              totalSimTime));
    rSimulator.AddSimulationModifier(p_summary);

    if (p.enableContinuousPvd)
    {
        boost::shared_ptr<ContinuousPvdModifier<DIM>> p_pvd(
            new ContinuousPvdModifier<DIM>(p.samplingMultiple));
        rSimulator.AddSimulationModifier(p_pvd);
    }
}

// ===================================================================
// Post-simulation output: profiler summary + CSV + zero-cell check
// ===================================================================

template<unsigned DIM>
void PrintSimulationSummary(AbstractCellPopulation<DIM>& rPopulation,
                             const std::string& outputDir)
{
    unsigned finalCells = rPopulation.GetNumRealCells();
    std::cout << "\nSIMULATION COMPLETE  |  Final cells: " << finalCells << std::endl;

    SimProfiler::Instance().PrintSummary();

    OutputFileHandler profHandler(outputDir, false);
    std::string profCsvPath = profHandler.GetOutputDirectoryFullPath() + "profiler.csv";
    SimProfiler::Instance().WriteSummaryCSV(profCsvPath);

    if (finalCells == 0)
        EXCEPTION("Simulation ended with zero cells");
}

// ===================================================================
// Wire lumen pressure force (shared between Node2d, Vertex2d, Node3d)
// ===================================================================

#include "LumenPressureForce.hpp"
#include "CellLumenForceWriter.hpp"

template<unsigned DIM>
void WireLumenPressure(OffLatticeSimulation<DIM>& rSimulator,
                        AbstractCellPopulation<DIM>& rPopulation,
                        const CryptBuddingParams& p)
{
    if (!p.enableLumenPressure) return;

    auto p_lumen = boost::make_shared<LumenPressureForce<DIM>>();
    p_lumen->SetPressure(p.lumenPressure);
    p_lumen->SetTrackCenter(true);
    p_lumen->SetWriteForce(true);
    rSimulator.AddForce(p_lumen);
    rPopulation.template AddCellWriter<CellLumenForceWriter>();
}

// ===================================================================
// Wire apical constriction force
// ===================================================================

#include "ApicalConstrictionForce.hpp"

template<unsigned DIM>
void WireApicalConstriction(OffLatticeSimulation<DIM>& rSimulator,
                             const CryptBuddingParams& p)
{
    if (!p.enableApicalConstriction) return;

    auto p_ac = boost::make_shared<ApicalConstrictionForce<DIM>>();
    p_ac->SetConstrictionStrength(p.apicalConstrictionStrength);
    rSimulator.AddForce(p_ac);
}

// ===================================================================
// Wire cell polarity force
// ===================================================================

#include "CellPolarityForce.hpp"

template<unsigned DIM>
void WireCellPolarity(OffLatticeSimulation<DIM>& rSimulator,
                       const CryptBuddingParams& p, double interactionCutoff)
{
    if (!p.enableCellPolarity) return;

    auto p_polarity = boost::make_shared<CellPolarityForce<DIM>>();
    p_polarity->SetBendingStrength(p.polarityBendingStrength);
    p_polarity->SetPolarityAlignmentStrength(p.polarityAlignmentStrength);
    p_polarity->SetInteractionCutoff(interactionCutoff);
    p_polarity->SetInitializeRadially(true);
    rSimulator.AddForce(p_polarity);
}

#endif // SIMULATIONHELPERS_HPP_
