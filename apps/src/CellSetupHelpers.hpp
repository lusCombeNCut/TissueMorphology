/*
 * CellSetupHelpers.hpp
 *
 * Shared helpers for cell cycle model creation, cell type assignment,
 * and generational cascade setup. Used by all Run*.hpp runners.
 */
#ifndef CELLSETUPHELPERS_HPP_
#define CELLSETUPHELPERS_HPP_

#include "AbstractCellCycleModel.hpp"
#include "UniformContactInhibitionCellCycleModel.hpp"
#include "UniformContactInhibitionGenerationalCellCycleModel.hpp"
#include "StochasticFourTypeCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"
#include "CryptBuddingParams.hpp"

/**
 * Create the appropriate cell cycle model based on feature flags.
 *
 * @param p           simulation parameters
 * @param dim         spatial dimension (2 or 3)
 * @param eqVolume    equilibrium volume for contact inhibition
 *                    (1.0 for node-based, target_area for 2D vertex, avgVol for 3D vertex)
 * @return raw pointer to the new model (ownership transferred to Cell constructor)
 */
inline AbstractCellCycleModel* CreateCellCycleModel(
    const CryptBuddingParams& p, unsigned dim, double eqVolume)
{
    if (p.enableStochasticFourType)
    {
        auto* pCycle = new StochasticFourTypeCellCycleModel();
        pCycle->SetDimension(dim);
        pCycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        pCycle->SetEquilibriumVolume(eqVolume);
        pCycle->SetTotalCycleMin(p.stemCycleMin);
        pCycle->SetTotalCycleMax(p.stemCycleMax);
        pCycle->SetTransitCycleRatio(p.taCycleRatio);
        pCycle->SetProbStemToStem(p.probStemToStem);
        pCycle->SetProbStemToPaneth(p.probStemToPaneth);
        pCycle->SetProbTaToTaEarly(p.probTaToTaEarly);
        pCycle->SetProbTaToTaLate(p.probTaToTaLate);
        pCycle->SetTransitionTime(p.transitionTime);
        return pCycle;
    }
    else if (p.enableGenerationalCascade)
    {
        auto* pCycle = new UniformContactInhibitionGenerationalCellCycleModel();
        pCycle->SetDimension(dim);
        pCycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        pCycle->SetEquilibriumVolume(eqVolume);
        pCycle->SetTotalCycleMin(p.stemCycleMin);
        pCycle->SetTotalCycleMax(p.stemCycleMax);
        pCycle->SetTransitCycleRatio(p.taCycleRatio);
        pCycle->SetMaxTransitGenerations(p.maxTransitGenerations);
        return pCycle;
    }
    else
    {
        auto* pCycle = new UniformContactInhibitionCellCycleModel();
        pCycle->SetDimension(dim);
        pCycle->SetQuiescentVolumeFraction(p.quiescentFraction);
        pCycle->SetEquilibriumVolume(eqVolume);
        pCycle->SetTotalCycleMin(p.stemCycleMin);
        pCycle->SetTotalCycleMax(p.stemCycleMax);
        pCycle->SetTransitCycleRatio(p.taCycleRatio);
        return pCycle;
    }
}

/**
 * Assign cell types using the stochastic 4-type model (stem/TA/Paneth/enterocyte).
 * Sets both mutation state and proliferative type + cell_type_id data.
 *
 * Only valid when p.enableStochasticFourType is true.
 */
inline void AssignStochasticFourType(
    CellPtr pCell, const CryptBuddingParams& p,
    boost::shared_ptr<AbstractCellMutationState> pWildType,
    boost::shared_ptr<AbstractCellMutationState> pTaMut,
    boost::shared_ptr<AbstractCellMutationState> pPanethMut,
    boost::shared_ptr<AbstractCellMutationState> pEcMut,
    boost::shared_ptr<AbstractCellProliferativeType> pTransit,
    boost::shared_ptr<AbstractCellProliferativeType> pDiff)
{
    double u = RandomNumberGenerator::Instance()->ranf();
    if (u < p.stemFraction)
    {
        pCell->SetCellProliferativeType(pTransit);
        pCell->SetMutationState(pWildType);         // WildType = Stem
        pCell->GetCellData()->SetItem("cell_type_id", 0.0);
    }
    else if (u < p.stemFraction + p.transitFraction)
    {
        pCell->SetCellProliferativeType(pTransit);
        pCell->SetMutationState(pTaMut);
        pCell->GetCellData()->SetItem("cell_type_id", 1.0);
    }
    else if (u < p.stemFraction + p.transitFraction + p.panethFraction)
    {
        pCell->SetCellProliferativeType(pDiff);
        pCell->SetMutationState(pPanethMut);
        pCell->GetCellData()->SetItem("cell_type_id", 2.0);
    }
    else
    {
        pCell->SetCellProliferativeType(pDiff);
        pCell->SetMutationState(pEcMut);
        pCell->GetCellData()->SetItem("cell_type_id", 3.0);
    }
}

/**
 * Set the initial generation on a generational cell cycle model
 * based on the cell's current proliferative type.
 *
 * Call after cell type assignment. No-op if generational cascade is disabled
 * or stochastic 4-type is enabled (they are mutually exclusive).
 */
inline void SetInitialGeneration(
    CellPtr pCell, AbstractCellCycleModel* pCycleBase, const CryptBuddingParams& p)
{
    if (!p.enableGenerationalCascade || p.enableStochasticFourType) return;

    auto* pGenCycle = static_cast<UniformContactInhibitionGenerationalCellCycleModel*>(pCycleBase);
    if (pCell->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
        pGenCycle->SetGeneration(0);
    else if (pCell->GetCellProliferativeType()->template IsType<TransitCellProliferativeType>())
        pGenCycle->SetGeneration(1);
    else
        pGenCycle->SetGeneration(p.maxTransitGenerations + 1);
}

#endif // CELLSETUPHELPERS_HPP_
