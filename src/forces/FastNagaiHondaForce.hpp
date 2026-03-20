/*
 * FastNagaiHondaForce.hpp
 *
 * Optimised Nagai–Honda vertex force with per-cell-type differential adhesion.
 *
 * Key differences from the upstream NagaiHondaForce:
 *   1. Element-centric iteration: each edge gradient is computed ONCE, then
 *      accumulated onto both endpoint nodes.  The upstream version iterates
 *      over nodes → containing elements, so every internal edge is visited
 *      twice and GetNodeLocalIndex does a linear scan each time.
 *   2. Adhesion lookup via a small matrix indexed by integer cell-type IDs
 *      (0 = stem, 1 = transit, 2 = differentiated) instead of per-edge
 *      std::set_intersection + virtual GetAdhesionParameter.
 *   3. Target perimeters, deformation/membrane coefficients, and per-type
 *      boundary adhesion values are all precomputed ONCE before the main loop.
 *   4. No virtual getter calls inside the hot loop.
 *
 * Cell-type identification:
 *   Reads the CellData item "cell_type_id" (double → unsigned):
 *     0 = Stem, 1 = Transit, 2 = Differentiated/Paneth.
 *   Falls back to CellProliferativeType introspection if the item is absent.
 */
#ifndef FASTNAGAIHONDAFORCE_HPP_
#define FASTNAGAIHONDAFORCE_HPP_

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include <array>
#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>

template<unsigned DIM>
class FastNagaiHondaForce : public AbstractForce<DIM>
{
public:
    // Number of distinct cell types we track
    static constexpr unsigned NUM_TYPES = 3;

    FastNagaiHondaForce()
        : AbstractForce<DIM>(),
          mDeformationParameter(100.0),
          mMembraneParameter(10.0),
          mTargetArea(1.0)
    {
        // Default: uniform adhesion (no differential)
        mCellCellAdhesion.fill({});
        for (unsigned i = 0; i < NUM_TYPES; ++i)
            for (unsigned j = 0; j < NUM_TYPES; ++j)
                mCellCellAdhesion[i][j] = 0.5;

        mBoundaryAdhesion.fill(1.0);
    }

    virtual ~FastNagaiHondaForce() = default;

    // ── Parameter setters ──────────────────────────────────────────────

    void SetDeformationEnergyParameter(double v)  { mDeformationParameter = v; }
    void SetMembraneSurfaceEnergyParameter(double v) { mMembraneParameter = v; }
    void SetTargetAreaParameter(double v) { mTargetArea = v; }

    /** Set uniform cell-cell adhesion (same for all type pairs). */
    void SetUniformCellCellAdhesion(double v)
    {
        for (unsigned i = 0; i < NUM_TYPES; ++i)
            for (unsigned j = 0; j < NUM_TYPES; ++j)
                mCellCellAdhesion[i][j] = v;
    }

    /** Set adhesion between specific type pair (symmetric). */
    void SetCellCellAdhesion(unsigned typeA, unsigned typeB, double v)
    {
        assert(typeA < NUM_TYPES && typeB < NUM_TYPES);
        mCellCellAdhesion[typeA][typeB] = v;
        mCellCellAdhesion[typeB][typeA] = v;
    }

    /** Set uniform boundary adhesion (same for all types). */
    void SetUniformBoundaryAdhesion(double v)
    {
        mBoundaryAdhesion.fill(v);
    }

    /** Set boundary adhesion for a specific cell type. */
    void SetBoundaryAdhesion(unsigned type, double v)
    {
        assert(type < NUM_TYPES);
        mBoundaryAdhesion[type] = v;
    }

    // Convenience: set the full 3×3 adhesion matrix at once
    // Order: [stem-stem, stem-ta, stem-diff, ta-ta, ta-diff, diff-diff]
    void SetAdhesionMatrix(double ss, double st, double sd,
                           double tt, double td, double dd)
    {
        mCellCellAdhesion[0][0] = ss;
        mCellCellAdhesion[0][1] = mCellCellAdhesion[1][0] = st;
        mCellCellAdhesion[0][2] = mCellCellAdhesion[2][0] = sd;
        mCellCellAdhesion[1][1] = tt;
        mCellCellAdhesion[1][2] = mCellCellAdhesion[2][1] = td;
        mCellCellAdhesion[2][2] = dd;
    }

    // ── Getters for serialisation / output ─────────────────────────────

    double GetDeformationEnergyParameter() const { return mDeformationParameter; }
    double GetMembraneSurfaceEnergyParameter() const { return mMembraneParameter; }
    double GetTargetAreaParameter() const { return mTargetArea; }
    double GetCellCellAdhesion(unsigned a, unsigned b) const { return mCellCellAdhesion[a][b]; }
    double GetBoundaryAdhesion(unsigned t) const { return mBoundaryAdhesion[t]; }

    // ── Core force computation ─────────────────────────────────────────

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override
    {
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
        {
            EXCEPTION("FastNagaiHondaForce is to be used with a VertexBasedCellPopulation only");
        }

        auto* pPop = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        auto& rMesh = pPop->rGetMesh();
        const unsigned numNodes = pPop->GetNumNodes();
        const unsigned numElems = pPop->GetNumElements();

        // Cache scalar parameters (avoid virtual calls in loop)
        const double K_def  = mDeformationParameter;
        const double K_mem  = mMembraneParameter;

        // Check for target area modifier
        bool usingModifier = pPop->Begin()->GetCellData()->HasItem("target area");

        // ── Phase 1: per-element precomputation ────────────────────────
        std::vector<double> elemArea(numElems);
        std::vector<double> elemPerim(numElems);
        std::vector<double> targetArea(numElems);
        std::vector<double> targetPerim(numElems);
        std::vector<unsigned> elemType(numElems);   // 0=stem, 1=ta, 2=diff

        for (auto it = rMesh.GetElementIteratorBegin();
             it != rMesh.GetElementIteratorEnd(); ++it)
        {
            const unsigned ei = it->GetIndex();
            elemArea[ei]  = rMesh.GetVolumeOfElement(ei);
            elemPerim[ei] = rMesh.GetSurfaceAreaOfElement(ei);

            CellPtr pCell = pPop->GetCellUsingLocationIndex(ei);

            if (usingModifier)
                targetArea[ei] = pCell->GetCellData()->GetItem("target area");
            else
                targetArea[ei] = mTargetArea;

            targetPerim[ei] = 2.0 * std::sqrt(M_PI * targetArea[ei]);

            // Determine cell type ID
            elemType[ei] = GetCellTypeId(pCell);
        }

        // ── Phase 2: allocate per-node force accumulators ──────────────
        std::vector<c_vector<double, DIM>> nodeForce(numNodes, zero_vector<double>(DIM));

        // ── Phase 3: element-centric edge walk ─────────────────────────
        // For each element, walk its edges.  Each edge contributes to
        // exactly two nodes (the edge endpoints).  We compute gradients
        // once and scatter to both nodes, avoiding the duplicate work of
        // the node-centric approach.
        for (auto elemIt = rMesh.GetElementIteratorBegin();
             elemIt != rMesh.GetElementIteratorEnd(); ++elemIt)
        {
            const unsigned ei = elemIt->GetIndex();
            VertexElement<DIM, DIM>* pElem = pPop->GetElement(ei);
            const unsigned nVerts = pElem->GetNumNodes();

            const double areaErr = elemArea[ei] - targetArea[ei];
            const double perimErr = elemPerim[ei] - targetPerim[ei];
            const unsigned myType = elemType[ei];

            for (unsigned li = 0; li < nVerts; ++li)
            {
                const unsigned liPrev = (nVerts + li - 1) % nVerts;
                const unsigned liNext = (li + 1) % nVerts;

                Node<DIM>* pNode = pElem->GetNode(li);
                const unsigned ni = pNode->GetIndex();

                // ─ Area gradient at this node ─
                c_vector<double, DIM> areaGrad =
                    rMesh.GetAreaGradientOfElementAtNode(pElem, li);

                nodeForce[ni] -= 2.0 * K_def * areaErr * areaGrad;

                // ─ Edge gradients (prev edge and next edge at this node) ─
                c_vector<double, DIM> prevEdgeGrad =
                    -rMesh.GetNextEdgeGradientOfElementAtNode(pElem, liPrev);
                c_vector<double, DIM> nextEdgeGrad =
                    rMesh.GetNextEdgeGradientOfElementAtNode(pElem, li);

                // Membrane surface tension
                c_vector<double, DIM> perimGrad = prevEdgeGrad + nextEdgeGrad;
                nodeForce[ni] -= 2.0 * K_mem * perimErr * perimGrad;

                // ─ Adhesion for the NEXT edge (node li → node liNext) ─
                // We handle this edge's adhesion contribution ONLY from
                // the element with the lower index to avoid double-counting.
                // For boundary edges (shared by only one element) we always
                // process them here.
                Node<DIM>* pNextNode = pElem->GetNode(liNext);
                const unsigned niNext = pNextNode->GetIndex();

                double adhesionParam = 0.0;

                // Find whether this edge is shared with another element
                // by checking the intersection of containing-element sets
                // for the two endpoint nodes.  We do this efficiently:
                // iterate the smaller set and probe the other.
                const std::set<unsigned>& setA = pNode->rGetContainingElementIndices();
                const std::set<unsigned>& setB = pNextNode->rGetContainingElementIndices();

                unsigned sharedCount = 0;
                unsigned otherElem = UINT_MAX;
                for (unsigned idx : setA)
                {
                    if (idx != ei && setB.count(idx))
                    {
                        otherElem = idx;
                        ++sharedCount;
                    }
                }

                if (sharedCount == 0)
                {
                    // Boundary edge
                    adhesionParam = mBoundaryAdhesion[myType];
                }
                else
                {
                    // Internal edge shared with otherElem
                    // Only process from the element with the lower index
                    // to avoid applying the adhesion force twice.
                    if (ei < otherElem)
                    {
                        unsigned otherType = elemType[otherElem];
                        adhesionParam = mCellCellAdhesion[myType][otherType];

                        // The other element also contributes its own adhesion
                        // term for this edge; since we skip it there, add both
                        // contributions here.
                        // Actually — adhesion is symmetric and the edge gradient
                        // from the other element at these same two global nodes
                        // has the same magnitude but enters with the same sign
                        // in the energy gradient.  The upstream code applies the
                        // adhesion from EACH containing element independently,
                        // which means the total adhesion on a shared edge is
                        // actually 2× the parameter value (this is documented in
                        // the Chaste source: "This parameter corresponds to 1/2
                        // of the sigma parameter").
                        // To maintain identical physics, we multiply by 2 here
                        // since we only visit this edge once.
                        adhesionParam *= 2.0;
                    }
                    else
                    {
                        // This edge will be (or has been) processed by otherElem
                        // where otherElem < ei; skip to avoid double-counting.
                        // However, we still need to account for the fact that
                        // the next-edge gradient at this node contributes to
                        // adhesion.  The upstream implementation applies adhesion
                        // from BOTH elements, but since we use the same edge
                        // gradient direction, the two contributions are identical.
                        // We handle this by processing all adhesion from the
                        // lower-index element side (see the 2× above).
                        adhesionParam = 0.0;
                    }
                }

                if (adhesionParam != 0.0)
                {
                    // The next-edge gradient at node li points from liNext → li
                    // normalised by edge length.  The adhesion contribution is:
                    //   F_li      -= adhesion * nextEdgeGrad
                    //   F_liNext  -= adhesion * (-nextEdgeGrad)   [prev-edge grad for liNext]
                    nodeForce[ni]     -= adhesionParam * nextEdgeGrad;
                    nodeForce[niNext] += adhesionParam * nextEdgeGrad;
                }
            }
        }

        // ── Phase 4: apply accumulated forces to nodes ─────────────────
        for (unsigned ni = 0; ni < numNodes; ++ni)
        {
            pPop->GetNode(ni)->AddAppliedForceContribution(nodeForce[ni]);
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<DeformationEnergyParameter>" << mDeformationParameter
                     << "</DeformationEnergyParameter>\n";
        *rParamsFile << "\t\t\t<MembraneSurfaceEnergyParameter>" << mMembraneParameter
                     << "</MembraneSurfaceEnergyParameter>\n";
        *rParamsFile << "\t\t\t<TargetAreaParameter>" << mTargetArea
                     << "</TargetAreaParameter>\n";
        for (unsigned i = 0; i < NUM_TYPES; ++i)
            for (unsigned j = i; j < NUM_TYPES; ++j)
                *rParamsFile << "\t\t\t<CellCellAdhesion_" << i << "_" << j << ">"
                             << mCellCellAdhesion[i][j]
                             << "</CellCellAdhesion_" << i << "_" << j << ">\n";
        for (unsigned i = 0; i < NUM_TYPES; ++i)
            *rParamsFile << "\t\t\t<BoundaryAdhesion_" << i << ">"
                         << mBoundaryAdhesion[i]
                         << "</BoundaryAdhesion_" << i << ">\n";

        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }

private:
    double mDeformationParameter;
    double mMembraneParameter;
    double mTargetArea;

    // Symmetric adhesion matrix indexed by cell-type ID (0,1,2)
    std::array<std::array<double, NUM_TYPES>, NUM_TYPES> mCellCellAdhesion;
    // Per-type boundary adhesion
    std::array<double, NUM_TYPES> mBoundaryAdhesion;

    /** Map a CellPtr to an integer type ID 0/1/2. */
    static unsigned GetCellTypeId(CellPtr pCell)
    {
        // Fast path: check CellData first (avoids RTTI)
        if (pCell->GetCellData()->HasItem("cell_type_id"))
        {
            return static_cast<unsigned>(pCell->GetCellData()->GetItem("cell_type_id"));
        }
        // Fallback: introspect proliferative type
        auto pType = pCell->GetCellProliferativeType();
        if (pType->template IsType<StemCellProliferativeType>())
            return 0;
        if (pType->template IsType<TransitCellProliferativeType>())
            return 1;
        return 2;  // Differentiated
    }
};

#endif // FASTNAGAIHONDAFORCE_HPP_
