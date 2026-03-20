/*
 * NonNeighbourCellRepulsionForce.hpp
 *
 * Short-range repulsion between cell centres for non-neighbouring cells
 * in a vertex-based population. This prevents different tissue folds from
 * passing through one another when they come into close proximity.
 */
#ifndef NONNEIGHBOURCELLREPULSIONFORCE_HPP_
#define NONNEIGHBOURCELLREPULSIONFORCE_HPP_

#include <set>
#include <vector>

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
class NonNeighbourCellRepulsionForce : public AbstractForce<DIM>
{
private:
    double mRepulsionStiffness;
    double mInteractionCutoff;

public:
    NonNeighbourCellRepulsionForce()
        : AbstractForce<DIM>(),
          mRepulsionStiffness(10.0),
          mInteractionCutoff(1.25)
    {
    }

    virtual ~NonNeighbourCellRepulsionForce() = default;

    void SetRepulsionStiffness(double stiffness)
    {
        mRepulsionStiffness = stiffness;
    }

    void SetInteractionCutoff(double cutoff)
    {
        mInteractionCutoff = cutoff;
    }

    double GetRepulsionStiffness() const
    {
        return mRepulsionStiffness;
    }

    double GetInteractionCutoff() const
    {
        return mInteractionCutoff;
    }

    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override
    {
        if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == nullptr)
        {
            EXCEPTION("NonNeighbourCellRepulsionForce is to be used with a VertexBasedCellPopulation only");
        }

        if (mRepulsionStiffness <= 0.0 || mInteractionCutoff <= 0.0)
        {
            return;
        }

        VertexBasedCellPopulation<DIM>* p_population =
            static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        std::vector<CellPtr> cells;
        std::vector<unsigned> location_indices;
        std::vector<c_vector<double, DIM>> centres;
        std::vector<std::set<unsigned>> neighbours;

        cells.reserve(p_population->GetNumRealCells());
        location_indices.reserve(p_population->GetNumRealCells());
        centres.reserve(p_population->GetNumRealCells());
        neighbours.reserve(p_population->GetNumRealCells());

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = p_population->Begin();
             cell_iter != p_population->End();
             ++cell_iter)
        {
            CellPtr p_cell = *cell_iter;
            unsigned location_index = p_population->GetLocationIndexUsingCell(p_cell);

            cells.push_back(p_cell);
            location_indices.push_back(location_index);
            centres.push_back(p_population->GetLocationOfCellCentre(p_cell));
            neighbours.push_back(p_population->GetNeighbouringLocationIndices(p_cell));
        }

        const double cutoff_sq = mInteractionCutoff * mInteractionCutoff;

        for (unsigned i = 0; i < cells.size(); ++i)
        {
            for (unsigned j = i + 1; j < cells.size(); ++j)
            {
                const unsigned loc_i = location_indices[i];
                const unsigned loc_j = location_indices[j];

                // Skip direct topological neighbours: adhesion/edge mechanics already handle them.
                if (neighbours[i].count(loc_j) > 0u)
                {
                    continue;
                }

                c_vector<double, DIM> separation = centres[i] - centres[j];
                const double distance_sq = inner_prod(separation, separation);

                if (distance_sq < 1e-12 || distance_sq >= cutoff_sq)
                {
                    continue;
                }

                const double distance = std::sqrt(distance_sq);
                c_vector<double, DIM> unit = separation / distance;
                const double overlap = mInteractionCutoff - distance;
                const double magnitude = mRepulsionStiffness * (overlap / mInteractionCutoff);

                c_vector<double, DIM> force_i = magnitude * unit;
                c_vector<double, DIM> force_j = -force_i;

                VertexElement<DIM, DIM>* p_elem_i = p_population->GetElement(loc_i);
                VertexElement<DIM, DIM>* p_elem_j = p_population->GetElement(loc_j);

                c_vector<double, DIM> per_node_i = force_i / static_cast<double>(p_elem_i->GetNumNodes());
                c_vector<double, DIM> per_node_j = force_j / static_cast<double>(p_elem_j->GetNumNodes());

                for (unsigned local_idx = 0; local_idx < p_elem_i->GetNumNodes(); ++local_idx)
                {
                    p_elem_i->GetNode(local_idx)->AddAppliedForceContribution(per_node_i);
                }
                for (unsigned local_idx = 0; local_idx < p_elem_j->GetNumNodes(); ++local_idx)
                {
                    p_elem_j->GetNode(local_idx)->AddAppliedForceContribution(per_node_j);
                }
            }
        }
    }

    void OutputForceParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<RepulsionStiffness>" << mRepulsionStiffness
                     << "</RepulsionStiffness>\n";
        *rParamsFile << "\t\t\t<InteractionCutoff>" << mInteractionCutoff
                     << "</InteractionCutoff>\n";

        AbstractForce<DIM>::OutputForceParameters(rParamsFile);
    }
};

#endif // NONNEIGHBOURCELLREPULSIONFORCE_HPP_