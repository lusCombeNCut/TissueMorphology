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
#include "MutableElement.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MonolayerVertexBasedCellPopulation.hpp"

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
        auto* p_vertex_pop = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);
        auto* p_monolayer_pop = dynamic_cast<MonolayerVertexBasedCellPopulation<DIM>*>(&rCellPopulation);

        if (!p_vertex_pop && !p_monolayer_pop)
        {
            EXCEPTION("NonNeighbourCellRepulsionForce is to be used with a VertexBasedCellPopulation or MonolayerVertexBasedCellPopulation only");
        }

        if (mRepulsionStiffness <= 0.0 || mInteractionCutoff <= 0.0)
        {
            return;
        }

        std::vector<CellPtr> cells;
        std::vector<unsigned> location_indices;
        std::vector<c_vector<double, DIM>> centres;
        std::vector<std::set<unsigned>> neighbours;

        const unsigned num_cells = rCellPopulation.GetNumRealCells();
        cells.reserve(num_cells);
        location_indices.reserve(num_cells);
        centres.reserve(num_cells);
        neighbours.reserve(num_cells);

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            CellPtr p_cell = *cell_iter;
            unsigned location_index = rCellPopulation.GetLocationIndexUsingCell(p_cell);

            cells.push_back(p_cell);
            location_indices.push_back(location_index);
            centres.push_back(rCellPopulation.GetLocationOfCellCentre(p_cell));
            neighbours.push_back(rCellPopulation.GetNeighbouringLocationIndices(p_cell));
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

                // Get elements via the appropriate population type (both return
                // subclasses of MutableElement, which provides GetNumNodes/GetNode).
                MutableElement<DIM, DIM>* p_elem_i;
                MutableElement<DIM, DIM>* p_elem_j;
                if (p_vertex_pop)
                {
                    p_elem_i = p_vertex_pop->GetElement(loc_i);
                    p_elem_j = p_vertex_pop->GetElement(loc_j);
                }
                else
                {
                    p_elem_i = p_monolayer_pop->GetElement(loc_i);
                    p_elem_j = p_monolayer_pop->GetElement(loc_j);
                }

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