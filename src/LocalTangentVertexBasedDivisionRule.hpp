/*
 * LocalTangentVertexBasedDivisionRule.hpp
 *
 * A cell division rule for 2D vertex-based populations that divides along
 * the LOCAL tangent direction computed from neighboring cells.
 *
 * Unlike RadialVertexBasedDivisionRule which uses the global centroid,
 * this rule computes the local surface direction from the cell's immediate
 * neighbors. This correctly handles:
 *   - Budding organoids where the global centroid doesn't reflect local geometry
 *   - Non-convex shapes
 *   - Any curved monolayer topology
 *
 * Algorithm:
 *   1. Find the dividing cell's neighboring elements (cells sharing edges)
 *   2. Compute the tangent as the vector between the two side neighbors
 *   3. Division axis = tangent direction (daughters placed side-by-side)
 *
 * This ensures daughters are placed tangentially along the monolayer
 * regardless of the global shape.
 */
#ifndef LOCALTANGENTVERTEXBASEDDIVISIONRULE_HPP_
#define LOCALTANGENTVERTEXBASEDDIVISIONRULE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"
#include <set>
#include <vector>
#include <algorithm>
#include <cmath>

template <unsigned SPACE_DIM>
class LocalTangentVertexBasedDivisionRule : public AbstractVertexBasedDivisionRule<SPACE_DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractVertexBasedDivisionRule<SPACE_DIM> >(*this);
    }

public:

    LocalTangentVertexBasedDivisionRule() {}

    virtual ~LocalTangentVertexBasedDivisionRule() {}

    /**
     * Compute local tangent direction from neighboring cells.
     * Returns a unit vector along the local monolayer surface.
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(
        CellPtr pParentCell,
        VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
    {
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
        c_vector<double, SPACE_DIM> cell_centre =
            rCellPopulation.rGetMesh().GetCentroidOfElement(elem_index);

        // Get neighboring elements (cells sharing an edge with this cell)
        std::set<unsigned> neighbor_indices = 
            rCellPopulation.rGetMesh().GetNeighbouringElementIndices(elem_index);

        if (neighbor_indices.size() < 2)
        {
            // Not enough neighbors - fall back to short axis
            return rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
        }

        // Collect neighbor centroids with their angular position relative to this cell
        std::vector<std::pair<double, c_vector<double, SPACE_DIM>>> angular_neighbors;
        
        for (unsigned nbr_idx : neighbor_indices)
        {
            c_vector<double, SPACE_DIM> nbr_centre =
                rCellPopulation.rGetMesh().GetCentroidOfElement(nbr_idx);
            
            c_vector<double, SPACE_DIM> delta = nbr_centre - cell_centre;
            double angle = std::atan2(delta[1], delta[0]);
            angular_neighbors.push_back(std::make_pair(angle, nbr_centre));
        }

        // Sort by angle
        std::sort(angular_neighbors.begin(), angular_neighbors.end(),
            [](const std::pair<double, c_vector<double, SPACE_DIM>>& a,
               const std::pair<double, c_vector<double, SPACE_DIM>>& b) {
                return a.first < b.first;
            });

        // For a 2D ring/annulus, we want the two "side" neighbors
        // These are the neighbors that form the local tangent direction
        // In an annulus, each cell typically has 2 side neighbors and possibly
        // inner/outer neighbors if there's multi-layering
        
        // Strategy: Use the two neighbors that are most opposite to each other
        // (largest angular separation) - these define the tangent
        double max_sep = 0.0;
        c_vector<double, SPACE_DIM> best_tangent = zero_vector<double>(SPACE_DIM);
        
        for (size_t i = 0; i < angular_neighbors.size(); ++i)
        {
            for (size_t j = i + 1; j < angular_neighbors.size(); ++j)
            {
                double sep = std::abs(angular_neighbors[j].first - angular_neighbors[i].first);
                // Normalize to [0, π]
                if (sep > M_PI) sep = 2.0 * M_PI - sep;
                
                if (sep > max_sep)
                {
                    max_sep = sep;
                    // Tangent direction connects these two neighbors
                    best_tangent = angular_neighbors[j].second - angular_neighbors[i].second;
                }
            }
        }

        double tangent_length = norm_2(best_tangent);
        if (tangent_length < 1e-10)
        {
            // Degenerate - fall back to short axis
            return rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
        }

        // Normalize tangent
        best_tangent /= tangent_length;
        
        // Return the LOCAL NORMAL (perpendicular to tangent)
        // In 2D, rotate tangent by 90 degrees: (x, y) -> (-y, x)
        // Division along the normal creates new edge perpendicular to it,
        // placing daughter cells side-by-side along the tangent (in the monolayer)
        c_vector<double, SPACE_DIM> local_normal = zero_vector<double>(SPACE_DIM);
        local_normal[0] = -best_tangent[1];
        local_normal[1] = best_tangent[0];
        
        return local_normal;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LocalTangentVertexBasedDivisionRule)

namespace boost {
namespace serialization {

template <class Archive, unsigned DIM>
inline void save_construct_data(
    Archive& ar, const LocalTangentVertexBasedDivisionRule<DIM>* t, const unsigned int file_version)
{
}

template <class Archive, unsigned DIM>
inline void load_construct_data(
    Archive& ar, LocalTangentVertexBasedDivisionRule<DIM>* t, const unsigned int file_version)
{
    ::new(t) LocalTangentVertexBasedDivisionRule<DIM>();
}

} // namespace serialization
} // namespace boost

#endif // LOCALTANGENTVERTEXBASEDDIVISIONRULE_HPP_
