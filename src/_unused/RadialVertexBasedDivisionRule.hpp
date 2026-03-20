/*
 * RadialVertexBasedDivisionRule.hpp
 *
 * A cell division rule for 2D vertex-based populations that always divides
 * along the radial direction (from the population centroid to the cell centre).
 *
 * This ensures that daughter cells are placed side-by-side within the
 * monolayer, rather than stacking inward/outward which would break
 * the single-layer topology.
 *
 * Rationale:
 *   In an annular vertex mesh, the default ShortAxisVertexBasedDivisionRule
 *   divides along the cell's shortest moment-of-inertia axis. When a cell
 *   becomes taller (radially) than it is wide (tangentially) — e.g. after
 *   a T1 swap — the short axis rotates to be tangential, producing a
 *   radial cut that stacks daughters. This rule prevents that by always
 *   using the geometric radial direction.
 *
 * The returned division vector is pointed radially from the population
 * centroid through the cell centre. DivideElementAlongGivenAxis() then
 * creates the new inter-daughter edge aligned with this vector, placing
 * the two daughter cells tangentially side-by-side.
 */
#ifndef RADIALVERTEXBASEDDIVISIONRULE_HPP_
#define RADIALVERTEXBASEDDIVISIONRULE_HPP_

#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "AbstractVertexBasedDivisionRule.hpp"

template <unsigned SPACE_DIM>
class RadialVertexBasedDivisionRule : public AbstractVertexBasedDivisionRule<SPACE_DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractVertexBasedDivisionRule<SPACE_DIM> >(*this);
    }

public:

    RadialVertexBasedDivisionRule() {}

    virtual ~RadialVertexBasedDivisionRule() {}

    /**
     * Return the radial direction from the population centroid to the cell
     * centre.  The element will be divided along this direction, so the two
     * daughter cells sit side-by-side tangentially on the monolayer.
     *
     * Falls back to the short-axis rule if the cell happens to sit exactly
     * at the centroid (degenerate case).
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(
        CellPtr pParentCell,
        VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation)
    {
        // Population centroid
        c_vector<double, SPACE_DIM> centroid = rCellPopulation.GetCentroidOfCellPopulation();

        // Cell centre
        unsigned elem_index = rCellPopulation.GetLocationIndexUsingCell(pParentCell);
        c_vector<double, SPACE_DIM> cell_centre =
            rCellPopulation.rGetMesh().GetCentroidOfElement(elem_index);

        // Radial direction
        c_vector<double, SPACE_DIM> radial = cell_centre - centroid;
        double r = norm_2(radial);

        if (r < 1e-10)
        {
            // Degenerate: cell sits at the centroid — fall back to short axis
            return rCellPopulation.rGetMesh().GetShortAxisOfElement(elem_index);
        }

        // Return normalised radial direction
        return radial / r;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RadialVertexBasedDivisionRule)

namespace boost {
namespace serialization {

template <class Archive, unsigned DIM>
inline void save_construct_data(
    Archive& ar, const RadialVertexBasedDivisionRule<DIM>* t, const unsigned int file_version)
{
}

template <class Archive, unsigned DIM>
inline void load_construct_data(
    Archive& ar, RadialVertexBasedDivisionRule<DIM>* t, const unsigned int file_version)
{
    ::new(t) RadialVertexBasedDivisionRule<DIM>();
}

} // namespace serialization
} // namespace boost

#endif // RADIALVERTEXBASEDDIVISIONRULE_HPP_
