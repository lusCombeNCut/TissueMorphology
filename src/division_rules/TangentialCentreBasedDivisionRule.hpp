/*
 * TangentialCentreBasedDivisionRule.hpp
 *
 * A cell division rule for centre-based (node-based) populations that places
 * daughter cells TANGENTIALLY — perpendicular to the cell's polarity vector.
 *
 * The cell polarity (stored as polarity_theta, polarity_phi in cell data)
 * defines the outward-facing normal. Daughters are displaced perpendicular
 * to this, keeping both on the monolayer surface.
 *
 * Falls back to centroid-based radial direction if polarity data is absent.
 */
#ifndef TANGENTIALCENTREBASEDDIVISIONRULE_HPP_
#define TANGENTIALCENTREBASEDDIVISIONRULE_HPP_

#include <cmath>
#include <boost/serialization/base_object.hpp>
#include "ChasteSerialization.hpp"
#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class TangentialCentreBasedDivisionRule : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    TangentialCentreBasedDivisionRule() {}

    virtual ~TangentialCentreBasedDivisionRule() {}

    /**
     * Return a pair of positions for the parent and daughter cells.
     *
     * Both are displaced tangentially (perpendicular to the radial direction
     * from the population centroid) by ±half the division separation distance.
     *
     * @param pParentCell  the dividing cell
     * @param rCellPopulation  the cell population
     * @return a pair of (daughter_position, parent_position)
     */
    virtual std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM>>
    CalculateCellDivisionVector(
        CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation)
    {
        // Get parent's current location
        c_vector<double, SPACE_DIM> parent_location =
            rCellPopulation.GetLocationOfCellCentre(pParentCell);

        // Get the division separation distance
        double separation = rCellPopulation.GetMeinekeDivisionSeparation();

        // Use cell polarity as the outward (radial) direction
        // Polarity is stored as spherical angles in cell data
        c_vector<double, SPACE_DIM> polarity = zero_vector<double>(SPACE_DIM);
        bool has_polarity = false;

        try
        {
            double theta = pParentCell->GetCellData()->GetItem("polarity_theta");

            if constexpr (SPACE_DIM == 3)
            {
                double phi = pParentCell->GetCellData()->GetItem("polarity_phi");
                polarity[0] = sin(theta) * cos(phi);
                polarity[1] = sin(theta) * sin(phi);
                polarity[2] = cos(theta);
            }
            else // SPACE_DIM == 2
            {
                polarity[0] = cos(theta);
                polarity[1] = sin(theta);
            }

            double p_norm = norm_2(polarity);
            if (p_norm > 1e-10)
            {
                polarity /= p_norm;
                has_polarity = true;
            }
        }
        catch (Exception&)
        {
            // Polarity data not available — fall back to centroid-based
        }

        // If no valid polarity, fall back to centroid-based radial direction
        if (!has_polarity)
        {
            c_vector<double, SPACE_DIM> centroid =
                rCellPopulation.GetCentroidOfCellPopulation();
            polarity = parent_location - centroid;
            double r = norm_2(polarity);
            if (r > 1e-10)
            {
                polarity /= r;
            }
            else
            {
                // Cell at centroid with no polarity — random direction
                RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                double angle = 2.0 * M_PI * p_gen->ranf();
                polarity[0] = cos(angle);
                polarity[1] = sin(angle);
                if constexpr (SPACE_DIM == 3)
                {
                    polarity[2] = 0.0;
                }
            }
        }

        // Compute tangential direction: perpendicular to polarity (the outward normal)
        c_vector<double, SPACE_DIM> tangential;

        if constexpr (SPACE_DIM == 2)
        {
            // 90° rotation of polarity: (-py, px)
            tangential[0] = -polarity[1];
            tangential[1] =  polarity[0];
        }
        else // SPACE_DIM == 3
        {
            // Project a random vector onto the plane perpendicular to polarity
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
            c_vector<double, 3> random_vec;
            random_vec[0] = p_gen->ranf() - 0.5;
            random_vec[1] = p_gen->ranf() - 0.5;
            random_vec[2] = p_gen->ranf() - 0.5;

            // t = v - (v·p) * p
            double dot = inner_prod(random_vec, polarity);
            tangential = random_vec - dot * polarity;

            double t_norm = norm_2(tangential);
            if (t_norm < 1e-10)
            {
                // Random vector was parallel to polarity — use a fixed fallback
                random_vec[0] = 1.0;
                random_vec[1] = 0.0;
                random_vec[2] = 0.0;
                dot = inner_prod(random_vec, polarity);
                tangential = random_vec - dot * polarity;
                t_norm = norm_2(tangential);
            }
            tangential /= t_norm;
        }

        // Place parent and daughter symmetrically along the tangential direction
        c_vector<double, SPACE_DIM> daughter_location =
            parent_location + 0.5 * separation * tangential;
        c_vector<double, SPACE_DIM> new_parent_location =
            parent_location - 0.5 * separation * tangential;

        return std::make_pair(daughter_location, new_parent_location);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(TangentialCentreBasedDivisionRule)

namespace boost {
namespace serialization {

template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive& ar,
    const TangentialCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>* t,
    const unsigned int file_version)
{
}

template <class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive& ar,
    TangentialCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>* t,
    const unsigned int file_version)
{
    ::new(t) TangentialCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>();
}

} // namespace serialization
} // namespace boost

#endif // TANGENTIALCENTREBASEDDIVISIONRULE_HPP_
