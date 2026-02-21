/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

*/

#ifndef CELLPOLARITYWRITER_HPP_
#define CELLPOLARITYWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"
#include <cmath>

/**
 * A cell writer that outputs cell polarity vectors for VTK visualization.
 *
 * The polarity is stored in cell data as spherical coordinates (theta, phi)
 * and is converted to a Cartesian unit vector for output. This allows
 * visualization in ParaView using Glyph filters.
 *
 * Cell data items used:
 *   - "polarity_theta": polar angle from z-axis (3D) or angle from x-axis (2D)
 *   - "polarity_phi": azimuthal angle (3D only)
 *
 * The output file is called cellpolarity.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Cell polarity" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellPolarityWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    CellPolarityWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellpolarity.dat")
    {
        this->mVtkVectorCellDataName = "Cell polarity";
        this->mOutputScalarData = false;
        this->mOutputVectorData = true;
    }

    /**
     * Overridden GetVectorCellDataForVtkOutput() method.
     *
     * Get the polarity vector associated with a cell, converted from
     * spherical coordinates (theta, phi) to Cartesian (px, py, pz).
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell.
     *
     * @return polarity vector associated with the cell
     */
    c_vector<double, SPACE_DIM> GetVectorCellDataForVtkOutput(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        c_vector<double, SPACE_DIM> polarity = zero_vector<double>(SPACE_DIM);

        try
        {
            double theta = pCell->GetCellData()->GetItem("polarity_theta");

            if constexpr (SPACE_DIM == 3)
            {
                double phi = pCell->GetCellData()->GetItem("polarity_phi");
                polarity[0] = sin(theta) * cos(phi);
                polarity[1] = sin(theta) * sin(phi);
                polarity[2] = cos(theta);
            }
            else if constexpr (SPACE_DIM == 2)
            {
                // In 2D, theta is the angle from the x-axis
                polarity[0] = cos(theta);
                polarity[1] = sin(theta);
            }
            else // SPACE_DIM == 1
            {
                polarity[0] = (theta >= 0) ? 1.0 : -1.0;
            }
        }
        catch (Exception&)
        {
            // If polarity data not found, return zero vector
        }

        return polarity;
    }

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its polarity.
     *
     * Outputs a line of space-separated values of the form:
     * [location index] [cell id] [x-pos] [y-pos] [z-pos] [px] [py] [pz] [theta] [phi]
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        unsigned cell_id = pCell->GetCellId();
        c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
        c_vector<double, SPACE_DIM> polarity = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);

        *this->mpOutStream << location_index << " " << cell_id << " ";
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            *this->mpOutStream << cell_location[i] << " ";
        }
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            *this->mpOutStream << polarity[i] << " ";
        }

        // Also output raw angles for debugging
        double theta = 0.0, phi = 0.0;
        try
        {
            theta = pCell->GetCellData()->GetItem("polarity_theta");
            phi = pCell->GetCellData()->GetItem("polarity_phi");
        }
        catch (Exception&) {}
        *this->mpOutStream << theta << " " << phi << " ";
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPolarityWriter)

#endif /* CELLPOLARITYWRITER_HPP_ */
