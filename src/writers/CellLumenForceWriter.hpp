/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

*/

#ifndef CELLLUMENFORCEWRITER_HPP_
#define CELLLUMENFORCEWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A cell writer that outputs lumen pressure force vectors for VTK visualization.
 *
 * The force components are read from cell data items set by LumenPressureForce
 * when its SetWriteForce(true) option is enabled.
 *
 * Cell data items used:
 *   - "lumen_force_x"
 *   - "lumen_force_y"
 *   - "lumen_force_z"
 *
 * In ParaView, apply a Glyph filter on the "Lumen force" vector to visualise
 * the direction and magnitude of the lumen pressure acting on each cell.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellLumenForceWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    CellLumenForceWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celllumenforce.dat")
    {
        this->mVtkVectorCellDataName = "Lumen force";
        this->mOutputScalarData = false;
        this->mOutputVectorData = true;
    }

    c_vector<double, SPACE_DIM> GetVectorCellDataForVtkOutput(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        c_vector<double, SPACE_DIM> force = zero_vector<double>(SPACE_DIM);

        try
        {
            force[0] = pCell->GetCellData()->GetItem("lumen_force_x");
            force[1] = pCell->GetCellData()->GetItem("lumen_force_y");
            if constexpr (SPACE_DIM >= 3)
            {
                force[2] = pCell->GetCellData()->GetItem("lumen_force_z");
            }
        }
        catch (Exception&)
        {
            // Force data not yet available — return zero vector
        }

        return force;
    }

    virtual void VisitCell(
        CellPtr pCell,
        AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
        unsigned cell_id = pCell->GetCellId();
        c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(pCell);
        c_vector<double, SPACE_DIM> force = GetVectorCellDataForVtkOutput(pCell, pCellPopulation);

        *this->mpOutStream << location_index << " " << cell_id << " ";
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            *this->mpOutStream << cell_location[i] << " ";
        }
        for (unsigned i = 0; i < SPACE_DIM; i++)
        {
            *this->mpOutStream << force[i] << " ";
        }
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellLumenForceWriter)

#endif /* CELLLUMENFORCEWRITER_HPP_ */
