/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

*/

#ifndef CELLCONTACTINHIBITIONSTATUSWRITER_HPP_
#define CELLCONTACTINHIBITIONSTATUSWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"
#include "CellLabel.hpp"

/**
 * A cell writer that outputs whether each cell is currently contact inhibited.
 *
 * Outputs 1.0 if the cell has a CellLabel (indicating contact inhibition is
 * preventing division), and 0.0 otherwise.
 *
 * The CellLabel is added by ContactInhibitionCellCycleModel and its subclasses
 * when a cell in G1 phase has volume below the quiescent threshold.
 *
 * The output file is called contactinhibition.dat by default. If VTK is switched on,
 * then the writer also specifies the VTK output for each cell, which is stored in
 * the VTK cell data "Contact inhibited" by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellContactInhibitionStatusWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
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
    CellContactInhibitionStatusWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("contactinhibition.dat")
    {
        this->mVtkCellDataName = "Contact inhibited";
    }

    /**
     * Overridden GetCellDataForVtkOutput() method.
     *
     * Get a double associated with a cell indicating whether it is contact inhibited.
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     *
     * @return 1.0 if cell has CellLabel (contact inhibited), 0.0 otherwise
     */
    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        return pCell->HasCellProperty<CellLabel>() ? 1.0 : 0.0;
    }

    /**
     * Overridden VisitCell() method.
     *
     * Visit a cell and write its contact inhibition status.
     *
     * Outputs a line of space-separated values of the form:
     * ...[contact inhibited status (0 or 1)] ...
     *
     * @param pCell a cell
     * @param pCellPopulation a pointer to the cell population owning the cell
     */
    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        unsigned status = pCell->HasCellProperty<CellLabel>() ? 1 : 0;
        *this->mpOutStream << status << " ";
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellContactInhibitionStatusWriter)

#endif /* CELLCONTACTINHIBITIONSTATUSWRITER_HPP_ */
