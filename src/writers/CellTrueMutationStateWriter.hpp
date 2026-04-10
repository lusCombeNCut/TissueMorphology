/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

This file is part of Chaste.

*/

#ifndef CELLTRUEMUTATIONSTATEWRITER_HPP_
#define CELLTRUEMUTATIONSTATEWRITER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractCellWriter.hpp"

/**
 * A cell writer that outputs the true mutation state colour for each cell.
 *
 * This replaces the default CellMutationStatesWriter which has a bug:
 * when a CellLabel is present (e.g. from ContactInhibitionCellCycleModel),
 * it substitutes CellLabel::GetColour() (default 5) for the actual
 * mutation state. Since colour 5 = PanethCellMutationState, every
 * contact-inhibited cell appears as Paneth in VTU output.
 *
 * This writer always outputs GetMutationState()->GetColour(), ignoring
 * any CellLabel.
 *
 * The output file is called results.vizmutationstates by default. The VTK
 * cell data name is "Mutation states" to match the original writer.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellTrueMutationStateWriter : public AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    CellTrueMutationStateWriter()
        : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.vizmutationstates")
    {
        this->mVtkCellDataName = "Mutation states";
    }

    double GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        return pCell->GetMutationState()->GetColour();
    }

    virtual void VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
    {
        *this->mpOutStream << pCell->GetMutationState()->GetColour() << " ";
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellTrueMutationStateWriter)

#endif /* CELLTRUEMUTATIONSTATEWRITER_HPP_ */
