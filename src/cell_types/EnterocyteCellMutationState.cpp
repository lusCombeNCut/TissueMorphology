/*
 * EnterocyteCellMutationState.cpp
 */
#include "EnterocyteCellMutationState.hpp"

EnterocyteCellMutationState::EnterocyteCellMutationState()
    : AbstractCellMutationState(6)
{
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(EnterocyteCellMutationState)
