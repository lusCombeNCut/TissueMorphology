/*
 * TACellMutationState.cpp
 */
#include "TACellMutationState.hpp"

TACellMutationState::TACellMutationState()
    : AbstractCellMutationState(4)
{
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(TACellMutationState)
