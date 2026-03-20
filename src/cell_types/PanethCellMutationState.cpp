/*
 * PanethCellMutationState.cpp
 */
#include "PanethCellMutationState.hpp"

PanethCellMutationState::PanethCellMutationState()
    : AbstractCellMutationState(5)
{
}

#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PanethCellMutationState)
