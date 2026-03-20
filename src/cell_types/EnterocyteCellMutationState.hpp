/*
 * EnterocyteCellMutationState.hpp
 *
 * Enterocyte (differentiated absorptive) cell mutation state for the
 * 4-cell-type stochastic crypt budding model (Montes-Olivas et al.).
 */
#ifndef ENTEROCYTECELLMUTATIONSTATE_HPP_
#define ENTEROCYTECELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class EnterocyteCellMutationState : public AbstractCellMutationState
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    EnterocyteCellMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(EnterocyteCellMutationState)

#endif /* ENTEROCYTECELLMUTATIONSTATE_HPP_ */
