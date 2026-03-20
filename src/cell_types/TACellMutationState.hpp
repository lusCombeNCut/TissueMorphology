/*
 * TACellMutationState.hpp
 *
 * Transit-amplifying cell mutation state for the 4-cell-type
 * stochastic crypt budding model (Montes-Olivas et al.).
 */
#ifndef TACELLMUTATIONSTATE_HPP_
#define TACELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class TACellMutationState : public AbstractCellMutationState
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    TACellMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(TACellMutationState)

#endif /* TACELLMUTATIONSTATE_HPP_ */
