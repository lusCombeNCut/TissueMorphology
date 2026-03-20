/*
 * PanethCellMutationState.hpp
 *
 * Paneth cell mutation state for the 4-cell-type
 * stochastic crypt budding model (Montes-Olivas et al.).
 */
#ifndef PANETHCELLMUTATIONSTATE_HPP_
#define PANETHCELLMUTATIONSTATE_HPP_

#include "AbstractCellMutationState.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class PanethCellMutationState : public AbstractCellMutationState
{
private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellMutationState>(*this);
    }

public:
    PanethCellMutationState();
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PanethCellMutationState)

#endif /* PANETHCELLMUTATIONSTATE_HPP_ */
