/*
 * UniformContactInhibitionCellCycleModel.hpp
 *
 * Contact-inhibition cell cycle model with uniform-random total cycle
 * duration.  Inherits all quiescent-volume logic from
 * ContactInhibitionCellCycleModel; the only change is a stochastic
 * SetG1Duration() that draws total cycle time from U(min, max) for
 * stem cells and scales by a ratio for transit cells.
 *
 * Phase durations S, G2 are set to 0 and M to 1 h so that
 * G1 ≈ totalCycle − 1 h.  When contact-inhibition triggers
 * quiescence the G1 phase is extended exactly as in the parent class.
 *
 *   Stem:         totalCycle ~ U(mTotalCycleMin, mTotalCycleMax)
 *   Transit:      totalCycle ~ U(mTotalCycleMin, mTotalCycleMax) × mTransitCycleRatio
 *   Differentiated: G1 = ∞  (never divides)
 */
#ifndef UNIFORMCONTACTINHIBITIONCELLCYCLEMODEL_HPP_
#define UNIFORMCONTACTINHIBITIONCELLCYCLEMODEL_HPP_

#include "ContactInhibitionCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>
#include <algorithm>
#include <cassert>

class UniformContactInhibitionCellCycleModel
    : public ContactInhibitionCellCycleModel
{
private:
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<ContactInhibitionCellCycleModel>(*this);
        archive & mTotalCycleMin;
        archive & mTotalCycleMax;
        archive & mTransitCycleRatio;
    }

    /** Minimum total cell-cycle duration for stem cells (hours). */
    double mTotalCycleMin;

    /** Maximum total cell-cycle duration for stem cells (hours). */
    double mTotalCycleMax;

    /**
     * Transit-amplifying cell cycle duration as a fraction of the
     * stem-cell cycle.  1.0 = same length, 0.5 = half.
     */
    double mTransitCycleRatio;

protected:

    /**
     * Override: draw total cycle from U(min, max), compute G1 = total − (S+G2+M).
     */
    void SetG1Duration() override
    {
        assert(mpCell != nullptr);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double sgm = GetSDuration() + GetG2Duration() + GetMDuration();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            double total = mTotalCycleMin
                         + p_gen->ranf() * (mTotalCycleMax - mTotalCycleMin);
            mG1Duration = std::max(mMinimumGapDuration, total - sgm);
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            double total = mTotalCycleMin
                         + p_gen->ranf() * (mTotalCycleMax - mTotalCycleMin);
            total *= mTransitCycleRatio;
            mG1Duration = std::max(mMinimumGapDuration, total - sgm);
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mG1Duration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }
    }

    /** Copy constructor (used by CreateCellCycleModel). */
    UniformContactInhibitionCellCycleModel(
        const UniformContactInhibitionCellCycleModel& rModel)
        : ContactInhibitionCellCycleModel(rModel),
          mTotalCycleMin(rModel.mTotalCycleMin),
          mTotalCycleMax(rModel.mTotalCycleMax),
          mTransitCycleRatio(rModel.mTransitCycleRatio)
    {
    }

public:
    UniformContactInhibitionCellCycleModel()
        : ContactInhibitionCellCycleModel(),
          mTotalCycleMin(12.0),
          mTotalCycleMax(14.0),
          mTransitCycleRatio(1.0)
    {
        // Collapse S and G2 so G1 ≈ totalCycle − M
        SetSDuration(0.0);
        SetG2Duration(0.0);
        SetMDuration(1.0);
    }

    AbstractCellCycleModel* CreateCellCycleModel() override
    {
        return new UniformContactInhibitionCellCycleModel(*this);
    }

    // ---- Setters / Getters -----------------------------------------

    void SetTotalCycleMin(double val)  { mTotalCycleMin = val; }
    void SetTotalCycleMax(double val)  { mTotalCycleMax = val; }
    void SetTransitCycleRatio(double r) { mTransitCycleRatio = r; }

    double GetTotalCycleMin()  const { return mTotalCycleMin; }
    double GetTotalCycleMax()  const { return mTotalCycleMax; }
    double GetTransitCycleRatio() const { return mTransitCycleRatio; }

    void OutputCellCycleModelParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<TotalCycleMin>"
                     << mTotalCycleMin << "</TotalCycleMin>\n";
        *rParamsFile << "\t\t\t<TotalCycleMax>"
                     << mTotalCycleMax << "</TotalCycleMax>\n";
        *rParamsFile << "\t\t\t<TransitCycleRatio>"
                     << mTransitCycleRatio << "</TransitCycleRatio>\n";
        ContactInhibitionCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(UniformContactInhibitionCellCycleModel)

#endif // UNIFORMCONTACTINHIBITIONCELLCYCLEMODEL_HPP_
