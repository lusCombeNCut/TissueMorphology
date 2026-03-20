/*
 * StochasticFourTypeCellCycleModel.hpp
 *
 * Contact-inhibition cell cycle model with stochastic, time-dependent
 * cell fate transitions for 4 intestinal cell types, following the
 * approach of Montes-Olivas et al. (2023).
 *
 * Cell types are distinguished by MutationState:
 *   - WildTypeCellMutationState  → Stem Cell (SC)
 *   - TACellMutationState        → Transit-Amplifying cell (TA)
 *   - PanethCellMutationState    → Paneth Cell (PC)  — non-proliferative
 *   - EnterocyteCellMutationState→ Enterocyte (EC)   — non-proliferative
 *
 * All proliferating cells use TransitCellProliferativeType (SC and TA);
 * PC and EC are set to DifferentiatedCellProliferativeType (G1 = ∞).
 *
 * Division rules:
 *   SC → SC  with probability  p_sc_sc
 *   SC → PC  with probability  p_sc_pc
 *   SC → TA  with probability  1 - p_sc_sc - p_sc_pc
 *
 *   TA → TA  with probability  p_ta_ta_early  (t < mTransitionTime)
 *   TA → EC  with probability  1 - p_ta_ta_early
 *          then after transition time:
 *   TA → TA  with probability  p_ta_ta_late   (t >= mTransitionTime)
 *   TA → EC  with probability  1 - p_ta_ta_late
 *
 *   PC and EC do not divide.
 *
 * Contact inhibition is inherited from ContactInhibitionCellCycleModel:
 *   - Cells in G1 are arrested (quiescent) when V_curr/V_eq < f_q
 *   - G1 duration is extended while quiescent
 *
 * Reference:
 *   Montes-Olivas, S., Marucci, L. & Homer, M. (2023). In-silico intestinal
 *   crypt–organoid model. Cell Proliferation, 56(1), e13370.
 */
#ifndef STOCHASTICFOURTYPECELLCYCLEMODEL_HPP_
#define STOCHASTICFOURTYPECELLCYCLEMODEL_HPP_

#include "ContactInhibitionCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "SimulationTime.hpp"

#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "WildTypeCellMutationState.hpp"
#include "TACellMutationState.hpp"
#include "PanethCellMutationState.hpp"
#include "EnterocyteCellMutationState.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include <cfloat>
#include <algorithm>
#include <cassert>

class StochasticFourTypeCellCycleModel
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
        archive & mProbStemToStem;
        archive & mProbStemToPaneth;
        archive & mProbTaToTaEarly;
        archive & mProbTaToTaLate;
        archive & mTransitionTime;
    }

    // ---- Cell cycle duration parameters ----

    /** Minimum total cell-cycle duration for stem cells (hours). */
    double mTotalCycleMin;

    /** Maximum total cell-cycle duration for stem cells (hours). */
    double mTotalCycleMax;

    /** TA cycle duration as a fraction of the stem-cell cycle. */
    double mTransitCycleRatio;

    // ---- Stochastic transition probabilities ----

    /**
     * Probability that a stem cell division produces another stem cell
     * (as the daughter). The remaining probability is split between
     * PC and TA.
     *   P(daughter = SC) = mProbStemToStem
     *   P(daughter = PC) = mProbStemToPaneth
     *   P(daughter = TA) = 1 - mProbStemToStem - mProbStemToPaneth
     */
    double mProbStemToStem;
    double mProbStemToPaneth;

    /**
     * Probability that a TA cell division produces another TA cell
     * (vs. an enterocyte). Changes at mTransitionTime.
     *   t < mTransitionTime:  P(daughter = TA) = mProbTaToTaEarly
     *   t >= mTransitionTime: P(daughter = TA) = mProbTaToTaLate
     */
    double mProbTaToTaEarly;
    double mProbTaToTaLate;

    /**
     * Simulation time (hours) at which TA transition probabilities
     * switch from early to late values. In Sandra's model this is
     * day 5 = 120 hours.
     */
    double mTransitionTime;

protected:

    /**
     * Override: draw total cycle from U(min, max), compute G1 = total − (S+G2+M).
     * PC and EC cells get G1 = ∞ (never divide).
     * SC and TA cells are both TransitCellProliferativeType, distinguished
     * by mutation state; both can proliferate.
     */
    void SetG1Duration() override
    {
        assert(mpCell != nullptr);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double sgm = GetSDuration() + GetG2Duration() + GetMDuration();

        // PC and EC are DifferentiatedCellProliferativeType → never divide
        if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mG1Duration = DBL_MAX;
            return;
        }

        // Check mutation state: PC and EC should never proliferate
        // (belt-and-suspenders: they should already be Differentiated,
        //  but guard against incorrect type assignment)
        if (mpCell->GetMutationState()->IsType<PanethCellMutationState>() ||
            mpCell->GetMutationState()->IsType<EnterocyteCellMutationState>())
        {
            mG1Duration = DBL_MAX;
            return;
        }

        // SC or TA — both proliferate
        if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            // Stem cell cycle duration
            double total = mTotalCycleMin
                         + p_gen->ranf() * (mTotalCycleMax - mTotalCycleMin);
            mG1Duration = std::max(mMinimumGapDuration, total - sgm);
        }
        else if (mpCell->GetMutationState()->IsType<TACellMutationState>())
        {
            // TA cell cycle duration (scaled)
            double total = mTotalCycleMin
                         + p_gen->ranf() * (mTotalCycleMax - mTotalCycleMin);
            total *= mTransitCycleRatio;
            mG1Duration = std::max(mMinimumGapDuration, total - sgm);
        }
        else
        {
            // Fallback for unknown mutation state: use stem cycle
            double total = mTotalCycleMin
                         + p_gen->ranf() * (mTotalCycleMax - mTotalCycleMin);
            mG1Duration = std::max(mMinimumGapDuration, total - sgm);
        }
    }

    /** Copy constructor (used by CreateCellCycleModel). */
    StochasticFourTypeCellCycleModel(
        const StochasticFourTypeCellCycleModel& rModel)
        : ContactInhibitionCellCycleModel(rModel),
          mTotalCycleMin(rModel.mTotalCycleMin),
          mTotalCycleMax(rModel.mTotalCycleMax),
          mTransitCycleRatio(rModel.mTransitCycleRatio),
          mProbStemToStem(rModel.mProbStemToStem),
          mProbStemToPaneth(rModel.mProbStemToPaneth),
          mProbTaToTaEarly(rModel.mProbTaToTaEarly),
          mProbTaToTaLate(rModel.mProbTaToTaLate),
          mTransitionTime(rModel.mTransitionTime)
    {
    }

public:
    StochasticFourTypeCellCycleModel()
        : ContactInhibitionCellCycleModel(),
          mTotalCycleMin(12.0),
          mTotalCycleMax(14.0),
          mTransitCycleRatio(1.0),
          mProbStemToStem(0.89),
          mProbStemToPaneth(0.09),
          mProbTaToTaEarly(0.9),
          mProbTaToTaLate(0.7),
          mTransitionTime(120.0)   // day 5 = 120 hours
    {
        // Collapse S and G2 so G1 ≈ totalCycle − M
        SetSDuration(0.0);
        SetG2Duration(0.0);
        SetMDuration(1.0);
    }

    /**
     * Override InitialiseDaughterCell to implement stochastic fate decisions.
     * Called on the DAUGHTER cell after division.
     *
     * The parent cell keeps its current type. The daughter's type is
     * determined stochastically based on the parent's mutation state
     * and the current simulation time.
     */
    void InitialiseDaughterCell() override
    {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        double u = p_gen->ranf();
        double sim_time = SimulationTime::Instance()->GetTime();

        // ---- Stem cell (WildType) divides ----
        if (mpCell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            // All daughters start as TransitCellProliferativeType
            boost::shared_ptr<AbstractCellProperty> p_transit =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_transit);

            if (u < mProbStemToPaneth)
            {
                // Daughter becomes Paneth cell (non-proliferative)
                boost::shared_ptr<AbstractCellProperty> p_paneth =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
                mpCell->SetMutationState(p_paneth);

                // PC doesn't proliferate
                boost::shared_ptr<AbstractCellProperty> p_diff =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
                mpCell->SetCellProliferativeType(p_diff);
            }
            else if (u < mProbStemToPaneth + mProbStemToStem)
            {
                // Daughter remains stem cell
                boost::shared_ptr<AbstractCellProperty> p_wt =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<WildTypeCellMutationState>();
                mpCell->SetMutationState(p_wt);
            }
            else
            {
                // Daughter becomes TA cell
                boost::shared_ptr<AbstractCellProperty> p_ta =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
                mpCell->SetMutationState(p_ta);
            }
        }
        // ---- TA cell divides ----
        else if (mpCell->GetMutationState()->IsType<TACellMutationState>())
        {
            boost::shared_ptr<AbstractCellProperty> p_transit =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_transit);

            // Time-dependent transition probabilities
            double p_ta_to_ta = (sim_time < mTransitionTime) ? mProbTaToTaEarly : mProbTaToTaLate;

            if (u < p_ta_to_ta)
            {
                // Daughter remains TA
                boost::shared_ptr<AbstractCellProperty> p_ta =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TACellMutationState>();
                mpCell->SetMutationState(p_ta);
            }
            else
            {
                // Daughter becomes Enterocyte (non-proliferative)
                boost::shared_ptr<AbstractCellProperty> p_ec =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EnterocyteCellMutationState>();
                mpCell->SetMutationState(p_ec);

                boost::shared_ptr<AbstractCellProperty> p_diff =
                    mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
                mpCell->SetCellProliferativeType(p_diff);
            }
        }
        // ---- PC divides (should not happen, but handle gracefully) ----
        else if (mpCell->GetMutationState()->IsType<PanethCellMutationState>())
        {
            boost::shared_ptr<AbstractCellProperty> p_paneth =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<PanethCellMutationState>();
            mpCell->SetMutationState(p_paneth);

            boost::shared_ptr<AbstractCellProperty> p_diff =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff);
        }
        // ---- EC divides (should not happen, but handle gracefully) ----
        else if (mpCell->GetMutationState()->IsType<EnterocyteCellMutationState>())
        {
            boost::shared_ptr<AbstractCellProperty> p_ec =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<EnterocyteCellMutationState>();
            mpCell->SetMutationState(p_ec);

            boost::shared_ptr<AbstractCellProperty> p_diff =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff);
        }

        // Call parent to set G1 duration, birth time, etc.
        // This will call SetG1Duration() which checks type and sets ∞ for PC/EC
        ContactInhibitionCellCycleModel::InitialiseDaughterCell();
    }

    AbstractCellCycleModel* CreateCellCycleModel() override
    {
        return new StochasticFourTypeCellCycleModel(*this);
    }

    // ---- Setters / Getters -----------------------------------------

    void SetTotalCycleMin(double val)  { mTotalCycleMin = val; }
    void SetTotalCycleMax(double val)  { mTotalCycleMax = val; }
    void SetTransitCycleRatio(double r) { mTransitCycleRatio = r; }

    double GetTotalCycleMin()  const { return mTotalCycleMin; }
    double GetTotalCycleMax()  const { return mTotalCycleMax; }
    double GetTransitCycleRatio() const { return mTransitCycleRatio; }

    void SetProbStemToStem(double p) { mProbStemToStem = p; }
    void SetProbStemToPaneth(double p) { mProbStemToPaneth = p; }
    void SetProbTaToTaEarly(double p) { mProbTaToTaEarly = p; }
    void SetProbTaToTaLate(double p) { mProbTaToTaLate = p; }
    void SetTransitionTime(double t) { mTransitionTime = t; }

    double GetProbStemToStem() const { return mProbStemToStem; }
    double GetProbStemToPaneth() const { return mProbStemToPaneth; }
    double GetProbTaToTaEarly() const { return mProbTaToTaEarly; }
    double GetProbTaToTaLate() const { return mProbTaToTaLate; }
    double GetTransitionTime() const { return mTransitionTime; }

    void OutputCellCycleModelParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<TotalCycleMin>"
                     << mTotalCycleMin << "</TotalCycleMin>\n";
        *rParamsFile << "\t\t\t<TotalCycleMax>"
                     << mTotalCycleMax << "</TotalCycleMax>\n";
        *rParamsFile << "\t\t\t<TransitCycleRatio>"
                     << mTransitCycleRatio << "</TransitCycleRatio>\n";
        *rParamsFile << "\t\t\t<ProbStemToStem>"
                     << mProbStemToStem << "</ProbStemToStem>\n";
        *rParamsFile << "\t\t\t<ProbStemToPaneth>"
                     << mProbStemToPaneth << "</ProbStemToPaneth>\n";
        *rParamsFile << "\t\t\t<ProbTaToTaEarly>"
                     << mProbTaToTaEarly << "</ProbTaToTaEarly>\n";
        *rParamsFile << "\t\t\t<ProbTaToTaLate>"
                     << mProbTaToTaLate << "</ProbTaToTaLate>\n";
        *rParamsFile << "\t\t\t<TransitionTime>"
                     << mTransitionTime << "</TransitionTime>\n";
        ContactInhibitionCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(StochasticFourTypeCellCycleModel)

#endif // STOCHASTICFOURTYPECELLCYCLEMODEL_HPP_
