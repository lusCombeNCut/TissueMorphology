/*
 * UniformContactInhibitionGenerationalCellCycleModel.hpp
 *
 * Contact-inhibition cell cycle model with:
 *   1. Uniform-random total cycle duration (as in UniformContactInhibitionCellCycleModel)
 *   2. Generation-based Stem → TA → Differentiated cascade (Meineke et al. 2001)
 *
 * Generation logic (from AbstractSimpleGenerationalCellCycleModel):
 *   - Stem cells (generation = 0) divide into 1 Stem + 1 TA (gen = 1)
 *   - TA cells increment generation at each division
 *   - When generation > mMaxTransitGenerations, cell becomes Differentiated
 *   - Differentiated cells are permanently arrested in G0
 *
 * Contact inhibition (from ContactInhibitionCellCycleModel):
 *   - Cells in G1 are arrested (quiescent) when compressed below threshold volume
 *   - G1 duration extended while quiescent
 *
 * Reference:
 *   Meineke, F.A., Potten, C.S. & Loeffler, M. (2001). Cell migration and
 *   organization in the intestinal crypt using a lattice-free model.
 *   Cell Proliferation, 34(4), 253–266. doi:10.1046/j.0960-7722.2001.00216.x
 */
#ifndef UNIFORMCONTACTINHIBITIONGENERATIONALCELLCYCLEMODEL_HPP_
#define UNIFORMCONTACTINHIBITIONGENERATIONALCELLCYCLEMODEL_HPP_

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

class UniformContactInhibitionGenerationalCellCycleModel
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
        archive & mGeneration;
        archive & mMaxTransitGenerations;
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

    /** The generation of this cell (stem cells have generation 0). */
    unsigned mGeneration;

    /** Maximum number of TA divisions before terminal differentiation. */
    unsigned mMaxTransitGenerations;

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
    UniformContactInhibitionGenerationalCellCycleModel(
        const UniformContactInhibitionGenerationalCellCycleModel& rModel)
        : ContactInhibitionCellCycleModel(rModel),
          mTotalCycleMin(rModel.mTotalCycleMin),
          mTotalCycleMax(rModel.mTotalCycleMax),
          mTransitCycleRatio(rModel.mTransitCycleRatio),
          mGeneration(rModel.mGeneration),
          mMaxTransitGenerations(rModel.mMaxTransitGenerations)
    {
    }

public:
    UniformContactInhibitionGenerationalCellCycleModel()
        : ContactInhibitionCellCycleModel(),
          mTotalCycleMin(12.0),
          mTotalCycleMax(14.0),
          mTransitCycleRatio(1.0),
          mGeneration(0),
          mMaxTransitGenerations(3)  // Meineke et al. 2001
    {
        // Collapse S and G2 so G1 ≈ totalCycle − M
        SetSDuration(0.0);
        SetG2Duration(0.0);
        SetMDuration(1.0);
    }

    /**
     * Override ResetForDivision to handle generation counting.
     * Called on the PARENT cell just before division.
     */
    void ResetForDivision() override
    {
        mGeneration++;

        // Stem cells reset their generation to 0 (they remain stem)
        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mGeneration = 0;
        }

        // If generation exceeds max, this cell becomes differentiated
        if (mGeneration > mMaxTransitGenerations)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
        }

        // Call parent to reset mReadyToDivide, mCurrentCellCyclePhase, etc.
        ContactInhibitionCellCycleModel::ResetForDivision();
    }

    /**
     * Override InitialiseDaughterCell to handle generation-based differentiation.
     * Called on the DAUGHTER cell after division.
     *
     * Key logic (Meineke model):
     *   - If parent was stem (gen=0), daughter gets gen=1 and becomes TA
     *   - Otherwise daughter inherits parent's generation (already incremented)
     *   - If generation > max, daughter becomes differentiated
     */
    void InitialiseDaughterCell() override
    {
        // If parent was stem cell, its generation was reset to 0.
        // Daughter must be TA with generation 1.
        if (mGeneration == 0)
        {
            mGeneration = 1;
        }

        // Daughter cells of stem cells become transit-amplifying
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);

        // Check if daughter should be differentiated
        if (mGeneration > mMaxTransitGenerations)
        {
            boost::shared_ptr<AbstractCellProperty> p_diff_type =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
            mpCell->SetCellProliferativeType(p_diff_type);
        }

        // Call parent to set G1 duration, birth time, etc.
        ContactInhibitionCellCycleModel::InitialiseDaughterCell();
    }

    AbstractCellCycleModel* CreateCellCycleModel() override
    {
        return new UniformContactInhibitionGenerationalCellCycleModel(*this);
    }

    // ---- Setters / Getters -----------------------------------------

    void SetTotalCycleMin(double val)  { mTotalCycleMin = val; }
    void SetTotalCycleMax(double val)  { mTotalCycleMax = val; }
    void SetTransitCycleRatio(double r) { mTransitCycleRatio = r; }

    double GetTotalCycleMin()  const { return mTotalCycleMin; }
    double GetTotalCycleMax()  const { return mTotalCycleMax; }
    double GetTransitCycleRatio() const { return mTransitCycleRatio; }

    void SetGeneration(unsigned generation) { mGeneration = generation; }
    unsigned GetGeneration() const { return mGeneration; }

    void SetMaxTransitGenerations(unsigned maxTransitGenerations)
    {
        mMaxTransitGenerations = maxTransitGenerations;
    }
    unsigned GetMaxTransitGenerations() const { return mMaxTransitGenerations; }

    void OutputCellCycleModelParameters(out_stream& rParamsFile) override
    {
        *rParamsFile << "\t\t\t<TotalCycleMin>"
                     << mTotalCycleMin << "</TotalCycleMin>\n";
        *rParamsFile << "\t\t\t<TotalCycleMax>"
                     << mTotalCycleMax << "</TotalCycleMax>\n";
        *rParamsFile << "\t\t\t<TransitCycleRatio>"
                     << mTransitCycleRatio << "</TransitCycleRatio>\n";
        *rParamsFile << "\t\t\t<MaxTransitGenerations>"
                     << mMaxTransitGenerations << "</MaxTransitGenerations>\n";
        ContactInhibitionCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(UniformContactInhibitionGenerationalCellCycleModel)

#endif // UNIFORMCONTACTINHIBITIONGENERATIONALCELLCYCLEMODEL_HPP_
