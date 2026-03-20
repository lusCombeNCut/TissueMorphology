/*

Copyright (c) 2005-2025, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef DIFFERENTIALADHESIONFORCE_HPP_
#define DIFFERENTIALADHESIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "GeneralisedLinearSpringForce.hpp"
#include "AbstractCellPopulation.hpp"

/**
 * A force class implementing differential adhesion between cell types.
 *
 * Extends GeneralisedLinearSpringForce to include cell-type specific
 * adhesion strengths. Different cell types (e.g., apical vs basal)
 * have different adhesion energies, driving cell sorting and layer formation.
 *
 * Adhesion matrix:
 * - Apical-Apical: strong
 * - Basal-Basal: strong  
 * - Apical-Basal: weak
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class DifferentialAdhesionForce : public GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>
{
    friend class boost::serialization::access;

private:

    /** Adhesion strength for apical-apical interactions */
    double mApicalApicalAdhesion;

    /** Adhesion strength for basal-basal interactions */
    double mBasalBasalAdhesion;

    /** Adhesion strength for apical-basal interactions (heterotypic) */
    double mApicalBasalAdhesion;

    /**
     * Boost Serialization method for archiving/checkpointing.
     *
     * @param archive  the archive
     * @param version  the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM> >(*this);
        archive & mApicalApicalAdhesion;
        archive & mBasalBasalAdhesion;
        archive & mApicalBasalAdhesion;
    }

public:

    /**
     * Constructor.
     */
    DifferentialAdhesionForce()
        : GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>(),
          mApicalApicalAdhesion(1.0),
          mBasalBasalAdhesion(1.0),
          mApicalBasalAdhesion(0.5)
    {
    }

    /**
     * Destructor.
     */
    virtual ~DifferentialAdhesionForce()
    {
    }

    /**
     * Get cell type from cell data.
     *
     * @param pCell pointer to the cell
     * @return cell type: 0 = basal, 1 = apical, -1 = undefined
     */
    int GetCellType(CellPtr pCell)
    {
        try
        {
            double is_apical = pCell->GetCellData()->GetItem("is_apical");
            return (is_apical > 0.5) ? 1 : 0;  // 1 = apical, 0 = basal
        }
        catch (Exception&)
        {
            return -1;  // Undefined
        }
    }

    /**
     * Overridden multiplication factor for spring constant.
     * 
     * Returns different spring constants based on cell-type pairing,
     * implementing differential adhesion.
     *
     * @param nodeAGlobalIndex index of one node
     * @param nodeBGlobalIndex index of the other node
     * @param rCellPopulation the cell population
     * @param isCloserThanRestLength whether the nodes are closer than rest length
     *
     * @return the multiplication factor
     */
    double VariableSpringConstantMultiplicationFactor(unsigned nodeAGlobalIndex,
                                                       unsigned nodeBGlobalIndex,
                                                       AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation,
                                                       bool isCloserThanRestLength)
    {
        // Get cells associated with nodes
        CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
        CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

        // Get cell types
        int type_A = GetCellType(p_cell_A);
        int type_B = GetCellType(p_cell_B);

        // If types undefined, use default spring constant
        if (type_A < 0 || type_B < 0)
        {
            return 1.0;
        }

        // Determine adhesion strength based on cell type pairing
        double adhesion_multiplier = 1.0;

        if (type_A == 1 && type_B == 1)
        {
            // Apical-Apical
            adhesion_multiplier = mApicalApicalAdhesion;
        }
        else if (type_A == 0 && type_B == 0)
        {
            // Basal-Basal
            adhesion_multiplier = mBasalBasalAdhesion;
        }
        else
        {
            // Apical-Basal (heterotypic)
            adhesion_multiplier = mApicalBasalAdhesion;
        }

        return adhesion_multiplier;
    }

    /**
     * Set apical-apical adhesion strength.
     *
     * @param adhesion adhesion multiplier
     */
    void SetApicalApicalAdhesion(double adhesion)
    {
        assert(adhesion >= 0.0);
        mApicalApicalAdhesion = adhesion;
    }

    /**
     * Get apical-apical adhesion strength.
     *
     * @return adhesion strength
     */
    double GetApicalApicalAdhesion() const
    {
        return mApicalApicalAdhesion;
    }

    /**
     * Set basal-basal adhesion strength.
     *
     * @param adhesion adhesion multiplier
     */
    void SetBasalBasalAdhesion(double adhesion)
    {
        assert(adhesion >= 0.0);
        mBasalBasalAdhesion = adhesion;
    }

    /**
     * Get basal-basal adhesion strength.
     *
     * @return adhesion strength
     */
    double GetBasalBasalAdhesion() const
    {
        return mBasalBasalAdhesion;
    }

    /**
     * Set apical-basal (heterotypic) adhesion strength.
     *
     * @param adhesion adhesion multiplier (typically < homotypic)
     */
    void SetApicalBasalAdhesion(double adhesion)
    {
        assert(adhesion >= 0.0);
        mApicalBasalAdhesion = adhesion;
    }

    /**
     * Get apical-basal adhesion strength.
     *
     * @return adhesion strength
     */
    double GetApicalBasalAdhesion() const
    {
        return mApicalBasalAdhesion;
    }

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<ApicalApicalAdhesion>" << mApicalApicalAdhesion << "</ApicalApicalAdhesion>\n";
        *rParamsFile << "\t\t\t<BasalBasalAdhesion>" << mBasalBasalAdhesion << "</BasalBasalAdhesion>\n";
        *rParamsFile << "\t\t\t<ApicalBasalAdhesion>" << mApicalBasalAdhesion << "</ApicalBasalAdhesion>\n";

        // Call method on direct parent class
        GeneralisedLinearSpringForce<ELEMENT_DIM, SPACE_DIM>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DifferentialAdhesionForce)

#endif /*DIFFERENTIALADHESIONFORCE_HPP_*/
