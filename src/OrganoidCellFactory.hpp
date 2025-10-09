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

#ifndef ORGANOIDCELLFACTORY_HPP_
#define ORGANOIDCELLFACTORY_HPP_

#include "UniformG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellData.hpp"
#include "SmartPointers.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractCellMutationState.hpp"
#include "AbstractCellProliferativeType.hpp"

/**
 * A utility class to create cells for organoid formation simulations.
 * 
 * This class creates cells with varying properties relevant to 
 * basement membrane stiffness effects on organoid morphology.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class OrganoidCellFactory
{
private:
    
    /** The stiffness of the basement membrane */
    double mBasementMembraneStiffness;
    
    /** Whether to add basement membrane interactions */
    bool mUseBasementMembrane;

public:

    /**
     * Constructor.
     * 
     * @param basementMembraneStiffness The stiffness parameter for basement membrane
     * @param useBasementMembrane Whether to include basement membrane effects
     */
    OrganoidCellFactory(double basementMembraneStiffness = 1.0, 
                       bool useBasementMembrane = true)
        : mBasementMembraneStiffness(basementMembraneStiffness),
          mUseBasementMembrane(useBasementMembrane)
    {
    }

    /**
     * Create cells for all nodes in a mesh.
     * 
     * @param rCells Reference to vector to populate with cells
     * @param numNodes Number of nodes to create cells for
     */
    void CreateCells(std::vector<CellPtr>& rCells, unsigned numNodes)
    {
        rCells.clear();
        rCells.reserve(numNodes);
        
        for (unsigned i = 0; i < numNodes; i++)
        {
            CellPtr p_cell = CreateCell(i);
            if (p_cell)  // Check that cell was created successfully
            {
                rCells.push_back(p_cell);
            }
        }
    }

    /**
     * Create a cell with the appropriate organoid properties.
     * 
     * @param nodeIndex Index of the node
     * @return a pointer to the created cell
     */
    CellPtr CreateCell(unsigned nodeIndex)
    {
        // Create a simple cell cycle model
        UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
        p_cycle_model->SetDimension(SPACE_DIM);
        
        // Create shared pointers for mutation state and proliferative type
        boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
        boost::shared_ptr<AbstractCellProliferativeType> p_type(new TransitCellProliferativeType);
        
        // Create the cell with mutation state and proliferative type
        CellPtr p_cell(new Cell(p_state, p_cycle_model));
        p_cell->SetCellProliferativeType(p_type);
        
        // Set birth time to random value to avoid synchrony
        double birth_time = -RandomNumberGenerator::Instance()->ranf() * 12.0;
        p_cell->SetBirthTime(birth_time);
        
        // Add basement membrane stiffness as cell data
        if (mUseBasementMembrane)
        {
            p_cell->GetCellData()->SetItem("basement_membrane_stiffness", mBasementMembraneStiffness);
        }
        
        // Add simple node-index based properties 
        // (In real implementation, would use actual node locations)
        if (nodeIndex < 5) // First few cells are stem-like
        {
            p_cell->GetCellData()->SetItem("cell_type", 0.0); // stem-like
        }
        else
        {
            p_cell->GetCellData()->SetItem("cell_type", 1.0); // differentiated-like
        }
        
        return p_cell;
    }
    
    /**
     * @return the basement membrane stiffness parameter
     */
    double GetBasementMembraneStiffness()
    {
        return mBasementMembraneStiffness;
    }
    
    /**
     * Set the basement membrane stiffness parameter
     * @param stiffness the new stiffness value
     */
    void SetBasementMembraneStiffness(double stiffness)
    {
        mBasementMembraneStiffness = stiffness;
    }
};

#endif /*ORGANOIDCELLFACTORY_HPP_*/