/*

Copyright (c) 2005-2021, University of Oxford.
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

#include "FixedBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "AbstractCellProperty.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FixedBoundaryCondition<DIM>::FixedBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mFixTopCells(false)
{
}

template<unsigned DIM>
void FixedBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    // Make sure that we are in the correct dimension - this code will be eliminated at compile time
    assert(DIM == 2); // LCOV_EXCL_LINE // this method only works in 2D at present

    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation) == nullptr)
    {
        EXCEPTION("FixedBoundaryCondition is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    {
        // Iterate over all nodes to update their positions
        for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
        {
            Node<DIM>* p_node = p_cell_population->GetNode(node_index);
            double node_height = rOldLocations.find(p_node)->second[DIM-1];

            if ((p_node->IsBoundaryNode()))
            {
                if (mFixTopCells)
                {
                    // Set node back to previous height
                    p_node->rGetModifiableLocation()[DIM-1] = node_height;
                }
                else if (node_height < 1.5)
                {
                    // Set node back to previous height
                    p_node->rGetModifiableLocation()[DIM-1] = node_height;
                }
                
            }
        }
    }
}

template<unsigned DIM>
bool FixedBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    bool boundary_condition_satisfied = true;

    /*
     * Here we verify that the boundary condition is still satisfied by 
     * checking that cells on the bottom boundary have not moved in the
     * y-direction.
     */
    // VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(this->mpCellPopulation);
    // {
    //     // Iterate over all nodes to update their positions
    //     for (unsigned node_index=0; node_index<p_cell_population->GetNumNodes(); node_index++)
    //     {
    //         Node<DIM>* p_node = p_cell_population->GetNode(node_index);
    //         double old_height = rOldLocations.find(p_node)->second[DIM-1];

    //         if ((p_node->IsBoundaryNode()))
    //         {
    //             if (mFixTopCells)
    //             {
    //                 // Check that the height has not changed
    //                 if (p_node->rGetModifiableLocation()[DIM-1] != old_height)
    //                 {
    //                     boundary_condition_satisfied = false;
    //                     break;
    //                 }
    //             }
    //             else if (old_height < 1.5)
    //             {
    //                 // Check that the height has not changed
    //                 if (p_node->rGetModifiableLocation()[DIM-1] != old_height)
    //                 {
    //                     boundary_condition_satisfied = false;
    //                     break;
    //                 }
    //             }
                
    //         }
    //     }
    // }

    return boundary_condition_satisfied;
}

template<unsigned DIM>
void FixedBoundaryCondition<DIM>::SetFixTopCells(bool fixTopCells)
{
    mFixTopCells = fixTopCells;
}

template<unsigned DIM>
bool FixedBoundaryCondition<DIM>::GetFixTopCells()
{
    return mFixTopCells;
}

template<unsigned DIM>
void FixedBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<FixTopCells>" << mFixTopCells << "</FixTopCells>\n";
    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class FixedBoundaryCondition<1>;
template class FixedBoundaryCondition<2>;
template class FixedBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedBoundaryCondition)

