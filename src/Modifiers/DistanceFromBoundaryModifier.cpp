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

#include "DistanceFromBoundaryModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MutableVertexMesh.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

template<unsigned DIM>
DistanceFromBoundaryModifier<DIM>::DistanceFromBoundaryModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
    mDepth(INT_UNSET)
{
}

template<unsigned DIM>
DistanceFromBoundaryModifier<DIM>::~DistanceFromBoundaryModifier()
{
}

template<unsigned DIM>
void DistanceFromBoundaryModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DistanceFromBoundaryModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void DistanceFromBoundaryModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    MutableVertexMesh<DIM,DIM>& r_mesh = static_cast<VertexBasedCellPopulation<DIM>* >(&rCellPopulation)->rGetMesh();

    // Iterate over elements
    for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = r_mesh.GetElementIteratorBegin();
         elem_iter != r_mesh.GetElementIteratorEnd();
         ++elem_iter)
    {
        unsigned elem_index = elem_iter->GetIndex();
        CellPtr cell_iter = rCellPopulation.GetCellUsingLocationIndex(elem_index);
        c_vector<double, 2> cell_location = rCellPopulation.GetLocationOfCellCentre(cell_iter);
        double y = cell_location[1];
        if ( (elem_iter->IsElementOnBoundary()) && (y > 1) )
        {
            // If element is on the boundary, distance to the boundary is zero
            rCellPopulation.GetCellUsingLocationIndex(elem_index)->GetCellData()->SetItem("boundary distance", 0);
        }
        else
        {
           // Else set the boundary distance to -1 (will be overridden below)
            rCellPopulation.GetCellUsingLocationIndex(elem_index)->GetCellData()->SetItem("boundary distance", -1);
        }
    }

    // todo: add a check for whether all cells have been labelled, regardless of depth

    for (unsigned bound_dist = 1; bound_dist <= mDepth; bound_dist ++)
    {
        // Find cells one step closer to boundary.
        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
            cell_iter != rCellPopulation.End();
            ++cell_iter)
        {
            unsigned cell_id = cell_iter->GetCellId();
            int dist_label = cell_iter->GetCellData()->GetItem("boundary distance");
            int bdi = bound_dist - 1;
            
            if ( dist_label == bdi )
            {
                // Find their neighbours
                std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);
                if (!neighbour_indices.empty())
                {
                    for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
                        iter != neighbour_indices.end();
                        ++iter)
                    {
                        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);
                        int nbr_dist_label = p_cell->GetCellData()->GetItem("boundary distance");
                        if (nbr_dist_label <  0)
                        {
                            // If the neighbour is unlabelled, assign boundary distance to label
                            p_cell->GetCellData()->SetItem("boundary distance", bound_dist);
                        }
                    }
                }
                
            }
        }
        
    }
}

template<unsigned DIM>
unsigned DistanceFromBoundaryModifier<DIM>::GetDepth()
{
    return mDepth;
}

template<unsigned DIM>
void DistanceFromBoundaryModifier<DIM>::SetDepth(unsigned depth)
{
    assert(depth >= 0.0);
    mDepth = depth;
}

template<unsigned DIM>
void DistanceFromBoundaryModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DistanceFromBoundaryModifier<1>;
template class DistanceFromBoundaryModifier<2>;
template class DistanceFromBoundaryModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DistanceFromBoundaryModifier)

