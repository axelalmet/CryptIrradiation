/*

Copyright (c) 2005-2020, University of Oxford.
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

#include "StemCellRetainerForce.hpp"
#include "StemCellProliferativeType.hpp"

template<unsigned DIM>
StemCellRetainerForce<DIM>::StemCellRetainerForce()
    : AbstractForce<DIM>(),
    mRetainerForceStrength(DOUBLE_UNSET)
{
}

template<unsigned DIM>
StemCellRetainerForce<DIM>::~StemCellRetainerForce()
{
}

// Method to get the strength of the retainer force
template<unsigned DIM>
double StemCellRetainerForce<DIM>::GetRetainerForceStrength()
{
    return mRetainerForceStrength;
}

// Method to set the strength of retainer force
template<unsigned DIM>
void StemCellRetainerForce<DIM>::SetRetainerForceStrength(double retainerForceStrength)
{
    mRetainerForceStrength = retainerForceStrength;
}

template<unsigned DIM>
void StemCellRetainerForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    // Get the force parameters
    double retainer_force_strength = GetRetainerForceStrength(); // Strength of migration force

    // Initialise the retainer force vector
    c_vector<double, DIM> retainer_force = zero_vector<double>(DIM);
    retainer_force[DIM - 1] = -1.0 * retainer_force_strength;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        // Only apply this migration force to fibroblasts
        boost::shared_ptr<AbstractCellProperty> p_cell_type = cell_iter->GetCellProliferativeType(); // Get the cell type

        if (p_cell_type->IsType<StemCellProliferativeType>())
        {
            unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
            rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(retainer_force); 
        }
    }
}

template<unsigned DIM>
void StemCellRetainerForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // No parameters to include
	*rParamsFile <<  "\t\t\t<RetainerForceStrength>"<<  mRetainerForceStrength << "</RetainerForceStrength> \n" ;

    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class StemCellRetainerForce<1>;
template class StemCellRetainerForce<2>;
template class StemCellRetainerForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StemCellRetainerForce)
