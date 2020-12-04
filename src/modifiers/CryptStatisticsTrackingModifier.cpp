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

#include "CryptStatisticsTrackingModifier.hpp"
#include <algorithm>
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "Debug.hpp"

template<unsigned DIM>
CryptStatisticsTrackingModifier<DIM>::CryptStatisticsTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mCryptTop(DOUBLE_UNSET)
{
}

template<unsigned DIM>
CryptStatisticsTrackingModifier<DIM>::~CryptStatisticsTrackingModifier()
{
}

template<unsigned DIM>
double CryptStatisticsTrackingModifier<DIM>::GetCryptTop()
{
    return mCryptTop;
}

template<unsigned DIM>
void CryptStatisticsTrackingModifier<DIM>::SetCryptTop(double cryptTop)
{
    mCryptTop = cryptTop;
}

template<unsigned DIM>
void CryptStatisticsTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CryptStatisticsTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

    // Set up the output file
    OutputFileHandler output_file_handler(outputDirectory + "/", false);
    mpCryptStatisticsFile = output_file_handler.OpenOutputFile("cryptstatistics.dat");

    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CryptStatisticsTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Initialise all the relevant data
    unsigned num_stem_cells = 0; // Number of stem cells
    unsigned num_transit_cells = 0; // Number of TA cells
    unsigned num_differentiated_cells = 0; // Number of differentiated cells 

    double min_height = DBL_MAX; // Initialise the candidate to track the minimum height

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the cell type
        boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = cell_iter->GetCellProliferativeType();

        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter); // Get node index
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index); // Get node
        double cell_height = p_node->rGetLocation()[DIM - 1]; // Get the cell location

        // Update the minimum height, if need be
        if (cell_height < min_height) 
        {
            min_height = cell_height;
        }

        if (p_cell_type->IsType<StemCellProliferativeType>()) // If the cell is a stem cell
        {
            num_stem_cells += 1;
        }
        else if (p_cell_type->IsType<TransitCellProliferativeType>()) // If cell is a TA cell
        {
            num_transit_cells += 1;
        }
        else // Should (!) be a differentiated cell
        {
            num_differentiated_cells += 1;
        }

    }

    // Update the data now
    double current_time = SimulationTime::Instance()->GetTime();

    // Write the data to file
    *mpCryptStatisticsFile << current_time << " ";
    *mpCryptStatisticsFile << mCryptTop << " "; // Maximum height
    *mpCryptStatisticsFile << min_height << " "; // Minimum height
    *mpCryptStatisticsFile << num_stem_cells << " "; // Number of stem cells
    *mpCryptStatisticsFile << num_transit_cells << " "; // Number of transit cells
    *mpCryptStatisticsFile << num_differentiated_cells << " "; // Number of differentiated cells
    *mpCryptStatisticsFile << "\n"; // New line

}

template<unsigned DIM>
void CryptStatisticsTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CryptStatisticsTrackingModifier<1>;
template class CryptStatisticsTrackingModifier<2>;
template class CryptStatisticsTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptStatisticsTrackingModifier)
