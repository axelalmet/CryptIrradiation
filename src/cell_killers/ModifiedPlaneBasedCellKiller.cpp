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

#include "ModifiedPlaneBasedCellKiller.hpp"
#include "SimulationTime.hpp"

template<unsigned DIM>
ModifiedPlaneBasedCellKiller<DIM>::ModifiedPlaneBasedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation,
                                                  c_vector<double, DIM> point,
                                                  c_vector<double, DIM> normal, 
                                                  std::string outputDirectory)
    : AbstractCellKiller<DIM>(pCellPopulation),
      mPointOnPlane(point),
      mOutputDirectory(outputDirectory)
{
    assert(norm_2(normal) > 0.0);
    mNormalToPlane = normal/norm_2(normal);

    // Sets up output file
	OutputFileHandler output_file_handler(mOutputDirectory + "/", false);
	mpSloughingDataFile = output_file_handler.OpenOutputFile("sloughingdata.dat");

}

template<unsigned DIM>
const c_vector<double, DIM>& ModifiedPlaneBasedCellKiller<DIM>::rGetPointOnPlane() const
{
    return mPointOnPlane;
}

template<unsigned DIM>
const c_vector<double, DIM>& ModifiedPlaneBasedCellKiller<DIM>::rGetNormalToPlane() const
{
    return mNormalToPlane;
}

template<unsigned DIM>
const std::string ModifiedPlaneBasedCellKiller<DIM>::rGetOutputDirectory() const
{
    return mOutputDirectory;
}

template<unsigned DIM>
void ModifiedPlaneBasedCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    // Get the current time
    double current_time = SimulationTime::Instance()->GetTime();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
        c_vector<double, DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

        if (inner_prod(cell_location - mPointOnPlane, mNormalToPlane) > 0.0)
        {
            // Record the current time, cell index, and position
            unsigned cell_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);

            *mpSloughingDataFile << current_time << "\t";

            *mpSloughingDataFile << cell_index << "\t";

            // Record the cell position
            for (unsigned i = 0; i != DIM-1U; i++)
            {
                *mpSloughingDataFile << cell_location[i] << "\t";
            }

             *mpSloughingDataFile << cell_location[DIM - 1] << "\t";

            // Line break
            *mpSloughingDataFile << "\n";

            cell_iter->Kill();
        }
    }
}

template<unsigned DIM>
void ModifiedPlaneBasedCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<PointOnPlane>";
    for (unsigned index=0; index != DIM-1U; index++) //Note: inequality avoids testing index < 0U when DIM=1
    {
        *rParamsFile << mPointOnPlane[index] << ",";
    }
    *rParamsFile << mPointOnPlane[DIM-1] << "</PointOnPlane>\n";

    *rParamsFile << "\t\t\t<NormalToPlane>";
     for (unsigned index=0; index != DIM-1U; index++) //Note: inequality avoids testing index < 0U when DIM=1
     {
         *rParamsFile << mNormalToPlane[index] << ",";
     }
     *rParamsFile << mNormalToPlane[DIM-1] << "</NormalToPlane>\n";

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class ModifiedPlaneBasedCellKiller<1>;
template class ModifiedPlaneBasedCellKiller<2>;
template class ModifiedPlaneBasedCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ModifiedPlaneBasedCellKiller)
