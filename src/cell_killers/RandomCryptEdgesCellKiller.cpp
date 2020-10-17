#include "RandomCryptEdgesCellKiller.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StromalCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

RandomCryptEdgesCellKiller::RandomCryptEdgesCellKiller(AbstractCellPopulation<2>* pCellPopulation)
    : AbstractCellKiller<2>(pCellPopulation),
    mApoptosisProbability(DOUBLE_UNSET)
{
}

// Method to get the probability of apoptosis
double RandomCryptEdgesCellKiller::GetApoptosisProbability()
{
	return mApoptosisProbability;
}

// Method to set the probability of apoptosis
void RandomCryptEdgesCellKiller::SetApoptosisProbability(double apoptosisProbability)
{
	mApoptosisProbability = apoptosisProbability;
}

// Get the minimal and maximal heights of the crypt epithelium.
c_vector<double, 2> RandomCryptEdgesCellKiller::GetCryptHeightExtremes()
{
	c_vector<double, 2> crypt_height_extremes;

    double min_height = DBL_MAX;
    double max_height = 0.0;

    for (AbstractCellPopulation<2>::Iterator cell_iter = (this->mpCellPopulation)->Begin();
        cell_iter != (this->mpCellPopulation)->End();
        ++cell_iter)
    {
        // Epithelial cells are defined
        if (!cell_iter->GetCellProliferativeType()->IsType<StromalCellProliferativeType>())
        {
            // Get the height
            double y = (this->mpCellPopulation)->GetLocationOfCellCentre(*cell_iter)[1]; 

            if (y < min_height)
            {
                min_height = y;
            }
            else if (y > max_height)
            {
                max_height = y;
            }
            
        }
    }
        
    crypt_height_extremes[0] = min_height; 
    crypt_height_extremes[1] = max_height;

    return crypt_height_extremes;
}

/*
 * Cell Killer that kills differentiated (maybe change this?) epithelial cells
 * at the crypt collar randomly. The probability calculation's based on a 
 * Binomial distribution with probability p. So it's easy to calculate what the 
 * probability within the interval [0, dt) should be.
 */
void RandomCryptEdgesCellKiller::CheckAndLabelCellsForApoptosisOrDeath()
{
    
    double apoptosis_probability = GetApoptosisProbability(); // Get probability of apoptosis
    c_vector<double, 2> crypt_heights = GetCryptHeightExtremes(); // Get min and max heights of crypt epithelium
    double min_height = crypt_heights[0];
    double max_height = crypt_heights[1];

    double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt

    double probability_of_dying = 1.0 - pow(1.0 - apoptosis_probability, dt); // We calculate the probability of dying now 

    for (AbstractCellPopulation<2>::Iterator cell_iter = (this->mpCellPopulation)->Begin();
        cell_iter != (this->mpCellPopulation)->End();
        ++cell_iter)
    {
        // Epithelial cells are defined
        if (cell_iter->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            // Get the height
            double y = (this->mpCellPopulation)->GetLocationOfCellCentre(*cell_iter)[1]; 

            if (y > min_height + 0.9 * (max_height - min_height)) // If the cell is in the upper 10%
            {
                double uniform_random_number = RandomNumberGenerator::Instance()->ranf();

                if ((uniform_random_number < probability_of_dying)&&(!cell_iter->HasApoptosisBegun()))
                {
                    cell_iter->StartApoptosis();
                }
            }
        }
    }

}

void RandomCryptEdgesCellKiller::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    AbstractCellKiller<2>::OutputCellKillerParameters(rParamsFile);
}




#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(RandomCryptEdgesCellKiller)
