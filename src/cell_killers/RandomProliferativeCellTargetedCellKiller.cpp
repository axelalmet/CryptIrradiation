#include "RandomProliferativeCellTargetedCellKiller.hpp"
#include "AbstractCellKiller.hpp"
#include "AbstractCellProperty.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
RandomProliferativeCellTargetedCellKiller<DIM>::RandomProliferativeCellTargetedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation),
    mApoptosisProbability(DOUBLE_UNSET)
{
}

// Method to get the probability of apoptosis
template<unsigned DIM>
double RandomProliferativeCellTargetedCellKiller<DIM>::GetApoptosisProbability()
{
	return mApoptosisProbability;
}

// Method to set the probability of apoptosis
template<unsigned DIM>
void RandomProliferativeCellTargetedCellKiller<DIM>::SetApoptosisProbability(double apoptosisProbability)
{
	mApoptosisProbability = apoptosisProbability;
}

/*
 * Cell Killer that kills stem cells (in the base) only. The probability calculation's based on a 
 * Binomial distribution with probability p. So it's easy to calculate what the 
 * probability within the interval [0, dt) should be.
 */
template<unsigned DIM>
void RandomProliferativeCellTargetedCellKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{
    
    double apoptosis_probability = GetApoptosisProbability(); // Get probability of apoptosis

    double dt = SimulationTime::Instance()->GetTimeStep(); // Get dt

    double probability_of_dying = 1.0 - pow(1.0 - apoptosis_probability, dt); // We calculate the probability of dying now 

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = (this->mpCellPopulation)->Begin();
        cell_iter != (this->mpCellPopulation)->End();
        ++cell_iter)
    {
        // Target stem cells only
        if (cell_iter->GetCellProliferativeType()->template IsType<StemCellProliferativeType>())
        {
            double uniform_random_number = RandomNumberGenerator::Instance()->ranf();

            if ((uniform_random_number < probability_of_dying)&&(!cell_iter->HasApoptosisBegun()))
            {
                cell_iter->StartApoptosis();
            }
        }
    }

}

template<unsigned DIM>
void RandomProliferativeCellTargetedCellKiller<DIM>::OutputCellKillerParameters(out_stream& rParamsFile)
{
    // Call direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class RandomProliferativeCellTargetedCellKiller<1>;
template class RandomProliferativeCellTargetedCellKiller<2>;
template class RandomProliferativeCellTargetedCellKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomProliferativeCellTargetedCellKiller)