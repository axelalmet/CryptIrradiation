#ifndef RANDOMPROLIFERATIVECELLTARGETEDCELLKILLER_HPP_
#define RANDOMPROLIFERATIVECELLTARGETEDCELLKILLER_HPP_
#include "CheckpointArchiveTypes.hpp"

#include "ChasteSerialization.hpp"
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"

/*
 * Cell killer to randomly remove differentiated cells at the crypt edge,
 * as modelled in Dunn et al. (2012)
 */
template<unsigned DIM>
class RandomProliferativeCellTargetedCellKiller : public AbstractCellKiller<DIM>
{
private:

    // Probability of cell at the edges being killed
    double mApoptosisProbability;

    friend class boost::serialization::access;
    
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
        archive & mApoptosisProbability;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
	RandomProliferativeCellTargetedCellKiller(AbstractCellPopulation<DIM>* pCellPopulation);

    /* 
     * Get the probability that a cell at the edges will be removed
     */
    double GetApoptosisProbability();

    /*
     * Set the probability of apoptosis
     */
    void SetApoptosisProbability(double apoptosisProbability);

    /* Get the crypt height extremes, where
     * [0] - height of crypt base
     * [1] - height at crypt top
     */
    c_vector<double, 2> GetCryptHeightExtremes();

    /**
     *  Loops over and kills cells by anoikis or at the orifice if instructed.
     */
    void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Outputs cell killer parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomProliferativeCellTargetedCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive, unsigned DIM>
        inline void save_construct_data(
            Archive & ar, const RandomProliferativeCellTargetedCellKiller<DIM> * t, const unsigned int file_version)
        {
            const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive, unsigned DIM>
        inline void load_construct_data(
            Archive & ar, RandomProliferativeCellTargetedCellKiller<DIM> * t, const unsigned int file_version)
        {
            AbstractCellPopulation<DIM>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)RandomProliferativeCellTargetedCellKiller<DIM>(p_cell_population);
        }
    }
}

#endif /* RANDOMPROLIFERATIVECELLTARGETEDCELLKILLER_HPP_ */
