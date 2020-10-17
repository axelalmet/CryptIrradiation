#ifndef RANDOMCRYPTEDGESCELLKILLER_HPP_
#define RANDOMCRYPTEDGESCELLKILLER_HPP_
#include "CheckpointArchiveTypes.hpp"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

#include "AbstractCellKiller.hpp"

/*
 * Cell killer to randomly remove differentiated cells at the crypt edge,
 * as modelled in Dunn et al. (2012)
 */

class RandomCryptEdgesCellKiller : public AbstractCellKiller<2>
{
private:

    // Probability of cell at the edges being killed
    double mApoptosisProbability;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<2> >(*this);
        archive & mApoptosisProbability;
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to a tissue
     * @param sloughOrifice whether to slough compressed cells at crypt orifice
     */
	RandomCryptEdgesCellKiller(AbstractCellPopulation<2>* pCellPopulation);

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
CHASTE_CLASS_EXPORT(RandomCryptEdgesCellKiller)

namespace boost
{
    namespace serialization
    {
        template<class Archive>
        inline void save_construct_data(
            Archive & ar, const RandomCryptEdgesCellKiller * t, const unsigned int file_version)
        {
            const AbstractCellPopulation<2>* const p_cell_population = t->GetCellPopulation();
            ar << p_cell_population;
        }

        template<class Archive>
        inline void load_construct_data(
            Archive & ar, RandomCryptEdgesCellKiller * t, const unsigned int file_version)
        {
            AbstractCellPopulation<2>* p_cell_population;
            ar >> p_cell_population;

            // Invoke inplace constructor to initialise instance
            ::new(t)RandomCryptEdgesCellKiller(p_cell_population);
        }
    }
}

#endif /* RANDOMCRYPTEDGESCELLKILLER_HPP_ */
