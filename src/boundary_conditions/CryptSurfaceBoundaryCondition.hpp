#ifndef CRYPTSURFACEBOUNDARYCONDITION_HPP_
#define CRYPTSURFACEBOUNDARYCONDITION_HPP_

#include "AbstractCellPopulationBoundaryCondition.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

/**
* A cell population boundary condition class, which restricts nodes to lie
* on the surface of a deformed test tube in the domain.
*/
template<unsigned DIM>
class CryptSurfaceBoundaryCondition : public AbstractCellPopulationBoundaryCondition<DIM>
{
private:

    /** The maximum distance from the surface that cells may be. */
    double mMaximumDistance;

    /* The target cell population number to shift the crypt length */
    double mTargetPopulation;

    /* The maximal crypt length, which is static */
    double mMaximumHeight;

    /* The minimum crypt height, which evolves according to the current crypt population */
    double mMinimumHeight;

    /* Rate at which crypt length remodels due to population size change */
    double mRemodellingRate;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationBoundaryCondition<DIM> >(*this);
        archive & mMaximumDistance;
        archive & mTargetPopulation;
        archive & mMaximumHeight;
        archive & mMinimumHeight;
        archive & mRemodellingRate;
    }

public:

    /**
     * Constructor.
     *
     * @param pCellPopulation pointer to the cell population
     * @param distance the maximum distance from the surface that cells may be (defaults to 1e-5)
     */
    CryptSurfaceBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation, double distance=1e-5);

    /**
     * Get the crypt target population
     */
    double GetTargetPopulation();

    /* 
     * Set the crypt target population (to model evolving crypt length)
     */
    void SetTargetPopulation(double TargetPopulation);

    /* 
     * Get the length remodelling rate
     */
    double GetRemodellingRate();

    /* 
     * Set the crypt length remodelling rate
     */
    void SetRemodellingRate(double RemodellingRate);

    /* 
     * Get the crypt maximum height
    */
    double GetMaximumCryptHeight();

    /*
      * Set the crypt maximum height
     */
    void SetMaximumCryptHeight(double maximumHeight);

    /* 
     * Update the minimum epithelium height
     */
    double UpdateMinimumCryptHeight();

    /**
     * Overridden ImposeBoundaryCondition() method.
     *
     * Apply the cell population boundary conditions.
     *
     * @param rOldLocations the node locations before any boundary conditions are applied
     */
    void ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations);

    /**
     * Overridden VerifyBoundaryCondition() method.
     * Verify the boundary conditions have been applied.
     * This is called after ImposeBoundaryCondition() to ensure the condition is still satisfied.
     *
     * @return whether the boundary conditions are satisfied.
     */
    bool VerifyBoundaryCondition();

    /**
     * Overridden OutputCellPopulationBoundaryConditionParameters() method.
     * Output cell population boundary condition parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptSurfaceBoundaryCondition)

namespace boost
{
namespace serialization
{
/**
* Serialize information required to construct a CryptSurfaceBoundaryCondition.
*/
template<class Archive, unsigned DIM>
inline void save_construct_data(Archive & ar, const CryptSurfaceBoundaryCondition<DIM>* t, unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
}

/**
 * De-serialize constructor parameters and initialize a CryptSurfaceBoundaryCondition.
*/
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, CryptSurfaceBoundaryCondition<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)CryptSurfaceBoundaryCondition<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*CryptSurfaceBoundaryCondition_HPP_*/    