
#include "CryptSurfaceBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SimulationTime.hpp"
#include "Debug.hpp"
/* This is a boundary condition which is pulled from Dunn et al. (2016)., where the authors fit a polynomial
 *  curve to imaging data from Paul Appleton in Inke Nathke's lab.  We will generalise some parameters so that
 * we can adapt the boundary condition to a colonic (mouse) crypt, if we want to do so.
*/


template<unsigned DIM>
CryptSurfaceBoundaryCondition<DIM>::CryptSurfaceBoundaryCondition(AbstractCellPopulation<DIM>* pCellPopulation,
                                                                      double distance)
    : AbstractCellPopulationBoundaryCondition<DIM>(pCellPopulation),
      mMaximumDistance(distance),
      mTargetPopulation(UNSIGNED_UNSET),
      mMaximumHeight(DOUBLE_UNSET),
      mMinimumHeight(0.0),
      mRemodellingRate(0.0)
{
    assert(mMaximumHeight > 0.0);
    assert(mMaximumDistance >= 0.0);

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(this->mpCellPopulation) == NULL)
    {
       EXCEPTION("A NodeBasedCellPopulation must be used with this boundary condition object.");
    }
    if (DIM < 3)
    {
        EXCEPTION("This boundary condition is only implemented in 3D.");
    }
}

template<unsigned DIM>
double CryptSurfaceBoundaryCondition<DIM>::GetTargetPopulation()
{
    return mTargetPopulation;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::SetTargetPopulation(double TargetPopulation)
{
    mTargetPopulation = TargetPopulation;
}

template<unsigned DIM>
double CryptSurfaceBoundaryCondition<DIM>::GetRemodellingRate()
{
    return mRemodellingRate;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::SetRemodellingRate(double RemodellingRate)
{
    mRemodellingRate = RemodellingRate;
}

template<unsigned DIM>
double CryptSurfaceBoundaryCondition<DIM>::GetMaximumCryptHeight()
{
    return mMaximumHeight;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::SetMaximumCryptHeight(double maximumHeight)
{
    mMaximumHeight = maximumHeight;
}

template<unsigned DIM>
double CryptSurfaceBoundaryCondition<DIM>::UpdateMinimumCryptHeight()
{
    double minimum_height = mMinimumHeight; // Get the current minimum height
    double remodelling_rate = GetRemodellingRate();
    double dt = SimulationTime::Instance()->GetTimeStep(); // Get the time step

    double target_population = GetTargetPopulation(); // Get the target population
    double current_population = (double)(this->mpCellPopulation)->rGetMesh().GetNumNodes();
    
    minimum_height += dt * remodelling_rate * (target_population - current_population);

    return minimum_height;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{

    // Get the crypt's minimum and maximum heights
    double min_height = UpdateMinimumCryptHeight();
    double max_height = GetMaximumCryptHeight();
    double scale_factor = 70.0/max_height;

    // Set the radius of the crypt base appropriately (an ellipse)
    double major_radius = 16.6968;
    double minor_radius = 15.3973 / scale_factor;

    // polynomial coefficients: p = p_1x^5 + p_2x^4 + p_3x^3 + p_4x^2 + p_5x + p_6
    // we scale the polynomial coefficients with the adjusted crypt heights
    double p_1 = 8.365027421137118e-07 * scale_factor * scale_factor * scale_factor * scale_factor * scale_factor;
    double p_2 = -1.612884286709950e-04 * scale_factor * scale_factor * scale_factor * scale_factor;
    double p_3 = 1.169561186127581e-02 * scale_factor * scale_factor * scale_factor;
    double p_4 = -3.912727534949710e-01 * scale_factor * scale_factor;
    double p_5 = 5.850759485536169e+00 * scale_factor;
    double p_6 = -1.497897551314556e+01;

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
        cell_iter != this->mpCellPopulation->End();
        ++cell_iter)
    {
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM> ellipse_centre = zero_vector<double>(DIM);
        ellipse_centre[DIM-1] = minor_radius + min_height;   // the shorter of the two radii in the ellipse

        // Where the cell originally sat along the z axis (will be used to put it back to its original height later)
        double original_height = cell_location[DIM-1];

        // The target radius depends on each node, and is defined separately for those cells in the crypt base,
        // and along the crypt walls
        double target_radius = 0.0;

        // Those cells that sit at the crypt base
        if (cell_location[DIM-1] < minor_radius + min_height)
        {
            double target_x_coordinate = sqrt((pow(major_radius, 2.0)*pow(minor_radius * scale_factor, 2.0)*pow(cell_location[DIM-3], 2.0)) /
                (pow(minor_radius * scale_factor, 2.0)*pow(cell_location[DIM-3], 2.0) + pow(major_radius, 2.0)*pow(cell_location[DIM-2], 2.0)));
            double target_y_coordinate = sqrt((pow(major_radius, 2.0)*pow(minor_radius * scale_factor, 2.0)*pow(cell_location[DIM-2], 2.0)) /
            (pow(minor_radius * scale_factor, 2.0)*pow(cell_location[DIM-3], 2.0) + pow(major_radius, 2.0)*pow(cell_location[DIM-2], 2.0)));

            target_radius = sqrt(pow(target_x_coordinate, 2.0) + pow(target_y_coordinate, 2.0));
        }

        // For those cells which sit above the elliptic base of the crypt
        if (cell_location[DIM-1] >= minor_radius + min_height)
        {
            double z = (max_height) / (max_height - min_height) * cell_location[DIM-1] - min_height; // To avoid the base of the function where it goes negative

            // Get the target radius for this cell (which sits higher in the crypt)
            target_radius = p_1*z*z*z*z*z + p_2*z*z*z*z + p_3*z*z*z + p_4*z*z + p_5*z + p_6; // interpolated crypt shape

            cell_location[DIM-1] = minor_radius + min_height;
        }

        // Find the radial distance between this cell and the surface
        double radius = norm_2(cell_location - ellipse_centre);

        // If the cell is too far from the surface of the ellipsoid / upper crypt
        if (fabs(radius - target_radius) > mMaximumDistance)
        {
            // ...move the cell back onto the surface

            c_vector<double, DIM> location_on_surface = target_radius*(cell_location - ellipse_centre)/radius + ellipse_centre;

            if (original_height >= minor_radius + min_height)
            {
                assert(fabs(location_on_surface[DIM-1] - minor_radius - min_height)<1e-8);
                location_on_surface[DIM-1] = original_height;
            }

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            if ( location_on_surface[DIM-1] < min_height )
            {
                location_on_surface[DIM-1] = min_height;
            }

            p_node->rGetModifiableLocation() = location_on_surface;

        }
    }
}

template<unsigned DIM>
bool CryptSurfaceBoundaryCondition<DIM>::VerifyBoundaryCondition()
{
    return true;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaximumDistance>" << mMaximumDistance << "</MaximumDistance>\n";
    *rParamsFile << "\t\t\t<MaximumHeight>" << mMaximumHeight << "</MaximumHeight>\n";
    *rParamsFile << "\t\t\t<MinimumHeight>" << mMinimumHeight << "</MinimumHeight>\n";
    *rParamsFile << "\t\t\t<TargetPopulation>" << mTargetPopulation << "</TargetPopulation>\n";
    *rParamsFile << "\t\t\t<RemodellingRate>" << mRemodellingRate << "</RemodellingRate>\n";

   // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class CryptSurfaceBoundaryCondition<1>;
template class CryptSurfaceBoundaryCondition<2>;
template class CryptSurfaceBoundaryCondition<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CryptSurfaceBoundaryCondition)