
#include "CryptSurfaceBoundaryCondition.hpp"
#include "NodeBasedCellPopulation.hpp"
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
      mTargetPopulation(DOUBLE_UNSET),
      mLengthRemodellingRate(DOUBLE_UNSET)
{
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
void CryptSurfaceBoundaryCondition<DIM>::SetTargetPopulation(double targetPopulation)
{
    mTargetPopulation = targetPopulation;
}

template<unsigned DIM>
double CryptSurfaceBoundaryCondition<DIM>::GetLengthRemodellingRate()
{
    return mLengthRemodellingRate;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::SetLengthRemodellingRate(double lengthRemodellingRate)
{
    mLengthRemodellingRate = lengthRemodellingRate;
}

template<unsigned DIM>
void CryptSurfaceBoundaryCondition<DIM>::ImposeBoundaryCondition(const std::map<Node<DIM>*, c_vector<double, DIM> >& rOldLocations)
{
    // Set the radius of the crypt base appropriately (an ellipse)
    double major_radius = 16.6968;
    double minor_radius = 15.3973;

    // polynomial coefficients: p = p_1x^5 + p_2x^4 + p_3x^3 + p_4x^2 + p_5x + p_6
    double p_1 = 8.365027421137118e-07;
    double p_2 = -1.612884286709950e-04;
    double p_3 = 1.169561186127581e-02;
    double p_4 = -3.912727534949710e-01;
    double p_5 = 5.850759485536169e+00;
    double p_6 = -1.497897551314556e+01;

    // Get the crypt target population and length remodelling rate
    double target_population = GetTargetPopulation();
    double length_remodelling_rate = GetLengthRemodellingRate();

    // Iterate over the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
        cell_iter != this->mpCellPopulation->End();
        ++cell_iter)
    {
        c_vector<double,DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

        c_vector<double,DIM> ellipse_centre = zero_vector<double>(DIM);
        ellipse_centre[DIM-1] = minor_radius;   // the shorter of the two radii in the ellipse

        // Where the cell originally sat along the z axis (will be used to put it back to its original height later)
        double original_height = cell_location[DIM-1];

        // The target radius depends on each node, and is defined separately for those cells in the crypt base,
        // and along the crypt walls
        double target_radius = 0.0;

        // Those cells that sit at the crypt base
        if (cell_location[DIM-1] < minor_radius)
        {
            double target_x_coordinate = sqrt((pow(major_radius, 2.0)*pow(minor_radius, 2.0)*pow(cell_location[DIM-3], 2.0)) /
                (pow(minor_radius,2)*pow(cell_location[DIM-3], 2.0) + pow(major_radius, 2.0)*pow(cell_location[DIM-2], 2.0)));
            double target_y_coordinate = sqrt((pow(major_radius, 2.0)*pow(minor_radius, 2.0)*pow(cell_location[DIM-2], 2.0)) /
            (pow(minor_radius,2)*pow(cell_location[DIM-3], 2.0) + pow(major_radius, 2.0)*pow(cell_location[DIM-2], 2.0)));

            target_radius = sqrt(pow(target_x_coordinate, 2.0) + pow(target_y_coordinate, 2.0));
        }

        // For those cells which sit above the elliptic base of the crypt
        if (cell_location[DIM-1] >= minor_radius)
        {
            double z = cell_location[DIM-1]; // To avoid the base of the function where it goes negative

            // Only bother translating if we actually have a non-zero remodelling rate
            if (target_population > 0.0)
            {
                unsigned num_cells = (this->mpCellPopulation)->rGetMesh().GetNumNodes();

                // Translate the height function
                z -= length_remodelling_rate * (target_population - (double)num_cells);
            }

            // Get the target radius for this cell (which sits higher in the crypt)
            target_radius = p_1*z*z*z*z*z + p_2*z*z*z*z + p_3*z*z*z + p_4*z*z + p_5*z + p_6; // interpolated crypt shape

            cell_location[DIM-1] = minor_radius;
        }

        // Find the radial distance between this cell and the surface
        double radius = norm_2(cell_location - ellipse_centre);

        // If the cell is too far from the surface of the ellipsoid / upper crypt
        if (fabs(radius - target_radius) > mMaximumDistance)
        {
            // ...move the cell back onto the surface

            c_vector<double, DIM> location_on_surface = target_radius*(cell_location - ellipse_centre)/radius + ellipse_centre;

            if (original_height >= minor_radius)
            {
                assert(fabs(location_on_surface[DIM-1]-minor_radius)<1e-8);
                location_on_surface[DIM-1] = original_height;
            }

            unsigned node_index = this->mpCellPopulation->GetLocationIndexUsingCell(*cell_iter);
            Node<DIM>* p_node = this->mpCellPopulation->GetNode(node_index);

            if ( location_on_surface[DIM-1] < 0 )
            {
                location_on_surface[DIM-1] = 0.0;
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
    *rParamsFile << "\t\t\t<TargetPopulation>" << mTargetPopulation << "</TargetPopulation>\n";
    *rParamsFile << "\t\t\t<LengthRemodellingRate>" << mLengthRemodellingRate << "</LengthRemodellingRate>\n";

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