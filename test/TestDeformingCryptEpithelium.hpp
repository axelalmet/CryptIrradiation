#ifndef TESTDEFORMINGCRYPTPITHELIUM_HPP_
#define TESTDEFORMINGCRYPTPITHELIUM_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "CylindricalHoneycombMeshGenerator.hpp" //Generates mesh
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "PlaneBasedDivisionRule.hpp" // Division rule where the plane of division is determined by the epithelial neighbours
#include "VoronoiDataWriter.hpp" //Allows us to visualise output in Paraview
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp" // Stops cells from proliferating
#include "StromalCellProliferativeType.hpp" // Cell type to model the underlying stroma
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp" // Uniform cell cycle model
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "GeneralisedLinearSpringForce.hpp" // The default from Chaste.
#include "NonPeriodicBasementMembraneForce.hpp" //BM force (based off Dunn et al. (2012))
#include "AnoikisCellKiller.hpp" // Anoikis-based cell killer
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "CrossSectionalCrypt";
static const std::string M_SS_OUTPUT_DIRECTORY = M_OUTPUT_DIRECTORY + "/DEFORMATION/";
static const double M_DT = 0.005;
static const double M_STEADY_STATE_TIME = 100.0;
static const double M_SS_SAMPLING_TIMESTEP = 0.1*M_STEADY_STATE_TIME/M_DT;

class TestDeformingCryptEpithelium : public AbstractCellBasedTestSuite
{
public:
	void TestCryptDeformation()
	{
		// To be careful, we reseed the random number generator 
		unsigned index = 28; // Just pick a random one (for some reason SJD picked this)
		RandomNumberGenerator::Instance()->Reseed(100*index);

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 15.0;
		double epithelial_stromal_stiffness = 15.0;
		double stromal_stromal_stiffness = 15.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 30;
		unsigned cells_up = 20;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 12.0;
		double target_curvature = 0.6;

        // Mechanical parameters for the spring force
        double spring_cutoff_length = 1.5; // Cut-off length

		// Generate the periodic mesh
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		//Get the initial non-ghost indices
		std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

		//Get the maximum width and height of the real nodes to define the monolayer
		double max_height = 0.0;

		// These parameters define the initial crypt boundary
		double crypt_edge_left = 0.4 * (double)cells_across;
		double crypt_edge_right = 0.6 * (double)cells_across;

        // Get the maximum height to set the initial epithelium
		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned node_index = location_indices[i];
			c_vector<double, 2> node_location = p_mesh->GetNode(node_index)->rGetLocation();
			double y = node_location[1];

			if (y > max_height)
			{
				max_height = y;
			}
		}


		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_stromal_type = CellPropertyRegistry::Instance()->Get<StromalCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_transit_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();
        boost::shared_ptr<AbstractCellProperty> p_cell_label = CellPropertyRegistry::Instance()->Get<CellLabel>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<location_indices.size(); i++)
		{

            unsigned node_index = location_indices[i];
			c_vector<double, 2> node_location = p_mesh->GetNode(node_index)->rGetLocation();
			double y = node_location[1];

            if (y == max_height)
            {
                // Initialise placeholder cell cycle
                UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel(); //Uniformly distributed cell cycle times
                p_cycle_model->SetDimension(2);

                // To avoid synchronisations in division
                double birth_time = 12.0 * RandomNumberGenerator::Instance()->ranf();
                p_cycle_model->SetBirthTime(-birth_time);

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->InitialiseCellCycleModel(); // For paranoia really.

                p_cell->SetCellProliferativeType(p_transit_type);

                cells.push_back(p_cell);
            }
            else // Set to be a non-dividing stromal cell
            {
                // Initialise placeholder cell cycle
                NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Uniformly distributed cell cycle times
                p_cycle_model->SetDimension(2);

                CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
                p_cell->InitialiseCellCycleModel(); // For paranoia really.

                p_cell->SetCellProliferativeType(p_stromal_type);

                // Also fix the bottom row of cells
                if (y == 0.0)
                {
                    p_cell->AddCellProperty(p_cell_label);
                }

                cells.push_back(p_cell);
            }

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		// Add our plane based division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new PlaneBasedDivisionRule<2,2>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		OffLatticeSimulation<2> simulator(cell_population);

		// Set steady state output directory
		simulator.SetOutputDirectory(M_SS_OUTPUT_DIRECTORY);
		simulator.SetDt(M_DT);
		simulator.SetSamplingTimestepMultiple(M_SS_SAMPLING_TIMESTEP); //Sample the simulation at every hour
		simulator.SetEndTime(M_STEADY_STATE_TIME); //Hopefully this is long enough for a steady state

		// Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
		// MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
		MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
		p_spring_force->SetCutOffLength(spring_cutoff_length);
		// p_spring_force->SetMeinekeSpringStiffness(epithelial_epithelial_stiffness);
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
		p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
		p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
		simulator.AddForce(p_spring_force);

		//Add basement membrane force
		MAKE_PTR(NonPeriodicBasementMembraneForce, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
		p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
        p_bm_force->SetLeftCryptBoundary(crypt_edge_left); // Left boundary that imposes non-zero target curvature
        p_bm_force->SetRightCryptBoundary(crypt_edge_right); // Right boundary that imposes non-zero target curvature
		p_bm_force->ApplyForceToCrypt(true);
		p_bm_force->ApplyVerticallyDependentTargetCurvature(false);
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		p_anoikis_killer->SetApoptosisProbability(0.0);
		simulator.AddCellKiller(p_anoikis_killer);

		//Fix the bottom row of cells
		c_vector<double, 2> point, normal;

		point(0) = 0.0;
		point(1) = 0.5;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		simulator.Solve();
	}
};

#endif /* TESTDEFORMINGCRYPTEPITHELIUM_HPP_ */
