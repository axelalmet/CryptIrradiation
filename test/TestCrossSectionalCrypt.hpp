#ifndef TESTCROSSSECTIONALCRYPT_HPP_
#define TESTCROSSSECTIONALCRYPT_HPP_

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
#include "SimpleWntContactInhibitionCellCycleModel.hpp" // Simple Wnt-based cell cycle model
#include "WntConcentration.hpp" // Singleton for Wnt concentration
#include "ContactInhibitionCellCycleModel.hpp" // Contact-inhibition-based cell cycle model (may be useful!)
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "VolumeTrackingModifier.hpp" // Modifier to track cell volumes
#include "FixedRegionPlaneBoundaryCondition.hpp" // Fixed-position boundary condition
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "LinearSpringForceWithVariableRestLength.hpp"
#include "GeneralisedLinearSpringForce.hpp" // The default from Chaste.
#include "NonPeriodicBasementMembraneForce.hpp" //BM force (based off Dunn et al. (2012))
#include "AnoikisCellKiller.hpp" // Anoikis-based cell killer
#include "RandomCryptEdgesCellKiller.hpp" // Random-apoptosis-based cell killer
#include "RandomProliferativeCellTargetedCellKiller.hpp" // Random-apoptosis-based cell killer
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "CrossSectionalCrypt";
static const std::string M_SS_OUTPUT_DIRECTORY = M_OUTPUT_DIRECTORY + "/SS/APOPTOSISATEDGES/";
static const std::string M_INJURY_OUTPUT_DIRECTORY = M_OUTPUT_DIRECTORY + "/INJURY";
static const double M_DT = 1.0/240.0; // Every 15 seconds
static const double M_STEADY_STATE_TIME = 50.0;
static const double M_INJURY_TIME = M_STEADY_STATE_TIME + 50.0;
static const double M_SS_SAMPLING_TIMESTEP = 0.1*M_STEADY_STATE_TIME/M_DT;
static const double M_INJURY_SAMPLING_TIMESTEP = 1.0/M_DT;

class TestCrossSectionalCrypt : public AbstractCellBasedTestSuite
{
public:
	void TestSteadyState()
	{
		// To be careful, we reseed the random number generator 
		unsigned index = 28; // Just pick a random one (for some reason SJD picked this)
		RandomNumberGenerator::Instance()->Reseed(100*index);

		//Set all the spring stiffness variables
		double epithelial_epithelial_stiffness = 45.0;
		double epithelial_stromal_stiffness = 45.0;
		double stromal_stromal_stiffness = 40.0;

		//Set the number of cells across and down for the array
		unsigned cells_across = 10;
		unsigned cells_up = 29;
		unsigned ghosts = 2; //Set the number of ghost node layers

		//Set the basement membrane force parameters
		double bm_stiffness = 12.0;
		double target_curvature = 0.3;

        // Mechanical parameters for the spring force
		// double resting_spring_length = 1.0; // Rest length
        // double spring_cutoff_length = 1.5; // Cut-off length
		double division_separation = 0.1;

        // Set the Wnt thresholds for differentiation
        double stem_threshold = 1.0; // Stem threshold before differentiating into transit cells
        double transit_threshold = 0.65; // Transit threshold before differentiating into terminally-differentiated cells
		double equilibrium_volume = 0.5*sqrt(3.0); // Equilibrium volume for contact inhibition
		double quiescent_fraction = 0.5; // Fraction to initiate contact inhibition

		// Set the probability of apoptosis at the edges
		double apoptosis_probability = 0.05;

		// Generate the periodic mesh
		CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, ghosts);
		Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

		//Get the maximum width and height of the real nodes to define the monolayer
		double max_height = 0.0;
		double max_width = 0.0;

		// Aim to get 20-25 cells down each side (See Van Leeuwen 2005)
		double width = (double) cells_across - 0.5;
		double height = ((double) cells_up - 1.0)*sqrt(3.0)*0.5;

		double centre_line = width*0.5;                         // The axis of symmetry down the centre of the crypt
		double radius = width*0.25;                             // Distance from centre line to edge of crypt

		// These parameters define the initial crypt boundary
		double crypt_edge_left = centre_line - radius;
		double crypt_edge_right = centre_line + radius;
		double crypt_base = height*0.25;

		c_vector<double,2> circle_centre;
		circle_centre[0] = centre_line;
		circle_centre[1] = crypt_base;

		// Re-defined the actual non-ghost indices that defines a crypt shape
		std::vector<unsigned> location_indices;
		std::vector<unsigned> ghost_node_indices; 

        // The initial mesh basically looks like a block-y |_|. The BM force
        // corrects the geometry to make it look more like a crypt shape, i.e. U
		for (unsigned i = 0; i < p_mesh->GetNumAllNodes(); i++)
		{
			c_vector<double, 2> node_location = p_mesh->GetNode(i)->rGetLocation();
			double x = node_location[0];
			double y = node_location[1];

			// Need to calculate the length of the vector between the node and the circle centre
			c_vector<double, 2> vector_to_circle_centre = node_location - circle_centre;

			double distance_to_centre = norm_2(vector_to_circle_centre);
			assert(distance_to_centre > 0);
			assert(!isnan(distance_to_centre));

			if ( ( (crypt_edge_left <= x) && (x <= crypt_edge_right)
					&& (y >= crypt_base) ) || (distance_to_centre <= radius) || (x < -1e-6) ||
					(x > width) || (y < -1e-6) || (y > height) )
			{
					ghost_node_indices.push_back(i);
			}
			else
			{
				location_indices.push_back(i);

				if (y > max_height)
				{
					max_height = y;
				}

				if (x > max_width)
				{
					max_width = x;
				}
			}
		}

		//Create shared pointers for cell and mutation states
		boost::shared_ptr<AbstractCellProperty> p_stromal_type = CellPropertyRegistry::Instance()->Get<StromalCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
		boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

		//Create tissue of cells. Initially we set them all to be differentiated
		std::vector<CellPtr> cells; //Create vector of cells

		for (unsigned i = 0; i<location_indices.size(); i++)
		{

			//Set stochastic duration based cell cycle
			NoCellCycleModel* p_cycle_model = new NoCellCycleModel(); //Uniformly distributed cell cycle times
			p_cycle_model->SetDimension(2);

			CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
			p_cell->InitialiseCellCycleModel(); // For paranoia really.

			p_cell->SetCellProliferativeType(p_stromal_type);

			cells.push_back(p_cell);

		}

		//Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
		cell_population.SetMeinekeDivisionSeparation(division_separation);

		// Add our plane based division rule
        boost::shared_ptr<AbstractCentreBasedDivisionRule<2,2> > p_division_rule_to_set(new PlaneBasedDivisionRule<2,2>());
        cell_population.SetCentreBasedDivisionRule(p_division_rule_to_set);

		//Create the epithelium of stem cells
		for (unsigned i = 0; i < location_indices.size(); i++)
		{
			unsigned cell_index = location_indices[i];
			CellPtr cell_iter = cell_population.GetCellUsingLocationIndex(cell_index);
			double x = cell_population.GetLocationOfCellCentre(cell_iter)[0];
			double y = cell_population.GetLocationOfCellCentre(cell_iter)[1];

			//"un-differentiate" the epithelium
			if (y == max_height)
			{
				Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);

				//Iterate over the elements (triangles) containing the nodes
				for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
						iter != p_node->ContainingElementsEnd();
						++iter)
				{
					bool element_contains_ghost_nodes = false;

					// Get a pointer to the element
					Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

					// Check whether its triangulation contains a ghost node
					for (unsigned local_index=0; local_index<3; local_index++)
					{	
						unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

						if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
						{
							element_contains_ghost_nodes = true;
							break; 				// This should break out of the inner for loop
						}
					}

					if(element_contains_ghost_nodes) // If the cell has ghost nodes in its neighbourhood, it's likely an epithelial cell
					{
                        // Set the cell cycle model to now depend on Wnt.
                        SimpleWntContactInhibitionCellCycleModel* p_cycle_model = new SimpleWntContactInhibitionCellCycleModel(); //Uniformly distributed cell cycle times
                        p_cycle_model->SetDimension(2);
                        p_cycle_model->SetWntStemThreshold(stem_threshold); // Set the Wnt threshold before stem cells differentiate
                        p_cycle_model->SetWntTransitThreshold(transit_threshold); // Set the Wnt threshold before transit cells differentiate
						p_cycle_model->SetEquilibriumVolume(equilibrium_volume); // Set equilibrium volume for contact inhibition
						p_cycle_model->SetQuiescentVolumeFraction(quiescent_fraction); // Volume fraction for contact inhibition

                        double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
                        cell_iter->SetBirthTime(-birth_time);
						cell_iter->SetCellProliferativeType(p_stem_type); // Initialise as stem cell
                        cell_iter->SetCellCycleModel(p_cycle_model); // Change to wnt-dependent cell cycle
					}
				}
			}
			else if ( (x >= crypt_edge_left - 1.0)&&(x <= crypt_edge_right + 1.0)&&(y > crypt_base - radius - 1.0) )
			{
				Node<2>* p_node = cell_population.GetNodeCorrespondingToCell(cell_iter);

				//Iterate over the elements (triangles) containing the nodes
				for (Node<2>::ContainingElementIterator iter = p_node->ContainingElementsBegin();
						iter != p_node->ContainingElementsEnd();
						++iter)
				{
					bool element_contains_ghost_nodes = false;

					// Get a pointer to the element
					Element<2,2>* p_element = cell_population.rGetMesh().GetElement(*iter);

					// Check whether its triangulation contains a ghost node
					for (unsigned local_index=0; local_index<3; local_index++)
					{
						unsigned nodeGlobalIndex = p_element->GetNodeGlobalIndex(local_index);

						if (cell_population.IsGhostNode(nodeGlobalIndex) == true)
						{
							element_contains_ghost_nodes = true;
							break; 				// This should break out of the inner for loop
						}
					}

					if(element_contains_ghost_nodes) // If the cell's attached to ghost nodes, it's likely epithelial.
					{
                        // Set the cell cycle model to now depend on Wnt.
                        SimpleWntContactInhibitionCellCycleModel* p_cycle_model = new SimpleWntContactInhibitionCellCycleModel(); //Uniformly distributed cell cycle times
                        p_cycle_model->SetDimension(2);
                        p_cycle_model->SetWntStemThreshold(stem_threshold); // Set the Wnt threshold before stem cells differentiate
                        p_cycle_model->SetWntTransitThreshold(transit_threshold); // Set the Wnt threshold before transit cells differentiate
						p_cycle_model->SetEquilibriumVolume(equilibrium_volume); // Set equilibrium volume for contact inhibition
						p_cycle_model->SetQuiescentVolumeFraction(quiescent_fraction); // Volume fraction for contact inhibition

                        double birth_time = 12.0*RandomNumberGenerator::Instance()->ranf(); //We would like the birth time to be ~U(0,12) and set in the past
                        cell_iter->SetBirthTime(-birth_time);
						cell_iter->SetCellProliferativeType(p_stem_type); // Initialise as stem cell
                        cell_iter->SetCellCycleModel(p_cycle_model); // Change to wnt-dependent cell cycle
					}
				}
			}
			else if (y == 0.0) // Add labels to the bottom row for the boundary conditoin
			{
				boost::shared_ptr<AbstractCellProperty> p_label = CellPropertyRegistry::Instance()->Get<CellLabel>();
				cell_iter->AddCellProperty(p_label);
			}
		}

		//Output data to vtk format so we can visualise it in Paraview
		cell_population.SetWriteVtkAsPoints(true);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(height + 3.0);

		OffLatticeSimulation<2> simulator(cell_population);

		// Set steady state output directory
		simulator.SetOutputDirectory(M_SS_OUTPUT_DIRECTORY);
		simulator.SetDt(M_DT);
		simulator.SetSamplingTimestepMultiple(M_SS_SAMPLING_TIMESTEP); //Sample the simulation at every hour
		simulator.SetEndTime(M_STEADY_STATE_TIME); //Hopefully this is long enough for a steady state

		// Add modifier to track cell volumes
		MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier); 

		// Add linear spring force (modified to have three different spring stiffnesses, depending on the type of pair)
		// MAKE_PTR(GeneralisedLinearSpringForce<2>, p_spring_force);
		MAKE_PTR(LinearSpringForceWithVariableRestLength<2>, p_spring_force);
		// p_spring_force->SetCutOffLength(spring_cutoff_length);
		// p_spring_force->SetMeinekeSpringStiffness(epithelial_epithelial_stiffness);
		p_spring_force->SetEpithelialEpithelialSpringStiffness(epithelial_epithelial_stiffness); //Default is 15
		p_spring_force->SetEpithelialStromalSpringStiffness(epithelial_stromal_stiffness); //Default is 15
		p_spring_force->SetStromalStromalSpringStiffness(stromal_stromal_stiffness); //Default is 15
		simulator.AddForce(p_spring_force);

		//Add basement membrane force
		MAKE_PTR(NonPeriodicBasementMembraneForce, p_bm_force);
		p_bm_force->SetBasementMembraneParameter(bm_stiffness); //Equivalent to beta in SJD's papers
		p_bm_force->SetTargetCurvature(target_curvature); //This is equivalent to 1/R in SJD's papers
		p_bm_force->ApplyForceToCrypt(true);
		p_bm_force->ApplyVerticallyDependentTargetCurvature(true);
		simulator.AddForce(p_bm_force);

		//Add anoikis-based cell killer
		MAKE_PTR_ARGS(AnoikisCellKiller, p_anoikis_killer, (&cell_population));
		p_anoikis_killer->SetApoptosisProbability(apoptosis_probability);
		simulator.AddCellKiller(p_anoikis_killer);

		// Add random-apoptosis-based cell killer
		// MAKE_PTR_ARGS(RandomCryptEdgesCellKiller, p_apoptosis_killer, (&cell_population));
		// p_apoptosis_killer->SetApoptosisProbability(apoptosis_probability);
		// simulator.AddCellKiller(p_apoptosis_killer);

		//Fix the bottom row of cells
		c_vector<double, 2> point, normal;

		point(0) = 0.0;
		point(1) = 0.5;
		normal(0) = 0.0;
		normal(1) = -1.0;
		MAKE_PTR_ARGS(FixedRegionPlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
		simulator.AddCellPopulationBoundaryCondition(p_bc1);

		simulator.Solve();

        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

		WntConcentration<2>::Destroy();
	}

	// void TestIrradiateCrypt()
	// {
	// 	double apoptosis_probability = 0.01; // Set the probability of killing the proliferative cells.

    //     // Load the steady state directory
	// 	OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(M_SS_OUTPUT_DIRECTORY, M_STEADY_STATE_TIME);

	// 	// Get the maximum height 
	// 	double max_height = 0.0;

	// 	for (AbstractCellPopulation<2>::Iterator cell_iter = p_simulator->rGetCellPopulation().Begin();
	// 			cell_iter != p_simulator->rGetCellPopulation().End();
	// 			++cell_iter)
	// 	{
	// 		double y = p_simulator->rGetCellPopulation().GetLocationOfCellCentre(*cell_iter)[1];

	// 		if (y > max_height)
	// 		{
	// 			max_height = y;
	// 		}
	// 	}

	// 	// Re-initialise the Wnt concentration
	// 	WntConcentration<2>::Instance()->SetType(LINEAR);
    //     WntConcentration<2>::Instance()->SetCellPopulation(p_simulator->rGetCellPopulation());
    //     WntConcentration<2>::Instance()->SetCryptLength(max_height);

	// 	// Set the new output directory, end time, and sampling timestep
	// 	p_simulator->SetOutputDirectory(M_INJURY_OUTPUT_DIRECTORY);
	// 	p_simulator->SetSamplingTimestepMultiple(M_INJURY_SAMPLING_TIMESTEP); //Sample the simulation at every hour
	// 	p_simulator->SetEndTime(M_INJURY_TIME); //Hopefully this is long enough for a steady state

	// 	// Add the cell killer to irradiate the crypt
	// 	MAKE_PTR_ARGS(RandomProliferativeCellTargetedCellKiller, p_irradiation_killer, (&p_simulator->rGetCellPopulation()));
	// 	p_irradiation_killer->SetApoptosisProbability(apoptosis_probability);
	// 	p_simulator->AddCellKiller(p_irradiation_killer);

	// 	p_simulator->Solve(); // Run the simulation again

	// 	WntConcentration<2>::Destroy();
	// }
};

#endif /* TESTCROSSSECTIONALCRYPT_HPP_ */
