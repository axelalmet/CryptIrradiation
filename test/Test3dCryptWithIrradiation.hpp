#ifndef TEST3DCRYPTWITHIRRADIATION_HPP_
#define TEST3DCRYPTWITHIRRADIATION_HPP_

#include <cxxtest/TestSuite.h> //Needed for all test files

#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellBasedEventHandler.hpp"

#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "SmartPointers.hpp" //Enables macros to save typing
#include "NodesOnlyMesh.hpp" // Nodes-only mesh, i.e. no Tessellation
#include "CellsGenerator.hpp" // Place-holder to generate cell cycles
#include "OffLatticeSimulation.hpp" //Simulates the evolution of the population
#include "NodeBasedCellPopulation.hpp"
#include "DifferentiatedCellProliferativeType.hpp" // Stops cells from proliferating
#include "NoCellCycleModel.hpp"
#include "SimpleWntContactInhibitionCellCycleModel.hpp" // Simple Wnt-based cell cycle model
#include "WntConcentration.hpp" // Singleton for Wnt concentration
#include "ContactInhibitionCellCycleModel.hpp" // Contact-inhibition-based cell cycle model (may be useful!)
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "VolumeTrackingModifier.hpp" // Modifier to track cell volumes
#include "FakePetscSetup.hpp" //Forbids tests running in parallel
#include "GeneralisedLinearSpringForce.hpp" // The default from Chaste.
#include "StemCellRetainerForce.hpp" // Force to keep stem cells in the base
#include "CryptSurfaceBoundaryCondition.hpp" // Boundary condition to keep the crypt cells on a test-tube-shaped surface
#include "PlaneBasedCellKiller.hpp" // Cell killer that removes cells when they're past a certain point
#include "RandomProliferativeCellTargetedCellKiller.hpp" // Cell killer that targets stem cells (used for irradiation)
#include "PetscSetupAndFinalize.hpp"

#include "Debug.hpp"

static const std::string M_OUTPUT_DIRECTORY = "3dCrypt";
static const std::string M_SS_OUTPUT_DIRECTORY = M_OUTPUT_DIRECTORY + "/SS/";
static const std::string M_INJURY_OUTPUT_DIRECTORY = M_OUTPUT_DIRECTORY + "/INJURY";
static const double M_DT = 0.005; // Every 15 seconds
static const double M_STEADY_STATE_TIME = 168.0;
static const double M_INJURY_TIME = M_STEADY_STATE_TIME + 50.0;
static const double M_SS_SAMPLING_TIMESTEP = 12.0/M_DT; // Let's see how it's doing every 12 hours
static const double M_INJURY_SAMPLING_TIMESTEP = 5.0/M_DT;

class Test3dCryptWithIrradiation : public AbstractCellBasedTestSuite
{
public:
	// void TestSteadyState()
	// {
	// 	// To be careful, we reseed the random number generator 
	// 	unsigned index = 1; // Just pick a random one 
	// 	RandomNumberGenerator::Instance()->Reseed(100*index);

    //     // Set the number of cells
    //     unsigned num_cells = 150; // This parameter really sets the number of stem cells in the base

    //     // Crypt Setup
    //     double cell_radius = 3.5;
    //     double crypt_length = 70.0; //70
    //     double crypt_radius = 8.0/M_PI*6.0; 
    //     // For this size domain there are about 75 cells in the bottom hemisphere 

    //     // Mechanical parameters for the forces
    //     double spring_cutoff_length = cell_radius * 3.0; // Cut-off length
	// 	double division_separation = 0.1;
	// 	double spring_stiffness = 30.0; // Spring stiffness
    //     double stem_retainer_force_strength = 75.0; // Force to keep the stem cells in the base

    //     double movement_threshold = 50.0;

	// 	// Set the probability of apoptosis at the edges
	// 	// double apoptosis_probability = 0.05;

    //     // Set the Wnt thresholds for differentiation and set the contact inhibition thresholds
    //     double stem_threshold = 0.9; // Stem threshold before differentiating into transit cells
    //     double transit_threshold = 0.6; // Transit threshold before differentiating into terminally-differentiated cells
	// 	double equilibrium_volume = 4.0 / 3.0 * M_PI * cell_radius * cell_radius * cell_radius; // Equilibrium volume for contact inhibition (we're dealing with spheres now)
	// 	double quiescent_fraction = 0.9; // Fraction to initiate contact inhibition

    //     // Parameters for boundary conditions
    //     double max_distance_from_surface = 0.0;
    //     double target_population = 0.0;
    //     double remodelling_rate = 0.0;

    //     // Create some starter nodes on the base of the crypt
    //     std::vector<Node<3>*> nodes;
    //     for(unsigned node_index= 0;  node_index<num_cells; node_index++)
    //     {
    //         double u = RandomNumberGenerator::Instance()->ranf();
    //         double v = RandomNumberGenerator::Instance()->ranf();

    //         double random_azimuth_angle = 2*M_PI*u;
    //         double random_zenith_angle = std::acos(2*v - 1);

    //         double x = 0.5*crypt_radius*cos(random_azimuth_angle)*sin(random_zenith_angle);
    //         double y = 0.5*crypt_radius*sin(random_azimuth_angle)*sin(random_zenith_angle);
    //         double z = 0.5*crypt_radius*(1.0 + cos(random_zenith_angle));

    //         // double x = crypt_radius/2.0 * sin(node_index*2.0*M_PI/num_cells);
    //         // double y = crypt_radius/2.0 * cos(node_index*2.0*M_PI/num_cells);
    //         // double z = 0.0;
    //         nodes.push_back(new Node<3>(node_index, false, x, y, z));
    //     }

    //     // Convert this to a NodesOnlyMesh
    //     NodesOnlyMesh<3> mesh;
    //     mesh.ConstructNodesWithoutMesh(nodes, cell_radius*3.0);

    //     // Create cells
    //     std::vector<CellPtr> cells;

    //     // Initialise the cell proliferative types and mutation states
    //     boost::shared_ptr<AbstractCellProperty> p_stem_type = CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>();
    //     boost::shared_ptr<AbstractCellProperty> p_transit_type = CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>();
	// 	boost::shared_ptr<AbstractCellProperty> p_wildtype_state = CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>();

    //     for (unsigned index = 0; index < nodes.size(); index++)
    //     {

    //         // Wnt-based contact inhibition cell cycle
	// 		SimpleWntContactInhibitionCellCycleModel* p_cycle_model = new SimpleWntContactInhibitionCellCycleModel();
    //         p_cycle_model->SetWntStemThreshold(stem_threshold);
    //         p_cycle_model->SetWntTransitThreshold(transit_threshold);
    //         p_cycle_model->SetQuiescentVolumeFraction(quiescent_fraction);
    //         p_cycle_model->SetEquilibriumVolume(equilibrium_volume);
	// 		p_cycle_model->SetDimension(3);

	// 		CellPtr p_cell(new Cell(p_wildtype_state, p_cycle_model));
	// 		// p_cell->InitialiseCellCycleModel(); // For paranoia really.

	// 		p_cell->SetCellProliferativeType(p_stem_type);

    //         // We set cells to be transit cells if they're in the upper hemisphere
    //         // of the initial node generation
    //         double z = nodes[index]->rGetLocation()[2]; // Get the node height

    //         if (z >= 0.5 * crypt_radius)
    //         {
    //             p_cell->SetCellProliferativeType(p_transit_type);
    //         }

    //         // Set the cell data
    //         p_cell->GetCellData()->SetItem("volume", equilibrium_volume);
    //         p_cell->GetCellData()->SetItem("Radius", cell_radius);

    //         cells.push_back(p_cell); // Add the cell to the vector
    //     }

    //     // Create cell population
    //     NodeBasedCellPopulation<3> cell_population(mesh, cells);
    //     cell_population.SetUseVariableRadii(true);
    //     cell_population.SetAbsoluteMovementThreshold(movement_threshold); // Flags if the cells move too far.
    //     cell_population.SetMeinekeDivisionSeparation(division_separation); // Set the initial division separation

    //     // Set the Wnt gradient
    //     WntConcentration<3>::Instance()->SetType(LINEAR);
    //     WntConcentration<3>::Instance()->SetCellPopulation(cell_population);
    //     WntConcentration<3>::Instance()->SetCryptLength(crypt_length);

	// 	OffLatticeSimulation<3> simulator(cell_population);

	// 	// Set steady state output directory
	// 	simulator.SetOutputDirectory(M_SS_OUTPUT_DIRECTORY);
	// 	simulator.SetDt(M_DT);
	// 	simulator.SetSamplingTimestepMultiple(M_SS_SAMPLING_TIMESTEP); //Sample the simulation at every hour
	// 	simulator.SetEndTime(M_STEADY_STATE_TIME); //Hopefully this is long enough for a steady state

	// 	// Add modifier to track cell volumes
	// 	MAKE_PTR(VolumeTrackingModifier<3>, p_volume_tracking_modifier);
	// 	simulator.AddSimulationModifier(p_volume_tracking_modifier); 

	// 	// Add generalised linear spring force
	// 	MAKE_PTR(GeneralisedLinearSpringForce<3>, p_spring_force);
	// 	p_spring_force->SetCutOffLength(spring_cutoff_length);
	// 	p_spring_force->SetMeinekeSpringStiffness(spring_stiffness);
	// 	simulator.AddForce(p_spring_force);

    //     // Add stem cell retainer force
    //     MAKE_PTR(StemCellRetainerForce<3>, p_retainer_force);
    //     p_retainer_force->SetRetainerForceStrength(stem_retainer_force_strength);
    //     simulator.AddForce(p_retainer_force);

    //     // Add the crypt surface boundary condition
    //     MAKE_PTR_ARGS(CryptSurfaceBoundaryCondition<3>, p_surface_bc, (&cell_population, max_distance_from_surface));
    //     p_surface_bc->SetTargetPopulation(target_population);
    //     p_surface_bc->SetLengthRemodellingRate(remodelling_rate);
    //     simulator.AddCellPopulationBoundaryCondition(p_surface_bc);

    //     // Add a plane-based cell killer to remove cells from the top
    //     c_vector<double, 3> crypt_top = crypt_length*unit_vector<double>(3,2);

    //     c_vector<double, 3> crypt_top_normal = unit_vector<double>(3,2);

    //     MAKE_PTR_ARGS(PlaneBasedCellKiller<3>, p_cell_killer, (&cell_population, crypt_top, crypt_top_normal));
    //     simulator.AddCellKiller(p_cell_killer);

	// 	simulator.Solve();

    //     CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Save(&simulator);

	// 	WntConcentration<3>::Destroy();
	// }

    // Test if the irradiation model works
    void TestInjury()
    {
		double apoptosis_probability = 0.5; // Set the probability of killing the proliferative cells.
        double crypt_length = 70.0; // Should probably set this parameter outside of both tests

        // Load the steady state directory
		OffLatticeSimulation<3>* p_simulator = CellBasedSimulationArchiver<3, OffLatticeSimulation<3> >::Load(M_SS_OUTPUT_DIRECTORY, M_STEADY_STATE_TIME);

		// Re-initialise the Wnt concentration
		WntConcentration<3>::Instance()->SetType(LINEAR);
        WntConcentration<3>::Instance()->SetCellPopulation(p_simulator->rGetCellPopulation());
        WntConcentration<3>::Instance()->SetCryptLength(crypt_length);

		// Set the new output directory, end time, and sampling timestep
		p_simulator->SetOutputDirectory(M_INJURY_OUTPUT_DIRECTORY);
		p_simulator->SetSamplingTimestepMultiple(M_INJURY_SAMPLING_TIMESTEP); //Sample the simulation at every hour
		p_simulator->SetEndTime(M_INJURY_TIME); //Hopefully this is long enough for a steady state
        
        // We will also remove the boundary condition and add it with new 
        p_simulator->RemoveAllCellPopulationBoundaryConditions();

        double target_population = (double)p_simulator->rGetCellPopulation().rGetMesh().GetNumNodes();
        double remodelling_rate = 0.1;
        double max_distance_from_surface = 0.0;

        // Add the crypt surface boundary condition
        MAKE_PTR_ARGS(CryptSurfaceBoundaryCondition<3>, p_surface_bc, (&p_simulator->rGetCellPopulation(), max_distance_from_surface));
        p_surface_bc->SetTargetPopulation(target_population);
        p_surface_bc->SetLengthRemodellingRate(remodelling_rate);
        p_simulator->AddCellPopulationBoundaryCondition(p_surface_bc);

        // We may need to turn the retainer force off via p_simulator->rGetForceCollection()
        boost::shared_ptr<StemCellRetainerForce<3> > p_retainer_force =
                boost::static_pointer_cast<StemCellRetainerForce<3> >(p_simulator->rGetForceCollection()[1]);
        p_retainer_force->SetRetainerForceStrength(0.0); // Turn off the retainer force (as we're killing stem cells)

		// Add the cell killer to irradiate the crypt
		MAKE_PTR_ARGS(RandomProliferativeCellTargetedCellKiller<3>, p_irradiation_killer, (&p_simulator->rGetCellPopulation()));
		p_irradiation_killer->SetApoptosisProbability(apoptosis_probability);
		p_simulator->AddCellKiller(p_irradiation_killer);

		p_simulator->Solve(); // Run the simulation again

		WntConcentration<3>::Destroy();

    }
};

#endif /* TEST3DCRYPTWITHIRRADIATION_HPP_ */
