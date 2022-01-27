#ifndef TESTFREETISSUE_HPP_
#define TESTFREETISSUE_HPP_

/* Begin by including the necessary header files. */
#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed if we would like to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed if we use GetIdentifier() method (which we do)
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp" //Proliferative cell type
#include "StemCellProliferativeType.hpp"
#include "SmartPointers.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "LeadingEdgeForce.hpp"
#include "SimpleNonUniformTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "FixedBoundaryCondition.hpp"
#include "NoCellCycleModel.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "VolumeTrackingModifier.hpp"
#include "CellVolumesWriter.hpp"
#include "CellLabel.hpp"
#include "BoundaryNodeVelocityWriter.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"
#include "FakePetscSetup.hpp"

class TestFreeTissue : public AbstractCellBasedTestSuite
{
private:

    void SetFilename(std::string& filename, double force)
    {
        std::stringstream forceValAsString;
        forceValAsString << force;
        
        filename = "WoundHealing/FreeTissue/LeadingEdgeForce_"
              + forceValAsString.str();
    };

public:
    /*
     * Generate the free tissue simulations to compare against full wound 
     * healing simulations. 
     */
    void TestGenerateFreeTissue()
    {
        // Simulation and sweep parameters
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_min"));
        double factive_min = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_min").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_max"));
        double factive_max = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_max").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_num_sweeps"));
        double factive_num_sweeps = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_num_sweeps").c_str());
        
        double end_time = 48.0;
        double num_sweeps = 10.0;
        /* Note that changing the seed will (slightly) alter the initial conditions,
        * and therefore the results, but the effect is negligible.
        */
        unsigned seed = 1;

        // Set up initial conditions
        RandomNumberGenerator::Instance()->Reseed(seed);
        // Generate mesh
        CylindricalHoneycombVertexMeshGenerator generator(10, 10);    // Parameters are: cells across, cells up
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Create vector of cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);

        // Loop over the nodes
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            // For each node we create a cell with no cell-cycle model
            NoCellCycleModel* p_model = new NoCellCycleModel();

            CellPtr p_cell(new Cell(p_state, p_model));

            //Randomly generate number
            double random_number = RandomNumberGenerator::Instance()->ranf();
            double cell_target_area = 0.9 + 0.2*random_number;
            // Set cell data
            p_cell->GetCellData()->SetItem("target area", cell_target_area);
            
            // Push the cell back into the vector of cellss
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        // Add data writers
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeVelocityWriter>();

        // Pass the cell population into an OffLatticeSimulation
        OffLatticeSimulation<2> simulator(cell_population);

        /// Update the volumes of the cells in CellData
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add Nagai Honda force
        MAKE_PTR(NagaiHondaForce<2>, p_nh_force);
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(50);
        p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1);
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1);
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1);
        simulator.AddForce(p_nh_force);

        /* We also make a pointer to the target area modifier and add it to the simulator.
        * Set the target area of individual cells to be drawn from U[0.9,1.1], so avoid
        * 'slanting' in steady state due to numerical issues.
        */
        MAKE_PTR(SimpleNonUniformTargetAreaModifier<2>, p_growth_modifier);
        p_growth_modifier->SetGrowthDuration(1);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Add boundary condition
        MAKE_PTR_ARGS(FixedBoundaryCondition<2>, p_bc, (&cell_population));
        p_bc->SetFixTopCells(false);
        simulator.AddCellPopulationBoundaryCondition(p_bc);

        // Add leading edge force
        MAKE_PTR(LeadingEdgeForce<2>, p_le_force);
        p_le_force->SetActiveForceParameter(0.0);        
        simulator.AddForce(p_le_force);

        // Set time step
        double dt = 0.005;

        // Create Numerical method
        MAKE_PTR(ForwardEulerNumericalMethod<2>, numerical_method);
        numerical_method->SetCellPopulation(&cell_population);
        
        // Output division locations
        simulator.SetOutputDivisionLocations(true);
        // Set output directory
        std::stringstream repAsString;
        repAsString << seed;
        std::string output_directory = "WoundHealing/FreeTissue/InitialConditions/Run_" + repAsString.str(); ;
        simulator.SetOutputDirectory(output_directory);
        simulator.SetNumericalMethod(numerical_method);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(5.0);

        // Run the simulation.
        simulator.Solve();

        // Save simulation
        CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Save(&simulator);

        //Tidying up
        SimulationTime::Instance()->Destroy();
        SimulationTime::Instance()->SetStartTime(0.0);
    
        for (unsigned index = 0; index<=factive_num_sweeps; index ++)  
        {
            Timer::Reset();
            RandomNumberGenerator::Instance()->Reseed(seed);
            double force_val;
            // Load simulation
            OffLatticeSimulation<2>* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulation<2> >::Load(output_directory,5.0);
            VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);

            // Remove the forces and boundaries - redefined here
            p_simulator->RemoveAllForces();
            p_simulator->RemoveAllCellPopulationBoundaryConditions();

            // Add Nagai Honda force
            MAKE_PTR(NagaiHondaForce<2>, p_force);
            p_force->SetNagaiHondaDeformationEnergyParameter(50);
            p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1);
            p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1);
            p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1);
            p_simulator->AddForce(p_force);

            /* Add a plane boundary condition so cells on the bottom boundary 
             * don't 'stretch out' into negative space
             */
            c_vector<double,2> point = zero_vector<double>(2);
            c_vector<double,2> normal = zero_vector<double>(2);
            normal(1) = -1.0;
            MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (p_cell_population, point, normal));
            p_simulator->AddCellPopulationBoundaryCondition(p_bc2);

            // Set time step and sampling time
            double dt = 0.0005;
            double sample_timestep = 0.5/dt;

            // Add leading edge force
            MAKE_PTR(LeadingEdgeForce<2>, p_motile_cell_force);
            force_val = factive_min + (double)index / double(factive_num_sweeps) * (factive_max-factive_min);
            p_motile_cell_force->SetActiveForceParameter(force_val);            
            p_simulator->AddForce(p_motile_cell_force);

            //Set output directory
            std::string filename;
            SetFilename(filename, force_val);
            std::cout << filename << std::endl;
            p_simulator->SetOutputDirectory(filename.c_str());
            p_simulator->SetNumericalMethod(numerical_method);
            p_simulator->SetDt(dt);
            p_simulator->SetSamplingTimestepMultiple(unsigned(sample_timestep));
            p_simulator->SetEndTime(end_time); //two days

            // Run the simulation.
            try
            {
                p_simulator->Solve();
            }
            catch (Exception& e)
            {
                // If it throws then we report the error message and go to the next simulation
                WARNING("Simulation didnt run" << filename << ".");
                WARNING(e.GetMessage());            
            }

            OutputFileHandler output_file_handler(filename, false);
            out_stream p_stream = output_file_handler.OpenOutputFile("timing.dat");
            *p_stream <<  Timer::GetElapsedTime();
            p_stream->close();

            PRINT_VARIABLE(Timer::GetElapsedTime());

            //Tidying up
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
            delete p_simulator;
        }
    }
};

#endif /* TESTFREETISSUE_HPP_ */
