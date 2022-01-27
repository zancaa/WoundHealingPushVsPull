#ifndef TESTBOUNDEDGROWTHMODEL_HPP_
#define TESTBOUNDEDGROWTHMODEL_HPP_

/* Begin by including the necessary header files. The first ones are common to all cell_based Chaste simulations */
#include <cxxtest/TestSuite.h> //Needed for all test files
#include "CellBasedSimulationArchiver.hpp" //Needed to save/load simulations
#include "AbstractCellBasedTestSuite.hpp" //Needed for cell-based tests: times simulations, generates random numbers and has cell properties
#include "CheckpointArchiveTypes.hpp" //Needed to use GetIdentifier() method
#include "CellBasedEventHandler.hpp"
#include "CellsGenerator.hpp"
#include "OffLatticeSimulation.hpp"
#include "OffLatticeSimulationWithLargeTissueStoppingEvent.hpp"
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
#include "DistanceFromBoundaryModifier.hpp"
#include "FixedBoundaryCondition.hpp"
#include "ContactInhibitionDistanceToWoundCellCycleModel.hpp"
#include "ShortAxisFixedHorizontalVertexBasedDivisionRule.hpp"
#include "VolumeTrackingModifier.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "CellVolumesWriter.hpp"
#include "BoundaryNodeVelocityWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"
#include "FakePetscSetup.hpp"

class TestBoundedGrowthModel : public AbstractCellBasedTestSuite
{
private:

    void SetFilename(std::string& filename, double force, double qvf, double d_max, int rep)
    {
        std::stringstream repAsString;
        repAsString << rep;
        std::stringstream forceValAsString;
        forceValAsString << force;
        std::stringstream qvfValAsString;
        qvfValAsString << qvf;
        std::stringstream dMaxValAsString;
        dMaxValAsString << d_max;
        
        filename = "WoundHealing/BoundedGrowth/dMax_"
              + dMaxValAsString.str()
              + "/LeadingEdgeForce_"
              + forceValAsString.str() 
              + "/MeanQuiescentVolumeFraction_"
              + qvfValAsString.str()      
              + "/Run_" + repAsString.str(); 
    };

public:
    /*
     * Wounded epidermal skin model with bounded growth of proliferative hub.
     */
    void TestPeriodicMonolayer()
    {
        // Simulation and sweep parameters
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_min"));
        double factive_min = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_min").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_max"));
        double factive_max = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_max").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mqvf_min"));
        double mqvf_min = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mqvf_min").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mqvf_max"));
        double mqvf_max = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mqvf_max").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-dmax_min"));
        double dmax_min = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-dmax_min").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-dmax_max"));
        double dmax_max = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-dmax_max").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-factive_num_sweeps"));
        double factive_num_sweeps = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-factive_num_sweeps").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-mqvf_num_sweeps"));
        double mqvf_num_sweeps = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mqvf_num_sweeps").c_str());

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-dmax_num_sweeps"));
        double dmax_num_sweeps = atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-dmax_num_sweeps").c_str());

        /* Reseed the random number generator via a command line argument */
        unsigned seed = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed");

        // Set up initial conditions
        RandomNumberGenerator::Instance()->Reseed(seed);
        // Generate mesh
        CylindricalHoneycombVertexMeshGenerator generator(10, 10);    // Parameters are: cells across, cells up
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Loop over the nodes. 
        for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
        {
            // For each node we create a cell with our cell-cycle model.
            ContactInhibitionDistanceToWoundCellCycleModel* p_model = new ContactInhibitionDistanceToWoundCellCycleModel();
            p_model->SetMinCellCycleDuration(6.0);
            p_model->SetMaxCellCycleDuration(16.0);
            p_model->SetMeanQuiescentVolumeFraction(0.05);
            p_model->SetIntervalWidthQuiescentVolumeFraction(0.01);
            p_model->SetDistanceToWoundMin(20);
            p_model->SetDistanceToWoundMax(20);

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_differentiated_type);
            p_cell->SetBirthTime(0.0);

            // Randomly generate number
            double random_number = RandomNumberGenerator::Instance()->ranf();
            double cell_target_area = 0.9 + 0.2*random_number;
            // Set cell data
            p_cell->GetCellData()->SetItem("target area", cell_target_area);
            
            // Push the cell back into the vector of cells.
            cells.push_back(p_cell);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        // Add data writers
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeVelocityWriter>();

        // Pass the cell population into an OffLatticeSimulationWithLargeTissueStoppingEvent
        OffLatticeSimulationWithLargeTissueStoppingEvent simulator(cell_population);

        // Update the volumes of the cells in CellData
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add Nagai Honda force
        MAKE_PTR(NagaiHondaForce<2>, p_nh_force);
        p_nh_force->SetNagaiHondaDeformationEnergyParameter(50);
        p_nh_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1);
        p_nh_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1);
        p_nh_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1);
        simulator.AddForce(p_nh_force);

        /* Make a pointer to the target area modifier and add it to the simulator.
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

        // Keep track of how far cells are from the wound
        MAKE_PTR(DistanceFromBoundaryModifier<2>, p_distance_modifier);
        p_distance_modifier->SetDepth(dmax_max+1.0);
        simulator.AddSimulationModifier(p_distance_modifier);

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
        std::string output_directory = "WoundHealing/BoundedGrowth/InitialConditions/Run_" + repAsString.str(); ;
        simulator.SetOutputDirectory(output_directory);
        simulator.SetNumericalMethod(numerical_method);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(5.0); 

        // Run the simulation
        simulator.Solve();

        // Save simulation
        CellBasedSimulationArchiver<2, OffLatticeSimulationWithLargeTissueStoppingEvent >::Save(&simulator);

        for (unsigned index = 0; index<=factive_num_sweeps; index ++)      
        {
            for (unsigned index2 = 0; index2<=mqvf_num_sweeps; index2 ++)
            {
                for (unsigned index3 = 0; index3<=dmax_num_sweeps; index3 ++)
                {
                    Timer::Reset();
                    RandomNumberGenerator::Instance()->Reseed(seed);
                    double force_val;
                    double qvf_val;
                    unsigned d_max;
                    // Load simulation
                    OffLatticeSimulationWithLargeTissueStoppingEvent* p_simulator = CellBasedSimulationArchiver<2, OffLatticeSimulationWithLargeTissueStoppingEvent >::Load(output_directory,5.0);
                    VertexBasedCellPopulation<2>* p_cell_population = static_cast<VertexBasedCellPopulation<2>*>(&(p_simulator->rGetCellPopulation()));

                    SimulationTime::Instance()->Destroy();
                    SimulationTime::Instance()->SetStartTime(0.0);

                    // Remove the forces and boundaries - redefined here
                    p_simulator->RemoveAllForces();

                    // Loop over cells to set up the leading edge and proliferative hub. 
                    for (typename VertexMesh<2,2>::VertexElementIterator elem_iter = p_cell_population->rGetMesh().GetElementIteratorBegin();
                        elem_iter != p_cell_population->rGetMesh().GetElementIteratorEnd();
                        ++elem_iter)
                    {
                        unsigned elem_index = elem_iter->GetIndex();
                        // Get cell associated with this element
                        CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(elem_index);

                        // For each node we create a cell with our cell-cycle model. 
                        ContactInhibitionDistanceToWoundCellCycleModel* p_cc_model = static_cast<ContactInhibitionDistanceToWoundCellCycleModel*>(p_cell->GetCellCycleModel());
                        p_cc_model->SetIntervalWidthQuiescentVolumeFraction(0.05);
                        qvf_val = mqvf_min + (double)index2 / (double(mqvf_num_sweeps)) * (mqvf_max-mqvf_min);
                        p_cc_model->SetMeanQuiescentVolumeFraction(qvf_val);
                        p_cc_model->SetDistanceToWoundMin(1);

                        // Set up proliferative hub
                        d_max = dmax_min + (double)index3 / double(dmax_num_sweeps) * (dmax_max-dmax_min);
                        p_cc_model->SetDistanceToWoundMax(d_max);

                        unsigned wound_dist = p_cell->GetCellData()->GetItem("boundary distance");
                        // Set up leading edge
                        if ( (wound_dist >= 1) && (wound_dist <= d_max))
                        {
                            p_cell->SetCellProliferativeType(p_transit_type);
                        }
                        p_cc_model->SetCellCycleDuration();

                        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
                        double mean = qvf_val;
                        double half_width = 0.025;

                        // Assign quiescent volume fraction
                        double qvf = (mean - half_width) + (2.0*half_width) * p_gen->ranf(); // U[MinCA,MaxCA]
                        p_cell->GetCellData()->SetItem("quiescent volume fraction", qvf);
                        p_cell->SetCellCycleModel(p_cc_model);

                        double birth_time = 11.0*RandomNumberGenerator::Instance()->ranf();
                        p_cell->SetBirthTime(-birth_time);
                    }

                    // Set the division rule
                    MAKE_PTR(ShortAxisFixedHorizontalVertexBasedDivisionRule<2>, p_division_rule_to_set);
                    p_cell_population->SetVertexBasedDivisionRule(p_division_rule_to_set);

                    // Add Nagai Honda force 
                    MAKE_PTR(NagaiHondaForce<2>, p_force);
                    p_force->SetNagaiHondaDeformationEnergyParameter(50);
                    p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1);
                    p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1);
                    p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1);
                    p_simulator->AddForce(p_force);

                    // Set time step and sampling time
                    double dt = 0.0005;
                    double sampling_multiple = 0.5/dt;

                    // Add leading edge force
                    MAKE_PTR(LeadingEdgeForce<2>, p_motile_cell_force);
                    force_val = factive_min + (double)index / double(factive_num_sweeps) * (factive_max-factive_min);
                    p_motile_cell_force->SetActiveForceParameter(force_val);        
                    p_simulator->AddForce(p_motile_cell_force);
                    
                    //Output division locations
                    p_simulator->SetOutputDivisionLocations(true);
                    //Set output directory
                    std::string filename;
                    SetFilename(filename, force_val, qvf_val, d_max, seed);
                    std::cout << filename << std::endl;
                    p_simulator->SetOutputDirectory(filename.c_str());
                    p_simulator->SetDt(dt);
                    p_simulator->SetSamplingTimestepMultiple(unsigned(sampling_multiple));
                    p_simulator->SetEndTime(48.0); //two days

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
                    // Tidy up
                    delete p_simulator;
                }
            }
        }
    }
};

#endif /* TESTBOUNDEDGROWTHMODEL_HPP_ */
