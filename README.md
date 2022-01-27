# WoundHealingProliferativeHubLeadingEdge

This project contains the code necesary for running the simulations presented in Zanca et al. "Push or pull? Cell proliferation and migration during wound healing".

If you are unfamiliar with Chaste, you may wish to look at some of the basic user tutorials for Chaste https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials.

## Installing
To install Chaste, see https://chaste.cs.ox.ac.uk/chaste/tutorials/release_2021.1/.

You will also need the source for the WoundHealingProliferativeHubLeadingEdge project.  This can be obtained by checking out the version from this github repository to the projects folder of the Chaste directory.

Navigate to the `projects` folder and clone this project into a new folder.
```
git clone https://github.com/zancaa/WoundHealingProliferativeHubLeadingEdge
```

## Documentation
There are four folders - `src`, `test`, `plot_scripts` and `build`.
 1. The `src` folder contains the classes necesary to run the simulation. These define the additional forces and boundary conditions not in the core chaste code.
 2. The `test` folder contains:
  * TestContinualGrowthModel - this file can be run to generate the results in Figures 2,3 and 4.
  * TestBoundedGrowthModel - this file can be run to generate the results in Figure 5.
  * run_[free/continual/bounded]_script.sh - Script to run multiple simultions
 3. The `plot_scripts` folder contains MATLAB scripts to generate the Figures.
 4. The `build` folder houses the compiled tests.

## Running simulations
The code for each set of simulations can be found in the `test` folder in the project directory. In order to run the tests, they must first be compiled, then executed with the minimum and maximum values, and the number of sweeps for the active force, mean quiescent volume fraction, and (for the bounded growth model) d_max and seed value as inputs. To compile the test, navigate to the main Chaste directory and run the following command:
```
scons b=GccOpt cl=0 co=1 projects/WoundHealingProliferativeHubLeadingEdge/test/TestContinualGrowthModel.hpp
```

Then to run the simulation:
```
cd projects/WoundHealingProliferativeHubLeadingEdge/test/
sh run_continual_script.sh
```
(For the bounded growth model, replace `continual` in the above code with `bounded`.)

## Running plot code

Simulation snapshots were generated with Paraview. All other plot code can be found in the project directory in the folder `plot_scripts`. All plotting code is in MATLAB. 

**NB**: the paper was developed with release version 2021.1.

For further information on using Chaste, see the extensive guide material (https://chaste.cs.ox.ac.uk/trac/wiki/ChasteGuides).
You may also wish to look at some of the basic user tutorials (https://chaste.cs.ox.ac.uk/trac/wiki/UserTutorials).
