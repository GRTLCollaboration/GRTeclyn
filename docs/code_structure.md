# Code structure

Yes, we know, GRTeclyn is a *big* code. At first the number of files will seem overwhelming, but with time you will start to learn where to find things and the structure will make (some) sense. 

On this page we provide some hints on how to find your way around the code, but in the end you just have to dive in and learn as you go.

Some useful references can be found in [Useful resources](useful_resources.md). One should look at the guides on C++ classes, inheritance and templating, which are used extensively in the code - some basic knowledge of these concepts is assumed below.

## Hierarchy of GRTeclyn

The code is designed to have 3 main levels in its hierarchy, as follows:

1. Specific Example related files, e.g. for BinaryBH - specific actions relevant to the BinaryBH example - key classes include BHAMR, BinaryBHLevel, SimulationParameters. Also important are the namespaces UserVariables and DiagnosticVariables in the BinaryBH examples folder and BinaryBHInitialData (the initial data).

   The functions that are specified at this level include things like setting initial data, calculating example-specific diagnostics, reading in example-specific parameters, specifying the tagging criterion.

This inherits most of the functionality from:

2. GRTeclyn - specific physics actions common to most   GR problems - key classes include GRAMR, GRAMRLevel, and SimulationParametersBase . See also the CCZ4UserVariables namespace (and most of the contents of the GRTeclyn Source folders).

   The functions that are specified at this level include things like performing the RHS calculation for CCZ4 and matter variables, calculating constraints, calculating the finite derivatives, checking for Nans, setting up and reading in the CCZ4 parameters.

This in turn inherits most of the functionality from:

3. AMReX - overall program flow relevant to any hyperbolic initial value problem with AMR - key classes: AMR, AMRLevel.

   The functions that are specified at this level include things like setting up the initial AMR hierarchy and performing the AMR regridding and the Runge-Kutta update.

## Where to find the files

Logically, all of the files related to level 1 : BinaryBH should be in the specific Example folder. However, there are a few exceptions:

- The initial data class BinaryBHInitialData is considered sufficiently general (ie, it will be used in many examples *without modification*) to be included in the `Source` code of GRTeclyn rather than in the Example folder itself, so it is in [Source/BlackHoles](https://github.com/GRTLCollaboration/GRTeclyn/tree/develop/Source/BlackHoles). For matter classes it is probable that you will want to put the initial condition code in the Example folder itself, as it is more likely to be problem specific (as in [InitialScalarData](https://github.com/GRTLCollaboration/GRTeclyn/blob/develop/Examples/ScalarField/InitialScalarData.hpp) ).

- Similarly many functions related to evolving the black holes, and tracking their punctures, are generally useful enough to be in `Source`. See in particular the [Source/BlackHoles folder](https://github.com/GRTLCollaboration/GRTeclyn/tree/develop/Source/BlackHoles).

- We often include tagging criteria in the Source code as they are reusable in many examples and provide a useful library of examples. See the [Source/TaggingCriteria folder](https://github.com/GRTLCollaboration/GRTeclyn/tree/develop/Source/Tagging).

Logically, all of the files related to level 2:GRTeclyn should be in the GRChombo/Source folder. However, there are exceptions:

- Some of the Example specific code is in here too if it is sufficiently general in use, as discussed above.
- Some functions that morally belong in AMReX are here too, as discussed below.

Logically, all of the files related to level 3:AMReX should be in the AMReX/src folder. However, there are exceptions:

- Some met up functions - tbc.

- The AMRInterpolator - tbc.

- The BoundaryConditions code - tbc.

Note that the main AMReX file is AMR.cpp, which controls the overall program flow. You don't really need to understand this file - you can just trust it will use the functions you give it in your Example Level file at the right moment, but if you want to really understand what is going on it is worth taking a quick look. If your question is "when and why is my code doing this step?" the answer probably lies in here somewhere. 

## Hooks and virtual functions

Each part of the AMReX/GRTeclyn/BinaryBH Hierarchy has some awareness of the part above and below it. This is provided by, for example, AMRLevel providing a virtual function that does nothing, but gets overwritten in GRAMRLevel. Functions may also provide hooks for certain actions, again defined via overwriting virtual functions, at certain points. Exciting detective work is often required for finding the connections between functions. The command `grep` is your friend here.

## A note on AMR (GRAMR/BHAMR) versus AMRLevel (GRAMRLevel/BinaryBHLevel)

Here we describe a key point which most users fail to grasp initially, and even experienced users have been known to get wrong - the difference between AMR and AMRLevel. It is always worth some extra thinking time, and probably also some outputting to check what is going on matches what you meant to do. (Don't ever feel ashamed to add a line `pout() << "I am here doing X on level " << m_level << endl;` to the code.)

AMR controls the program flow for the entire hierarchy - it knows that, for example, 6 levels of refinement exist, with the coarsest level having a certain value of `coarsest_dx`, etc.

AMRLevel is then a class for which an instance is created for each of these six levels. So there are 6 copies of it that get called in turn, in an order that is determined by the AMR class (as described in the previous section). Each instance has its own value for the level specific parameters like grid spacing `m_dx` and `m_level`. Any instructions in an AMRLevel class will happen on each level in turn, and won't affect the other levels unless you explicitly ask them to talk to each other. Usually they only know about and can access data on the levels above and below them, but they can appeal to the AMR class for wider control (this is required to use the AMRInterpolator, for example). 

So, for example, if you write in the postTimeStep() function a command to write out `"hello world"`, in the pout files, you will get this output on every level of the hierarchy, after each of its timesteps conclude. In one coarsest time-step level 0 will write out once, level 1 will write out 2^1 times, level 2 will write out 2^2 times, etc. This will be a lot of output.

If instead you want something to *only* happen once every coarse time-step, you will need to bracket it with an if statement that requires `if(m_level == 0)` so only the level 0 instance of the class takes the action. An example of this is something like writing out a global diagnostic. You probably don't want this to output at all the intermediate times, and even on the coarsest timestep, you probably only want it written out once and not by all 6 levels.

If you want something to happen on *every* level but *only at a time which is a multiple of the coarsest timestep*, you need to have an if statement that requires this, i.e. ``if(at_level_timestep_multiple(0))`` as in the BinaryBHLevel. An example of this is something like calculating the values of a diagnostic variable across the whole grid (ie, we need the calculation done for data on every level), where we plan to output only on the coarsest timestep. It would of course not crash the code to call this on every level postTimeStep, but it would be a big waste of computing time since the diagnostic would only be output once, whereas on the 6th level there will be 2^6 timesteps (and therefore 2^6 -1 redundant calculations of the diagnostic) in between the outputs. 

Getting this wrong can significantly slow down the code, and can be a source of incorrect results where, for example, finest level data is not updated before output. 

If you use the code a lot, at some point you will get it wrong, despite having now been warned about it. But at least it will make you feel less bad to know that others have done the same.