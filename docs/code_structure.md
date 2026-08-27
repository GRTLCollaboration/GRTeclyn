# Code structure

Yes, we know, GRTeclyn is a *big* code. At first the number of files will seem overwhelming, but with time you will start to learn where to find things and the structure will make (some) sense.

On this page we provide some hints on how to find your way around the code, but in the end you just have to dive in and learn as you go.

Some useful references can be found in [Useful resources](useful_resources.md). One should look at the guides on C++ classes, inheritance and templating, which are used extensively in the code - some basic knowledge of these concepts is assumed below. The [Matter classes](matter_classes.md) page also contains a short introduction to how templates are used in the RHS code.

## Hierarchy of GRTeclyn

The code has three main levels in its hierarchy:

1. **Example-specific code**, such as `BinaryBHLevel` or `ScalarFieldLevel`. This level sets the initial data, registers state variables and [diagnostics](diagnostics.md), reads example-specific parameters, selects the RHS calculations and specifies the tagging criterion. The files normally live in the corresponding folder under `Examples`.

2. **GRTeclyn**, which provides physics and infrastructure shared by several examples. The main classes are `GRAMR` and `GRAMRLevel`; most of the other code at this level lives under `Source`. It includes the CCZ4 equations, matter classes, finite derivatives, boundary conditions, parameter handling, interpolation and extraction.

3. **AMReX**, which controls the overall program flow for a block-structured AMR evolution. The important base classes are `amrex::Amr` and `amrex::AmrLevel`. AMReX constructs the grid hierarchy, regrids it, fills boundary and ghost cells, advances each level and writes plot files and checkpoints.

The inheritance chain for a typical level is therefore

```text
amrex::AmrLevel -> GRAMRLevel -> ScalarFieldLevel (or another example level)
```

and for the object which manages the whole hierarchy it is

```text
amrex::Amr -> GRAMR -> an optional example-specific AMR class
```

For example, `ScalarFieldAMR` adds the particle interpolators used by the Scalar Field example, while `BHAMR` adds puncture tracking for the Binary Black Hole example. An example which needs no extra hierarchy-wide data can use `GRAMR` directly, as Klein Gordon does.

## Where to find the files

The top-level folders divide the code by purpose:

- `Examples` contains complete programs which can be built and run. Each example supplies a main file, a level class, `SimulationParameters.hpp`, `StateVariables.hpp`, a parameter file and any problem-specific initial data or physics.

- `Source` contains reusable GRTeclyn code. Its current subfolders are `BlackHoles`, `CCZ4`, `GRTeclynCore`, `Grids`, `IO`, `Maths`, `Matter`, `ParticleInterpolator` and `Tagging`.

- `Tests` contains unit and regression tests. If you change reusable code in `Source`, this is a good place to look for both examples of its intended use and tests which should be extended.

- `Tools` contains development and build tools, while `docs` contains the pages you are reading.

Code should normally be placed at the most general level at which it can be reused. Problem-specific initial data belongs in its example folder; for example, the oscillaton initial data is in [`Examples/ScalarField/OscillatonInitialData.hpp`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/ScalarField/OscillatonInitialData.hpp). General black hole initial data and puncture tracking live in [`Source/BlackHoles`](https://github.com/GRTLCollaboration/GRTeclyn/tree/main/Source/BlackHoles), because several examples may use them. Reusable tagging criteria similarly live in [`Source/Tagging`](https://github.com/GRTLCollaboration/GRTeclyn/tree/main/Source/Tagging).

Header files commonly contain small or templated implementations. Where that would make a header unwieldy, GRTeclyn puts the definitions in a matching `.impl.hpp` file and includes it at the end of the main header. Ordinary non-templated functions may instead be defined in a `.cpp` file.

## Important files in an example

The level class, such as `ScalarFieldLevel`, is the best place to start. It inherits `GRAMRLevel` and overrides hooks for the work specific to that example:

- `variableSetUp()` registers the evolved state and the available diagnostics;

- `initData()` fills the initial state;

- `specificEvalRHS()` calculates the RHS during each Runge-Kutta stage;

- `specificAdvance()` and `specificUpdateODE()` perform any extra work around an update;

- `specificPostTimeStep()` performs work after a completed level time step; and

- `tag_cells()` selects cells for refinement.

`StateVariables.hpp` defines the component indices, output names, boundary parities and asymptotic values for everything carried in the evolution state. `SimulationParameters.hpp` gathers the parameter checks needed by that example. The main file creates the level factory and AMR object, initializes the hierarchy and runs the coarse time-step loop.

## Hooks and virtual functions

Each part of the AMReX/GRTeclyn/example hierarchy has some awareness of the part above and below it. AMReX calls virtual functions on `GRAMRLevel`; the GRTeclyn implementation performs the common work and calls hooks such as `specificEvalRHS()` or `specificPostTimeStep()` which the example level overrides.

The [`DefaultLevelFactory`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/GRTeclynCore/DefaultLevelFactory.hpp) tells AMReX which example level class to construct on every refinement level. This is another use of templating: `DefaultLevelFactory<ScalarFieldLevel>` creates `ScalarFieldLevel` objects, while the same factory code can be reused for other examples.

If your question is "when and why is my code doing this step?", begin at the corresponding function in the example level, then follow the override into `GRAMRLevel` and finally into AMReX. The command `rg` (or `grep`) is your friend here.

## A note on AMR versus AMRLevel

Here we describe a key point which most users fail to grasp initially, and even experienced users have been known to get wrong - the difference between AMR and AMRLevel. It is always worth some extra thinking time, and probably also some outputting to check what is going on matches what you meant to do. (Don't ever feel ashamed to add a line which prints "I am here doing X on level Y" to the code.)

The AMR object controls the program flow for the entire hierarchy. It knows, for example, how many levels of refinement exist, the grid spacing and time step on each level, and when the levels need to be advanced or regridded.

An AMRLevel object represents just one level of that hierarchy. If six refinement levels currently exist, there are six example-level objects, each with its own values of `Geom().CellSize()`, `Level()` and the state on that level. Instructions in a level class happen on each level in an order determined by the AMR object. A level can access its parent hierarchy through `get_gramr_ptr()` when wider coordination is needed.

So, for example, if you write a command in `specificPostTimeStep()` which outputs `"hello world"`, you will get this output after every time step on every level. With a refinement ratio of two, during one level 0 time step, level 1 normally takes two steps, level 2 takes four, and so on. This will be a lot of output.

If instead you want something to happen only once per coarse time step, it will usually be initiated by the level 0 object, so bracket it with `if (Level() == 0)`. Writing a global diagnostic or performing an extraction are common examples.

Sometimes work must first be done on every level at a time shared with level 0, before level 0 combines or outputs the result. In that case, test whether the current time is a multiple of the level 0 time step and make sure all levels have current data before doing the level 0 operation. Calculating a diagnostic on every fine-level substep would not usually be wrong, but most of those results would never be output and the extra work can significantly slow the code.

Getting this wrong can be a source of both poor performance and incorrect results, for example if the finest-level data has not been updated before a reduction or extraction. If you use the code a lot, at some point you will get it wrong, despite having now been warned about it. But at least it will make you feel less bad to know that others have done the same.
