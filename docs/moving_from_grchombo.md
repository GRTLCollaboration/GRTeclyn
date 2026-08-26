# Moving from GRChombo

This page is intended to help those coming from using the numerical relativity code [GRChombo](https://github.com/GRTLCollaboration/GRChombo) that was the predecessor code to GRTeclyn. If you are new to GRChombo and GRTeclyn, you can skip it!

## Why port to GRTeclyn?

GRTeclyn is a new and improved version of the GRChombo code. It supports GPUs, meaning your code can run as much as 10-20x faster if you have access to these resources. Even if you only have access to CPUs, GRTeclyn should be faster to run and has a few useful improvements and new features. Most importantly, GRTeclyn is now the default code for our collaboration and so all new features and updates will be focused here.

## Will it be painful?

No - it might even be fun! Porting your code will help you to understand both codes better, and most of the changes are cosmetic rather than fundamental, since Chombo and AMReX share common ancestry in BoxLibs. You should find most of GRTeclyn looks pretty familiar. In this page we will provide some detailed tips to get you started, but the basic idea is: look at how the binary or scalar field examples have changed from GRChombo to GRTeclyn, and then copy the same approach for your own examples.

## Key changes

- **No data_t**

In news that no one will be sad about, there is no longer any explicit vectorisation of the right hand side, so no more `data_t`s. Variables that live on the grid, and coordinates, are now all `amrex::Real`s, which is a type that can be compiled as either a double or single precision number. For this reason, you can just use normal `if` statements to test if some quantity is larger or smaller than another etc (but note that branching generally degrades GPU performance).

- **Accessing variables**

We now load and store variables directly to the grid, rather than creating a local copy of the variables at each point and then using those (i.e. the thing that was usually called `vars`). There is therefore no `enum_mapping_function` nonsense. We still have a way to access the variables in a readable form, but it uses pointers to the grid values, and so is much more efficient. The new way means you have to do `vars.chi()` to get a scalar variable, `vars.shift(i)` to get the i-th component of a vector and `vars.h(i,j)` to get the i,j component of a tensor. If you take a look at the new `CCZ4RHS` class you will get the idea. Note the use of curved (not square) brackets for accessing components in a tensor. (This arises because these are now all AMReX arrays under the hood).

- **Diagnostics are calculated when needed**

GRTeclyn no longer carries a separate diagnostic state alongside the evolution state. Diagnostics are registered as AMReX derived variables and calculated from the current evolution state only when they are requested for a plot file or extraction. See [**Diagnostics**](diagnostics.md) for how to register an existing diagnostic or implement your own.

- **Derivatives and symmetries**

We make use of symmetries where possible to save fewer values. This is particularly the case for symmetric tensors (`h` and `A`) and second derivatives. Consider the quantity `d2.h`. It is a rank 4 tensor with symmetries in the first and second tensorial (i,j) and third and fourth derivative (k,l) indices. You need to declare this as the right kind of object - `Tensor::Sym12Sym34Rank4` (this tells you that indices 1 and 2 are symmetric, as are indices 3 and 4, and that the overall rank is 4). Having done this, you can then index into it as d2.h(i,j,k,l) in the expected way. (Note again that all brackets are now curved and not square). As before derivative indices come last.

- **`Boxloops()` is now `ParallelFor`**

Previously when you wanted to iterate over the grid points on each level, you called `BoxLoops`. The equivalent in AMReX is `ParallelFor` and this is where the GPU magic happens. `Parallel4` runs the operations on each cell on the grid on GPUs simulataneously for each box. The function it calls is by default is no longer called `compute`, instead it is the`operator()` method of the class, so if we first declare the class
```
TraceARemoval trace_A_removal;
```
We then access the function we want in the ParallelFor loop as
```
    amrex::ParallelFor(state_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz)
                       {
                           trace_A_removal(ix, iy, iz, state_arrays[box_no]);
                       });
```
where we are passing the cell indices, and the state array (the variables) in that box.
Looking at the trace removal function in one of the examples will be the best way to see what is going on.

- **Telling it what to run on GPUs**

Certain functions, e.g. those responsible for the RHS calculations, must be instrumented by GPU macros, `AMREX_GPU_DEVICE` if they are to be run on the GPU. This is true regardless of the GPU architecture (the AMReX backend will expand these as appropriate for CUDA, SYCL etc.) They may also need `AMREX_FORCE_INLINE` as an additional function qualifier. If you start with an existing diagnostic or RHS file and copy this, it should have everything you need.

- **Particle interpolation**

All extraction functions are implemented using particles, which know their location and can extract variable values with 4th order accuracy from the grid. See the sections on the [**Particle Interpolator**](particle_interpolator.md) and [**Extraction**](extraction.md) for more details.

- **Parameters**

The parameters are now read in where used in each class, rather than passing around a long `m_p.params` struct. Parameter checking is preserved in the `SimulationParameters.hpp` class for key parameters, but this is optional if you add your own classes and parameters. See [**Parameters**](parameters.md) for more details.

Two key changes in the parameters that may trip you up are that to not regrid on a level you need to set the regrid interval entry to `-1` rather than zero. Checkpoints are also not viewable in AMReX, so you should output the variables you want to view to plot files, and this is done slightly differently to in GRChombo.
