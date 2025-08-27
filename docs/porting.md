This page is intended to help those coming from using the numerical relativity code [GRChombo](https://github.com/GRTLCollaboration/GRChombo) that was the predecessor code to GRTeclyn. If you are new to GRChombo and GRTeclyn, you can skip it!

## Why port to GRTelcyn?

GRTeclyn is a new and improved version of the GRChombo code. It supports GPUs, meaning your code can run as much as 10x faster if you have access to these resources. Even if you only have access to CPUs, GRTeclyn should be faster to run and has a few useful improvements and new features. Most importantly, at some point GRTeclyn will be the default code for our collaboration and so all new features and updates will be focused here. 

## Will it be painful?

No - it might even be fun! Porting your code will help you to understand both codes better, and most of the changes are cosmetic rather than fundamental, since Chombo and AMReX share common ancestry in BoxLibs. You should find most of GRTelcyn looks pretty familiar. In this page we will provide some detailed tips to get you started, but the basic idea is: look at how the binary or scalar field examples have changed from GRChombo to GRTeclyn, and then copy the same approach for your own examples.

## Key changes

1. No `data_t` ...

2. `Boxloops()` is now `ParallelFor`.

3. Storing and loading of variables is replaced by direct access to the grid ... wrapper access to specific variables for readability.

4. Certain functions, e.g. those responsible for the RHS calculations, must be instrumented by GPU macros, `AMREX_GPU_DEVICE` if they are to be run on the GPU. This is true regardless of the GPU architecture (the AMReX backend will expand these as appropriate for CUDA, SYCL etc.) They may also need `AMREX_FORCE_INLINE` as an additional function qualifier. 

5. Puncture tracking is done with spectator particles.
 