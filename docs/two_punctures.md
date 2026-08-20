# TwoPunctures initial data

The code creates Bowen-York initial data for two puncture black holes using a single domain spectral method. The method is described in Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy, "A single-domain spectral method for black hole puncture data" PRD 70, 064011 (2004) [arXiv:gr-qc/0404056](https://arxiv.org/abs/gr-qc/0404056).

## Build instructions

To build the BinaryBH example with TwoPunctures support, you need to do the following:

1. Clone the TwoPunctures repo next to GRTeclyn.

By default, the build system will assume the TwoPunctures repository is cloned next to GRTeclyn i.e. `/path/to/GRTeclyn/../TwoPunctures` (like AMReX) but this can be overridden by setting the `TWOPUNCTURES_HOME` variable, either

- on the command line during the build e.g.
```
make -j 4 USE_TWOPUNCTURES=TRUE TWOPUNCTURES_HOME=/path/to/TwoPunctures
```

- in the Make.local-pre or Make.local files located in amrex/Tools/GNUMake

- by setting it as an environment variable in your shell (also by adding it to your `.bashrc` file)
```
export TWOPUNCTURES_HOME=/mypath/to/twopunctures
```

2. Make sure you have GSL installed/available in your environment (e.g. by loading a gsl module). If it's installed in a non-standard place or the module doesn't set `CPATH` and `LIBRARY_PATH` environment variables, you can set an environment variable `GSL_HOME` to its base installation directory.

3. To build the BinaryBH GRTeclyn example with TwoPunctures initial data, set the makefile variable `USE_TWOPUNCTURES = TRUE`.

```
cd /path/to/GRTeclyn/Examples/BinaryBH
```
then
```
make -j 4 USE_TWOPUNCTURES=TRUE
```

Note that TwoPunctures will be appended to the filename stem/build config, that is the executable will look something like `BinaryBH3d.gnu.MPI.TwoPunctures.ex`.

## Parameters

We provide an example parameter file for TwoPunctures runs in the `BinaryBH` example folder.

They are all scoped by the prefix `two_punctures`. They should be self explanatory and are mostly the same as in GRChombo and other TwoPunctures implementations. Commented out values have their default value assigned to them.

See the [Parameters guide](parameters.md) for more info.

A few tips:

- the offsets are by default in the x direction, which means for circular inspirals the main momenta should be in the y direction to keep them in the x-y plane. For head on mergers it is conventional to switch to a collision along the z axis (meaning that the weyl scalars only have m=0 components for our chosen tetrad).

- Since TwoPunctures does not support interpolation on the GPU, the interpolation of data onto the grid is done on the CPU into a MultiFab which has a host memory arena. The data is then copied asynchronously to the usual state MultiFab. In the case spectral interpolation is used (`two_punctures.use_spectral_interpolation = true`), the initial data interpolation is very slow and this is further exacerbated by the larger number of cells per MPI process one typically has when running on GPUs. It may take several hours. The easiest workaround is to compute the initial data in a CPU-only run with more MPI processes, write a checkpoint file and then restart on GPUs from that.

!!! warning "The `swap_xz` convention"

    This option is useful when doing head on mergers, as it makes sense for the merger axis to be aligned with the one used for the extraction of Weyl scalars (since then only $m=0$ modes are active). However, the convention here is a bit confusing. TwoPunctures always solves with the punctures on its internal x axis. Setting `two_punctures.swap_xz = true` exchanges the x and z axes during interpolation, so `offset_plus` and `offset_minus` describe physical z offsets instead of x offsets. This is an axis exchange, not a proper rotation, and the momentum and spin inputs use the internal TwoPunctures convention. In particular, with `swap_xz` enabled, a physical spin in the positive z direction is entered as a negative x component. Check vector components and signs carefully when enabling this option.
