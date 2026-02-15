# Building on CPUs

Unlike in the predecessor code GRChombo, you don't have to separately build AMReX and GRTeclyn. Building an example in GRTeclyn will compile the relevant parts of AMReX too.

But you should have checked that you have the set up described in [Prerequisites](prerequisites.md) before trying to compile an example.

## Building an example

Navigate to the example directory, e.g. for the BinaryBH example
```bash
cd GRTeclyn/Examples/BinaryBH
```

The default build options are set in [Make.defs](https://github.com/AMReX-Codes/amrex/blob/development/Tools/GNUMake/Make.defs). They are very similar to Chombo's, and, like Chombo, can be overriden on the command line, e.g. by doing something like
```
make -j 4 USE_OMP=TRUE
```

If want to consistently override these, it is best to create a `Make.local-pre` file at
```
amrex/Tools/GNUMake/Make.local-pre
```
which sets the build configuration variables as you wish. There is an example in `${AMREX_HOME}/Tools/GNUMake/Make.local.template`. You might want to change some of the following: 

* Set `USE_OMP = TRUE` to use OpenMP
* If you don't want to use the GNU compiler `g++` , change `COMP = gnu` to 
  `COMP = intel-llvm` (for the Intel LLVM-based C++ compiler `icpx`) or 
  `COMP = llvm` (for LLVM `clang++`).
* If you don't have an MPI implementation available, set `USE_MPI = FALSE`.

See the separate page [Building and Running on GPUs](building_gpus.md) for more details on the options required for using GPUs.

For more detailed build configuration information, consult the [AMReX documentation](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html)

Now start building AMReX and the example. For example, to build with 4 processes (which will speed it up), do
```bash
make -j 4
```
A new `tmp_build_dir` directory will be created to store the compiled object and auxiliary files. Assuming all is well, you should have an executable in the current directory of the form `main<config>.ex` e.g. `main3d.gnu.MPI.OMP.ex`.

If you want to clean the current build configuration you can do
```bash
make cleanconfig
```
> Note that this is similar to `make clean` for GRChombo examples.

If you want to clean all build configuration, you can do
```bash
make clean
```
> Note that this is similar to `make realclean` for GRChombo examples.

Tip: To speed up subsequent compilations, you can install
[Ccache](https://ccache.dev/) and then add
```
USE_CCACHE = TRUE
```
to your `Make.local-pre` configuration file. This is particularly useful for GPU compilation which tends to be a bit slower. Alternatively just go and make a cup of tea and it should be done when you get back. It's good to take a break.

## Running the code

You can run the example by passing the parameter file as the first argument. For example, for launching with MPI, you might do something like

```bash
mpiexec -n 4 ./main3d.gnu.MPI.OMP.ex params.txt
```

On a cluster you will probably need to submit a jobscript for the run. You should consult the system documentation, but some example jobscripts can be found in [Example configs for specific HPC systems](example_configs.md) along with the `make.local-pre` files and module choices if relevant.