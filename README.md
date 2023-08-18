# GRTeclyn

GRTeclyn (previously referred to as GRAMReX) is a new numerical relativity code
that is currently under development.  It is a port of the [GRChombo
code](https://github.com/GRChombo/GRChombo) (based on the Chombo libraries) to
the [AMReX](https://amrex-codes.github.io/) library in order to leverage AMReX's
GPU offload capabilities and ongoing active development.

The AMReX documentation can be found [here](https://amrex-codes.github.io/amrex/docs_html).

### Feature status

Currently only the BinaryBH example has been [partially] ported.

Please see the table below for a list of features/capabilities which have been
ported/implemented for this example. Note that this should be correct at the
time this table was created but may not be correct at the time you read this.

| Feature | Ported/Implemented | Notes |
| --- | --- | --- |
| CCZ4 Evolution | :heavy_check_mark:  | Only fourth order spatial derivatives |
| Boundary Conditions | :heavy_check_mark: :grey_question: | All except mixed BCs have been implemented |
| Plot/checkpoint files + restart | :heavy_check_mark: | AMReX-specific, not HDF5 |
| AMR/tagging criterion | :heavy_check_mark: :grey_question: | Very basic tagging without puncture tracking |
| Diagnostic variables | :heavy_check_mark: :grey_question: | (e.g. Constraints, $\Psi_4$) Implemented but needs refactoring to remove hardcoding and allow generalisation |
| AMR Interpolator | :x: | |
| GW extraction | :x: | Requires AMR Interpolator |
| Puncture tracking | :x: | Requires AMR Interpolator |
| TwoPunctures initial data | :x: | |


## Obtaining and building the code

### Prerequisites

You will need the following software
* Git
* GNU Make >= 3.81
* Python >= 2.7
* A Unix-like environment with `perl` and `sed` commands
* C++ compiler with C++17 support (e.g. GCC >= 8, Clang >= 6, Intel Classic >= 19.14)
* MPI implementation (optional)

Note that the C++17 requirement means that older compilers you have used to
build [GR]Chombo may not work with GRTeclyn but so long as you have a more recent
compiler available, if you can build GRChombo, you can probably build GRTeclyn
too.

### Obtaining the code

First `cd` into a directory you are happy to clone the code into. For
simplicity, we shall assume that is your home directory (so adjust any commands 
below accordingly if not).

The AMReX source code is hosted on
[GitHub](https://github.com/AMReX-Codes/amrex). Clone this with a command such
as

```bash
git clone https://github.com/AMReX-Codes/amrex.git
```
It will be cloned to the `amrex` directory.

Next, clone the GRTeclyn repository using a command such as
```bash
git clone https://github.com/GRChombo/GRTeclyn.git
```

> **Note**
> We have assumed that you have cloned both of
> these repositories to the same directory so that the `amrex` and `GRTeclyn`
> directories share the same parent directory. If you want to clone AMReX
> elsewhere, make sure to set the `AMREX_HOME` environment variable
> appropriately e.g. 
> ```bash
> export AMREX_HOME=/path/to/amrex
> ```

### Building the BinaryBH example

> **Warning** 
> See [this
> page](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#gcc-on-macos)
> if you want to build on macOS using Homebrew GCC.

If you are on a cluster make sure you have loaded modules which provide you with
the software listed under [prerequisites](#prerequisites) (e.g. you might need a
compiler newer than the default).

Now, navigate to the BinaryBH example directory
```bash
cd ~/GRTeclyn/Examples/BinaryBH
```

The default build options are set in
[Make.defaults](Tools/GNUMake/Make.defaults). They are very
similar to Chombo's, and, like Chombo, can be overriden on the command line.
If want to consistently override these, it is best to create a `Make.local-pre`
file at
```
amrex/Tools/GNUMake/Make.local-pre
```
which sets the build configuration variables as you wish. You might want to
change some of the following: 

* Set `USE_OMP = TRUE` to use OpenMP or `USE_CUDA = TRUE` to use CUDA GPU
  offload support.
* If you don't want to use the GNU compiler `g++` , change `COMP = gnu` to 
  `COMP = intel-llvm` (for the Intel LLVM-based C++ compiler `icpx`) or 
  `COMP = llvm` (for LLVM `clang++`).
* If you don't have an MPI implementation available, set `USE_MPI = FALSE`.

For more detailed build configuration information, consult the [AMReX
documentation](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html)

Now start building AMReX and the BinaryBH example. For example, to build with 4
processes, do
```bash
make -j 4
```
A new `tmp_build_dir` directory will be created under `~/GRTeclyn` to store the
compiled object and auxiliary files. Assuming all is well, you should have an
executable in the current directory of the form `main<config>.ex` e.g.
`main3d.gnu.MPI.OMP.ex`.

## Running the code

You can run the example as for GRChombo by passing the parameter file as the
first argument. For example, for launching with MPI, you might do something like

```bash
mpiexec -n 4 ./main3d.gnu.MPI.OMP.ex ./<params_file>.txt
```

## Unit tests

There are some tests in the [Tests](Tests) directory that are automatically
run using GitHub actions when new commits are pushed to this repository. Please
consult the [Tests README](Tests/README.md) on how to build and run them.

Note that only some GRChombo tests have been ported.

## License

GRTeclyn is licensed under the BSD 3-Clause License. Please see [LICENSE](LICENSE) for details.