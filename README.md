# Building and running the GRAMReX BinaryBH example with SYCL

## Obtaining and building the code

### Prerequisites

You will need the following software
* Git
* GNU Make >= 3.81
* Python >= 2.7
* A Unix-like environment with `perl` and `sed` commands
* Intel oneAPI DPC++ compiler
* Intel OpenCL Offline Compiler (OCLOC)
* Intel oneMKL
* MPI implementation (optional)

### Obtaining the code

First `cd` into a directory you are happy to clone the code into. For
simplicity, I will assume that is your home directory (so adjust any commands 
below accordingly if not).

The upstream AMReX source code is hosted on
[GitHub](https://github.com/AMReX-Codes/amrex). However, we will be using a
branch on [my fork](https://github.com/mirenradia/amrex). Clone this with a command such
as

```bash
git clone -b enhancement/sycl_aot_build_speedup https://github.com/AMReX-Codes/amrex.git
```
It will be cloned to the `amrex` directory with the
`enhancement/sycl_aot_build_speedup` branch checked out.

The GRAMReX repository is currently private so first check you have access by
navigating to https://github.com/GRChombo/GRAMReX. If you're reading this on
your own device, then you have the necessary permissions. If you don't have
access, please let me (@mirenradia) know so I can give them to you. Clone the
repository with a command such as

```bash
git clone -b test/sycl_aot_pvc_notes https://github.com/GRChombo/GRAMReX.git
```
It will cloned to the `GRAMReX` directory with the `test/sycl_aot_pvc_notes`
branch checked out. 
> :information_source: **Note**
> I have assumed that you have cloned both of
> these repositories to the same directory so that the `amrex` and `GRAMReX`
> directories share the same parent directory. If you want to clone AMReX
> elsewhere, make sure to set the `AMREX_HOME` environment variable
> appropriately e.g. 
> ```bash
> export AMREX_HOME=/path/to/amrex
> ```

### Building the code

First make sure you have set up your environment to ensure you have access to
the software listed under [prerequisites](#prerequisites). You may need to load
appropriate modules or `source` the oneapi `setvars.sh` script. The latter is
usually done using the command
```bash
source /opt/intel/oneapi/setvars.sh
```

The default build options can be found in
[Make.defaults](./Tools/GNUMake/Make.defaults). They can be modified either on
the command line or by creating a file with the overridden variables at
```
$HOME/amrex/Tools/GNUMake/Make.local-pre
```

Create this file with the following contents
```
DEBUG = TRUE
USE_SYCL = TRUE
SYCL_AOT = TRUE
AMREX_INTEL_ARCH = pvc
SYCL_MAX_PARALLEL_LINK_JOBS = 16
```
Set the last variable to an appropriate number (e.g. the number of CPU cores you
have access to). You can set `AMREX_INTEL_ARCH` to any value accepted by the
`-device` argument to `ocloc`. To check what these are, use the command
```
ocloc compile --help
```
If you don't want to use MPI, add
```
USE_MPI = FALSE
```

If you want to add additional compiler/linker flags, you can create a file at
```
$HOME/amrex/Tools/GNUMake/Make.local
```
with the following contents
```
CXXFLAGS += <extra compiler flags>

LDFLAGS += <extra linker flags>
```

A new `tmp_build_dir` directory will be created under `$HOME/GRAMReX` to store
the compiled object and auxiliary files. Assuming all is well, you should have
an executable in the current directory of the form `main<config>.ex` e.g.
`main3d.sycl.DEBUG.MPI.ex`

## Running the code

To run the code, simply execute it (using `mpiexec` if built with MPI) with the
parameter file e.g. [inputs.test](./Examples/BinaryBH/inputs.test) as the only
argument e.g. (with MPI)

```bash
mpiexec -n 2 ./main3d.sycl.DEBUG.MPI.ex ./inputs.test
```