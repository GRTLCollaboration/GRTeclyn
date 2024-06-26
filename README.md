# GRTeclyn Scalar Fields example hands-on activity

## Introduction

This hands-on exercise involves porting the RHS for Scalar Field
example. An empty function `specificEvalRHS` has already been provided for you in
`ScalarFieldLevel.cpp`.

Other functions in `ScalarFieldLevel` have already been defined:
* `initData`: setup the initial conditions (scalar bubble)
* `specificAdvance` : impose trace free `Aij`, positive `chi` and `alpha`
* `specificUpdateODE`: trace free `Aij`
* `errorEst`: tag cells for refinement
* `derive`: calculate derived quantities e.g. `Ham` and `Mom`





### Status
| Feature | Ported/Implemented | Notes |
| --- | --- | --- |
| CCZ4 Evolution | :heavy_check_mark:  | Only fourth order spatial derivatives |
| Boundary Conditions | :heavy_check_mark: :grey_question: | All except mixed BCs have been implemented |
| Plot/checkpoint files + restart | :heavy_check_mark: | AMReX-specific, not HDF5 |
| AMR/tagging criterion | :heavy_check_mark: | Very basic tagging without puncture tracking |
| Diagnostic variables | :heavy_check_mark: | (e.g. Constraints, $\Psi_4$) Implemented but needs refactoring to remove hardcoding and allow generalisation|
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
* C compiler with C99 support
* C++ compiler with C++17 support (e.g. GCC >= 8, Clang >= 6, Intel Classic >= 19.14)
* MPI implementation (optional)

Note that the C++17 requirement means that older compilers you have used to
build [GR]Chombo may not work with GRTeclyn but so long as you have a more recent
compiler available, if you can build GRChombo, you can probably build GRTeclyn
too. If you the accounts provided for you on Dawn, these requirements won't be a problem. 

### Logging into CSD3 and Dawn 

Research Computing Services has kindly provided us with training accounts on Dawn but 
feel free to try this example locally on your own system (the example parameter
file is relatively small so should run fine on a laptop) if you have the above
software.

> [!IMPORTANT]
> Refer to Kacper's email for first time login instructions

To access a login node use:

```ssh <your_user_id>@login.hpc.cam.ac.uk```

Next download the code. 

If the CSD3 login nodes are slow, try using the [web interface](https://login-web.hpc.cam.ac.uk/pun/sys/dashboard/) or use the `login-p-X` nodes. 

### Obtaining the code

First `cd` into a directory you are happy to clone the code into. For
simplicity, I will assume that is your home directory (so adjust any commands 
below accordingly if not).

The AMReX source code is hosted on
[GitHub](https://github.com/AMReX-Codes/amrex). Clone this with a command such
as

```bash
git clone https://github.com/AMReX-Codes/amrex.git
```
It will be cloned to the `amrex` directory.

Then clone the GRTeclyn repository with a command such as

```bash
git clone -b training/202406_grtl_meeting https://github.com/GRTLCollaboration/GRTeclyn.git
```
Note that the above command automatically checks out the
`training/202406_grtl_meeting` branch but if you missed out this flag make
sure to check it out using e.g.
```bash
git checkout training/202406_grtl_meeting
```

The solution is in the branch `enhancement/scalar_fields`. 

> **Note**
> I have assumed that you have cloned both of
> these repositories to the same directory so that the `amrex` and `GRTeclyn`
> directories share the same parent directory. If you want to clone AMReX
> elsewhere, make sure to set the `AMREX_HOME` environment variable
> appropriately e.g. 
> ```bash
> export AMREX_HOME=/path/to/amrex
> ```

### Building the Scalar Field example

> **Warning** 
> If you are on macOS and want to use GCC built with Homebrew (which
> has the major version number appended to the executable e.g. `g++-9`), you
> will need to create a `Make.local` (like Chombo's `Make.defs.local` file) in 
> ```
> amrex/Tools/GNUMake/Make.local
> ```
> with something like the following lines
> ```makefile
> CXX = g++-9
> CC  = gcc-9
> FC  = gfortran-9
> F90 = gfortran-9
> 
> INCLUDE_LOCATIONS += /usr/local/include
> ```
> 
> See [this
> page](https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#gcc-on-macos)
> for more details.  

On a CSD3 login node, request an interactive session using:

```
salloc -A training-dawn-gpu -p pvc --nodes=1 --ntasks=<#ntasks> --time=01:00:00
```

This is example is quite small so 1 or 2 tasks (usually MPI ranks) is
sufficient. If you want to use the GPUs add `--gres=gpu:<#gpu>` which
will add #gpu to your allocation. Note that you can only request a
maximum of 4 per node. Please see the Dawn notes
[here](https://docs.hpc.cam.ac.uk/hpc/user-guide/pvc.html) and also
Miren's guide for running on the PVCs
[here](https://github.com/GRTLCollaboration/GRTeclyn/issues/67)

Some modules will be loaded by default, but just to be sure: 
 
```
module purge
module load default-dawn
module load intel-oneapi-compilers/2024.1.0/gcc/wadpqv2p
module load intel-oneapi-mpi/2021.12.0/oneapi/nbxgtwyb
module load intel-oneapi-mkl/2024.1.0/oneapi/xps7uyz6
```

which provide you with the recommended versions of Intel compilers and
MPI.


The main build options are located in
```
~/GRTeclyn/Tools/GNUMake/Make.defaults
```

In this file, you can adjust several variables such as:

* Set `USE_OMP = TRUE` to use OpenMP.
* If you don't want to use the GNU compiler `g++` , change `COMP = gnu` to 
  `COMP = intel` (for the Intel classic compiler `icpc`) or `COMP = llvm` (for 
  LLVM `clang++`). 
* If you don't have an MPI implementation available, set `USE_MPI = FALSE`.
* Optionally, set `USE_SYCL = TRUE` to compile for the PVCs. You can also just use the CPUs though. 

With the above modules loaded, use `COMP=intel-llvm` for Dawn. 


Like Chombo, there is also a file controlling the AMReX build options:

```
~/amrex/Tools/GNUMake/Make.defs
```

but you shouldn't have to set these for this exercise. If you have
some environment variables that you would preferrentially like to
define, you can do so in `~/amrex/Tools/GNUMake/Make.local-pre` and
this will be sourced by `Make.defs`.


Now, navigate to the Scalar Field example directory
```bash
cd ~/GRTeclyn/Examples/ScalarField
```

Source directories and includes are defined in the `GNUMakefile`. You
shouldn't have to modify this file for this exercise.



Now start building AMReX and the Scalar Field example with the command
```bash
make -j 4 
```
Note that you can add build options here too e.g. `make -j 4 USE_MPI=FALSE` to turn off MPI.

A new `tmp_build_dir` directory will be created to store the compiled object
and auxiliary files. Assuming all is well, you should have an executable in the
current directory of the form `main<config>.ex` e.g. `main3d.intel-llvm.MPI.OMP.ex`.

## Running and visualizing the example

### Running the example

There is an example parameter file called: `params.txt`. You must fill
in the output directories with your own! If you open the file, you'll
find that many of the parameters have been inherited from GRChombo and
thus their names and functions are the same. Note that not all
parameters have been implemented (e.g. `output_path`). Consult the
table in the [AMReX renamed parameters](#amrex-renamed-parameters)
section below.

Running the example is the same as for GRChombo examples e.g. (with MPI):
```bash
export OMP_NUM_THREADS=2
mpiexec -n 4 ./main3d.intel-llvm.MPI.OMP.ex ./params.txt
```

There's no need to resubmit this as another job if you're already on Dawn compute node for compilation.
Otherwise, you can use the following command: 
> ```bash
> srun -p <partition> -t 5 -n 4 -c 2 -A <account> --export=ALL,OMP_NUM_THREADS=2 ./main3d.intel-llvm.MPI.OMP.ex ./params.txt
> ```

With the provided parameters, it should evolve 3 timesteps and write 1 plot
file (actually directories/folders rather than single files, each about 26 MB):
```
plt00000  plt00001 ...
```
This should not take longer than a few minutes on a laptop.

### AMReX new/renamed parameters
Below are some parameters you will find in the Scalar Field `params.txt` that
are either new or differ in name to the GRChombo equivalent. Most of their names
are self-explanatory.
| GRChombo parameter name | GRTeclyn parameter name | Notes |
| --- | --- | --- |
| `chk_prefix` | `amr.check_file` | Still a prefix in GRTeclyn |
| `plot_prefix` | `amr.plot_file` | Still a prefix in GRTeclyn |
| `restart_file` | `amr.restart` | |
| - | `amr.file_name_digits` | Number of timestep digits in the filenames. |
| `checkpoint_interval` | `amr.check_int` | |
| `plot_interval` | `amr.plot_int` | |

### Visualizing the plot files
The easiest way to visualize the plot files on Dawn is to use the `fsnapshot`
tool suppiled by AMReX. Go to `~/amrex/Tools/Plotfile` and type `make`.
You can make a plots of slices through the simulation volume using:
```
fsnapshot.intel-llvm.ex -v phi -p Palette plt_xxxxxx
```
where `phi` is the value of the field being plotted. 
Please refer to the guide [here](https://amrex-codes.github.io/amrex/docs_html/Post_Processing.html#fsnapshot). 

The `pltxxxxx` directories can be opened using ParaView (>=5.7), VisIt or yt.
Note that, unlike GRChombo, checkpoints cannot be visualized.

If you are using ParaView, open the group of `plt...` directories and then
select "AMReX/BoxLib Grid Reader". 

If you are using VisIt, open the `Header` file in one of the `pltxxxxx`
directories.

## Code structure

### Class structure

The structure of the code is very similar to GRChombo. Here is a diagram showing
the inheritance/aggregation relationship between the main classes for this
example in GRTeclyn:

```mermaid
graph BT
    B[AmrCore] --> C[AmrMesh]
    A[Amr] --> B[AmrCore]
    L1[AmrLevel] -.-> A[Amr]
    L2[AmrLevel] -.-> A
    L3[AmrLevel] -.-> A
    LG1[GRAMRLevel] --> L1
    LG2[GRAMRLevel] --> L2
    LG3[GRAMRLevel] --> L3
    LS1[ScalarFieldLevel] --> LG1
    LS2[ScalarFieldLevel] --> LG2
    LS3[ScalarFieldLevel] --> LG3
    SD1[StateData] -.-> LS1
    SD2[StateData] -.-> LS2
    SD3[StateData] -.-> LS3
    MF1["MultiFab (new data)"] -.-> SD1
    MF2["MultiFab (old data)"] -.-> SD1
```

### [GR]Chombo to GRTeclyn/AMReX class dictionary

Here are a list of some of the main classes in [GR]Chombo and their equivalents
in GRTeclyn/AMReX:

| [GR]Chombo class | GRTeclyn/AMReX class  | Description |
| ---           | ---           | ---         |
| `AMR`         | `Amr`         | Class which manages the whole AMR hierarchy and knows about all of the levels |
| `AMRLevel`    | `AmrLevel`    | Class which manages things on a single level |
| `GRLevelData` | `MultiFab`    | Class which stores the actual array data on a single level |
| `FArrayBox`   | `FArrayBox`   | Class which stores the actual array data in a single box |
| - | `StateData` | Class which stores the new and old `MultiFab`, boundary conditions, etc. In Chombo, the `GRLevelData` objects are just stored directly in `GRAMRLevel` |


### Further information

More information about AMReX can be found in its
[documentation](https://amrex-codes.github.io/amrex/docs_html/).


