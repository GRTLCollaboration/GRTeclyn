# Example configs for specific HPC systems

In general you should start with an example hybrid MPI/OpenMP or GPU jobscript from your cluster documentation, and adapt it for your GRChombo run.

However, some example jobs for systems we often use are collected here, along with the `make.local-pre` setup and some helpful tips - to be copied at your own risk!

---

## Mac laptop

This is really for only for debugging. 

Install command line tools on your Mac using the command `xcode-select --install` in the terminal.

In the `make.local-pre` file:
```
COMP = llvm
USE_OMP = FALSE
USE_MPI = FALSE
```

then to run `./main3d.llvm.ex params.txt`.

---

## Cosma8 (CPUs)

The modules used
```
module load gnu_comp openmpi
```
Or you can also use the Intel modules for better performance: 

```
module load intel_comp/2025.1.1 compiler-rt/latest tbb/latest  compiler/latest mpi/latest
```
NB: You will need to set the option `--mpi=pmi2` in your Slurm script to launch Intel MPI correctly. 

In the `make.local-pre` file we have
```
COMP = gnu
USE_OMP = TRUE

# Options you might want to turn on
#TINY_PROFILE=TRUE
#PROFILE=FALSE
#USE_MPI = FALSE
#USE_HDF5 = TRUE

```
Note that `TINY_PROFILE=TRUE` and `PROFILE=FALSE` are a _pair_ of options, that is to use `TINY_PROFILER`, `PROFILE` must be set to `FALSE`. 

And the jobscript:

```
#!/bin/bash -l

# number of nodes
#SBATCH --nodes 4
#! How many tasks
#SBATCH --ntasks-per-node=64
#! How many CPUs per tasks should be allocated (for OpenMP threading)
#! The product of this and ntasks-per-node should be 128
#SBATCH --cpus-per-task=2
#SBATCH -J GRTLTest1
#SBATCH -o output.out
#SBATCH -e error.err
#SBATCH -p cosma8
#SBATCH -A dp0xx
#SBATCH --exclusive
#SBATCH -t 6:00:00
#SBATCH --mail-type=ALL                          # notifications for job done & fail
#SBATCH --mail-user=xxx@university.ac.uk

# load the required modules
module load gnu_comp openmpi

#!Print info
module list
pwd
date

#! Are you using OpenMP?
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#! Number of nodes and tasks per node allocated by SLURM (do not change):
mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

#! Number of nodes and tasks per node allocated by SLURM (do not change):
numnodes=$SLURM_JOB_NUM_NODES
numtasks=$SLURM_NTASKS

#! Number of MPI tasks to be started by the application per node and in total (do not change):
np=$[${numnodes}*${mpi_tasks_per_node}]

#! Full path to application executable: 
application="./main3d.gnu.MPI.OMP.ex"

#! Run options for the application:
options="params.txt"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"

# Run the program

mpirun --map-by ppr:$mpi_tasks_per_node:node -np $SLURM_NTASKS $application $options
```

---

## Cosma8 (GPUs)

### AMD MI300X partition

The modules used
```
module load hipcc/6.3amd openmpi/5.0.3
```

In the `make.local-pre` file we have
```
COMP = gnu
AMREX_USE_GPU=TRUE
USE_HIP=TRUE         
# for MI300                                                 
AMREX_AMD_ARCH=gfx942
# Option for MI200 - need to switch if using this
# AMREX_AMD_ARCH=gfx90a
```

And the jobscript (uses all 8 GPUs to avoid overlap with other users, but be considerate!):

```
#!/bin/bash
#SBATCH -J grteclyn-bbh
#SBATCH -A do018
#SBATCH -p mi300x
#SBATCH -N 1 # there is only one CPU node running the GPUs
#SBATCH -n 8 # one MPI process per GPU (8 GPUs)
#SBATCH -c 12
#SBATCH --exclusive
#SBATCH -t 2:00:00
                                                                
module load hipcc/6.3amd openmpi/5.0.3

exec="./main3d.hip.MPI.HIP.ex"
args="params.txt"

mpirun -n ${SLURM_NTASKS} ${exec} ${args}
```

---

## Tursa

Tursa is primarily a GPU system with Nvidia A100s as the GPU but there are CPU queues available as well. We will focus on the GPU part for this wiki. 

No modules are loaded by default and you'll notice that the system wide modules available are a bit sparse. These are the modules that Juliana uses for GRTeclyn/AMReX workloads:

```
module load gcc/9.3.0 ucx/1.15.0-cuda12.3 openmpi/4.1.5-gcc9-cuda12 ucx/1.15.0-gcc9-cuda12 cuda/12.3
```

You can save these modules for later reuse with `module save <name of module environment>`
Then every time you log back into Tursa, you can load the same environment using `module restore <name of module environment>`

### Building GRTeclyn
Once the above modules have been loaded, you can compile the code in your working directory (for example, `Examples/BinaryBH`) using:
`make -j 4 USE_CUDA=TRUE CUDA_ARCH=80`

`CUDA_ARCH` is necessary because the code will generate a binary only valid for that type of Nvidia architecture, which are grouped by generation. For example: A100s are 80, V100s are 70 and Grace Hoppers are 90. Tursa exclusively contains A100s. 

### Running GRTeclyn
An example of a single GPU run: 
```
srun -A <your project code> -p gpu --qos=<priority> --nodes=1 --ntasks=1 --gres=gpu:1 --time=00:10:00 ./main3d.gnu.MPI.CUDA.ex params_test.txt
```
where `<your project code>` is the name of your allocation and `<priority>` sets your priority level. Choose from `standard`, `dev`, `low`, or `high`. Read more about the possible options [here](https://epcced.github.io/dirac-docs/tursa-user-guide/scheduler/)

NB: The `gpu` queue above is general. There are two types of A100s on Tursa, those with 40GB memory per card and those with 80GB memory per card. If you need a particular amount of memory, you can choose between them by setting either `-p gpu-a100-40` or `-p gpu-a100-80`. Using `-p gpu` will give you whichever is available first and so increases your chances of the job starting earlier. 

---

## Swirles

Swirles is a local cluster hosted by DAMTP, funded with the support of the ExCALIBUR programme and so open to all UK researchers. It contains a mix of Nvidia and Intel GPUs (codenamed Ponte Vecchio or PVCs). 

To get started, git clone the repo located here: 
```
git clone https://github.com/COSMOS-CTC-Cambridge/swirles-training.git
```
In it, you will find introductory notes on the system and most importantly example scripts to set up an environment for using the GPUs and for submitting to the Slurm queue.

For example, to build GRTeclyn for the Intel GPUs: 

1. Start an interactive session from the head node:
```
srun -p pvc --nodes=1 --ntasks=28 --gres=gpu:1 --time=1:00:00 --pty bash
``` 
This will give you a quarter of the node and 1 GPU for 1 hour (each node has 4 GPUs and 112 cores spread over 2 sockets). If you want your own node add `--exclusive` and change to `--gres=gpu:4`. Do not build GRTeclyn on the head node if you want to run on the Intel GPUs - it has a different architecture. 

2. Source the environment from the git repo `swirles-training`:
```
source ~/swirles-training/Environments/pvc-env.sh
``` 
This will set up the relevant modules as well as some environment variables necessary for using Intel MPI on the system. It will also expose the Intel GPUs as separate tiles. 

3. Build the Binary BH example:
```
make -j 28 USE_SYCL=TRUE
```
You may also like to set some ancillary AMReX variables to help Intel GPU performance, such as `SYCL_SUB_GROUB_SIZE=16`. For more information, see [this](https://amrex-codes.github.io/amrex/docs_html/GPU.html#sycl-configuration-variables) part of the AMReX documentation and [this](https://github.com/GRTLCollaboration/GRTeclyn/issues/67) GitHub issue. 

4. Run the example:
```
mpirun -np 2 ./main3d.sycl.MPI.ex ./params_test.txt
```
Notice that I am only running with 2 MPI ranks even though I have access to 28 cores (if submitted exactly as above). This is because AMReX is heavily optimized for 1 MPI rank to 1 GPU. In this case, it is running 2 MPI ranks for the 2 tiles on a single Intel GPU card (again only possible if the environment variables have been set up as above). 

Example Slurm submission scripts are also available in the `swirles-training` repository. 

NB: To install the pre-commit hook for GitHub, you will need to load the newer version of Python on the system: 
```
module load python/3.10.13/gcc/t6nyisbb
```
then setup a virtual environment
```
python -m venv .precommit
source ~/.precommit/bin/activate
```
then `pip install` the package:
```
python -m pip install pre-commit
```

You can then move it to a more convenient location once installed e.g. `~/.local/bin`

## CSD3

...