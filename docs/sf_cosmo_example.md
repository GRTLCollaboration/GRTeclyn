# Running the scalar field example

---

NOTE: THE FOLLOWING EXAMPLE HAS NOT YET BEEN PORTED... UNDER CONSTRUCTION...

---

Here we explain the cosmological scalar field example found in the code.

## Physical scenario

This page describes running the Cosmo example using the parameters found in [this parameter file](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/params_cheap.txt).

The example is not very physically motivated but provides a good starting point for developing cosmological examples with some useful examples of how to implement diagnostics. It simulates an expanding spacetime with a scalar field that has a sinusoidal dependence in the x direction and a quadratic potential. The mass of scalar field is chosen to be equal to the wave number such that the constraints are satisfied initially (in a more general case one would need to solve the constraints numerically).

### Initial data

In this example, we are choosing the scalar field profile (see, [`InitialScalarData.hpp`](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/InitialScalarData.hpp)) as

$$ A \sin(2 n \pi x/L) $$

where $A$, $n$ and $L$ are the amplitude, mode and length of the simulation box, respectively.

With the quadratic potential (see [`Potential.hpp`](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/Potential.hpp))

$$ V(\phi) = \frac{1}{2} m^2 \phi^2 $$

and the choice where $n = m$, the energy density,

$$ \rho = \frac{1}{2} (\partial_x \phi)^2 + V(\phi) = \frac{A^2 m^2}{2} $$

This is spatially constant, which makes it easier to make sure that the constraints are initially satisfied.

The example is set up by solving $A$ from Friedmann's equation in an FLRW spacetime

$$ H^2 = \frac{8  \pi G \rho}{3} = \frac{8  \pi G}{3} \left(\frac{A^2 m^2}{2}\right)  $$

where $n=1$, the Hubble length $L_H \sim 1/H = 1$.

The initial data for extrinsic curvature $K$ that satisfies the Hamiltonian constraint is calculated from this (note that $K=-3H$ and stored on the grid as can be seen in `InitialK.hpp`. Note that we also assign values of initial `Pi = 0`, `lapse = chi = h11 = h22 = h33 = 1` in `InitialScalarData.hpp`.

### Diagnostics

All diagnostics for this particular cosmological simulation are specified in [`CosmoDiagnostics.hpp`](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/CosmoDiagnostics.hpp). As a guide of how a lineout could be done (extracting data on points along a line), the example also shows how to extract values of energy density along the x-axis of the grid.

A class called `CustomExtraction` has been added to do the extraction of the values of components at specified points using the `AMRInterpolator` which interpolates data at any coordinates on the grid, using data from the finest available AMR level (see [`CustomExtraction.hpp`](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/CustomExtraction.hpp)). The given class is written to extract the data along x-axis and output it to `lineout.dat`. One could do a fancier lineout rather than just along an axis by modifying the class. It is possible to specify a starting point of the extraction in `CosmoLevel.cpp` at `specific_post_timestep()`. The number of extracted points, `lineout_num_points`, can be specified in `params.txt`.

### Regridding

Note that with the default parameters regriding should not happen, but an example of a possible tagging criteria is already provided with this example at [`Source/TaggingCriteria/CosmoHamTaggingCriterion.hpp`](https://github.com/GRTLCollaboration/GRChombo/blob/main/Source/TaggingCriteria/CosmoHamTaggingCriterion.hpp). This tagging criterion tags the value of the sum of the magnitude of every term in the Hamiltonian constraint - this is usually a good choice to focus regridding in areas of strong curvature or high density. If one wants to turn the regridding on, one could change relevant parameter(s) under the Grid parameters section in params file, or modify this tagging criterion. In this example, the density and curvature are quite uniform, and so it doesn't make sense to regrid any particular areas.

## Running the example

Go to `Examples/ScalarFieldCosmo` and type:

`make all`

If the example compiles with no errors, there should be an executable will have been produced with a name which depends on your compiler options but will look something like:

`Main_ScalarFieldCosmo3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.COSMA8.Intel2022.ex`

Normally you will not want to output the resulting hdf5 files to your home directory (where you would normally store the code), so you should uncomment and amend the following line in the params.txt file to give the (full, not relative) path to your data directory, where large output files are stored, e.g.

```
output_path = "/full/path/to/your/data/directory/"
```

To run the example, you are required to submit a batch job so that the simulation is queued and run on the compute nodes. The jobscript should be written and submitted to the queue. There should be examples for MPI and hybrid OpenMP/MPI jobs in your cluster documentation. Note that you may also need to reload the module files within the jobscript, depending on how your cluster works. Some examples of systems we use can be found [here](https://github.com/GRChombo/GRChombo/wiki/Jobscript-tips-and-examples). The number of cores you need will depend on the job and the memory per core on your cluster. If you run out of memory, you may want to switch to using several OpenMP threads per task, rather than running a pure MPI job.

Here is an example jobscript for the Cosma8 machine (see your cluster docs for a template for your own system):

```bash
#!/bin/bash -l

#! How many nodes
#SBATCH --nodes 1
### NB cosma8 has 128 cores per node so product of these = 128
#SBATCH --ntasks-per-node=64 ## Total tasks will be this * number of nodes
#SBATCH --cpus-per-task=2 ## Equal to number of OMP threads per task (set below)
#! How many CPUs per tasks should be allocated (for OpenMP threading)
#SBATCH -J CosmoEx_cheap
#SBATCH -o out_%J.out
#SBATCH -e err_%J.err
#SBATCH -p cosma8
#SBATCH -A <account name>
#SBATCH --exclusive
#SBATCH -t 01:00:00
#SBATCH --mail-type=NONE                          # notifications for job done & fail
#SBATCH --mail-user=NONE

module purge
## Cosma 8
module load gnu_comp/7.3.0
module load intel_comp/2022.1.2
module load compiler/2022.0.2 mkl/2022.0.2 mpi/2021.5.1
module load parallel_hdf5/1.12.0

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
application="Main_ScalarFieldCosmo3d_ch.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.COSMA8.Intel2022.ex"
#! Run options for the application:
options="params_cheap.txt"

#! Work directory (i.e. where the job will run):
workdir="$SLURM_SUBMIT_DIR"

# Run the program

mpirun -ppn $mpi_tasks_per_node -np $SLURM_NTASKS $application $options
```

The job can be submitted with a command which depends on the scheduler, again for example:

`sbatch jobscript.sh`

## Checking the outputs

### Viewing data in VisIt

The example should produce a series of hdf5 files, which you can then view to check what has been simulated (see [[Visualising outputs|Visualising outputs]] for further details).

One can check the evolution of chi, the conformal factor of the metric, which can be converted to e-folding number $\text{ln}(a) = \text{ln}\left(\frac{1}{\sqrt{\chi}}\right)$.

### Plots from data files

Some simple plots can also be generated using Python script provided in the directory via the command

`$ python3 plot.py `

If one runs the example with `parmas_cheap.txt`, the resulting plots should be as shown below.

The first graph ([plot_efold.pdf](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/plot_efold.pdf)) is a plot of the e-folding number over time in the code unit. The plot shows that the example simulates a decelerated expanding universe.

![plot_efold](https://github.com/user-attachments/assets/b49d27e4-5403-40ba-a558-15d31139bbd7)

The second plot ([plot_rho_mean.pdf](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/plot_rho_mean.pdf)) is the mean energy density across the grid as a function of e-folding number.

![plot_rho_mean](https://github.com/user-attachments/assets/4973281f-4b70-4238-a767-cdae01622473)

The last figure ([plot_lineout.pdf](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/plot_lineout.pdf)) is a lineout extraction of energy density along x-axis of the simulation grid. The plot extracts 10 data points using AMRInterpolator to interpolate the data.

![plot_lineout](https://github.com/user-attachments/assets/9f1a382e-7b57-4dfa-9547-faa15a6ff5fd)

## Load balancing
Note that your job will output files like pout.0 for each processor. Looking at these files during the run allows you to see how fast it is going. More output can be produced by setting

    verbose = 1

in the params.txt file. These files will also tell you how many boxes are being processed by each processor - for good load balancing and to avoid using more processors than you need, ideally all the processors should have one (or several) boxes each, and there should be at least one per processor.

**Note that Chombo does not adjust the box sizes automatically based on the number of processes you choose so you almost certainly need to adjust them for your specific application and job configuration to avoid wasting resources and/or achieve good speeds.**

More information on load balancing is given in the section on [performance optimisation](https://github.com/GRChombo/GRChombo/wiki/Performance-optimisation) in the main GRChombo wiki. This is highly recommended reading.

### Expected performance
The following examples are taken from the jobs running on the [Cosma8](https://dirac.ac.uk/memory-intensive-durham) system.

- Running with [params_cheap.txt](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/params_cheap.txt)

This params file specifies the grid size as $16^3$ with `max_` and `min_box_size = 8`. Therefore, the grid is divided into $(16/8)^3 = 8$ boxes. Running this on 8 processes we will have one box per rank, and if we check `pout.0`, we should have a running speed of around 1000 M/hr. With this speed, the job will finish within a couple of minutes - this is not a realistic job size but provides a quick test that everything is working ok.

- Running with [params.txt](https://github.com/GRTLCollaboration/GRChombo/blob/main/Examples/ScalarFieldCosmo/params.txt)

This params file specifies the grid size as $128^3$ with `max_` and `min_box_size = 16`. Therefore, the grid is divided into $(128/16)^3 = 512$ boxes. Running this on 512 processes will assign one box to each rank. If we check `pout.0`, we should find a running speed of around 100 M/hr. With this speed, the job will finish within an hour. Using half the number of processes would mean there would be twice as many boxes per rank, and so we would expect the speed to roughly halve too (although this is only approximate).

It is well worth checking that you can get these expected speeds on any new system in order to check that things are working properly. Cosma8 is a very fast system, but most systems should give the same order of magnitude. If you are getting speeds of 1-10M/hr with a similar number of processes, something is wrong (probably a memory bottleneck, or a bad compiler option).

***
## Restarting from GRTresna - running/restarting from the initial data solver

If one wants to use [`GRTresna`](https://github.com/GRTLCollaboration/GRTresna) to solve the initial data, please see GRTresna's [wiki page](https://github.com/GRTLCollaboration/GRTresna/wiki). Once the solver successfully runs and outputs the final hdf5 file (usually named `InitialDataFinal.3d.hdf5`), the example can then be run starting from this initial data by using `InitialDataFinal.3d.hdf5` as a checkpoint file to restart the simulation from. To do so, make sure that all relevant parameters are consistent between GRChombo and GRTresna, uncomment and specify the path to the checkpoint file in the params file

```
# If restarting from the IC Solver or checkpoint, uncomment this and put the
# file into the hdf5 folder (need to create this if not there)
restart_file = /path/to/your/InitialDataFinal.3d.hdf5/
```

Params files and results of the example using GRTresna can be found in [`/ScalarFieldCosmo/GRTresna_restart_setup`](https://github.com/GRTLCollaboration/GRChombo/tree/main/Examples/ScalarFieldCosmo/GRTresna_restart_setup).