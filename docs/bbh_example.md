# Running the binary black hole example

The following gives some instructions on running the binary black hole example.

---

NOTE: THE FOLLOWING EXAMPLE HAS NOT YET BEEN PORTED... UNDER CONSTRUCTION...

---

## Physical scenario

This page describes running the Binary BH example using the parameters found in [this parameter file](https://github.com/GRTLCollaboration/GRTeclyn/blob/develop/Examples/BinaryBH/params.txt).

The initial conditions consist of a superposition of two equal mass black holes, with initial momenta and positions chosen to give approximately 10 circular orbits before merger.

Note that finding the right momenta for an initial circular trajectory is in general a hard problem in NR - simply setting the momentum to the Newtonian approximation will result in elliptic orbits even for well separated initial BHs. One must use the Post-Newtonian approximations, and then adjust the momentum manually over several iterations to achieve an initial eccentricity of less than 1%. We are grateful to Dr Sebastian Khan at Cardiff University for providing us with these initial data values.

Note that this is a realistic example but not production quality for generating BBH waveforms (one should probably have a larger wave zone, and better resolution at the horizon).

The initial data solves the Momentum constraint exactly but uses an approximation for low boosts to solve the Hamiltonian constraint (see [Baumgarte and Shapiro](https://doi.org/10.1017/CBO9781139193344), pg 73-74) which is accurate up to order P⁴, where P is the momentum of the BH. One could also use our [[Two Punctures Initial Data | TwoPunctures Initial Data]] extension (best for vacuum binaries) or [GRTresna](https://github.com/GRTLCollaboration/GRTresna) (for binaries with additional field content) to solve the Hamiltonian constraint if higher accuracy is required.

## Computational set up

See the general comments on compiling and running examples in [[Running examples | Running examples]].

For this specific example, with `max_box_size = min_box_size = 16` (set in `params.txt`) there will be approximately 800 boxes on the finest levels (which are the most important for the load balancing - see [[Performance optimisation | Performance optimisation]]), so for good load balancing, one could run the example with 200 ranks on 200 cores with 4 boxes per core. However, normally the memory per rank will not be sufficient, and the job will either run out of memory and crash, or (worse) just run *very* slowly. In this case better performance can be obtained by using 200 ranks x 4 openMP threads = 800 cores. The OpenMP threads then share the work in each box, providing a further speed up, whilst having more cores per rank provides more memory. On the [CSD3 HPC system in Cambridge](https://www.hpc.cam.ac.uk/high-performance-computing) one achieves roughly 50 M/hr at t = 10M in this configuration, setting `OPT = HIGH` (which turns off checking of assertion statements), `DEBUG = FALSE` and using AVX512 vectorisation. Since the stop time is t = 2200M, this implies a run time of just less than 2 days.
(If fewer cores are available, an alternative is to set `max_box_size = min_box_size = 32` and run on 96 ranks x 4 openMP threads = 384 cores, which gives roughly two boxes per rank, and a speed of approx 25 M/hr.)

Here is an example jobscript for the CSD3 machine (see your cluster docs for a template for your own system):

```bash

    #! Name of the job:
    #SBATCH -J BBH
    #! Which project should be charged:
    #SBATCH -A <account name>
    #! How many whole nodes should be allocated?
    #SBATCH --nodes=25
    #! How many (MPI) tasks will there be in total? (<= nodes*32)
    #! The skylake nodes have 32 CPUs (cores) each.
    #SBATCH --ntasks=200
    #! How many CPUs per tasks should be allocated (for OpenMP threading)
    #SBATCH --cpus-per-task=4
    #! How much wallclock time will be required?
    #SBATCH --time=1-12:00:00
    #! For 6GB per CPU, set "-p skylake"
    #SBATCH -p skylake

    #! Number of nodes and tasks per node allocated by SLURM (do not change):
    mpi_tasks_per_node=$(echo "$SLURM_TASKS_PER_NODE" | sed -e  's/^\([0-9][0-9]*\).*$/\1/')

    #! Full path to application executable:
    application="Main_BinaryBH3d.Linux.64.mpicxx.ifort.OPTHIGH.MPI.OPENMPCC.ex"
    #! Run options for the application:
    options="params.txt"

    #! Are you using OpenMP?
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

    #! The following variables define a sensible pinning strategy for Intel MPI tasks -
    #! this should be suitable for both pure MPI and hybrid MPI/OpenMP jobs:
    export I_MPI_PIN_DOMAIN=omp:compact # Domains are $OMP_NUM_THREADS cores in size
    export I_MPI_PIN_ORDER=scatter # Adjacent domains have minimal sharing of caches/sockets

    CMD="mpirun -ppn $mpi_tasks_per_node -np $SLURM_NTASKS $application $options"

    eval $CMD
```


## Checking the outputs

### Viewing data in VisIt

See [[Visualising Ouputs | Visualising Outputs]] for further details on visualising with VisIt.

The most enlightening variable to look at in the outputs is chi, the conformal factor of the metric, which goes to zero at the centres of the BHs. When the Pseudocolour plot is viewed on a slice through the centre of the grid, normal to the z axis, the black holes should be seen to inspiral and merge.

Note that you can look at these files even when the job is still running. It is always a good idea to check the files after a few timesteps.

### Puncture plots

The example also outputs some ASCII datafiles for post processing. Firstly, it gives a file `BinaryBHChk_Punctures.dat` which gives the location at each timestep of the two BHs.

This can be plotted using a `gnuplot` command:

`f='BinaryBHChk_Punctures.dat'; L=256; plot f using ($2-L):($3-L) with lines, f using ($5-L):($6-L) with lines`

or a `python` script such as:

```python

# A simple python script to plot the puncture
# tracks over time. This helps with setting up
# circular orbits. The params.txt file should
# give around 8-9 orbits before merger.

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# location of center of grid
L = 256
# half the separation of punctures
r = 6
# output data from running merger
data = np.loadtxt("BinaryBHChk_Punctures.dat")

# make the plot
fig = plt.figure()

# first puncture
x1 = data[:,1]-L
y1 = data[:,2]-L
plt.plot(x1,y1)

# second puncture
x2 = data[:,4]-L
y2 = data[:,5]-L
plt.plot(x2,y2)

# make the plot look nice
plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.ylim(-r-1,r+1)

# save as png image
plt.savefig("PunctureTracks.png")
```

The resulting image should look like this one:

![](https://github.com/GRTLCollaboration/GRChombo/wiki/files/BBHPunctures.png)

### Gravitational Wave plots

The other data files starting with `Weyl_integral_` give the spin weight -2 spherical harmonic components of the Weyl Scalar Ψ₄ as calculated at each of the selected extraction radii. The aim of these files is to check that the signal has converged and so is similar at the different extraction radii, and also that one obtains a correct signal for a binary merger.

Again, these can be plotted using a `gnuplot` command, for example:

`f='Weyl_integral_22.dat'; plot f using ($1-50):2 with lines, f using ($1-100):4 with lines`

or a `python` script such as:

```python

# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

# coord locations of extraction radii
R0 = 50
R1 = 100
# Total ADM mass
M = 1.0
# The mode, as text
mode = "22"
# output data from running merger
data = np.loadtxt("Weyl_integral_" + mode + ".dat")

# make the plot
fig = plt.figure()

# first radius
r0 = R0 + M*np.log(R0/(2.0*M) - 1.0)
timedata0 = (data[:,0] - r0) / M
fluxdata0 = data[:,1]
plt.plot(timedata0, fluxdata0, ':', lw = 0.5, label="r0")

# second radius
r1 = R1 + M*np.log(R1/(2.0*M) - 1.0)
timedata1 = (data[:,0] - r1) / M
fluxdata1 = data[:,3]
plt.plot(timedata1, fluxdata1, '-', lw = 0.75, label="r1")

# make the plot look nice
plt.xlabel("time t / M")
plt.ylabel("Re (Psi4) el, em = " + mode)
plt.xlim(0, 2200)
plt.ylim(-0.061, 0.061)

# save as png image
filename = "Weyl_" + mode + ".png"
plt.savefig(filename)
```

The resulting image should look like this one:

![](https://github.com/GRTLCollaboration/GRChombo/wiki/files/WeylScalar22.png)

### What next?

Congratulations - you have successfully run a binary black hole simulation with GRTeclyn! We suggest that you try amending some of the variables in params.txt and rerun the code to understand what they do (see the [[Guide to parameters | Guide to parameters]] for helpful info).