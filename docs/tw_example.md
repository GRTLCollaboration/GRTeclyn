# Running the Teukolsky wave example

This page provides some background information for the Teukolsky wave example found in the code, as well as some example output. The example is relatively cheap, and uses only 1 coarse grid without refinement.

TO DO: provide a measure for this! Is this runnable on a laptop?

## Physical Scenario

This example evolves Teukolsky waves [Teukolsky:1982nz](https://inspirehep.net/literature/182784); a linearized quadrupolar wave in the TT gauge. In particular, the initial data follows the Eppley construction [Eppley:1979](https://ui.adsabs.harvard.edu/abs/1979sgrr.work..275E/abstract) which constructs a superposition of an in- and outgoing wave such that the initial data is at a moment of time-symmetry and the momentum constraint is trivially satisfied. More precisely, we employ the seed function in [Hilditch:2013cba](https://inspirehep.net/literature/1254684), which generalizes this construction to start the wave packets at a distance $r_0$ away from the origin.

## Computational set up

### Problem-specific parameters

The seed function is controlled by an amplitude $A$ and a width $\sigma$.

The parameters in `params.txt` are such that `regularize_r` $\ll \sigma \ll r_0$. Note that `regularize_r` is a small parameter that regularizes the initial data near the origin (see [Hilditch:2015aba](https://inspirehep.net/literature/1362166) for a comment on the seed function used here). The example has not been tested for parameters that do not respect the above order.

### Grid symmetries

Currently, the example has implemented three kinds of Teukolsky wave: even parity with $M=0,2$ or odd parity with $M=2$. The former can be run with octant symmetry, but the latter can be evolved with quadrant symmetry at best (reflection symmetry around the $xy$-plane is broken).

## Output

### Viewing 3D data and 2D slices

See [**Visualising Ouputs**](visualising_outputs.md) for details on visualising. Note that checkpoint files are not viewable with AMReX, only plot files are.

You probably want to look at the real and imaginary parts of the Weyl4 scalar, `Re_Weyl4` and `Im_Weyl4`, which should display wavefronts.

The Hamiltonian constraint `Ham` is useful to visualise as well, especially at $t=0$: if the amplitude $A$ is sufficiently large, `Ham` scales quadratically with $A$ as it is dominated by non-linear contributions that do not solve Einstein's equations. For small $A$, discretization errors ruin this scaling.

### Plots from data files

The example incorporates the extraction of the Weyl scalar (following e.g. the BinaryBH example), and hence gravitational waves can be visualised. For the $M=0$ ($M=2$) example, the Weyl components $\psi_{20}$ ($\psi_{22}$ and $\psi_{2-2}$) show two bursts, separated by a time delay $\sim 2r_0$, representing the outgoing wave followed by the ingoing wave, after it has passed through the origin.

TO DO: add plot