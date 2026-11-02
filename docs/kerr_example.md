# Running the Kerr black hole example

Here we explain the Kerr black hole example found in the code, which evolves a single rotating black hole.

The example provides a good starting point for simulations of spinning black holes. The Kerr spacetime is stationary, so nothing physical should change during the evolution: the black hole should stay at the centre of the grid, and the evolution comes from the gauge adjusting (e.g. the lapse starts from a pre-collapsed value rather than its stationary one) and from numerical error. This also makes it a cheaper test of black hole evolution than the [**Binary BH example**](bbh_example.md).

## Physical scenario

This page describes running the KerrBH example using the parameters found in [this parameter file](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/KerrBH/params.txt).

### Initial data

The initial data is set in [`KerrBHInitialData.hpp`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/BlackHoles/KerrBHInitialData.hpp). Note that this lives in `Source/BlackHoles` rather than in the example folder, so that it can be reused by other examples.

We use the Kerr metric in puncture form with the radial coordinate $r$ introduced by [Liu, Etienne and Shapiro](https://arxiv.org/abs/1001.4077), which is related to the Boyer-Lindquist radius $r_{BL}$ by

$$ r_{BL} = r \left( 1 + \frac{r_+}{4r} \right)^2, \quad r_\pm = M \pm \sqrt{M^2 - a^2} $$

where $M$ and $a$ are the mass and spin parameter of the black hole, set by `kerr.mass` and `kerr.spin` (with $|a| \le M$). The horizon is then at $r = r_+/4$. Note that unlike the quasi-isotropic radial coordinate, where the horizon radius shrinks to zero as $|a| \to M$, here it stays finite ($r = M/4$ in the extreme Kerr limit), which makes it more accurate for rapidly spinning black holes.

The direction of the spin is set by `kerr.spin_direction`, which can be any non-zero vector. The solution is calculated in coordinates rotated so that the spin lies along the z-axis, and the tensors are then rotated back to the grid coordinates. By default the black hole is placed at `geometry.center`, and one can move it with `kerr.offset`.

The CCZ4 variables are `chi`$= det(\gamma)^{-1/3}$ and `h`$= \chi \gamma_{ij}$. In the coordinates where the spin is along the z-axis, only the $r\phi$ and $\theta\phi$ components of the extrinsic curvature $K_{ij}$ are non-zero, so its trace `K` is zero and `A`$=\tilde A_{ij}$ is only non-zero because of the rotation. In these coordinates the shift `shift`$=\beta^i$ also only has a $\phi$ component. Note that we do not use the analytic Kerr lapse, but a pre-collapsed lapse `lapse`$= \sqrt{\chi}$ instead.

The conformal connection `Gamma`$=\tilde{\Gamma}^i$ is non-zero since the conformal metric is not flat. This is calculated numerically from `h` using [`GammaCalculator`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Source/CCZ4/GammaCalculator.hpp) once the rest of the initial data has been set (see `initData()` in [`KerrBHLevel.cpp`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/KerrBH/KerrBHLevel.cpp)).

We currently impose a fixed hierarchy of 2:1 refinement grids centred on the black hole using the `FixedGridsTagger` (see [Tagging criteria](tagging.md)). Since the black hole should not move, regridding is turned off with `amr.regrid_int = -1`.

### Diagnostics

The Hamiltonian and momentum constraints, `Ham` and `Mom1`, `Mom2`, `Mom3`, are output in the plot files by setting `amr.derive_plot_vars = constraints`. These are diagnostics, so they are only calculated when a plot file is written (see [Diagnostics](diagnostics.md) for more details).

## Computational set up

Please read the [Performance Optimisation](performance_optimisation.md) guide to understand how to divide the domain into boxes that are shared over the MPI ranks, as this is crucial for obtaining good performance on HPC systems. When GRTeclyn regrids (with these parameters only at the first step), it will output the grid setup to the output file. This tells you how many boxes are on each level, and what their sizes are.

The parameters should be mostly self explanatory if you are familiar with NR, but you can look at the [**Parameters**](parameters.md) guide for more details.

The coarsest level has $64^3$ cells and the domain is $128M$ on each side. There are 4 refined levels, so the finest resolution is $\Delta x = M/8$. With `amr.max_grid_size = 32`, the coarsest level is divided into $(64/32)^3 = 8$ boxes and each refined level, which has $96^3$ cells, into $(96/32)^3 = 27$ boxes, so there are 116 boxes in total.

Note that Sommerfeld boundary conditions are used on all faces. This is because a black hole whose spin is tilted away from the z-axis is not symmetric under $z \to -z$, so one cannot use a reflective boundary. The smaller [`params_test.txt`](https://github.com/GRTLCollaboration/GRTeclyn/blob/main/Examples/KerrBH/params_test.txt), which is used for the regression test, has the spin along the z-axis, so there it is fine to use a reflective boundary on the lower z face and only half of the domain.

On an M1 MacBook Pro with no MPI, we obtain a speed of around 12 M/hr with the provided parameters. With this speed, the run to `evolution.stop_time = 20` takes about two hours. If you are getting speeds much lower than this, something is wrong. Each plot file is around 175 MB, and each checkpoint file is several times larger since it holds all the evolved variables, so make sure you have a few GB of disk space free.

## Checking the outputs

### Viewing 3D data and 2D slices

See [**Visualising Outputs**](visualising_outputs.md) for details on visualising. Note that checkpoint files are not viewable with AMReX, only plot files are.

One can check the conformal factor `chi`, which goes to zero at the centre of the black hole, and the `lapse`. These are written to the plot files together with the constraints.
