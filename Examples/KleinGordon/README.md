# Klein-Gordon Solver

This example solves the Klein-Gordon equation in the absense of GR. Please note that this example does not have a dynamical metric, and simply assumes a Minkowski-flat background metric. It is intended to provide a simple example that implements a real scalar wave equation, showing how different initial data, diagnostics and parameters can be set.
```math
\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi + V(\phi),
```

In the example parameter file, `params_test.ini`, you will find a section of parameters prefaced by `klein_gordon` - these control the behaviour of the Klein-Gordon simulation. 
The most important one of these is `klein_gordon.model`, for which you will need to enter a string to select the particular model being simulated. The model string must be one of the following:  
#### Wave  
This will solve the wave equation:
```math
\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi + \frac{1}{2} m^2 \phi^2
```
Without the potential term, the analytic solution is:
```math
$$ \phi (\bf{r},t) = \exp\left[i(\bf{k_r \cdot \bf{r}} - \omega t)\right] $$
```
We also assume $\omega = k_r$ since $c=1$.

#### SineGordon1D 
This will solve the 1D Sine Gordon equation:
```math
\frac{\partial^2 \phi}{\partial t^2} = \frac{\partial^2 \phi}{\partial x^2} - \sin \phi
```
The analytic breather solution is:
```math
\phi(x,t) = 4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta x)} \right)
```
where $\beta = \sqrt{1-\alpha^2}$
Note that it is a 1D solution embedded in a 3D simulation volume so `AMREX_SPACEDIM` is still set to 3.
#### SineGordon3D
This will solve the 3D Sine Gordon  equation:
```math
\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi - \sin \phi
```
The analytic pseudo-breather solution is
```math
\phi(x,y,z,t) = 4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta x)} \right) \times
                4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta y)} \right) \times
                4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta z)} \right)

```

Other settings in the parameter file:
* `klein_gordon.alpha` - (only for `SineGordon` models) controls the frequency of the breather mode and must be less than 1.
* `klein_gordon.scalar_mass` - (only for `Wave` model) mass of the scalar field in the potential
* `klein_gordon.wave_vector` - (only for the `Wave` model) this is the characteristic wave number, $$k_r$$, in the above equations. 
* `klein_gordon.initial_time` - starting time for the initial conditions, this is important for the analytic solution to match up to numerical solution.

There are three derived parameters available: 
* `phi_analytic` - the analytic solution to the field value at each cell
* `Pi_analytic` - the analytic solution to the first derivative of the field value at each cell
* `rho` - the value of the energy density at each cell.

Note that the analytic solution to the wave equation assumes no potential, i.e. $$m=0$$.

These will be outputted to the plot files if `amr.derive_plot_vars` is set with the name of the variable. Use `amr.derive_plot_vars = ALL` to print all the derived variables. 
