This solves the Klein-Gordon equation
```math
\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi + V(\phi),
```
based on the wave example in the [AMReX guided tutorials](https://github.com/AMReX-Codes/amrex-tutorials).

Choose from three possible models in the input parameter file:
- "Wave" for the wave equation:
```math
\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi + \frac{1}{2} m^2 \phi^2
```
Without the potential term, the analytic solution is:
```math
$$ \phi (\bf{r},t) = \exp\[i(\bf{k_r \cdot \bf{r}} - \omega t)\] $$
```
and $\omega = k_r$ since $c=1$.

- "SineGordon1D" for the 1D Sine Gordon breather solution to the equation:
```math
\frac{\partial^2 \phi}{\partial t^2} = \frac{\partial^2 \phi}{\partial x^2} - \sin \phi
```
The analytic solution is:
```math
\phi(x,t) = 4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta x)} \right)
```
where $\beta = \sqrt(1-\alpha^2)$
Note that it is a 1D solution embedded in a 3D simulation volume so `AMREX_SPACEDIM` is still set to 3.
- "SineGordon3D" for the 3D Sine Gordon pseudo-breather solution to the equation:
$$\frac{\partial^2 \phi}{\partial t^2} = \nabla^2 \phi - \sin \phi$$
The analytic solution is
```math
\phi(x,y,z,t) = 4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta x)} \right) \times
                4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta y)} \right) \times
                4 \arctan \left(\frac{\beta \cos(\alpha t)}{\alpha \cosh(\beta z)} \right)

```


Other settings in the parameter file:
* `alpha` - (only for SineGordon models) controls the frequency of the breather mode and must be less than 1.
* `scalar_mass` - (only for Wave model) mass of the scalar field in the potential
* `initial_time` - starting time for the initial conditions, this is important for the analytic solution to match up