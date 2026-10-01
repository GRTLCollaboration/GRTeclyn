# The apparent-horizon finder: method

A summary of how `AHFinder` locates a marginally outer-trapped surface. The
companion document `ah_finder_implicit_ptc.md` is the investigation record
behind these choices.

## Surface representation

The surface is star-shaped about a chosen centre $c$ and parameterised by a
radius $h$ on a ring (latitude $\times$ longitude) grid, so grid point $ip$ sits
at

$$ x^k_{ip} = c^k + h_{ip}\, \hat n^k_{ip}, $$

with $\hat n_{ip}$ a fixed unit direction. Angular derivatives of $h$ come from
a finite-difference stencil on that grid. Each grid point carries one particle,
and the CCZ4 metric is interpolated onto the particles.

## The expansion

Writing the surface as the level set $F = r - h$, with

$$ F_i = \hat n_i - (\nabla h)_i, \qquad
   \lambda = \sqrt{\gamma^{ij} F_i F_j}, \qquad
   s_i = F_i / \lambda, $$

the outgoing null expansion is

$$ \Theta = \nabla_i s^i - K + s^i s^j K_{ij}, $$

evaluated from the interpolated $\chi$, $h_{ij}$, $K$, $A_{ij}$ and the first
derivatives of $\chi$ and $h_{ij}$, with $\gamma_{ij} = h_{ij}/\chi$ and
$K_{ij} = (A_{ij} + \tfrac13 h_{ij} K)/\chi$. A horizon is a surface with
$\Theta = 0$.

## Pseudo-transient continuation

$\Theta(h) = 0$ is solved by PTC. Each step solves

$$ \left( \frac{1}{\delta t} I + J \right) \delta h = -\Theta(h_n),
   \qquad J = \frac{\partial \Theta}{\partial h}, $$

with `amrex::GMRESMLMG`, matrix-free: $J$ is applied as a Brown--Saad
finite-difference directional derivative,

$$ J v \simeq \frac{\Theta(h_n + \varepsilon v) - \Theta(h_n)}{\varepsilon}. $$

$v$ is the direction GMRES asks the operator to act on, one component per grid
point in the same indexing as $h$; only $\delta h$ ever moves the surface. The
step length balances truncation against cancellation near
$\sqrt{\epsilon_{\text{mach}}}$, scaled by the two vectors involved:

$$ \varepsilon = \sigma\,\sqrt{\epsilon_{\text{mach}}}\;
   \frac{1 + \|h_n\|_2}{\|v\|_2}, $$

with $\sigma$ = `ah_finder.fd_eps_scale`, normally $1$. Dividing by $\|v\|_2$
keeps the operator linear whatever normalisation GMRES uses; the
$1 + \|h_n\|_2$ scales the perturbation to the surface while degrading to an
absolute step as $h_n \to 0$.

Small $\delta t$ gives a heavily damped, robustly globalised step; large
$\delta t$ approaches Newton. $\delta t$ is grown by an SER rule,
$\delta t \leftarrow \delta t \cdot r\, \|\Theta_n\| / \|\Theta_{n+1}\|$, and a
backtracking line search on $\|\Theta\|_\infty$ accepts or shortens each step.

## The frozen Jacobian and its shift

$\Theta$ depends on $h$ both through the surface geometry and through the metric
evaluated at the moving sample point $x(h)$. The mat-vec holds the metric fixed
at $h_n$, so the second dependence is dropped and each step costs one metric
interpolation rather than one per Krylov iteration.

The dropped term carries no derivative of $\delta h$ -- the metric at particle
$ip$ depends only on $h_{ip}$ -- so it is *diagonal*:

$$ J_{\text{exact}} = J_{\text{frozen}} + \mathrm{diag}\big(c(x)\big),
   \qquad c(x_{ip}) = \frac{\partial \Theta_{ip}}{\partial \Phi}\,
   \hat n^k \partial_k \Phi . $$

Knowing $c$ recovers the exact Jacobian at frozen cost. The mat-vec applies

$$ \big(J_{\text{frozen}} + \mathrm{diag}(\max(1/\delta t,\ c(x)))\big)\,\delta h, $$

so once the SER ramp has grown past the local $c$ the operator is exactly
$J_{\text{exact}}$ there, and below it the step retains its PTC damping. This
also keeps $1/\delta t$ clear of the pole of $(I/\delta t + J_{\text{frozen}})$.
For an isolated horizon of mass $m$, $c \to 3/m^2$.

## Measuring $c$

$c$ is obtained by transporting the frozen fields to the displaced radius.
Moving particle $ip$ outwards by $\delta h$ moves its sample point along
$\hat n_{ip}$, so every field $\Phi$ entering $\Theta$ changes by

$$ \delta \Phi = \delta h\; \hat n^k \partial_k \Phi . $$

Applying that shift to the frozen arrays and re-evaluating $\Theta$ gives the
diagonal directly, as a central difference in field space:

$$ c(x) = \frac{\Theta(\Phi + \epsilon\, \hat n^k\partial_k\Phi)
               - \Theta(\Phi - \epsilon\, \hat n^k\partial_k\Phi)}{2\epsilon}. $$

Both the field values and their first derivatives must be transported, so this
needs $\partial_k$ of all interpolated components and $\partial_j\partial_k$ of
$\chi$ and $h_{ij}$. The surface does not move, so no particle placement is
involved.

The measurement is gated: $c$ only influences the step once $1/\delta t$ has
fallen to it, so it is re-measured when
$1/\delta t < 2\,c_{\text{last}}$ and the previous value is reused otherwise.

## Convergence and cost

The iteration stops when $\|\Theta\|_\infty$ falls below `ah_finder.tolerance`.
A typical isolated horizon converges in 14 PTC steps to
$\Theta \sim 10^{-13}$; the common horizon of a close binary takes around 21.
The converged surface's area and irreducible mass are reported.

The attainable rate is set by the resolution of the grid the metric is
interpolated from, not by the solver.

## Limits

The star-shaped parameterisation has a finite basin. For a single mass $m$ the
expansion of a coordinate sphere is

$$ \Theta(r) = \frac{8r(2r-m)}{(2r+m)^3}, $$

which vanishes at the horizon $r = m/2$ and at infinity, with a maximum at

$$ r_* = m\left(1 + \tfrac{\sqrt3}{2}\right) \approx 1.87\,m . $$

Beyond $r_*$, $\mathrm{d}\Theta/\mathrm{d}r < 0$ and no step can recover, so the
initial guess must satisfy $r < r_*$. Guessing small is safe; guessing large is
not.

## Principal parameters

| parameter | meaning |
|---|---|
| `ah_finder.tolerance` | convergence threshold on $\|\Theta\|_\infty$ |
| `ah_finder.max_iter` | PTC iteration cap |
| `ah_finder.r` | SER growth factor for $\delta t$ |
| `ah_finder.newton_shift_auto` | measure $c$ each step |
| `ah_finder.newton_shift_diagonal` | apply $c$ pointwise rather than as its mean |
| `ah_finder.newton_shift_analytic` | obtain $c$ by field transport |

## Key files

- `Source/AHFinder/AHFinder.hpp`, `AHFinder.impl.hpp` -- PTC loop, mat-vec,
  interpolation, $\Theta$, shift measurement
- `Source/AHFinder/AHGeometry.hpp` -- ring grid, directions, derivatives of $h$,
  surface diagnostics
- `Source/AHFinder/AHJacobianOp.hpp` -- matrix-free operator for GMRES
- `Source/AHFinder/AHFinderParameters.hpp` -- input parameters
