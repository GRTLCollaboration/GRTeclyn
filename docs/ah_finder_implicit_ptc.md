# Apparent horizon finder: implicit pseudo-transient continuation

This document summarises the implicit solver used by `AHFinder` to locate an
apparent horizon, and the solver experiments that shaped its current design.

## Problem

The horizon is represented as a level set `r = h(theta, phi)` over a sphere,
discretised on a latitude x longitude "ring" grid (`AHGeometry`, flat index
`idx = i * ring_size + j`). The finder drives the surface radius `h` to the
place where the expansion `Theta(h) = 0`.

`Theta` depends on `h` two ways:

1. **directly**, through the surface derivatives `grad_h` / `hess_h`
   (`AHGeometry`, where the antipodal pole coupling enters); and
2. **through the metric** `gamma_ij(x)` sampled at the surface point
   `x = center + h * direction`, which is interpolated from the Cartesian
   AMReX grid onto particles placed on the surface.

## Method: pseudo-transient continuation (PTC)

Each step is an implicit-Euler / inexact-Newton update of the pseudo-time flow
`dh/dt = -Theta(h)`:

```
(I/dt + J) delta_h = -Theta(h_n),   h_{n+1} = h_n + delta_h,   J = dTheta/dh
```

- The `I/dt` term is Levenberg-Marquardt/Tikhonov damping: small `dt` gives a
  heavily damped, robust step far from the solution; `dt -> infinity` recovers
  Newton near it.
- `dt` is grown by a **Switched Evolution Relaxation (SER)** rule as the
  residual falls, so the method approaches Newton as it converges.

This replaced the previous explicit second-order damped-wave relaxation
(`h_dot = v - eta*h`, `v_dot = -c^2 Theta`) stepped by `amrex::TimeIntegrator`.
The velocity field `v` and the damping coefficient `eta` are gone; the implicit
step needs neither.

### Jacobian: matrix-free, frozen-metric JFNK

`J` is applied matrix-free as a finite-difference directional derivative
(Jacobian-free Newton-Krylov):

```
J v ~= (Theta(h_n + eps*v) - Theta(h_n)) / eps
```

with `eps` set by the Brown-Saad rule `eps = sqrt(macheps) * (1 + ||h_n||) / ||v||`.
The operator apply is `Fapply(out, in) = in/dt + J in`.

The metric is **frozen**: `interpolate_metric(h_n)` runs once per PTC step, and
the mat-vec perturbs only the cheap, purely local tensor algebra
(`theta_from_metric`), reading the already-interpolated metric. This neglects
the `d(metric)/dh` term, giving an *approximate* Jacobian, which PTC /
inexact-Newton tolerate. It is also what keeps the mat-vec cheap and clean (see
"Unfreezing the Jacobian" below).

### Linear solver

The linear system is solved with AMReX's MLMG framework used as a
**single-level, matrix-free BiCGStab**:

- coarsening disabled (`LPInfo().setMaxCoarseningLevel(0)`);
- `BottomSolver::bicgstab` (BiCGStab handles the non-symmetric operator and
  needs only a mat-vec);
- bottom/final smoothing disabled (no `Fsmooth`/relaxation is implemented for
  this matrix-free operator; there are no geometric multigrid levels);
- `setThrowException(true)` with a `try/catch` that keeps whatever increment the
  solver reached -- an inexact linear solve is acceptable for PTC.

The pole topology (`AHGeometry::neighbours()` antipodal coupling at the
pole-adjacent rings) is handled **inside the operator** by reusing the existing
ring-grid stencil, so the solver never has to express it. Stock multigrid
(MLABecLaplacian etc.) cannot represent this anisotropic, non-symmetric,
antipodally-coupled, wider-than-5-point operator, which is why the framework is
used only as a Krylov driver.

### SER pseudo-timestep control

Implemented in `update_dt`:

```
ratio = clamp(r * theta_old / theta_new, [m_dt_shrink, m_dt_grow])
dt   *= ratio
dt    = max(dt, m_min_dt)
```

- There is **no upper cap** on `dt`. The implicit solve is unconditionally
  stable, so a CFL-style cap is meaningless; robustness comes from step
  rejection instead. (An earlier CFL cap `max_dt = cfl_factor * min_ring_spacing`
  scaled as ~1/num_particles and pinned `dt` far too small at high particle
  counts, stalling convergence -- it was removed, along with `min_ring_spacing`
  and the `cfl_factor` parameter.)
- **Step rejection**: a step is accepted only if it *strictly* reduces the
  residual (`theta_new < theta_old`). A non-improving or zero increment (e.g. a
  BiCGStab breakdown at large `dt`) is rejected; `h_n` and its frozen metric are
  restored and `dt` is shrunk. This lets `dt` self-limit at whatever value the
  linear solver can still make progress with.

### Parameters

Runtime (`ah_finder.*`, see `Tests/AHFinderUnitTest/params_test.txt`):

| parameter         | value  | meaning                                    |
|-------------------|--------|--------------------------------------------|
| `tolerance`       | 1e-4   | convergence threshold on inf-norm of Theta |
| `r`               | 1.15   | SER target growth factor                   |
| `max_iter`        | 200    | max PTC iterations                         |
| `linear_rel_tol`  | 1e-6   | per-step BiCGStab relative tolerance       |
| `linear_abs_tol`  | 0.0    | per-step BiCGStab absolute tolerance       |

Hard-coded (in `AHFinder::init`): initial `dt = 1e-2`, `m_min_dt = 1e-4`,
`m_dt_shrink = 0.8`, `m_dt_grow = 1.25`, `m_theta_floor = 1e-12`.

## Result

On the binary-BH test (`Tests/AHFinderUnitTest`, `num_particles = 1024`) the
finder converges in **21 iterations** to inf-norm `Theta = 1.7e-5`, surface area
`15.907`, irreducible mass `0.5625`. The full unit-test suite passes (11/11),
and results are bit-identical on 1, 2 and 4 MPI ranks.

## Solver experiments (what did not work, and why)

These were tried and rejected; recorded so they are not re-attempted.

### GMRES instead of BiCGStab (frozen metric)

GMRES is not available as an MLMG `BottomSolver` in this AMReX version; it
exists only as the standalone `amrex::GMRESMLMG` driver (run unpreconditioned,
since the operator has no smoother/coarsening to precondition with). With the
frozen (approximate) Jacobian, GMRES converged **worse** than BiCGStab: it
reached a higher `dt` early but then **stalled near convergence** at
`Theta ~ 5e-3` and drove `dt` to the floor, hitting `max_iter`. Accurately
solving an *inaccurate* Jacobian produces steps that are correct for the wrong
linear model and overshoot; BiCGStab's inexact/partial solves were effectively
better damped there. Tuning GMRES (restart length) would not help -- a better
linear solve makes the overshoot worse.

### Unfreezing the Jacobian (exact Newton Jacobian)

Making the mat-vec re-interpolate the metric at the trial surface
(`interpolate_metric(h_n + eps*v)`), so `J` includes the `d(metric)/dh` term
and becomes the exact Newton Jacobian:

- **Robustness up**: `dt` climbed to ~0.34 (3x the frozen BiCGStab ceiling of
  ~0.11); BiCGStab no longer broke down there.
- **But it hit a hard residual floor at `Theta ~ 2.8e-3`** and could not
  converge (reject-and-shrink until `dt` collapsed to the floor, then
  `max_iter`). Cause: the metric is only known via **piecewise-polynomial grid
  interpolation**, so `Theta . interp(h)` is only piecewise-smooth. The FD
  directional derivative through it carries interpolation-derivative error
  (largest near the punctures, where `gamma` is steep), which biases `J`. Once
  `||Theta||` nears that bias level, no Newton step reduces it.
- It is also ~20x more expensive per step (every mat-vec becomes a particle
  re-query + two `interp()` + an MPI reduce).

### GMRES with the unfrozen Jacobian

Traced the **same path** as unfrozen BiCGStab and stalled at the **same floor**
(`Theta ~ 3e-3`, first plateau at 5e-3 then 3e-3). This confirmed the floor is a
property of the operator's interpolation noise, not the Krylov method: no
solver, BiCGStab or GMRES, can solve below the noise level baked into the
operator.

### Is the floor a grid-resolution problem?

Partly, but not the *surface* grid: the frozen Jacobian reaches `Theta ~ 2e-5`
on the same ring grid, so 1024 surface points are not the limit. The unfrozen
floor is set by the **Cartesian metric-grid spacing and interpolation order**,
because the unfrozen mat-vec differentiates through that interpolation.
Refining the Cartesian grid (or using smoother/higher-order interpolation)
would lower the unfrozen floor -- but it would not be worthwhile: the frozen
Jacobian already reaches `~2e-5` at the current resolution, cheaply, precisely
because it never differentiates the interpolation.

## Conclusion

The frozen-metric, matrix-free BiCGStab PTC step is the chosen design: it is
the cheapest per step, avoids the interpolation-noise floor entirely, self-limits
`dt` via step rejection, and converges to well below tolerance in ~21 iterations.

## Key files

- `Source/AHFinder/AHFinder.impl.hpp` -- PTC loop (`find()`), the mat-vec
  closure, `interpolate_metric` / `theta_from_metric` split, `update_dt`.
- `Source/AHFinder/AHFinder.hpp` -- members and SER clamp declarations.
- `Source/AHFinder/AHJacobianOp.{hpp,cpp}` -- custom single-level matrix-free
  `MLCellLinOp` (the `Fapply` mat-vec, ring-grid <-> MultiFab scatter/gather).
- `Source/AHFinder/AHFinderState.hpp` -- `AHState` (now `h` only).
- `Source/AHFinder/AHFinderParameters.hpp` -- runtime parameters.
- `Source/AHFinder/AHGeometry.{hpp,impl.hpp}` -- ring grid, stencil, surface
  derivatives, area diagnostics (reused by the mat-vec, incl. antipodal
  `neighbours()`).
- `Tests/AHFinderUnitTest/` -- regression test and parameters.
