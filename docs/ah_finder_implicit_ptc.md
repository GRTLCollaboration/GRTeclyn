# Apparent horizon finder: implicit pseudo-transient continuation

This document summarises the implicit solver used by `AHFinder` to locate an
apparent horizon, and the solver experiments that shaped its current design.

## Problem

The horizon is represented as a level set $r = h(\vartheta, \varphi)$ over a
sphere, discretised on a latitude $\times$ longitude "ring" grid (`AHGeometry`,
flat index `idx = i * ring_size + j`). The finder drives the surface radius $h$
to the place where the expansion $\Theta(h) = 0$.

$\Theta$ depends on $h$ two ways:

1. **directly**, through the surface derivatives $\nabla h$ / $\nabla\nabla h$
   (`AHGeometry`, where the antipodal pole coupling enters); and
2. **through the metric** $\gamma_{ij}(x)$ sampled at the surface point
   $x = x_c + h\,\hat{n}$, which is interpolated from the Cartesian AMReX grid
   onto particles placed on the surface.

## Method: pseudo-transient continuation (PTC)

Each step is an implicit-Euler / inexact-Newton update of the pseudo-time flow
$\dot h = -\Theta(h)$:

$$
\left(\frac{I}{\mathrm{d}t} + J\right)\delta h = -\Theta(h_n),
\qquad h_{n+1} = h_n + \delta h,
\qquad J = \frac{\partial\Theta}{\partial h}
$$

- The $I/\mathrm{d}t$ term is Levenberg-Marquardt/Tikhonov damping: small
  $\mathrm{d}t$ gives a heavily damped, robust step far from the solution;
  $\mathrm{d}t \to \infty$ recovers Newton near it.
- $\mathrm{d}t$ is grown by a **Switched Evolution Relaxation (SER)** rule as the
  residual falls, so the method *tries* to approach Newton as it converges. It
  does not actually get there -- the frozen Jacobian imposes a finite ceiling
  on $\mathrm{d}t$; see "dt does not run away" below.

This replaced the previous explicit second-order damped-wave relaxation
($\dot h = v - \eta h$, $\dot v = -c^2\Theta$) stepped by
`amrex::TimeIntegrator`. The velocity field $v$ and the damping coefficient
$\eta$ are gone; the implicit step needs neither.

### Jacobian: matrix-free, frozen-metric JFNK

$J$ is applied matrix-free as a finite-difference directional derivative
(Jacobian-free Newton-Krylov):

$$
J v \approx \frac{\Theta(h_n + \varepsilon v) - \Theta(h_n)}{\varepsilon}
$$

with $\varepsilon$ set by the Brown-Saad rule
$\varepsilon = \sqrt{\epsilon_{\mathrm{mach}}}\,(1 + \lVert h_n\rVert) / \lVert v\rVert$.
The operator apply is $F_{\mathrm{apply}}(v) = v/\mathrm{d}t + J v$.

The metric is **frozen**: `interpolate_metric(h_n)` runs once per PTC step, and
the mat-vec perturbs only the cheap, purely local tensor algebra
(`theta_from_metric`), reading the already-interpolated metric. This neglects
the $\partial\gamma/\partial h$ term, giving an *approximate* Jacobian, which
PTC / inexact-Newton tolerate. It is also what keeps the mat-vec cheap and clean
(see "Unfreezing the Jacobian" below).

### Linear solver

The linear system is solved with **single-level, matrix-free, unpreconditioned
GMRES** (`amrex::GMRESMLMG`):

- coarsening disabled (`LPInfo().setMaxCoarseningLevel(0)`);
- GMRES is not available as an MLMG `BottomSolver` in this AMReX version, so
  `amrex::GMRESMLMG` is used instead. It wraps an `amrex::MLMG` built on the
  operator, but MLMG is only the *operator host*: GMRES drives the Krylov
  iteration and reaches the mat-vec through `MLMG::applyPrecond`;
- `usePrecond(false)`. There are no geometric multigrid levels to precondition
  with and no `Fsmooth`/relaxation is implemented for this matrix-free
  operator, so with the preconditioner off MLMG's V-cycle and bottom solver are
  never entered;
- the iteration count is capped (`linear_max_iter`) and whatever increment
  GMRES reached is kept -- an inexact linear solve is acceptable for PTC, and
  unlike `MLMG::solve()` GMRES returns a failure status rather than throwing,
  so no exception handling is needed.

Because GMRES is driven directly, `linear_rel_tol` / `linear_abs_tol` are the
Krylov solver's own convergence tolerances. (Passed to `MLMG::solve()` they
would have set only the *outer* MLMG target, leaving the Krylov solve at
`bottom_reltol`.)

The pole topology (`AHGeometry::neighbours()` antipodal coupling at the
pole-adjacent rings) is handled **inside the operator** by reusing the existing
ring-grid stencil, so the solver never has to express it. Stock multigrid
(MLABecLaplacian etc.) cannot represent this anisotropic, non-symmetric,
antipodally-coupled, wider-than-5-point operator, which is why the framework is
used only as a Krylov driver.

### SER pseudo-timestep control

Implemented in `update_dt`:

$$
\rho = \mathrm{clamp}\!\left(r\,\frac{\Theta_{\mathrm{old}}}{\Theta_{\mathrm{new}}},\;
  \left[\sigma_{\mathrm{shrink}},\, \sigma_{\mathrm{grow}}\right]\right),
\qquad
\mathrm{d}t \leftarrow \max\!\left(\rho\,\mathrm{d}t,\; \mathrm{d}t_{\min}\right)
$$

where $\sigma_{\mathrm{shrink}}$, $\sigma_{\mathrm{grow}}$ and
$\mathrm{d}t_{\min}$ are `m_dt_shrink`, `m_dt_grow` and `m_min_dt`.

- There is **no upper cap** on $\mathrm{d}t$. The implicit solve is
  unconditionally stable, so a CFL-style cap is meaningless; robustness comes
  from step rejection instead. (An earlier CFL cap
  $\mathrm{d}t_{\max} = \mathrm{cfl} \cdot \Delta_{\mathrm{ring}}$ scaled as
  $\sim 1/N_{\mathrm{part}}$ and pinned $\mathrm{d}t$ far too small at high
  particle counts, stalling convergence -- it was removed, along with
  `min_ring_spacing` and the `cfl_factor` parameter.)
- **Step rejection**: a step is accepted only if it *strictly* reduces the
  residual ($\Theta_{\mathrm{new}} < \Theta_{\mathrm{old}}$). A non-improving or
  near-zero increment (e.g. GMRES stagnating on the ill-conditioned operator at
  large $\mathrm{d}t$) is rejected; $h_n$ and its frozen metric are restored and
  $\mathrm{d}t$ is cut by $\sigma_{\mathrm{shrink}}$. This lets $\mathrm{d}t$
  self-limit at whatever value the frozen Jacobian can still produce a usable
  direction at. The cut is deliberately aggressive (0.5) -- see "What a
  rejection actually means" below.

### Globalisation: backtracking line search

A rejected step throws away a GMRES solve, which is by far the most expensive
part of the iteration. Before rejecting, the same direction is retried at
shorter step lengths:

$$
h_{\mathrm{trial}} = h_n + \alpha\,\delta h,
\qquad \alpha = 1,\ \tfrac12,\ \tfrac14,\ \ldots
$$

halving $\alpha$ (by `m_backtrack_factor`) up to `max_backtracks` times and
accepting the first trial that reduces the residual. Because $\delta h$ solves
an *approximate* (frozen-metric) Jacobian system, the direction is usually good
even when the full step overshoots, so this converts most would-be rejections
into progress at the cost of one `interpolate_metric` + `theta_from_metric` per
backtrack.

The step length also feeds back into the controller: if $\alpha < 1$ the linear
model over-predicted at this $\mathrm{d}t$, so $\mathrm{d}t$ is scaled by
$\alpha$ (the measured over-prediction factor) rather than grown by SER. Only
full steps grow $\mathrm{d}t$. Setting `max_backtracks = 0` disables the line
search and recovers the plain accept-or-reject behaviour exactly.

### Parameters

Runtime (`ah_finder.*`, see `Tests/AHFinderUnitTest/params_test.txt`):

| parameter               | value | meaning                                    |
|-------------------------|-------|--------------------------------------------|
| `tolerance`             | 1e-4  | convergence threshold on inf-norm of $\Theta$ |
| `r`                     | 1.15  | SER target growth factor                   |
| `max_iter`              | 200   | max PTC iterations                         |
| `linear_rel_tol`        | 1e-6  | per-step GMRES relative tolerance          |
| `linear_abs_tol`        | 0.0   | per-step GMRES absolute tolerance          |
| `linear_max_iter`       | 1000  | max GMRES iterations per PTC step          |
| `linear_restart_length` | 1000  | GMRES restart length                       |
| `max_backtracks`        | 4     | max line-search halvings before rejecting  |
| `fd_eps_scale`          | 1.0   | multiplier on the Brown-Saad FD step $\varepsilon$ |
| `unfreeze_jacobian`     | 0     | 1 = exact Newton Jacobian on every step    |
| `unfreeze_every`        | 0     | exact Jacobian every $N$th step (0 = never) |
| `unfreeze_on_reject`    | 0     | exact Jacobian on the step after a reject  |
| `jacobian_diagnostic`   | 0     | 1 = print the frozen/exact Rayleigh probe  |
| `newton_shift`          | 0.0   | cap $\mathrm{d}t$ at $1/c$ (0 = no cap)    |
| `newton_shift_auto`     | 0     | 1 = measure $c$ each step instead          |

Hard-coded (in `AHFinder::init`): initial $\mathrm{d}t = 10^{-2}$,
$\mathrm{d}t_{\min} = 10^{-4}$, $\sigma_{\mathrm{shrink}} = 0.5$,
$\sigma_{\mathrm{grow}} = 1.25$, $\Theta_{\mathrm{floor}} = 10^{-12}$,
`m_backtrack_factor` $= 0.5$.

$\sigma_{\mathrm{shrink}}$ nominally does double duty as the lower clamp on the
SER ratio in `update_dt()`, but that clamp cannot bind for $r \geq 1$:
`update_dt()` is only reached on an accepted *full* step, where
$\Theta_{\mathrm{new}} < \Theta_{\mathrm{old}}$ and so the ratio exceeds $r$. In
practice it is purely the rejection response.

## Result

On the binary-BH test (`Tests/AHFinderUnitTest`, `num_particles = 1024`,
`n_cell = 128`) the finder converges in **20 iterations** to inf-norm
$\Theta = 3.2\times10^{-5}$, surface area $15.929$, irreducible mass $0.5629$.
The full unit-test suite passes (11/11).

Two changes got it there from the original 35 iterations:

| configuration                              | iterations | $\Theta$ | area        |
|--------------------------------------------|------------|----------|-------------|
| no line search (`max_backtracks = 0`)      | 35         | 1.88e-05 | 15.92948823 |
| + line search, $\sigma_{\mathrm{shrink}} = 0.8$ | 22    | 3.46e-05 | 15.92948832 |
| + $\sigma_{\mathrm{shrink}} = 0.5$ (current)    | 20    | 3.19e-05 | 15.92948829 |

Raising `max_backtracks` from 4 to 8 changes nothing, so 4 is already past the
point of diminishing returns. The surface area agrees to 7 significant figures
across every variant, so the horizon found is unchanged -- only the path to it
is shorter.

## Known issue: pseudo-timestep limit cycle

Before the line search was added, the run settled from about iteration 12 into
a strict accept/reject alternation, wasting every second iteration (a rejected
step still costs a full GMRES solve plus two `interpolate_metric` sweeps).
Roughly 40% of iterations were rejections. Two independent causes:

**1. $\sigma_{\mathrm{grow}}$ and $\sigma_{\mathrm{shrink}}$ were exact
reciprocals** ($1.25 \times 0.8 = 1$), so the controller could not settle. Once
$\mathrm{d}t$ reached the largest value $\mathrm{d}t^*$ the frozen Jacobian
supports: succeed at $\mathrm{d}t^*$ $\to$ grow to $1.25\,\mathrm{d}t^*$ $\to$
fail $\to$ shrink back to *exactly* $\mathrm{d}t^*$ $\to$ succeed $\to \ldots$
The logged $\mathrm{d}t$ values alternated between two numbers whose ratio is
1.25 to 10 significant figures.

Worse, the SER ratio $r\,\Theta_{\mathrm{old}}/\Theta_{\mathrm{new}}$ was
**saturated at $\sigma_{\mathrm{grow}}$ on every accepted step of the whole run**
(with $r = 1.2$ the residual ratio never drops below $\sim 1.04$), so SER was not
adapting at all -- the controller was pure bang-bang, $\times 1.25$ on accept and
$\times 0.8$ on reject.

**2. The linear solve was being truncated on every step.** With the original
`linear_max_iter = 100` / `linear_restart_length = 30`, every PTC step used
exactly 104 `theta_from_metric` calls -- i.e. GMRES hit its iteration cap every
single time and never reached `linear_rel_tol`. See "The linear budget was the
real bottleneck" below; raising it removes essentially all the rejections.

### What a rejection actually means: an ascent direction

Adding the line search removed most of the alternation, but left a cluster of
three consecutive *full* rejections (iterations 18-20 of the 22-iteration run),
with $\Theta$ pinned at $3.37\times10^{-4}$ while $\mathrm{d}t$ crawled
$0.178 \to 0.142 \to 0.114 \to 0.091$. Even $\alpha = 1/16$ failed.
Instrumenting the line search settled what was going on:

```
it 20  alpha 1       inf 3.3663e-4 -> 3.5306e-4    l2 1.6746e-4 -> 1.6991e-4
       alpha 0.5                  -> 3.4270e-4                 -> 1.6861e-4
       alpha 0.25                 -> 3.3888e-4                 -> 1.6801e-4
       alpha 0.125                -> 3.3776e-4                 -> 1.6773e-4
       alpha 0.0625               -> 3.3719e-4                 -> 1.6759e-4
```

The residual approaches its old value *from above*, and the excess halves as
$\alpha$ halves. That is a positive directional derivative: $\delta h$ is an
**ascent** direction, not an overshooting descent direction. No step length and
no choice of merit norm can rescue it. (In particular the 2-norm rises in
lockstep with the inf-norm, which rules out the earlier guess that the
non-smooth inf-norm merit function was to blame.)

The ascent diagnosis above is a direct measurement and stands on its own. What
*caused* the direction to be an ascent direction was, however, misattributed at
the time to the frozen Jacobian: it was the truncated linear solve. (An earlier
note here claimed a 3000-iteration GMRES run reproduced every $\delta h$ bit for
bit. That experiment was invalid -- the parameter override never reached the
solver, because `AHFinderUnitTest.cpp` inserts its own `params_test.txt` at
`argv[1]`, so both runs used identical inputs. See below for the corrected
measurement.)

Independently of the cause, the controller response below is still the right
one, and is still in place. A rejection is a
qualitatively worse failure than an overshoot -- it says the direction is
unusable, so $\mathrm{d}t$ must fall far enough to change the direction in *one*
step. $\sigma_{\mathrm{shrink}}$ was raised from 0.8 to 0.5, and the
three-iteration ratchet collapses to a single rejection:

```
iter 17: theta = 3.366332057e-04, dt = 0.1776356839
iter 18: theta = 3.366332057e-04, dt = 0.08881784197   <- one rejection
iter 19: theta = 1.286079133e-04, dt = 0.1110223025
iter 20: theta = 3.187803896e-05, dt = 0.1387778781    <- converged
```

Cause 1 above is still present in principle, but $\sigma_{\mathrm{grow}}$ and
$\sigma_{\mathrm{shrink}}$ are no longer reciprocal ($1.25 \times 0.5 = 0.625$),
so the exact 2-cycle is gone. The remaining inefficiency is that SER still grows
$\mathrm{d}t$ blindly by $\times 1.25$ on every accepted full step and so
periodically walks back into the ascent regime.

### dt does not run away: the Jacobian sets a hard ceiling

In textbook PTC $\mathrm{d}t \to \infty$ as the residual falls, recovering
Newton and its quadratic convergence. **That does not happen here.** Running to
`tolerance = 1e-12`:

```
iter 21: theta = 2.462e-05, dt = 0.1735
iter 22: theta = 2.462e-05, dt = 0.0867   <- reject
iter 23: theta = 4.806e-06, dt = 0.1084
iter 24: theta = 1.287e-06, dt = 0.1355
iter 25: theta = 9.895e-07, dt = 0.1694
iter 26: theta = 9.895e-07, dt = 0.0847   <- reject
...
iter 36: theta = 6.187e-11, dt = 0.1262   (converged, 36 iterations)
```

$\Theta$ falls six further orders of magnitude and $\mathrm{d}t$ never leaves
the $0.08$-$0.17$ band. The $\mathrm{d}t$ at which the step is rejected is
essentially constant across the whole run ($0.1776$, $0.1735$, $0.1694$,
$0.1654$, $0.1616$).

> The account below is qualitatively right but was later made precise: the
> ceiling is the pole of $(I/\mathrm{d}t + J_{\mathrm{frozen}})$ at
> $\mathrm{d}t = -1/\lambda_{\min} = 0.1247$, measured directly. See "Why the
> exact Jacobian is worse" below.

The cause is the frozen metric. As $\mathrm{d}t \to \infty$ the increment tends
to $J_{\mathrm{approx}}^{-1}(-\Theta)$, i.e. a **chord / modified-Newton** step
rather than a Newton step. That iteration converges only if the neglected
$\partial\gamma/\partial h$ term is small enough; here it is not, so the
$I/\mathrm{d}t$ regularisation is permanently load-bearing and
$\mathrm{d}t_{\max} \sim 0.17$ is a property of the Jacobian error, not of
how close the surface is to the horizon. The observed convergence is linear at
roughly $0.3$ per accepted iteration with no acceleration, exactly as a chord
method predicts -- not the quadratic rate a true Newton limit would give.

Two consequences:

- Raising $\mathrm{d}t$ further is not a tuning problem. It needs a better
  Jacobian, and unfreezing the metric was tried and rejected (see below).
- The controller now runs a **4-iteration cycle**: grow $\times 1.25$ three
  times, reject, halve. The net factor is $1.25^3 \times 0.5 = 0.9766 < 1$, so
  $\mathrm{d}t$ leaks slowly downward forever. Making the factors non-reciprocal
  removed the exact 2-cycle but only lengthened the period; no choice of constant
  factors avoids this. The proper fix is a memory of the smallest $\mathrm{d}t$
  known to have failed, with growth clamped just below it, so the controller
  parks at $\mathrm{d}t_{\max}$ instead of rediscovering it every fourth
  iteration. That would remove the $\sim 10\%$ of iterations currently spent on
  rejections (4 of 36 above). It has not been implemented.

## The linear budget was the real bottleneck

All of the above was measured with `linear_max_iter = 100`,
`linear_restart_length = 30`. Those defaults were chosen on the reasoning that
PTC only needs an *inexact* solve, so a small Krylov budget is a saving. The
instrumented `theta_evals` counter showed the flaw: every PTC step consumed
exactly 104 `theta_from_metric` calls, i.e. GMRES hit the cap every step and
never once converged to `linear_rel_tol`.

Re-running with `linear_max_iter = linear_restart_length = 1000` (8 MPI ranks,
`tolerance = 1e-10`, 1024 particles):

| configuration            | iters | rejections | metric interps | solve time |
|--------------------------|-------|------------|----------------|------------|
| `100` / `30` (old)       | 33    | 4          | 55             | 8.80 s     |
| `1000` / `1000`          | 35    | 1          | 47             | 7.53 s     |

GMRES actually needs **75 to $\sim 300$ Krylov iterations** per step, rising with
$\mathrm{d}t$ as $(I/\mathrm{d}t + J)$ becomes ill-conditioned. Giving it that
budget is *faster* in wall clock despite roughly tripling the mat-vecs, because a
`theta_from_metric` call is $\sim 2$ ms while the `interpolate_metric` sweeps a
rejected step forces are $\sim 60$ ms each. Trading cheap Krylov iterations for
expensive rejections is the wrong way round.

More importantly, with a converged linear solve the **only** rejection in the
whole run is the last one, at the $10^{-10}$ tolerance floor. The accept/reject
limit cycle documented above was largely an artefact of the truncated solve, not
an intrinsic property of the frozen Jacobian. (The $\alpha < 1$ line-search
damping steps remain, at iterations 12, 19 and 26.)

This also revives the "less accurate linear solve is better" note recorded under
BiCGStab below: the opposite is true here.

## Hybrid Jacobian policies (`unfreeze_every`, `unfreeze_on_reject`)

Since the exact Jacobian is $\sim 50\times$ more expensive per step, the obvious
idea is to use the cheap frozen one normally and pay for the exact one only
sometimes: either on a fixed cadence (`unfreeze_every = N`) or only once the
cheap one has demonstrably failed (`unfreeze_on_reject = 1`).

The switch is necessarily per **PTC step**, not per mat-vec: GMRES builds its
Arnoldi basis assuming one fixed linear operator, so alternating frozen and
unfrozen applies *inside* a single solve would break the Krylov relation. Each
solve therefore uses one operator throughout; only which operator varies from
step to step. Both policies are implemented and both are **off by default**,
because the measurement below says the exact Jacobian is not worth buying at any
cadence.

### The exact Jacobian gives a worse step at every residual level

With `unfreeze_every = 5` and a converged linear solve, comparing each exact
step against the frozen run at the same iteration:

| iter | frozen $\Theta$ | exact $\Theta$ |
|------|----------------|---------------|
| 1    | 0.13946217     | 0.14042014    |
| 6    | 0.07130136     | 0.07640638    |
| 11   | 0.00449907     | 0.00523318    |
| 16   | 0.00094307     | 0.00130943    |
| 21   | 7.632e-06      | 2.152e-05     |
| 26   | 1.951e-08      | 7.497e-08     |

Every exact step is worse, including step 1, where the residual is $0.14$ --
nowhere near any roundoff floor -- and where GMRES converged on *both* systems
(71 Krylov iterations for the exact operator, 75 for the frozen one). The run
then stalls at $\Theta \sim 1.2\times10^{-9}$ instead of converging.

`unfreeze_on_reject = 1` behaves even worse, for a structural reason: with a
converged linear solve the frozen Jacobian is not rejected until the very last
step, at the tolerance floor, where the exact Jacobian is rejected too. Since
the trigger is "the previous step was rejected", it then unfreezes on every
subsequent step, burning $\sim 180$ metric interpolations each, until `max_iter`.
Wall clock 55 s and climbing, versus 7.53 s frozen.

### It is not the finite-difference step

The natural suspicion is the FD step: the per-particle displacement is
$\varepsilon v_i \sim 1.7\times10^{-8}$, and unlike `theta_from_metric`,
`interpolate_metric` is only *piecewise* smooth in the particle positions
(interpolation stencil and AMR level selection change across cell boundaries).
That would make the difference quotient noise-dominated. Sweeping `fd_eps_scale`
(the multiplier on the Brown-Saad $\varepsilon$) rules it out:

| `fd_eps_scale` | $\Theta$ after exact step 1 |
|----------------|----------------------------|
| $1$            | 0.1404201364               |
| $10^{2}$       | 0.1404201421               |
| $10^{4}$       | 0.1404186871               |
| $10^{6}$       | 0.1489250606               |

Scales $1$ and $10^{2}$ agree to 7 significant figures, so the directional
derivative is well resolved and stable over two decades of $\varepsilon$ --
exactly what a *noise-free* difference quotient looks like. ($10^{4}$ and
$10^{6}$ degrade as expected from FD truncation error, and $10^{6}$ destroys the
solve outright.) The unfrozen mat-vec is computing a genuine, converged
derivative of the composite map $h \mapsto \Theta(h, \gamma(h))$; that derivative
simply yields a worse Newton direction than the frozen approximation does.

### Why the exact Jacobian is worse: it is the frozen one, shifted

`ah_finder.jacobian_diagnostic = 1` probes both operators at the *same* surface
with the same test modes, so they are compared as operators rather than through
the trajectories they produce. It reports the Rayleigh quotient
$q = \langle v, J v\rangle / \langle v, v\rangle$ for three low-order modes --
the uniform $\ell = 0$ breathing mode and $\ell = 1$, $\ell = 2$ modes aligned
with the binary axis, where the surface comes closest to the punctures. ($J$ is
not symmetric, so $q$ is an indicative eigenvalue estimate, not an eigenvalue;
the amplification law below justifies trusting it here.) It costs a few extra
metric interpolations per step and does not alter the trajectory.

At the initial surface ($h = 1.2$), `n_cell = 128`:

| mode        | $q_{\mathrm{frozen}}$ | $q_{\mathrm{exact}}$ | difference |
|-------------|-----------------------|----------------------|------------|
| $\ell = 0$  | $-7.0884$             | $+3.4656$            | 10.554     |
| $\ell = 1x$ | $-0.1289$             | $+10.4358$           | 10.565     |
| $\ell = 2x$ | $+13.4041$            | $+23.9522$           | 10.548     |

Two things fall out of this, and between them they explain everything above.

**1. The unfrozen term is a pure spectral shift.** The difference
$q_{\mathrm{exact}} - q_{\mathrm{frozen}}$ is the same for all three modes to
0.2%, even though the modes' own eigenvalues span a factor of 200. The
$\partial\gamma/\partial h$ contribution is therefore, to excellent accuracy, a
multiple of the identity:

$$
J_{\mathrm{exact}} \approx J_{\mathrm{frozen}} + c\,I,
\qquad c \approx 10.55\ \text{(initial surface)},\quad 11.96\ \text{(converged)}
$$

which is what a *zeroth-order* term
$(\partial\Theta/\partial\gamma)(\hat{n}^k \partial_k \gamma)$ should look like:
it multiplies $h$ pointwise rather than differentiating it, and
$\lvert\partial_r \ln\gamma\rvert$ varies little over the surface. So the frozen
Jacobian is not a noisy or partial approximation to the exact one. It is the
exact one with a known constant subtracted.

That makes the two PTC systems algebraically identical up to a change of
timestep:

$$
\frac{I}{\mathrm{d}t} + J_{\mathrm{exact}}
  = \frac{I}{\mathrm{d}t} + c\,I + J_{\mathrm{frozen}}
  = \frac{I}{\mathrm{d}t_{\mathrm{eff}}} + J_{\mathrm{frozen}},
\qquad
\mathrm{d}t_{\mathrm{eff}} = \frac{\mathrm{d}t}{1 + c\,\mathrm{d}t}
  \xrightarrow[\mathrm{d}t \to \infty]{} \frac{1}{c} \approx 0.084
$$

**Using the exact Jacobian is the same as using the frozen one with
$\mathrm{d}t$ capped at $1/c \approx 0.084$.** The SER controller keeps growing
$\mathrm{d}t$ believing it is approaching Newton, but the effective timestep
saturates, so the exact run is permanently over-damped -- uniformly short steps
at every residual level, never catastrophically wrong. That is exactly the
signature in the table above.

**2. The frozen Jacobian is indefinite, and that is why it is fast.**
$q_{\mathrm{frozen}}$ is *negative* for $\ell = 0$ ($-7.09$ at the initial
surface, drifting to $-8.02$ as the surface settles) and slightly negative for
$\ell = 1x$. So $(I/\mathrm{d}t + J_{\mathrm{frozen}})$ has a genuine pole at
$\mathrm{d}t = -1/\lambda_{\min} \approx 1/8.022 = 0.1247$, and below it the near
singularity acts as a step amplifier. The `ampl` column
($\lVert\delta h\rVert / (\mathrm{d}t\,\lVert\Theta\rVert)$, $1$ if
$I/\mathrm{d}t$ dominated) matches the $\ell = 0$ prediction
$1/(1 + q\,\mathrm{d}t)$ to within 0.5% for the first eleven iterations:

| iter | $\mathrm{d}t$ | $1 + q\,\mathrm{d}t$ | predicted `ampl` | measured `ampl` |
|------|--------|------------|------------------|-----------------|
| 1    | 0.0100 | 0.9291     | 1.076            | 1.077           |
| 5    | 0.0244 | 0.8221     | 1.216            | 1.215           |
| 9    | 0.0596 | 0.5413     | 1.847            | 1.845           |
| 11   | 0.0931 | 0.2582     | 3.874            | 3.855           |
| 12   | 0.1164 | 0.0632     | 15.82            | 15.62  (reject) |

The step is essentially *entirely* the $\ell = 0$ mode of the frozen Jacobian,
and the frozen solve is an accidental over-relaxation: its missing $-c\,I$ partly
cancels the $I/\mathrm{d}t$ regularisation and recovers a step several times
longer than $-\mathrm{d}t\,\Theta$. The exact operator is positive definite on
all three modes, so $\mathrm{ampl} = 1/(1 + q\,\mathrm{d}t) < 1$ always -- it can
never amplify at all.

So the answer to the open question is that the exact Jacobian is not a worse
Jacobian; it is a *better conditioned* one, and PTC with an $I/\mathrm{d}t$
regulariser rewards the ill-conditioned one. The frozen operator's defect happens
to point in the useful direction.

This also supersedes the explanation in "dt does not run away" above. The
ceiling on $\mathrm{d}t$ is not a vague "property of the Jacobian error": it is
the pole of $(I/\mathrm{d}t + J_{\mathrm{frozen}})$ at $\mathrm{d}t = 0.1247$.
Every rejection in the converged-solve run occurs at $\mathrm{d}t$ in
$[0.106, 0.126]$, i.e. immediately below or on top of it, and the final one
(iteration 34, $\mathrm{d}t = 0.1262$, just past the pole) records
$\mathrm{ampl} = 247$. The accept/reject cycle is SER walking $\mathrm{d}t$ into
a pole and being thrown back.

### Where the shift comes from

The shift is not a numerical accident; it is forced by the differential order of
the term the frozen Jacobian drops.

`theta_from_metric` takes two kinds of input: the frozen fields
$\Phi = (\chi, h_{ij}, A_{ij}, K, \partial_k\chi, \partial_k h_{ij})$, sampled at
the particle positions by the last `interpolate_metric`, and $h$ itself, which
enters as the radius $r$ and through the $\nabla h$, $\nabla\nabla h$ built by
`set_h_derivatives`. Perturbing $h$ moves the sample point radially, so

$$
\underbrace{\frac{\mathrm{d}\Theta}{\mathrm{d}h}}_{J_{\mathrm{exact}}}
= \underbrace{\left.\frac{\partial\Theta}{\partial h}\right|_{\Phi}}_{J_{\mathrm{frozen}}}
+ \underbrace{\frac{\partial\Theta}{\partial\Phi}\,\hat{n}^k\partial_k\Phi}_{\text{the shift}}
$$

$J_{\mathrm{frozen}}$ contains $\nabla$ and $\nabla\nabla$ acting on $\delta h$:
it is a second-order elliptic operator, the principal part. The dropped term
carries **no derivative of $\delta h$ at all** -- it multiplies $\delta h$
pointwise by the radial derivative of the metric. A zeroth-order term is a
multiplication operator, i.e. diagonal. So
$J_{\mathrm{exact}} - J_{\mathrm{frozen}} = \mathrm{diag}(c(x))$ was structurally
guaranteed. The empirical content of the measurement is only that $c(x)$ is
near-*constant*, collapsing the diagonal to $c\,I$.

### What sets the magnitude

For conformally flat, time-symmetric data ($\gamma_{ij} = \psi^4\delta_{ij}$,
$K_{ij} = 0$) the expansion of a coordinate sphere of radius $r$ is

$$
\Theta = \frac{2}{\psi^2 r} + \frac{4\psi'}{\psi^3}
$$

Differentiating at fixed $(\psi, \psi')$ versus totally, and using the horizon
condition $\psi'/\psi = -1/(2r)$ to eliminate $\psi'$:

$$
\lambda_{\mathrm{frozen}} = -\frac{2}{\psi^2 r^2},
\qquad
c = -\frac{1}{\psi^2 r^2} + \frac{4\psi''}{\psi^3}
$$

The shift is therefore dominated by $4\psi''/\psi^3$ -- the **radial curvature of
the conformal factor**. It is large here because the punctures make $\psi$ steep.

This is directly applicable here, because the test surface is **not** the common
horizon: `AHFinderUnitTest` centres the finder on `bh1` with
`guess_radius = m/2 = 0.25`, and the surface converges to a nearly round sphere
of coordinate radius $h \approx 0.2216$ around that one puncture, with the
companion a perturbation $2.0$ away. It is a slightly tidally distorted
Schwarzschild horizon.

Including the companion's contribution to $\psi$ ($\psi = 1 + m/2\rho + m/2\rho_2$,
so $\psi \approx 2.125$ at the initial surface and $2.253$ at the converged one):

| quantity | model | measured |
|---|---|---|
| $\lambda_{\mathrm{frozen}}$ ($\ell = 0$), initial surface | $-7.087$ | $-7.0884$ |
| $\lambda_{\mathrm{frozen}}$ ($\ell = 0$), converged        | $-8.023$ | $-8.02$   |
| $c$, converged                                             | $+12.06$ | $+11.98$  |

Better than 1% on all three, with no fitted quantities. The isolated
($m = 0.5$, $\psi = 2$) limit gives the round numbers $\lambda = -2/m^2 = -8$ and
$c = 3/m^2 = 12$.

This also explains why the cap $1/c$ lands *below* the pole
$-1/\lambda_{\min}$ rather than by luck: for the isolated profile the pole is at
$m^2/2$ and the cap at $m^2/3$, a fixed ratio of $3/2$ independent of $m$.
Measured, $0.1247 / 0.0835 = 1.494$.

### The diagonal is measured, and it really is constant

$\Theta$ at a ring-grid point depends on the metric **only at that same point**
(`theta_from_metric` is a per-particle loop over the frozen arrays), so
$J_{\mathrm{exact}} - J_{\mathrm{frozen}}$ is exactly diagonal, and the
per-particle difference of the two applies on the uniform mode *is* that
diagonal:

$$
c(x_{ip}) = \texttt{jv\_unfrozen}[ip] - \texttt{jv\_frozen}[ip]
$$

`measure_newton_shift()` already computes both vectors and keeps only the mean.
Dumping the difference (under `jacobian_diagnostic`, to
`shift_profile_<iter>.csv`) costs nothing and does not perturb the trajectory --
the run is bit-identical at 17 iterations and area $15.92948831$. Over the 1024
surface points:

| surface | mean $c$ | s.d. | s.d./mean | full spread |
|---|---|---|---|---|
| initial ($h = 0.25$)   | 10.554 | 0.030 | 0.3% | 1.4% |
| converged ($h \approx 0.2216$) | 11.982 | 0.068 | 0.6% | 2.8% |

**So $c(x)$ is genuinely flat, not merely flat-looking to the probe modes.** The
variation that does exist is structured rather than noise -- a clean dipole along
the binary axis, tracking the companion's tidal field:

| $\cos$(angle to companion) | distance to companion | $c$ | deviation |
|---|---|---|---|
| $-0.816$ | 2.186 | 11.9143 | $-0.57\%$ |
| $-0.246$ | 2.066 | 11.9583 | $-0.20\%$ |
| $+0.246$ | 1.957 | 11.9992 | $+0.14\%$ |
| $+0.816$ | 1.825 | 12.0718 | $+0.75\%$ |

which is exactly what $c \simeq 4\psi''/\psi^3$ predicts: the near side sits in a
slightly deeper companion potential.

Two consequences:

- **A diagonal shift buys very little.** Replacing the scalar $c$ with the
  measured $\mathrm{diag}(c(x))$ corrects a 0.6% error in a term of size 12,
  against a principal part of size 8 -- under 1% of the operator. It cannot
  explain a 0.1 contraction rate. (An earlier draft of this document said it
  would buy *nothing*; the shift sweep below shows a small irreducible floor at
  the optimal scalar that the diagonal is the right size to account for, so the
  honest statement is "very little, and only after the scalar is already
  optimal".)
- **The residual factor of 10 is therefore not the shift approximation.** With
  $\mathrm{d}t = 1/c$ the solved operator is $J_{\mathrm{exact}}$ to better than
  1%, so the capped iteration *should* be near-quadratic, and it is instead
  cleanly linear at $0.1$. The cause is found in the next section: it is grid
  resolution.

The measurement also redirects the validation question. Unequal masses would
change $c \sim 3/m^2$ but leave it uniform over each *individual* horizon, so it
is not the sensitive case. What makes $c(x)$ uniform here is that the surface is
a near-sphere at essentially constant $\rho$ from a single puncture. The test
that breaks that assumption is a genuinely **non-spherical** surface -- the
common horizon shortly after merger, or a highly spinning puncture -- where
$\psi''$ really does vary over the surface.

## Newton at frozen cost (`newton_shift`, `newton_shift_auto`)

The shift result is directly exploitable. Since
$I/\mathrm{d}t + J_{\mathrm{frozen}} = J_{\mathrm{frozen}} + c\,I
= J_{\mathrm{exact}}$ when $1/\mathrm{d}t = c$, **capping $\mathrm{d}t$ at $1/c$
makes the cheap frozen operator an exact Newton Jacobian** -- with no metric
re-interpolation inside the solve at all. Below the cap the step is a
Levenberg-Marquardt damped Newton step, which is the globalisation one wants
anyway, so the cap composes with the existing SER ramp and line search: SER
climbs, parks at $1/c$, and stays there.

The cap also keeps $\mathrm{d}t$ off the pole at $0.125$, since
$1/c \approx 0.084$ sits below it. That alone removes every rejection.

`newton_shift` sets $c$ directly; `newton_shift_auto = 1` measures it at each
step from the $\ell = 0$ mode (two extra metric interpolations -- about the cost
of one rejected step). On the equal-mass binary test, `tolerance = 1e-10`,
8 ranks:

| configuration            | iters | rejects | metric interps | $\Theta$ evals | solve  |
|--------------------------|-------|---------|----------------|---------------|--------|
| uncapped (SER only)      | 35    | 3       | 187            | 6918          | 13.2 s |
| `newton_shift = 11.98`   | 18    | 0       | **19**         | 3420          | 6.0 s  |
| `newton_shift_auto = 1`  | 17    | 0       | 52             | 3184          | 8.3 s  |

(Times are medians of three runs taken back to back; run-to-run spread on this
machine is $\sim 30\%$, so only the ratio is meaningful. Iteration and
interpolation counts are exact.) The surface area agrees to all 10 printed digits
with every other variant, $15.92948831$.

Once the cap binds, convergence is clean and rejection-free -- one decade per
iteration, exactly:

```
iter 11: theta = 1.431e-03, dt = 0.10434
iter 12: theta = 1.218e-05      <- cap reached, dt constant from here
iter 13: theta = 1.115e-06
iter 14: theta = 1.126e-07
iter 15: theta = 1.121e-08
iter 16: theta = 1.122e-09
iter 17: theta = 1.123e-10
iter 18: theta = 1.128e-11
```

**This is not quite Newton.** A true Newton iteration would be quadratic; what is
observed is *linear* convergence at a rate of almost exactly $0.1$. The factor of
10 per step is the residual error in
$J_{\mathrm{exact}} \approx J_{\mathrm{frozen}} + c\,I$ -- the shift is constant
to 0.2% across the three probed modes but evidently not across the full spectrum,
leaving roughly a 10% Jacobian error and hence a $0.1$ contraction. That is still
far better than the uncapped run, which mixes $\sim 0.3$ contraction with
periodic rejections.

### The shift has to be right

Sweeping `newton_shift` with everything else fixed shows a sharp optimum at the
measured $c$, which is itself the strongest confirmation that the mechanism is
$I/\mathrm{d}t + J_{\mathrm{frozen}} = J_{\mathrm{exact}}$ rather than merely
"a smaller $\mathrm{d}t$ is more stable":

| `newton_shift` | $\mathrm{d}t$ cap | iters |
|----------------|----------|-------|
| 6.0            | 0.167    | 35 (cap sits *above* the pole, so never binds) |
| 9.0            | 0.111    | 31    |
| 11.98          | 0.0835   | **18** |
| 16.0           | 0.0625   | 38    |
| 24.0           | 0.0417   | 78    |
| 48.0           | 0.0208   | 200 (did not converge) |

Too small and the cap is inoperative -- worse, a cap between $1/0.125$ and $c$
parks $\mathrm{d}t$ near the pole. Too large and the iteration is over-damped; at
$c = 48$ it fails to reach $10^{-10}$ within `max_iter` at all. Useful accuracy is
roughly $\pm 20\%$, which is why $c$ should be measured rather than guessed: it is
a property of the data and of where the surface sits ($10.55$ at the initial
guess, $11.98$ at the horizon here), so a hand-tuned value will not transfer to
another spacetime or another initial radius. Prefer `newton_shift_auto = 1`.

The $c = 48$ row also re-exposes the missing stagnation exit: the solver runs to
`max_iter` and then prints "converged" at $1.23\times10^{-10}$, above the
$10^{-10}$ tolerance. See open question 4.

### Where the factor of 10 comes from: it is the grid

With $\mathrm{d}t$ capped at the measured $1/c$ the iteration should be
near-quadratic and is instead cleanly linear at $0.1$. Four hypotheses were
tested; three are refuted and the fourth is confirmed.

**The merit function is not it.** Several apparent-horizon papers report on a
combination of an $\infty$-norm and a 2-norm, on the grounds that an
$\infty$-norm is set by one particle and can hide (or fake) the behaviour of the
rest of the surface. Printing both alongside each other (the
`jacobian_diagnostic` rms/argmax output) settles it here: from iteration 12 on,

| | contraction per step | ratio to $\|\Theta\|_\infty$ |
|---|---|---|
| $\|\Theta\|_\infty$ | 0.0995 | 1 |
| $\|\Theta\|_2/\sqrt{N}$ | 0.0990 | 0.345, constant |

The two norms contract at the same rate and their ratio does not drift, i.e. the
residual field decays *self-similarly* -- it is not one stubborn point holding
back a fast bulk. The `argmax` particle wanders (512, 448, 256, 801, 668, 691,
580, 703) before settling, which is the signature of a field decaying uniformly
rather than of a localised defect. A combined $\infty$/2-norm merit would
therefore change nothing about the convergence *rate* here. It may still be worth
having for **line-search robustness** -- accepting a step that reduces the bulk
residual while one point transiently worsens -- but that is a separate question
from the rate, and the accept/reject record does not currently show that failure
mode.

**The Krylov tolerance is not it.** `linear_rel_tol` of $10^{-6}$, $10^{-10}$ and
$10^{-14}$ give residuals agreeing to three significant figures at every
iteration ($9.4259$, $9.4121$, $9.4169 \times 10^{-7}$) and identical rates. Only
cost moves: cumulative $\Theta$ evaluations $3286 \to 4317 \to 5048$. The linear
solves are already converged far beyond what the outer iteration can use.

**A better scalar shift helps, but only so far.** Sweeping `newton_shift` finely
around the measured value, the asymptotic rate falls monotonically and bottoms
out *above* the measured mean:

| `newton_shift` | 11.0 | 11.4 | 11.8 | 11.98 | 12.05 | **12.1** | 12.15 | 12.6 |
|---|---|---|---|---|---|---|---|---|
| iters | 42 | 26 | 20 | 18 | 17 | **16** | 17 | 20 |
| rate | 0.75 | 0.35 | 0.175 | 0.100 | 0.06 | **0.03** | 0.03 | 0.12 |

so a $1\%$ error in $c$ costs a factor of 3 in rate -- consistent with the error
propagation above, where $M = cI + J_{\mathrm{frozen}}$ has dominant eigenvalue
only $11.98 - 8.02 \approx 3.96$ and so amplifies shift errors threefold. But
even at the optimum the rate is $0.03$, not quadratic.

**It is the grid resolution.** Refining the Cartesian grid with
`newton_shift_auto = 1`:

| `n_cell` | iters | measured $c$ | asymptotic rate |
|---|---|---|---|
| 96  | 200 (diverged) | 5780 | -- |
| 128 | 17 | 11.98200 | 0.0995 |
| 160 | 16 | 11.98110 | $\sim 0.05$ |
| 192 | 15 | 11.99227 | $\sim 0.015$ |

This is the result that reorganises the picture. The measured $c$ is flat to
$0.1\%$ across the three working resolutions -- the shift is fully converged --
while the rate improves by nearly an order of magnitude. Rate $\propto \Delta
x^{3}$ to $\Delta x^{4}$, which is the order of the interpolation stencil.

The interpretation: the function being rooted is not the continuum $\Theta(h)$
but the *discrete* one, which includes `interpolate_metric`. Its interpolation
error is oscillatory on the cell scale, so while the error itself is
$O(\Delta x^4)$ and invisible in $c$ (a smooth, mode-averaged quantity), its
derivative with respect to particle position -- which is precisely what the true
Jacobian contains and what the frozen operator omits -- is far rougher. The
capped operator is an excellent model of the smooth part of $\mathrm{d}\Theta/
\mathrm{d}h$ and no model at all of the grid-scale part, and the latter is what
sets the $0.1$ floor. At `n_cell = 96` the puncture is under-resolved badly
enough that this term dominates outright and the measured $c$ is meaningless.

Consequences:

- The $0.1$ rate is a discretisation artefact, not a defect of the shift model.
  At production resolutions it will improve on its own.
- The $\sim 12.1$ optimum in the scalar sweep is not the "true" $c$; it is the
  scalar that best *compensates* the grid term at `n_cell = 128`. It will not
  transfer, which is a further argument for `newton_shift_auto` over a tuned
  constant.
- There is no point chasing the diagonal $\mathrm{diag}(c(x))$ until this floor
  is below the $0.6\%$ the diagonal represents -- i.e. not at 128, possibly at
  192 and above.

### Still open

`newton_shift_auto` is enabled in the unit test but the code default is `0`
(uncapped), because the shift's near-constancy has only been verified on this one
configuration. Before making it the default it is worth checking that $c$ is
still mode-independent for a spinning or unequal-mass binary, and for a surface
started further from the horizon.

### The measured shift itself is resolution-converged

(Distinct from the section above, which is about the convergence *rate*. This one
is about the *value* of $c$.) One candidate explanation was that the unfrozen finite difference is a
well-converged derivative of the *interpolant* rather than of the metric -- a
bias, worst near the punctures where $\partial_k \gamma$ is largest, and one that
a `fd_eps_scale` sweep cannot detect. Refining the Cartesian grid tests it
directly. $q$ at the initial surface:

| `n_cell` | $q_{\mathrm{frozen}}$ ($\ell = 0$) | $q_{\mathrm{exact}}$ ($\ell = 0$) | shift $c$ |
|----------|--------------------|-------------------|-----------|
| 64       | $-7.0786$          | $2.7962$          | 9.875     |
| 128      | $-7.0884$          | $3.4656$          | 10.554    |
| 192      | $-7.0884$          | $3.4904$          | 10.579    |

$q_{\mathrm{frozen}}$ is resolution-independent to 0.14% throughout, as it must
be -- it never touches $\partial\gamma/\partial h$. $q_{\mathrm{exact}}$ *is*
badly under-resolved at `n_cell = 64` (24% low), but 128 and 192 agree to 0.7%,
so the unfrozen term is converged at the production resolution. Interpolation
bias is real at coarse resolution and not the explanation here.

The practical conclusion is unchanged: keep the metric frozen.

## Solver experiments (what did not work, and why)

These were tried and rejected; recorded so they are not re-attempted.

### BiCGStab instead of GMRES (frozen metric)

The original implementation used MLMG's `BottomSolver::bicgstab` as the Krylov
driver. It was replaced by GMRES; `BottomSolver::bicgstab` remains available if
it needs to be compared against.

Note the earlier caveat, recorded when GMRES was first trialled: with the frozen
(approximate) Jacobian, GMRES reached a higher $\mathrm{d}t$ early but then
**stalled near convergence** at $\Theta \sim 5\times10^{-3}$ and drove
$\mathrm{d}t$ to the floor, hitting `max_iter`, whereas BiCGStab converged. The
suspected cause is that accurately solving an *inaccurate* Jacobian produces
steps that are correct for the wrong linear model and overshoot, and BiCGStab's
inexact/partial solves were effectively better damped.

**That suspicion did not survive measurement.** The stall was the *truncated*
GMRES solve, and the fix was to iterate the linear solve harder, not less hard
(see "The linear budget was the real bottleneck"). Do not reach for a looser
`linear_rel_tol` or a lower `linear_max_iter` if the stall reappears.

### Unfreezing the Jacobian (exact Newton Jacobian)

> Superseded by "Hybrid Jacobian policies" above, which repeats these
> experiments with a converged linear solve and controls for the FD step. The
> conclusion (keep the metric frozen) is unchanged, but the *reason* recorded
> below -- interpolation noise in the FD derivative -- was measured and ruled
> out. Retained for the $\mathrm{d}t$-ceiling and cost numbers.

Making the mat-vec re-interpolate the metric at the trial surface
(`interpolate_metric(h_n + eps*v)`), so $J$ includes the
$\partial\gamma/\partial h$ term and becomes the exact Newton Jacobian:

- **Robustness up**: $\mathrm{d}t$ climbed to $\sim 0.34$ ($3\times$ the
  frozen-Jacobian ceiling of $\sim 0.11$); the Krylov solve no longer broke down
  there.
- **But it hit a hard residual floor at $\Theta \sim 2.8\times10^{-3}$** and could
  not converge (reject-and-shrink until $\mathrm{d}t$ collapsed to the floor, then
  `max_iter`). Cause: the metric is only known via **piecewise-polynomial grid
  interpolation**, so $\Theta \circ \mathrm{interp}(h)$ is only piecewise-smooth.
  The FD directional derivative through it carries interpolation-derivative error
  (largest near the punctures, where $\gamma$ is steep), which biases $J$. Once
  $\lVert\Theta\rVert$ nears that bias level, no Newton step reduces it.
- It is also far more expensive per step (every mat-vec becomes a particle
  re-query + two `interp()` + an MPI reduce). Timed on 8 ranks with
  `tolerance = 1e-10`: `interpolate_metric` goes from 1 call per PTC step to one
  per Krylov iteration, so a fully unfrozen run reached only iteration 10 in
  422 s where the frozen one reached iteration 10 in $\sim 8$ s -- and with a
  *worse* residual there ($3.4\times10^{-2}$ unfrozen vs $7.1\times10^{-3}$
  frozen), on an identical $\mathrm{d}t$ schedule with no rejections in either.

### GMRES with the unfrozen Jacobian

Traced the **same path** as unfrozen BiCGStab and stalled at the **same floor**
($\Theta \sim 3\times10^{-3}$, first plateau at $5\times10^{-3}$ then
$3\times10^{-3}$). This confirmed the floor is a property of the operator's
interpolation noise, not the Krylov method: no solver, BiCGStab or GMRES, can
solve below the noise level baked into the operator.

### Is the floor a grid-resolution problem?

Partly, but not the *surface* grid: the frozen Jacobian reaches
$\Theta \sim 2\times10^{-5}$ on the same ring grid, so 1024 surface points are not
the limit. The unfrozen floor is set by the **Cartesian metric-grid spacing and
interpolation order**, because the unfrozen mat-vec differentiates through that
interpolation. Refining the Cartesian grid (or using smoother/higher-order
interpolation) would lower the unfrozen floor -- but it would not be worthwhile:
the frozen Jacobian already reaches $\sim 2\times10^{-5}$ at the current
resolution, cheaply, precisely because it never differentiates the interpolation.

**Caveat:** every part of that account has since been superseded. The
`fd_eps_scale` sweep shows the unfrozen difference quotient is not
noise-dominated, and **the floor no longer reproduces at all**: with the
backtracking line search and a converged linear budget in place, the unfrozen
run descends smoothly past $3\times10^{-3}$ to $3.2\times10^{-9}$ by iteration
25. The reason to keep the metric frozen is not that unfreezing stalls -- it is
cost and rate, quantified in the next section.

### Does unfreezing pay off at higher resolution? No, and the gap widens

The natural expectation is that at production resolutions the exact Jacobian
should eventually earn its cost. It does not. Running both operators at two
resolutions to $10^{-10}$ (or `max_iter = 25`, whichever comes first):

| | $\Theta$ reached | iters | metric interps | wall (8 ranks) | asymptotic rate |
|---|---|---|---|---|---|
| `n_cell = 128`, frozen + cap | $9.7\times10^{-11}$ | 17 | **52** | **5.2 s** | 0.0995 |
| `n_cell = 128`, unfrozen | $3.2\times10^{-9}$ | 25 (hit cap) | 4088 | 264 s | 0.152, 0.122, 0.099 |
| `n_cell = 192`, frozen + cap | $2.3\times10^{-11}$ | 15 | **46** | **13.5 s** | $\sim 0.015$ |
| `n_cell = 192`, unfrozen | $2.0\times10^{-9}$ | 25 (hit cap) | 3669 | 898 s | 0.159, 0.131, 0.108 |

The unfrozen rate is **the same at both resolutions**, while frozen + cap
improves by $6\times$. That resolution-independence is the tell: the unfrozen
iteration's rate is not set by its operator at all, it is set by the **SER ramp**.
Look at its $\mathrm{d}t$ column -- still climbing geometrically at iteration 25
($0.145 \to 2.65$, a factor 1.25 per step), with $I/\mathrm{d}t \approx 0.38$
still not small against $\lVert J\rVert \approx 8$. The exact-Jacobian run never
actually reaches Newton; it is watching its own damping decay, and the residual
falls at whatever rate $\mathrm{d}t$ grows. Refining the grid cannot speed that
up, because the controller does not know about the grid.

Frozen + cap, by contrast, parks $\mathrm{d}t$ at $1/c$ on the *first* step and
is exactly-Newton from then on, so its rate is set by the only error left --
the grid -- and improves as $\Delta x^{3\text{--}4}$.

So the ordering does not cross over; it diverges. At 128 the frozen solve is
$51\times$ faster, at 192 it is $66\times$ faster, and it reaches a residual two
orders of magnitude lower. **There is no resolution at which unfreezing becomes
the right choice.**

## Close binaries: where the scalar shift degrades

The open question this document has carried is whether $c$ stays near-constant on
a genuinely **non-spherical** surface. It does not. To test it, the unit test's
finder target was made configurable (`test.ah_offset`, measured from the domain
centre like `bh1.offset`, and `test.ah_guess_radius`; defaults reproduce the
individual-horizon run bit-for-bit) and the punctures were moved together to
search for the *common* horizon, with $m = 0.5$ each.

| puncture positions | result | mean $c$ | s.d. | full spread | surface distortion |
|---|---|---|---|---|---|
| $\pm 1.0$ (individual horizon) | 17 iters | 11.98 | 0.6% | 2.8% | ~0 |
| $\pm 0.2$ | **84 iters** | 2.699 | 23% | 86% | 18% |
| $\pm 0.3$ | stalls at $3.1\times10^{-3}$ | 2.147 | 56% | 216% | 43% |
| $\pm 0.4$ | stalls at $6.0\times10^{-2}$ | 1.071 | 161% | 686% | 88% |
| $\pm 0.5$ | stalls at $1.9\times10^{-1}$ | 0.227 | 774% | 3667% | 130% |

Distortion is $(h_{\max} - h_{\min})/\bar h$. The pattern is exactly what
$c \simeq 4\psi''/\psi^3$ predicts: on a peanut-shaped surface $\psi''$ is
completely different on the axis (pointing at a puncture) and at the equator
(pointing into empty space between them). At $\pm 0.4$, for instance,
$c_{\mathrm{axis}} = 3.58$ against $c_{\mathrm{equator}} = 0.07$; by $\pm 0.5$ the
equatorial value has gone **negative** ($-0.48$), which a positive scalar cannot
represent even in principle.

Two consequences for the $\mathrm{d}t$ cap specifically. The mean becomes
meaningless as a summary once the s.d. exceeds it, and -- worse -- the mean
*falls*, so the cap $1/c$ rises: at $\pm 0.5$ it is $1/0.227 = 4.4$, far above the
pole, so the cap silently stops protecting against the pole at all. That is the
mechanism behind the stalls, which all show $\mathrm{d}t$ collapsed to the
$10^{-4}$ floor.

**The exact Jacobian stalls in the same places.** Re-running $\pm 0.3$ and
$\pm 0.4$ with `unfreeze_jacobian = 1` stalls at $7.5\times10^{-3}$ and
$4.8\times10^{-2}$, against $3.1\times10^{-3}$ and $6.0\times10^{-2}$ frozen,
also with $\mathrm{d}t$ on the floor. On that evidence alone it looked as though
the stalls were a property of the problem -- no common horizon yet, or a surface
not representable as a star-shaped $r = h(\theta,\phi)$ -- rather than of the
shift. **That reading was wrong**, and the diagonal below disproves it: $\pm 0.3$
converges fine once the shift is applied pointwise. The exact Jacobian fails
there for a different reason (see below), not because the horizon is absent.

A pleasing confirmation in passing: the common horizon's $c = 2.70$ against the
$3/M^2 = 3$ predicted for a single hole of the total mass $M = 1$, while the
individual horizon gives $11.98$ against $3/m^2 = 12$ for $m = 0.5$. The same
formula covers both, four-fold apart in magnitude.

### The diagonal shift (`newton_shift_diagonal`)

The measurement already produces the full $\mathrm{diag}(c(x))$ -- the scalar is
only its mean -- so applying the diagonal costs no extra interpolations. It is
implemented in the mat-vec, not in `AHJacobianOp`: the operator takes its
mat-vec as a callback from `AHFinder`, so nothing about the MLMG plumbing needs
to change.

The right generalisation is not to replace $I/\mathrm{d}t$ by
$\mathrm{diag}(c)$ outright. Capping $\mathrm{d}t$ at $1/c$ is the same thing as
using a per-step shift of $\max(1/\mathrm{d}t_{\mathrm{SER}}, c)$, so the
pointwise version is

$$
\text{shift}_{ip} = \max\!\left(1/\mathrm{d}t_{\mathrm{SER}},\; c(x_{ip})\right)
$$

which reduces exactly to the current behaviour when $c(x)$ is constant. Where the
SER ramp has grown past the local $c$ the operator becomes
$J_{\mathrm{frozen}} + \mathrm{diag}(c) = J_{\mathrm{exact}}$ *there*, and
elsewhere it keeps the heavier PTC damping. Taking the max also keeps the shift
positive where $c(x)$ is not, which is what happens on strongly distorted
surfaces. The scalar $\mathrm{d}t$ cap is dropped when the diagonal is on;
applying both would only re-impose the mean on every particle.

Measured, it is better everywhere -- including on the near-spherical surface
where the diagonal was predicted to be irrelevant:

| case | scalar | diagonal |
|---|---|---|
| individual horizon ($c$ spread 0.6%) | 17 iters, $\Theta = 9.7\times10^{-11}$ | **14 iters**, $\Theta = 6.1\times10^{-13}$ |
| common horizon, $\pm 0.2$ (spread 86%) | 84 iters, 18685 $\Theta$-evals, 34.1 s | **22 iters**, 3171 evals, 7.8 s |
| common horizon, $\pm 0.3$ (spread 216%) | stalls at $3.1\times10^{-3}$ | **24 iters**, $\Theta = 9.0\times10^{-11}$ |
| common horizon, $\pm 0.4$ | stalls at $6.0\times10^{-2}$ | stalls at $5.7\times10^{-2}$ |

Three things worth drawing out.

It helps on the *individual* horizon too, which was not expected: with a 0.6%
spread the diagonal itself is nearly the scalar, and the gain comes from
somewhere else -- dropping the $\mathrm{d}t$ cap. Because $\max(1/\mathrm{d}t,
c_{ip})$ is non-singular pointwise however large $\mathrm{d}t$ becomes,
$\mathrm{d}t$ is free to keep growing (to $0.218$, against the cap's $0.104$)
without ever walking into the pole. The cap was protecting against the pole by
refusing to grow; the diagonal protects against it directly, so it does not have
to.

The $\pm 0.3$ row is the one that corrects the earlier reading above. A common
horizon *does* exist there and *is* representable as a star-shaped surface; both
the scalar shift and the exact unfrozen Jacobian simply failed to find it. That
the exact Jacobian fails where the diagonal succeeds is worth stating plainly:
the diagonal is not merely a cheap approximation to Newton, it is a
better-globalised iteration than exact Newton, because $\max(1/\mathrm{d}t, c)$
supplies a per-point floor that keeps the operator conditioned while
$(I/\mathrm{d}t + J_{\mathrm{exact}})$ has no such guard and collapses
$\mathrm{d}t$ to the floor.

$\pm 0.4$ and beyond still fail, and there the "no horizon / not star-shaped"
explanation does look right: the standard deviation of $c(x)$ ($1.73$) exceeds
its mean ($1.07$) and $c$ goes negative at the equator, which is the signature of
a surface that is not a single smooth horizon.

`newton_shift_diagonal` defaults to `0`, and requires `newton_shift_auto = 1`
(which is what measures the diagonal); `check_params` enforces that. With it off,
the scalar path is bit-identical to before. On the evidence above the default
should probably be flipped, but the scalar remains the default for now.

### Diagnostics

**The mean $c$ is a good health indicator and should be checked.** It should come
out $O(1/M^2)$ for the mass enclosed. The failure modes seen here are all visible
in it: $5780$ when the puncture is under-resolved at `n_cell = 96`, and a
standard deviation exceeding the mean when the surface is too distorted for a
single smooth horizon. Both warrant a warning. The per-step print reports the
spread as `shift = <mean> +- <s.d.> (diag)` when the diagonal is enabled.

### What the shift costs, and why the analytic route was dropped

`measure_newton_shift()` finite-differences the metric, so it costs two extra
`interpolate_metric()` calls per PTC step on top of the one the step needs
anyway. On the individual-horizon baseline (`n_cell = 128`, 1024 particles,
8 ranks) that is:

| | interps | interp time | total |
|---|---|---|---|
| `newton_shift_auto = 1` | 52 (3/step) | 3.20 s (60%) | 5.35 s |
| fixed `newton_shift` | 19 (1/step) | 1.16 s (31%) | 3.74 s |

so the measurement is $\approx 38\%$ of runtime. The obvious fix is to compute
$c$ analytically from interpolated derivative data, collapsing three calls into
one. That was measured and **is not worth it.** Two things kill it:

1. `theta_from_metric()` reads 14 state values *and* 21 first-derivative values
   (`d1_gamma` feeds the Christoffels). Its radial derivative therefore needs
   $\partial_k$ of all 14 (42 values) *plus* $\partial_j\partial_k$ of the 7
   metric components (42 values). Payload per particle goes $35 \to 98$, so
   three calls $\times\,35$ becomes one call $\times\,98$ — a $7\%$ reduction in
   interpolated values, not the $3\times$ a naive count suggests.
2. Query cost is payload-dominated, not overhead-dominated. Timing the two
   existing queries separately, and adding a third probe query of the same
   shape, gives state (14 values, includes particle placement) $2.02$ s, deriv
   (21 values) $1.04$ s, probe deriv (21 values) $1.05$ s over 52 calls. A
   second 21-value query costs exactly what the first did, i.e. per-query cost
   is additive in the payload with little fixed overhead to amortise.

Putting those together, the analytic step would cost roughly $39 + 40 + 40
\approx 119$ ms per step against $3 \times 59 \approx 176$ ms now: a $\sim 12\%$
saving on total runtime, for a second derivative query, a hand-derived
$\partial_r\Theta$ and a new validation path. The cheaper win is the gate below,
which recovers most of the same $38\%$ with no new derivation.

**That cost estimate was too pessimistic, and the feature was built anyway** —
see "The analytic shift" below. The error was assuming the analytic route needs
a full `interpolate_metric()` call; it does not, because the particles are
already placed at $h_n$, so it needs only the three extra *queries* with no
placement and no refresh. Measured, it halves the cost of a measurement rather
than shaving 32% off it.

### Gating the measurement on whether $c$ can matter

Both uses of $c$ — the scalar cap $dt \leftarrow \min(dt, 1/c)$ and the diagonal
$\max(1/dt, c(x))$ — only bite once $1/dt$ has fallen to $c$. While $1/dt$ is
above it the measured value is computed and discarded. On the individual horizon
$1/dt$ runs $80 \to 13.4$ over iterations 1–9 against $c \approx 10.5$–$11.5$,
so $c$ first binds at iteration 10: nine of fourteen steps measured for nothing.

A blind "every $N$ steps" cadence would be wrong, because $c$ drifts $13.5\%$
over the solve, monotonically upwards and several percent per step early on —
and the shift sweep showed a $1\%$ error in $c$ costs a factor of $\sim 3$ in
rate. But **the drift and the sensitivity never overlap**: all the drift is in
the masked phase, and by the time $c$ binds it has converged, the last three
values agreeing to 1 part in $10^5$. So the rule is a gate, not a cadence:

> re-measure when $1/dt < \texttt{m\_shift\_gate\_margin} \cdot c_{\text{last}}$,
> or when no valid $c$ exists yet (first step, or a non-positive reading);
> otherwise reuse the last value.

`m_shift_gate_margin = 2.0` is an internal guard rail, not an input parameter:
it buys one factor of two of warning before $c$ can influence the step, and $c$
rises monotonically in every case measured, so it only has to cover one step's
worth of drift. Nothing has to be tuned per problem — where $c$ is large the
gate simply fires earlier. At guess radius $0.05$, $c = 114$ against $1/dt = 80$,
so it measures from iteration 1, which is correct because $c$ is doing work
there from the start.

Measured on the individual horizon (`n_cell = 128`, 1024 particles, 8 ranks),
three runs each:

| | steps measured | interps | interp time | solve |
|---|---|---|---|---|
| ungated, diagonal | 14/14 | 43 | 2.60–3.13 s | 4.53–5.79 s |
| gated, diagonal | 8/14 | 31 | 1.99–2.04 s | 4.24–4.39 s |
| ungated, scalar | 17/17 | 52 | 3.33 s | 6.82 s |
| gated, scalar | 11/17 | 40 | 2.98 s | 7.15 s |

**The converged result is bit-identical in every case** — same iteration count,
same final $\Theta$ (6.075140391e-13 diagonal, 9.659784084e-11 scalar), same
final $c$ to all digits printed — which is the point: the skipped values were
masked by the `max` anyway. The $\pm 0.2$ common horizon is likewise unchanged
at 21 iterations and $\Theta = 8.447198496 \times 10^{-11}$, with interpolations
falling $64 \to 38$ and only 8 of 21 steps measuring.

Interpolation call count falls $28\%$ (diagonal) and $23\%$ (scalar); solve time
follows for the diagonal at roughly $21\%$ off the median. The single scalar-mode
timing pair above shows the gated run *slower*, which is machine noise — the
repeated diagonal runs show the ungated spread ($4.53$–$5.79$ s) is wider than
the difference being measured. Call count is the reliable number on this box.

### The analytic shift (`newton_shift_analytic`)

`measure_newton_shift()` gets $c$ by re-interpolating the metric at a perturbed
surface and differencing two Jacobian applies -- a difference of two finite
differences, both in $h$. The analytic route never moves the surface.

Displacing particle $ip$ by $\delta h$ moves its sample point along
$\hat n = \texttt{direction(ip)}$, so every field $\Phi$ that
`theta_from_metric()` reads changes by $\delta h\, \hat n^k \partial_k \Phi$.
Transporting the frozen arrays by that amount and re-evaluating $\Theta$ applies
the chain rule exactly: there is no finite difference in $h$ at all. What
remains is a central difference in *field* space, where $\Theta$ is a smooth
algebraic function, taken at the $\sqrt[3]{\epsilon_{\text{mach}}}$ step optimal
for it.

Because `theta_from_metric()` reads 14 state values **and** 21 first-derivative
values, both have to be transported: the values need $\partial_k \Phi$ and the
derivatives need $\partial_j\partial_k \Phi$. So the extra data is

- $\partial_k$ of $K$ and $A_{ij}$ (21 values; the base query already carries
  $\partial_k$ of $\chi$ and $h_{ij}$),
- $\partial_j\partial_k$ of $\chi$ and $h_{ij}$ (42 values, split across two
  queries to stay inside the `num_components` entry limit),

and `gamma_ij = h_ij/chi` is rebuilt from the transported values rather than
transported itself, since it is derived rather than interpolated. `transport(0)`
restores the frozen state exactly, so the surface is left as it was found.

**The key cost point:** these are three *queries*, not an `interpolate_metric()`
call. The particles are already at $h_n$, so there is no `set_particle_positions()`
and no refresh -- which is what makes the route pay, and what the estimate above
got wrong.

**Validation.** Against the interpolating measurement, pointwise over all 1024
particles at iteration 0 (`jacobian_diagnostic = 1` dumps both):

| | mean | min | max |
|---|---|---|---|
| finite difference | 10.553985547 | 10.502600 | 10.648500 |
| analytic | 10.553985156 | 10.502600 | 10.648500 |

Mean relative difference $7.4\times 10^{-8}$, max $9.5\times 10^{-6}$ -- and that
max is one ulp of the CSV's six-significant-figure output, so the true agreement
is better than the comparison can resolve.

**Cost**, individual horizon, diagonal shift with the gate on, three runs each:

| | interps | interp time | solve |
|---|---|---|---|
| `newton_shift_analytic = 0` | 31 | 2.01–2.28 s | 4.29–5.02 s |
| `newton_shift_analytic = 1` | 23 | 1.59–1.79 s | 3.92–4.33 s |

Roughly $24\%$ off interpolation and $16\%$ off the solve, on top of the gate.
Against the original ungated interpolating measurement (median $5.60$ s) the two
changes together take the solve to $3.98$ s, about $29\%$. Scalar mode benefits
equally: 17 iterations either way, interpolations $40 \to 29$.

**It buys no accuracy, and that is itself the useful result.** The final shift
comes out $11.98199595$ analytic against $11.98199574$ interpolated, agreeing to
$2\times 10^{-8}$. The shift sweep put the optimum near $12.1$, about $1\%$
above both. That $1\%$ gap is therefore **not** finite-difference error in the
measurement -- it survives removing the FD in $h$ entirely. It is the grid
discretisation, consistent with the resolution study above.

Off by default. Iteration counts and converged surfaces are unchanged when it is
switched on ($\Theta = 4.9\times10^{-13}$ against $6.1\times10^{-13}$, both far
below tolerance).

## What a user should actually set

Most of this document is an investigation record. The operational summary is much
shorter, and the intent is that a user sets **two** parameters:

```
ah_finder.tolerance = 1e-8     # what "found" means, in |Theta|_inf
ah_finder.max_iter  = 100      # give-up budget
```

Everything else should be a default that is right without tuning. That is not yet
true of the code -- `ah_finder` currently registers 15 parameters and
`newton_shift_auto` defaults to `0` -- so the recommended cleanup is:

- **`newton_shift_auto = 1` should become the default.** It is the single change
  that makes the solver behave: it removes the $\mathrm{d}t$ ceiling, the
  accept/reject cycling and the $\mathrm{d}t$-leak, roughly halves the iteration
  count, and costs two extra metric interpolations per step against the $\sim 50$
  the cap saves. It is also self-calibrating -- $c$ is measured from the data
  each step -- which matters because $c$ is a property of the surface and the
  spacetime ($10.55$ at the initial guess, $11.98$ at the horizon here) and a
  hand-set value does not transfer. No production input file sets any `ah_finder`
  parameter today, so nothing depends on the current default.
- **`newton_shift` (the hand-set scalar) should go**, or be demoted to a
  diagnostic. Its only use was reaching the behaviour `newton_shift_auto` now
  reaches automatically, and the sweep above shows the best fixed value is
  resolution-dependent, so exposing it invites mis-tuning.
- **The `unfreeze_*` family should be demoted to diagnostics** (or removed). The
  measurements above show unfreezing is $51$--$66\times$ slower and converges to a
  worse residual, at every resolution tested, with the gap widening as the grid
  refines. There is no setting of these a user should reach for.
- **`linear_*`, `fd_eps_scale`, `r`, `max_backtracks`** are already fine at their
  defaults and are tuning knobs for this document's kind of work, not for users.
  The one thing worth keeping documented is that the linear budget must stay
  generous (`linear_max_iter` and `linear_restart_length` at 1000): truncating
  the Krylov solve to "save time" costs far more in rejected steps.

That leaves `tolerance` and `max_iter` as the user-facing surface, which is the
right size for this component.

On resolution specifically: the $0.1$ contraction rate that this document spends
a lot of effort explaining is a `n_cell = 128` artefact and **gets better on its
own** at production resolutions ($\sim 0.015$ at 192). A user running a
well-resolved binary should see the frozen + cap solver converge in $\sim 15$
iterations and a few tens of metric interpolations, without touching anything.
The one resolution caveat that *does* bite is the other end: at `n_cell = 96` the
puncture is under-resolved, the measured $c$ comes back as $\sim 5780$, and the
solve diverges. A sanity check on the reported shift (it should be $O(1/m^2)$)
would catch that and is worth adding.

## Conclusion

The frozen-metric, matrix-free PTC step is the chosen design: it is the cheapest
per step, produces a better Newton direction than the exact Jacobian at every
residual level measured, and self-limits $\mathrm{d}t$ via step rejection, with a
backtracking line search recovering most of the steps that rejection would
otherwise discard. The Krylov driver is unpreconditioned GMRES, and it must be
given a budget large enough to actually converge (`linear_max_iter` and
`linear_restart_length` both 1000): the cheap thing is a Krylov iteration, the
expensive thing is a rejected step.

The single most useful thing measured here is that the frozen and exact Jacobians
differ by a constant, $J_{\mathrm{exact}} \approx J_{\mathrm{frozen}} + c\,I$.
That one fact explains why the exact Jacobian looked worse, why $\mathrm{d}t$ had
a ceiling, and why the solver cycled through rejections -- and it turns the frozen
operator into a near-Newton one for free, by capping $\mathrm{d}t$ at $1/c$
(`newton_shift_auto`). Do not chase the exact Jacobian; measure the shift instead.

Open questions, in rough order of value:

1. Flip `newton_shift_diagonal` to `1` by default (and `newton_shift_auto` with
   it). Now implemented and measured better in every case tried: 14 iterations
   vs 17 on the individual horizon, 22 vs 84 on a distorted common horizon, and
   convergence rather than stagnation at puncture separation $0.6$. Held at `0`
   pending validation on a spinning puncture and on a genuinely post-merger
   surface.
2. Warn when the measured $c$ is unphysical -- it should be $O(1/M^2)$ for the
   enclosed mass. Both known failure modes show up there: $5780$ for an
   under-resolved puncture, and a standard deviation exceeding the mean for a
   surface that is not a single smooth horizon.
3. Investigate the star-shaped surface parameterisation as the real limit on
   close binaries. With the diagonal shift the solver now reaches puncture
   separation $0.6$; at $0.8$ and beyond everything stalls, and there the
   sign-indefinite $c(x)$ suggests a surface a single-centre
   $r = h(\theta,\phi)$ genuinely cannot represent.
4. There is no stagnation exit: when the residual stops moving, the solver runs
   to `max_iter` rather than giving up, and `find()` then prints "converged"
   regardless -- as the `newton_shift = 48` row above shows, it will report
   success at a residual above `tolerance`.
5. $\mathrm{d}t$ leaking downward is fixed as a side effect of the cap (the
   controller parks at $1/c$), but only when the cap is enabled. Uncapped, the
   growth/shrink factors still cannot settle.

Question 1 in earlier revisions of this document -- why the exact Jacobian gives
a worse direction -- is answered above: it does not give a worse direction, it
gives a shorter one, because it lacks the frozen operator's negative $\ell = 0$
eigenvalue and so cannot cancel any of the $I/\mathrm{d}t$ damping.

## Key files

- `Source/AHFinder/AHFinder.impl.hpp` -- PTC loop (`find()`), the mat-vec
  closure, `interpolate_metric` / `theta_from_metric` split, `update_dt`,
  backtracking line search.
- `Source/AHFinder/AHFinder.hpp` -- members and SER clamp declarations.
- `Source/AHFinder/AHJacobianOp.{hpp,cpp}` -- custom single-level matrix-free
  `MLCellLinOp` (the `Fapply` mat-vec, ring-grid <-> MultiFab scatter/gather).
- `Source/AHFinder/AHFinderState.hpp` -- `AHState` (now `h` only).
- `Source/AHFinder/AHFinderParameters.hpp` -- runtime parameters.
- `Source/AHFinder/AHGeometry.{hpp,impl.hpp}` -- ring grid, stencil, surface
  derivatives, area diagnostics (reused by the mat-vec, incl. antipodal
  `neighbours()`).
- `Tests/AHFinderUnitTest/` -- regression test and parameters.
