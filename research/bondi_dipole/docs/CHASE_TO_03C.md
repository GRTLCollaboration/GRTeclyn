# Reaching v = 0.3 c — realisation plan for the Bondi-dipole chase cell

**What this is.** A build-and-run plan for the one cell the campaign does not
have and cannot get from the current code: a mixed canonical/phantom pair ridden
long enough for its speed to become a *relativistic* measurement rather than a
Newtonian one. `GPU_RUN_PAPER.md` section 4 sizes the problem and section 5
sketches the recentring box; this document is the worked version — the arithmetic
checked against the archive, the code named file by file, the seams costed, the
validation ladder given pass gates that can actually be met, and the physics risk
that dominates the whole exercise stated up front rather than at the end.

**Status.** Nothing here is launched. No code exists yet. This is the plan that
has to be agreed before either is true.

**The one-line summary.** The static 64-box cannot report a speed above about
`0.06 c`, because the leading star's envelope reaches the sponge at `t ≈ 400–500`
and the run stops meaning anything after that. Everything below is about removing
that ceiling for about twelve GPU-hours and roughly six hundred lines of new,
default-off code.

---

## 1. What "0.3 c" costs, exactly

A pair holding a constant proper acceleration `a` follows hyperbolic motion:

```
v(t) = a t / sqrt(1 + (a t)^2)          x(t) = ( sqrt(1 + (a t)^2) - 1 ) / a
```

so the target speed fixes `a·t` and the distance, and nothing else:

| target `v` | `a·t` needed | `γ` | distance × `a` | Newton over-predicts `v` by |
|---|---|---|---|---|
| 0.05 c | 0.05006 | 1.00125 | 0.00125 | 0.13 % |
| 0.10 c | 0.10050 | 1.00504 | 0.00504 | 0.50 % |
| 0.20 c | 0.20412 | 1.02062 | 0.02062 | 2.06 % |
| **0.30 c** | **0.31449** | **1.04828** | **0.04828** | **4.83 %** |

The last column is the point of the run. Newton says the speed is `a·t`;
relativity says it is `a·t/γ`. At `0.3 c` those differ by 4.8 %, which is far
outside every error bar the campaign carries, and the discrepancy grows smoothly
from zero — so the whole late trajectory is the measurement, not just its
endpoint.

### Which acceleration to use

The campaign has three defensible values, and they set the bill:

| configuration | measured `a` | `t` to 0.3 c | distance travelled | GPU-hours at 5.5 h / 1000 t |
|---|---|---|---|---|
| `d = 8`, N = 128 | 2.32e-04 | 1,360 | 208 | 7.5 |
| `d = 10`, N = 128 | 1.463e-04 | 2,150 | 330 | 11.8 |
| `d = 10`, converged (N = 256) | 1.596e-04 | 1,970 | 302 | 10.8 |

The N = 128 cost figure is a measurement, not an estimate: ten archived cells ran
at 5.5 GPU-hours per 1000 units of `t` and did not vary by more than 2 %.

**Recommendation: run the chase at `d = 10`, N = 128, and budget 12 GPU-hours.**
It costs 45 % more than `d = 8` and it is worth it. `d = 10` is the configuration
the resolution ladder converged on, it is the one every control in the paper is
read against, and its stars (rms radius ≈ 5.3 each) are not sitting inside one
another the way they are at `d = 8`. Buying a cheaper cell in a geometry the
paper does not otherwise use would mean defending a second configuration for no
scientific return.

An N = 128 run will follow its own `a`, so plan on `t = 2,200` and stop on the
measured speed rather than on the clock. If the true acceleration is the
converged 1.596e-04, the target arrives at `t ≈ 1,970` and the cell finishes an
hour early.

### The milestones on the way

Using the converged `a = 1.6e-04`, the run passes through:

| `v` | reached at `t ≈` | drift `x ≈` | what it is worth |
|---|---|---|---|
| 0.058 c | 400 | 11.5 | the edge of what the static box can already do |
| 0.10 c | 630 | 31 | first speed no existing cell can report |
| 0.20 c | 1,280 | 129 | Newtonian deviation becomes larger than the resolution spread |
| 0.30 c | 1,970 | 302 | the headline |

This matters operationally: the cell is *readable at every stage*. It can be
stopped at any milestone and still have produced a result the campaign did not
have. There is no all-or-nothing point in it.

---

## 2. Why the present setup stops at about 0.06 c

Not because of the wall clock. Because of the box.

The leading (canonical) star drifts toward the `+x` face. The sponge — the band
of extra dissipation that is the only reason mixed cells survive past `t = 60` —
starts at radius 24 and saturates at 32, measured from the box centre. The star
is not a point: its rms radius is 5.3 and grows, its 99 %-mass radius is 9.0, and
it sits half a separation ahead of the pair's midpoint.

Envelope reaches the sponge when `drift + d/2 + r₉₉ > 24`, i.e. at `drift ≈ 10`.

That is not a projection. It is what the existing `t = 400` cell did:

| `t` | canonical core `x` | phantom core `x` | distance of leading core from centre |
|---|---|---|---|
| 0 | 37.00 | 27.00 | 5.0 |
| 200 | 39.97 | 29.86 | 8.0 |
| 300 | ≈ 44.0 | ≈ 33.2 | 12.0 |
| 400 | 49.55 | 37.70 | **17.5** |

At `t = 400` the canonical core is 17.5 from the centre and its 99 % envelope
reaches radius 26.5 — inside the sponge band. The ramp is quartic, so the
dissipation it is feeling there is still weak (about 1.6 % of full strength), but
it is no longer zero, it grows as the fourth power of how far in the star goes,
and by `drift ≈ 18` the envelope is in the saturated region.

**So `longrun_pair_d10_t400_L64_N128_lev0` is already at the edge of its box**,
and its final quarter is the first data in the campaign taken with the leading
star's halo inside the sponge. That is worth recording in `Findings.md`
independently of this plan; here it fixes the ceiling:

> Without recentring, no cell in a 64-box can report a trajectory past
> `drift ≈ 10–18`, which is `t ≈ 370–500`, which is `v ≈ 0.053–0.072 c`.

A fifth of the way to the target, and the missing four-fifths are exactly the
part where relativity and Newton stop agreeing.

---

## 3. The four ways out, and which one to build

| route | what it changes | cost to 0.3 c | new code | risk |
|---|---|---|---|---|
| **A. Long box** | domain stretched along `x` to 352, transverse side kept at 64 | 1.7 GPU-days, 33 GB | none — per-direction cell counts already exist | none new, but the cost grows with every further target |
| **B. Recentring box** | data periodically carried back to the box centre; box stays 64 | **~12 GPU-hours, 6 GB** | ~600 lines, default-off | two fabricated/discarded shells per shift |
| **C. Heavier stars** | a denser star family, so `a` is larger everywhere | scales as `1/M²` | none — a solver knob | the mass is capped by stability, see below |
| **D. Initial boost** | start the pair already moving | skips the slow early phase | none — `grtresna_boost_lumps` exists | that code path has never been validated |

**Build B. Use C as a lever only if B's cost becomes the binding constraint, and
treat D as unusable until someone validates it.**

### Why not A

The long box is the honest fallback and it should stay on the books, but it is
the wrong default for three reasons. It costs three times as much. Its memory
(33 GB) leaves no room on the card for anything else. And its cost keeps growing:
0.5 c needs an 800-long box, 0.9 c needs 5,728 and 537 GB, which is off the node.
Route B's cost is linear in run time and its memory never changes — that is the
whole point of it.

The one thing A has that B does not is the wake. See §4.7.

### What C is actually worth — the outlook's "3× mass" is too optimistic

`GPU_RUN_PAPER.md` section 4 suggests taking a star three times heavier and
banking a 9× saving. The star family table
(`results/bondi-dipole-runaway/stars/star_family.csv`) does not support 3×:

| ω | canonical `M` | ratio to ω = 0.75 | 90 %-mass radius | on the stable branch? |
|---|---|---|---|---|
| 0.750 | 0.014350 | 1.00 | 5.24 | yes (surveyed) |
| 0.700 | 0.019172 | 1.34 | 5.59 | **unmeasured** |
| 0.670 | 0.023321 | 1.62 | 5.77 | **unmeasured — claimed branch edge** |
| 0.650 | 0.026859 | 1.87 | 5.98 | below the claimed edge |
| 0.600 | 0.039892 | 2.78 | 6.62 | unstable |

A 3× mass means ω ≈ 0.59, which is below the stable branch. At the branch edge
the gain is **1.62×, so the saving is 2.6×, not 9×**. And the branch edge itself
is not measured: the stability survey covered ω = 0.75, 0.80, 0.85, 0.90, all
*above* the design point. The "branch runs from ω ≈ 0.67 up" is inherited, not
tested here.

So route C costs, before it saves anything: a stability cell at ω = 0.67 (1.1 h),
a mass-matched phantom partner retuned on CPU (the family table interpolates to
**ω ≈ 0.682** for `|M₋| = 0.02332`, to be confirmed by re-running
`star_family_scan.py`), and a fresh `a·d² = GM` check to show the point-mass law
still holds for the heavier star. That is two GPU-hours and a day of care to save
about seven GPU-hours — not worth it for 0.3 c. It becomes worth it at 0.5 c and
essential at 0.9 c, which is the follow-up paper's problem.

### Why D is out

`grtresna_boost_lumps` is set to 0 in every cell the campaign has ever run. An
unvalidated code path that writes the initial data is the single worst place to
take a shortcut, because its failure mode is a plausible-looking trajectory. If
anyone wants it, it gets a validation cell of its own first: boost to a speed the
unboosted run also reaches, and check the two agree from there on. That is a
separate piece of work, not a lever for this one.

---

## 4. The recentring box, in detail

### 4.1 The idea, and the invariant

The domain never changes: 64 on a side, 128 cells, centre at 32, sponge at 24/32,
Sommerfeld boundaries — byte-identical to every mixed cell in the archive. What
changes is that when the pair has drifted `Δ` toward the leading face, **every
field is copied back by exactly `Δ/dx` cells** and `Δ` is added to a running
odometer. The physical trajectory is then

```
true displacement(t) = position on the grid(t) + odometer(t)
```

and the pair is back where the sponge and the outer boundary expect it to be.

The shift is always a whole number of cells, so it is a pure relabelling of which
cell holds which value — **no interpolation, exact to the last bit**. This is the
same discipline that fixed the initial-data alignment (artefact rule 1): shifts
that land on the mesh are free, shifts that land between cells are not, so only
the former are ever taken.

### 4.2 The shift primitive

AMReX provides both halves directly, and the whole operation is about fifteen
lines:

```
MultiFab tmp(state.boxArray(), state.DistributionMap(), NUM_VARS, 0);
MultiFab::Copy(tmp, state, 0, 0, NUM_VARS, 0);   // exact copy
tmp.shift(IntVect{-n, 0, 0});                     // relabels boxes AND fabs
state.ParallelCopy(tmp, 0, 0, NUM_VARS, 0, 0, Periodicity::NonPeriodic());
fill_front_sliver(state, n);                      // the only invented data
```

`FabArray::shift` (`Src/Base/AMReX_FabArray.H:3056`) shifts the BoxArray and every
FAB's index space together, which is precisely the relabelling wanted.
`ParallelCopy` then moves each cell to its new home across ranks and boxes.

**The trap.** `ParallelCopy` writes only the intersection. The `n` cell-layers at
the leading face have no source and are left holding their *pre-shift* values —
stale data that looks entirely plausible. They must be filled explicitly, and a
test must prove they were. This is the single most likely bug in the module.

### 4.3 What has to be shifted, and what has to be reset

| thing | why |
|---|---|
| the new state data | the obvious one |
| **the old state data** | `AmrLevel::RK` starts each step from `get_old_data`; leaving it unshifted mixes two frames and the run dies (or worse, does not) |
| ghost cells | not shifted — regenerated by the next `FillPatch` from the shifted interior and the boundary condition |
| the FillPatcher cache | `resetFillPatcher()` after any shift. It is only used on levels > 1, which the module forbids, but the reset is one line and removes the question |

### 4.4 The two seams, and why neither reaches the stars

With a threshold of 2 units the shift is 4 cells (2 units at `dx = 0.5`):

| face | what happens | where it lives | sponge strength there |
|---|---|---|---|
| leading (`+x`) | 4 cell-layers have no source and are invented | `x` from 62 to 64, i.e. radius ≥ 30 | 32 % to 100 % of full |
| trailing (`−x`) | 4 cell-layers of already-departed wake are discarded | radius ≥ 30 | same |

**Nothing inside radius 30 is ever fabricated or discarded.** The pair sits within
2 of the centre and its 99 % envelope reaches radius 16 — the seams are 14 units
of vacuum away from the nearest matter, inside a band whose entire purpose is to
absorb whatever arrives there.

The transverse faces are untouched: the shift is along one axis only.

**Filling the leading sliver.** Two options, and the module takes the first:

1. *Asymptotic falloff* (default). Each variable relaxes toward its asymptotic
   value as `1/r`: `v_new = v_∞ + (v_src − v_∞)·r_src/r_new`, with `v_∞` read
   from the same `nonzero_asymptotic_vars` / `nonzero_asymptotic_values` the
   Sommerfeld boundary already uses (χ, h11, h22, h33, lapse → 1; everything else
   → 0). This is exactly the model the outer boundary condition already assumes,
   so the fabricated cells are continuous with what the boundary would have
   produced anyway.
2. *Copy the outermost surviving layer* (a switch, for comparison). Cruder, and
   useful precisely because it is different: if the trajectory depends on which
   fill is used, the seam is doing something and the result is not believable.

### 4.5 Deciding when to shift

After every `treadmill_check_interval` steps (default 20 — the pair moves at most
0.06 units in that time even at 0.3 c), one reduction over the grid returns each
sector's core position, and the pair's midpoint is their mean.

**The weight must be the field squared, not the mass.** The canonical star carries
`+M` and the phantom `−M`, so a mass-weighted centroid of the pair divides by a
total that is approximately zero. The weight is `|Φ|²` per sector, restricted to
a search ball seeded at the last known core — the same construction
`RLLumpTracker.hpp` uses to steer the pump's spotlight, and the same core
definition `sector_dynamics.dat` uses so the two agree.

The module gets **its own copy** of that reduction rather than calling into the
reinforcement-learning path, so the two stay uncoupled.

**Threshold: 2 units, not the 8 sketched in `GPU_RUN_PAPER.md` section 5.** Four
reasons, all of them free:

- the fabricated sliver is 4 cells instead of 16;
- the discarded wake shell is 4 cells instead of 16;
- the pair stays within 2 of the centre instead of 8, so the sponge's geometric
  relationship to the stars is restored almost continuously instead of drifting a
  quarter of the way to the boundary between shifts;
- and the cost is nil: 165 shifts over the whole `d = 10` run instead of 41.

**There is no chatter risk.** A shift moves the midpoint from `+2` to `0` and the
pair only ever travels forward, so the threshold cannot be crossed twice in
quick succession. The trigger is monotone by construction.

**Shift accounting for the recommended cell** (`d = 10`, N = 128, 330 units
travelled, threshold 2):

| quantity | value |
|---|---|
| shifts over the whole run | 165 |
| first shift at | `t ≈ 165` |
| interval between shifts near 0.3 c | 6.7 units of `t` ≈ 670 steps |
| data moved per shift | 2 × 0.69 GB (new + old state) |
| time per shift | milliseconds — device-to-device |
| total cost of every shift in the run | well under a second |

### 4.6 The odometer, and the one thing that must survive a restart

The odometer is a single integer — cumulative cells shifted — and it is the only
genuinely new piece of state in the run.

**It has to go into the checkpoint.** Without it a restart resumes with the
trajectory silently reset to zero, and the resulting plot looks completely
plausible while being wrong. That is the worst failure mode available here, so it
gets explicit treatment: written by `specific_pre_checkpoint`, read back by
`specific_post_restart`, and *verified* by checking that the first `treadmill.dat`
row after a restart continues the last row before it rather than starting again.

### 4.7 What this cell cannot answer

The trailing wake is discarded, so this run yields **no radiated momentum and no
wave zone**. It answers how fast the pair goes, and nothing else. The wave-zone
cell stays in the paper matrix as the only source for the radiation field, and the
long box (route A) stays on the books as the only way to have both at once.

Say this in the paper text rather than leaving a referee to notice it.

### 4.8 Why not simply move the coordinates

A comoving gauge would be less code. It is rejected on presentation grounds, not
technical ones. The entire claim is *the pair accelerates*. Measuring that
acceleration in a gauge purpose-built to follow the pair is a far harder result to
defend, and it invites exactly the objection the whole control matrix exists to
kill. Translating the data leaves the evolution equations, the gauge and the
boundary condition untouched; only the labels change.

### 4.9 The honest statement of what changes

It is worth being precise, because the tempting claim ("it is an exact copy, so
nothing changes") is not quite true and a referee will find the gap.

Translating the data *is* an exact symmetry of the evolution equations, which are
translation invariant. It is **not** a symmetry of the sponge and the outer
boundary, which are anchored to the grid. So the recentring box is not the same
simulation as the static box; it is a simulation in which the truncation of the
domain follows the source instead of being left behind by it.

That difference is in the right direction, and for a stated reason: the exact
spacetime is unbounded, every finite box is a truncation, and a truncation
centred on the source is a better approximation than one the source is walking
out of. What the recentring box additionally does is discard outgoing wake before
it can be reflected back inward — which removes an error the static box has,
rather than adding one.

The practical consequence is that **validation cell 3 must not be given a
round-off pass gate**. See §7.

---

## 5. Code plan

The module rule applies without exception: new behaviour lives in its own file,
behind its own switch, writing its own output file, and no existing column
contract is widened.

### 5.1 Files

| piece | file | lines | notes |
|---|---|---|---|
| the shift, the sliver fill, the odometer | **new** `Source/Grids/GridTreadmill.hpp` | ~200 | knows nothing about this campaign — takes a MultiFab, an axis, a cell count and a table of asymptotic values |
| the core tracker and the trigger | **new** `Examples/RadialRecipe/RecentringBox.hpp` | ~180 | sector-aware; its own copy of the ball-restricted `|Φ|²` reduction, not a call into `RLLumpTracker.hpp` |
| parameters | `Examples/RadialRecipe/SimulationParameters.hpp` | ~15 | one block, all defaults off |
| the call site | `Examples/RadialRecipe/RadialRecipeLevel.cpp` | ~20 | end of `specificPostTimeStep()`, level 0 only |
| checkpoint / restart of the odometer | `Examples/RadialRecipe/RadialRecipeLevel.cpp` | ~40 | overrides `specific_pre_checkpoint` and `specific_post_restart` |
| unit tests | **new** `Tests/GridTreadmillTest/` | ~150 | registered in `Tests/TestCases.hpp` like every other test |

Roughly 600 lines, and about 400 of them are the two new headers.

### 5.2 Parameters — every one default-off or default-inert

| parameter | default | meaning |
|---|---|---|
| `treadmill_enabled` | `0` | the switch; every archived cell is unaffected |
| `treadmill_threshold` | `2.0` | shift when the pair's midpoint is this far from the box centre |
| `treadmill_axis` | `0` | which axis the pair runs along |
| `treadmill_check_interval` | `20` | steps between centroid reductions |
| `treadmill_ball_radius` | `8.0` | search-ball radius for the core tracker |
| `treadmill_fill_mode` | `0` | `0` = asymptotic `1/r`, `1` = copy outermost layer |

### 5.3 Output — `small_data/treadmill.dat`, its own file

| column | why it is there |
|---|---|
| `time`, `step` | the join key for everything else |
| `core_x_canon_grid`, `core_x_phantom_grid` | where the tracker thinks the stars are, in grid coordinates |
| `midpoint_grid` | what the trigger actually compares against the threshold |
| `cells_shifted` | 0 on most rows; the shift size on the rows where one happened |
| `offset_cells`, `offset_length` | the odometer, in cells and in code units |
| `midpoint_true` | `midpoint_grid + offset_length` — the trajectory the paper plots |

This deliberately duplicates information that `sector_dynamics.dat` also carries.
That is the point: `sector_dynamics.dat` is computed after the fact by the
plotfile consumer at whatever cadence plotfiles are written, and `treadmill.dat`
is computed in-situ every 20 steps. On a 215,000-step run that is the difference
between a few hundred trajectory samples and ten thousand. **The two must agree
where they overlap, and that is a pass gate (§7), not an assumption.**

### 5.4 Restrictions — refuse rather than guess

The module aborts at startup, with a message naming the reason, if any of these
is true. Every one of them is a case where doing something reasonable-looking
would be worse than stopping:

- `max_level > 0` — the shift is single-level by construction. Every cell in this
  campaign is single-level by decision, so this costs nothing in practice.
- the shift axis is periodic — the whole seam argument assumes it is not.
- the boundary on the shift axis is not Sommerfeld — the sliver fill is built on
  what that boundary assumes.
- the threshold is not a whole number of cells — the exactness claim dies here,
  and this is the one that will actually catch someone.
- `treadmill_threshold + d/2 + r₉₉ > sponge_inner_radius` — the configuration
  would sponge the star. Computable at startup from the parameters, and worth a
  hard stop rather than a slowly poisoned run.

### 5.5 Tests

Unit tests, run on the CPU in seconds, in `Tests/GridTreadmillTest/`:

1. **The shift is exact.** Fill a MultiFab with `f(i,j,k) = i + 1000·j + 10⁶·k`,
   shift by `n`, and require every cell in the overlap to hold the value that was
   `n` cells ahead of it — bit-for-bit, not to a tolerance.
2. **The sliver is written.** Poison the leading `n` layers with a sentinel before
   the shift and require none of it to survive. This is the `ParallelCopy` trap in
   §4.2 and it is why the test exists.
3. **Round trip.** Shift `+n` then `−n` on data whose leading and trailing shells
   are the asymptotic values; require the interior to be bit-identical to the
   original.
4. **The odometer survives a checkpoint.** Write, read back, compare.

Test 3 is the unit-level version of validation cell 1, and it should be green
before any GPU time is spent.

---

## 6. Running it — wrapper, restart, analysis

### 6.1 New launcher knobs

`run_pair_selfgrav.sh` gains four passthroughs, in the same style as the sponge
block it already has:

```
BONDI_TREADMILL=1                # -> --extra-override treadmill_enabled=1
BONDI_TREADMILL_THRESHOLD=2.0    # -> treadmill_threshold
BONDI_TREADMILL_CHECK=20         # -> treadmill_check_interval
BONDI_TREADMILL_FILL=0           # -> treadmill_fill_mode
```

Nothing else in the launcher changes, and a cell that does not set
`BONDI_TREADMILL` is byte-identical to what it is today.

### 6.2 Checkpoints — two things the current setup gets wrong for a long cell

Every cell in this campaign runs with `checkpoint_interval = -1`. That is right
for a 1.1-hour cell and wrong for a twelve-hour one. Switch it on at about
20,000 steps (≈ 200 units of `t`, ≈ 1.1 h of wall clock), which gives eleven
checkpoints over the run at ~1.4 GB each.

**And move them off the scratch disk.** `amr.check_file` currently points into
`/tmp` node-local scratch, and the wrapper's end-of-episode cleanup
(`core/scratch.py::purge_plotfile_scratch`) deletes checkpoint directories
unconditionally — they are treated as transients, with an explicit comment saying
the consumer never reads them. A finished chase cell would therefore lose exactly
the artefact needed to extend it. Override `amr.check_file` to a path inside the
cell's own run directory and prune old checkpoints by hand.

### 6.3 Restarting

The restart key is **`amr.restart`**, not `restart_file`. (`restart_file` is read
in `AMReXParameters.hpp:121` to set a flag that nothing uses; the file-existence
check next to it is inside `#if 0`. Do not reach for it.) `Amr::init` restarts if
and only if `amr.restart` is a non-empty checkpoint path.

Two ways to use it, and the second is the one to use:

1. Through the wrapper with `--extra-override amr.restart=<chkdir>`. It works —
   unknown keys are passed straight through to `params.txt` — but the wrapper
   also insists on a fresh run directory and re-runs the four-hour elliptic solve
   whose output the restart then ignores.
2. Re-launch the binary directly against the finished cell's own `params.txt`
   with `amr.restart` and a raised `stop_time` appended. No solve, no new
   directory, and the run continues into the same output files. This is the right
   shape for a one-off long cell and it needs no wrapper change at all.

Whichever is used, the restart is not believed until the odometer check in §4.6
has been run on `treadmill.dat`.

### 6.4 Analysis

The consumer needs no change. `sector_dynamics.dat` keeps reporting grid
positions, which is correct and should stay that way — a column that silently
switched between grid and true coordinates would be a trap. The join happens in
analysis:

```
true_x(t) = core_x_grid(t) + offset_length(t)
```

with `offset_length` taken from `treadmill.dat` at the same time. Because the
offset is piecewise constant and the shift is exact, the joined trajectory is
continuous across every shift — and **checking that it is continuous, to
round-off, at all 165 shifts is the cheapest and sharpest health check the run
has.** It should be run automatically, not by eye.

Add to `results/bondi-dipole-runaway/analysis/` one new script that does the join,
fits `x(t)` to the hyperbolic form, and reports `v(t)` against both the
relativistic and the Newtonian prediction. It does not touch the existing scripts.

### 6.5 The movie

With recentring switched on, the pair sits still in the middle of every frame and
the runaway becomes invisible — the box is moving with it. Any movie from this
cell must either be captioned with the odometer or rendered with the frame origin
offset by it. Worth deciding before launch, because
[artefact rule 6](../../../grteclyn-wrapper/README.md) applies: frames rendered
without cached slices cannot be redrawn afterwards.

### 6.6 Plotfile cadence

At `plot_interval = 80` steps a 215,000-step run writes 2,690 plotfiles of
0.66 GB each — 1.8 TB through the scratch disk, for a trajectory that
`treadmill.dat` is already recording ten times more densely. Set
`BONDI_PLOT_INTERVAL=400` (one plotfile every 4 units of `t`, 538 of them). That
still leaves the consumer's `sector_dynamics.dat` with more than enough samples to
cross-check the in-situ tracker.

---

## 7. Validation — nothing is believed before these

Each cell is 1.1 GPU-hours unless stated. Each has a gate, and the gate is
checked before the next cell starts. Total before the chase: about 5.5 GPU-hours.

### 0. Unit tests (free, CPU, seconds)

`Tests/GridTreadmillTest/` green, all four cases of §5.5. If test 2 (the sliver
poison test) is not present, the ladder does not start.

### 1. Seam damage — forced shift with zero net displacement (2 cells, 2.2 h)

The headline cell (`d = 10`, N = 128, `t = 200`) with the treadmill forced onto a
fixed schedule that nets to zero: 4 cells forward, 4 cells back, forever. The pair
ends where it started on the grid, so **every difference from the archive is seam
damage and nothing else**.

Run it twice, at two schedules:

| cell | forced shift every | seam pairs by `t = 200` |
|---|---|---|
| 1a | 2,000 steps | 10 |
| 1b | 500 steps | 40 |

*Gate.* Both reproduce the archived trajectory. Any residual must **not** grow
with the seam count — if 1b's departure is about four times 1a's, the seams are
injecting something and the design is wrong, however small the number looks at
`t = 200`. That scaling test is the point of running two; a single cell can hide
a systematic inside a tolerance.

### 2. Physics at low speed — recentred against the archive (1 cell, 1.1 h)

The headline cell with recentring genuinely armed and `treadmill_threshold = 2`,
so about two shifts land inside `t ≤ 200`.

*Gate — and read §4.9 before setting it.* **Not** round-off. The recentring box
suppresses returning wake that the static box keeps, so the two are not the same
run and demanding agreement to `1.4e-04` would fail a correct implementation. The
defensible gate is that recentring must not be the dominant error:

- drift and fitted acceleration agree with the archive to **better than 1 %**
  (for scale: doubling the box moved the drift 4 %, and N = 128 → N = 192 moved
  the acceleration 10 %);
- the in-situ `treadmill.dat` midpoint and the consumer's `sector_dynamics.dat`
  midpoint agree wherever both exist — this is what licenses using the dense
  in-situ trajectory for the long run;
- the joined trajectory is continuous, to round-off, across every shift.

### 3. Schedule independence (1 cell, 1.1 h)

The same cell at `treadmill_threshold = 4`, which halves the number of shifts and
doubles the size of every seam.

*Gate.* Agrees with cell 2 to better than 1 %. **This is the sharpest test in the
ladder** and it is internal — it needs no archive, no other resolution and no
assumption about what the right answer is. If the physics is invariant, when and
how often the box is re-centred cannot matter. If it does matter, the seam is
doing something and the result is not usable at any speed.

### 4. Fill independence (1 cell, 1.1 h)

Cell 2 repeated with `treadmill_fill_mode = 1` — the leading sliver filled by
copying the outermost surviving layer instead of by the asymptotic law.

*Gate.* Agrees with cell 2 to better than 1 %. Two deliberately different wrong
answers in a region that should not matter, giving the same trajectory, is good
evidence that the region does not matter.

### 5. Boundary sensitivity **at speed** (1 cell, 1.1 h, run from a checkpoint)

Cells 1–4 all run at `t ≤ 200`, where the pair is barely moving; none of them can
detect the thing that actually worries a referee, which is that the outer boundary
assumes a static centre and a pair at a third of light speed is not static.

The cheap version of that test uses the chase cell's own checkpoints: from the
checkpoint at `t ≈ 1,200` (`v ≈ 0.19 c`), restart a 200-unit branch with the
sponge's inner radius moved from 24 to 20 and compare against the main run over
the same window.

*Gate.* The two branches agree on velocity to better than 1 %. If moving the outer
treatment inward by four units changes the speed at 0.2 c, the outer region is
participating and the result needs the bigger box.

### 6. Box size at speed — the expensive one, and what it really costs

The proper version is a recentring cell at `L = 96` / `N = 192` (same `dx`, so the
same physics at a larger radius) compared against `L = 64` / `N = 128` over a
window at speed. At 3.4× the cells that is about **19 GPU-hours per 1,000 units of
`t`**, so a comparison window ending at `t = 700` costs ≈ 13 GPU-hours.

The ideal version — comparing against the long box, which has no seams at all —
is better science and worse engineering: it needs per-direction cell counts
(`N1_full`/`N2_full`/`N3_full`, which GRTeclyn supports and the wrapper does not
currently emit, and which are mutually exclusive with the `N_full` the template
always writes) and an elliptic solve box large enough to cover a 160-long domain.
Budget it as a stretch goal, not a gate.

**Recommendation: gate on cell 5, book cell 6 only if cell 5 fails or if a referee
asks.**

---

## 8. The risk that dominates everything — is `a` actually constant?

Every cost and every milestone above assumes the pair holds its acceleration. The
campaign's own data does not yet establish that, and this is a larger threat to
the chase cell than anything in the code.

### What the archive shows

Coordinate separation at `t = 200`, `d₀ = 10`:

| rung | separation at `t = 200` | opening |
|---|---|---|
| N = 128 | 10.107 | +0.107 |
| N = 192 | 10.029 | +0.029 |
| N = 256 | 10.013 | +0.013 |

The gap opening shrinks by a factor of eight across the ladder, which says most of
it at N = 128 is discretisation, not physics — good news, and consistent with the
rigid pair the paper claims.

But at N = 128 it **grows very fast with time**. In the `t = 400` cell the
separation runs 10.000 → 10.107 (`t = 200`) → 10.773 (`t = 320`) → 11.844
(`t = 400`): the opening is 17 times larger at `t = 400` than at `t = 200`, which
is roughly a fourth power of time over that window.

`a ∝ 1/d²`, so an opening to 11.84 should have cut the acceleration by 29 %. The
reported acceleration over the same window did not move — it sat at 1.418e-04 over
the last two thirds and 1.444e-04 in the final quarter. **Those two statements
cannot both be clean.** Either the separation measure is contaminated at late
times, or the acceleration fit is dominated by the earlier part of its window, or
the pair is not obeying `a = GM/d²` any more.

And the cell that would settle it is the same cell that was, by then, drifting its
leading star into the sponge (§2). It cannot settle it.

### What this means for the chase

- The `t` and cost figures in §1 are **upper-bound-shaped if `a` holds and
  meaningless if it does not.** If the separation really opens, `a` falls, and
  0.3 c arrives late or not at all.
- The chase cell is not a bet on the answer; it is the measurement of it. It logs
  separation continuously and is readable at every milestone, so a decaying
  acceleration shows up early — around `t ≈ 500`, long before the budget is spent.
- The right response to a decaying `a` is not to abandon the run. A pair whose gap
  opens while it accelerates is itself a result, and the recentring box is the only
  instrument that can watch it happen.

### The cell that de-risks it, and it is worth buying

**`chase_pair_d10_recentred_L64_N192_lev0`, to `t = 400`, ≈ 11 GPU-hours.**

At N = 192 with recentring, the run reaches `t = 400` with the pair still at the
box centre and the sponge nowhere near it — the first uncontaminated look at the
late-time separation the campaign has ever had. It answers, in one cell, both
open questions: whether the gap opening converges away at length as it does at
`t = 200`, and whether the treadmill's trajectory converges at the same rate the
static ladder did.

Run it **before** the 12-hour chase, not after.

---

## 9. Cost and order of work

| phase | what | GPU-hours | CPU |
|---|---|---|---|
| 1 | write `GridTreadmill.hpp`, `RecentringBox.hpp`, params, call site, checkpoint hooks, unit tests | 0 | — |
| 2 | validation cells 1a, 1b, 2, 3, 4 | 5.5 | 5 solves, reusable from the archive if the config is identical |
| 3 | `chase_pair_d10_recentred_L64_N192_lev0` to `t = 400` — the acceleration-law cell | 11 | 1 solve (~4 h, N = 384) |
| 4 | `chase_pair_d10_v03c_recentred_L64_N128_lev0` to 0.3 c | 12 | 1 solve (~4 h) |
| 5 | validation cell 5, from the chase's own checkpoint | 1.1 | none |
| | **total** | **≈ 30** | ≈ 8 h, overlappable |

Under two days of card time in total, and the phases are sequential by design —
the gate between each is the entire value of doing it this way. Phase 2's five
cells are independent of one another and can share the cards; phases 3 and 4 are
one cell each.

Launch manually, one cell at a time, as every other phase of this campaign has
been. Do not build a queue for this.

### If it works, what the next target costs

The recentring box's cost is linear in run time and its memory never changes, so
the follow-up targets are simple arithmetic rather than a new design:

| target | `t` needed | distance | GPU-hours | memory |
|---|---|---|---|---|
| 0.30 c | 1,970 | 302 | 11 | 6 GB |
| 0.50 c | 3,620 | 969 | 20 | 6 GB |

(at `d = 10` with the converged `a = 1.596e-04`). Past about 0.6 c the
Lorentz-contracted star needs a finer cell and the arithmetic changes; at 0.3 c
the contraction is 4.8 % and no refinement is needed, which is one more reason
0.3 c is the right first target.

---

## 10. Failure modes and what each one looks like

| symptom | almost certainly | how it is caught |
|---|---|---|
| trajectory jumps by exactly the shift size | odometer sign or units wrong | the continuity check at the first shift, `t ≈ 165` |
| NaN on the step after the first shift | the old state data was not shifted (§4.3) | immediately, by the existing NaN check |
| plausible trajectory, wrong slope, no jumps | the leading sliver kept its stale values (§4.2) | unit test 2 — which is why it exists |
| slow drift that scales with the number of shifts | seam damage | validation cell 1b against 1a |
| constraint norm steps up at every shift | the fill disagrees with what the boundary expects | validation cell 4; if both fill modes step, the sliver is too far in and the threshold must come down |
| restart resumes with the trajectory reset to zero | the odometer was not checkpointed | the first-row-continues-last-row check (§4.6) |
| checkpoints missing after the run finishes | the wrapper purged them from scratch (§6.2) | prevented, by putting `amr.check_file` in the run directory |
| separation opens, acceleration decays | physics, not code | §8 — report it; it is a result |

---

## 11. What needs deciding before any of this starts

1. **`d = 8` or `d = 10`?** `d = 8` is 7.5 GPU-hours, `d = 10` is 12. The
   recommendation is `d = 10`, because it is the configuration the ladder
   converged on and every control is read against, and because at `d = 8` the two
   stars' envelopes substantially overlap.
2. **Buy the acceleration-law cell first (11 GPU-hours) or go straight to the
   chase?** The recommendation is to buy it. It is the only clean look at the
   late-time separation the campaign can get, and without it §8 stays open no
   matter how the chase comes out.
3. **Gate on the cheap boundary check (cell 5) or the expensive one (cell 6)?**
   The recommendation is cell 5, with cell 6 held in reserve.
4. **Stop at 0.3 c, or keep the card and continue to 0.5 c for another nine
   hours?** This can be decided while the run is in flight; the checkpoints make
   it a free option rather than a commitment.

Nothing in this document is on the required matrix for the current paper. It is
the follow-up paper's headline, and the recentring box is the instrument it needs.
