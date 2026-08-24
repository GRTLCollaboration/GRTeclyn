# The paper campaign — plan and run-book

**Goal.** Produce every number, figure and table in
`research/bondi_dipole/bondi_dipole.tex` from runs on the corrected code path,
so that nothing in the paper rests on the old campaign — which ran with the
double-subtracted potential (canonical stars were collapsing) and the displaced
initial-data metric (stars born off the centre of their own wells). The old
paper's "gap closing", "contact contamination" and "halo bias" storyline was
built on those artefacts; the corrected result is cleaner and stronger: **both
stars accelerate together at a = GM/d², the separation stays put, total
momentum stays zero.**

**Method.** One cell at a time. Each phase has a launch command, a cost, and a
pass gate; nothing in the next phase starts until the gate is green. No
pipeline script — the checking between steps is the point.

The eight artefact rules at the top of `grteclyn-wrapper/README.md` are the
constitution of this campaign; this plan just applies them cell by cell.

---

## The scoreboard — every run the paper needs

Tick a box only when the cell has passed its gate and been moved into
`runs/bondi/staging/archive/`. Details for each are in the phase it belongs to.

**Already in the archive (corrected path, gates passed):**

- [x] `runaway_pair_d08_L64_N128_lev0` — separation scan point, strong-signal pictures
- [x] `runaway_pair_d10_L64_N128_lev0` — the headline cell; base rung of the ladder
- [x] `runaway_pair_d12_L64_N128_lev0` — separation scan point
- [x] `runaway_pair_d16_L64_N128_lev0` — separation scan point
- [x] `runaway_pair_d20_L64_N128_lev0` — **evolved to t=200 on 2026-08-24**; the widest separation, run to decide whether the close-pair excess in a·d² is physics or a systematic floor. **It is physics**: a·d² = 0.01436 on the t ≥ 5 convention, **0.07% from GM**, and the five-point power law tightens to d^−2.028 (from −2.046 on four points). Drift +0.7106, separation 20.000 → 20.013, px_total 4.9e-06
- [x] `control_lone_canonical_L64_N128_lev0` — single-star null (box centre)
- [x] `control_lone_phantom_L64_N128_lev0` — single-star null (off-centre, the sharper test)
- [x] `stability/canonical_w{075,080,085,090}_L64_N128_lev0` — stability survey, t=120

**Required — the paper is not submittable without these:**

- [x] `smoke_mpi_evo` — phase 0, **passed 2026-08-22**: 50/50 steps on cards 0+1, no segfault, exit 0. Two-GPU evolution is available at `max_level=0`. Cell deleted.
- [x] `control_pair_pp_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; barycentre 32.00073, drift +0.00073, min χ = 0.97947. Gap unmeasurable (see phase 1)
- [x] `control_pair_mm_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; barycentre 31.99972, drift −0.00028, min χ = 1.00000 — which means only that no *well* forms; the phantom stars make hills the diagnostic cannot see. Gap unmeasurable from the tracker (see phase 1); the `_frames` re-run supplies it
- [x] `control_mirror_mp_d10_L64_N128_lev0` — phase 1, **evolved to t=200 on 2026-08-22**; runaway reverses exactly: displacement and acceleration both −1.00002× the archived cell. Frameless by decision
- [x] `runaway_pair_d10_L64_N192_lev0` — phase 2, **evolved to t=200 on 2026-08-22**; drift +3.0139, separation 10.000 → 9.930, a = 1.611e-04. Middle rung of the ladder
- [x] `runaway_pair_d10_L64_N256_lev0` — phase 2, **evolved to t=200 on 2026-08-22**; drift +3.0016, separation 10.000 → 9.915, a = 1.596e-04. **The ladder converges**: N192 and N256 agree to 0.4% on drift and 0.9% on acceleration, N128 sits 4% low. Frames + slice cache present (headline movie)
- [x] `control_pair_pp_d10_L64_N192_lev0` — phase 2, **evolved to t=200 on 2026-08-22**; solved by the `converged` door at the flat 0.002 gate. min χ = 0.97920. The two stars merge at t ≈ 35 and the box activity ×7 is merger ejecta, not an artefact — see "the same-sign pairs merge" below; centroid still pinned to ±7e-4 over the full run
- [x] `control_pair_pp_d10_L64_N256_lev0` — phase 2, **evolved to t=200 on 2026-08-23**; solved by the `converged` door in 13 of 50 NL iterations (Ham 0.00086%, Mom 0.00102% against the 0.002 gate), 28 min. min χ = 0.97914. Frames + slice cache present — **the frames are the per-star measurement**: the two wells close 8.78 → merged at t = 33.6, then the remnant rings down (the chi movie's collide-and-bounce). See "the same-sign pairs merge" below
- [x] `massscale_pair_d10_w0804_L64_N128_lev0` — phase 3, **evolved to t=200 on 2026-08-22**; lighter phantom (M = −0.011472, 79.95% of matched). Separation-corrected fit gives the canonical star's pull at 0.809 of the matched cell against 0.7995 predicted (1.2%) and the phantom's at 0.973 against 1.000 (2.7%) — the pull follows the partner's mass. The pair is no longer rigid: the gap closes 10.000 → 9.408. Ratios only, never the constants — see phase 3
- [x] `wavezone_pair_d10_L128_N256_lev0` — phase 4, **evolved to t=200 on 2026-08-22**; drift +2.7606 (4% from the L=64 N128 cell, so the box is not driving the result). **The r·ψ₄ gate fails, informatively**: ψ₄(l=2) falls as r^−4.0, not r^−1, across R = 16/24/32/40 — all four shells are in the near zone and no radiative tail is measurable. Reported as a null, see phase 4

**Optional — only if the paper wants the figure:**

- [x] `longrun_pair_d10_t400_L64_N128_lev0` — **evolved to t=400 on 2026-08-22**; drift +11.5177, a = 1.418e-04 over the last two thirds and 1.444e-04 in the final quarter — **the acceleration is steady, not growing**, which retires the late-time uptick seen at t=200. Separation opens 10.000 → 11.741. Frames + slice cache present
- [x] `control_pair_mm_d10_L64_N192_lev0` — **evolved to t=200 on 2026-08-23**; completes the null ladder. Centroid drift −0.00035 over the full run against the PP rung's +0.00048 at the same grid — both four orders of magnitude below the runaway's +2.88, and flat rather than shrinking across resolution, which marks them as the measurement floor. min χ exactly 1.00000 — the diagnostic's blindness, not flat geometry: phantom stars make *hills* (χ up to 1.011 on the R=8 shell) and min χ only sees wells. Boundary growth 7.0×, matching the PP rung digit for digit. The `_frames` re-run answered the per-star question on 2026-08-24: **the phantom pair merges like the PP pair** (single lump from t = 32.8, against PP's t = 33.6 — a 2.4% difference across a grid change) — the field-overlap force is sign-blind and ~35× gravity; see "the two interaction channels"
- [x] `control_pair_mm_d10_L64_N128_lev0_frames` — **evolved to t=200 on 2026-08-24**; the archived MM cell's physics re-run with frames on, because frames are the *only* per-star measurement a same-sign pair has. Centroid drift −0.00022 over the full run (in line with −0.00026/−0.00035 on the frameless rungs). **The phantom pair merges**: the two |χ−1| hills close 8.57 → 7.56 (t=20) → 4.54 (t=32) and are a single lump from t = 32.8, against the PP pair's t = 33.6. Gravity's sign is opposite between those two cells, so the driving force is sign-blind — this is the campaign's sharpest mechanism control. Tracking in `well_tracking.dat` via `analysis/track_wells.py`; 19 movies rendered at fixed colour scale
- [~] `control_pair_mm_d10_L64_N256_lev0` — ~17 h; MM null, finest rung. **Deliberately not run** (2026-08-23): 17 GPU-hours and 58 GB for a third point on a null ladder whose first two rungs already agree, when its first two rungs already agree digit for digit. Launch script is staged if the decision is revisited
- [x] `massratio_heavyphantom_d10_L64_N128_lev0` — phase 3b, **evolved to t=200 on 2026-08-23**; reversed mass ordering (canonical ω = 0.81, M+ = 0.010721, |M−|/M+ = 1.333). **Gate PASSED on the sign**: the gap OPENS, 10.000 → 10.603, against 10.000 → 10.003 matched and 10.000 → 9.408 for the lighter phantom. See phase 3b result
- [x] `massratio_w088_r060_d10_L64_N128_lev0` — phase 3b, **evolved to t=200 on 2026-08-23**; phantom ω = 0.88, ratio 0.597. Canonical pull ratio 0.597 against 0.597 predicted (−0.1%, the tightest in the campaign); gap closes 10.000 → 8.618. See phase 3b result
- [x] ~~`massratio_*_r040_d10_L64_N128_lev0`~~ — phase 3b; **cannot be built.** The phantom branch has a floor at |M−| ≈ 0.00791 (0.55 of matched, near ω = 0.94) and a 0.40 rung would need 0.00574. The bound is itself a result — see phase 3b
- [x] `amrcheck_pair_d10_L64_N128_lev1` — **evolved to t=200 on 2026-08-23**; referee-proofing. **Level 1 was never created** — the tagger never fired, exactly as predicted — and drift and acceleration match the uniform cell to **0.001% and 0.002%**. Mesh refinement changes nothing here
- [x] `deepsolve_pair_d10_L64_N128_lev0` — **evolved to t=200 on 2026-08-24**; referee-proofing the other way. The headline cell with the elliptic tolerance alone tightened 8.6e-04% → 2e-04% (exit at pass 16, Ham 1.94e-04%, a 4.4× deeper solve). The drift is unchanged to four decimals over the whole run (2.88122 vs 2.88117 at t = 200, 0.0015%; a = 1.4634e-04 in both), while N128 → N192 changes it by 4.6%. **The base rung is grid-limited, not solve-limited** — and the t = 0 constraint violation barely moved (2.022e-04 → 1.999e-04), because it is grid-transfer noise. See "the deep-solve twin"
- [ ] `chase_pair_d08_v03c_Lx352_L64_N128_lev0` — ~1.7 GPU-days; ride the runaway to 0.3c in a long box. **Superseded by the recentring box — implemented and implementation-tested 2026-08-23; plan and build record in `CHASE_TO_03C.md`.** Follow-up paper material either way

---

## Findings for the paper — what each cell established, verified against the packed data

Every number below was recomputed on 2026-08-23 from
`results/bondi-dipole-runaway/campaign/` (the paper-facing extracts), not
copied from earlier notes. `a` is the quadratic coefficient of the pair
midpoint over the last two thirds unless a window is named.

### The one-line verdicts

| cell | for the paper |
|---|---|
| `runaway_pair_d10_L64_N128` | the headline at base resolution: drift +2.8815, a = 1.4634e-04, px_total 3.7e-05 |
| `runaway_pair_d10_L64_N192` | drift +3.0139, a = 1.6106e-04 — brackets the finest rung from above |
| `runaway_pair_d10_L64_N256` | **the quotable cell**: drift +3.0016, a = 1.5958e-04, px_total 3.8e-06 |
| `runaway_pair_d{08,12,16,20}` | with d10: a ∝ d^−2.028 over five separations; a·d² → GM to 0.1% from d = 12 out, so the close-pair excess is finite-size, not a floor |
| `control_mirror_mp` | swap the sectors, the runaway reverses: drift ratio −1.000022, a ratio −1.000025 |
| `control_lone_{canonical,phantom}` | a star alone moves ≤ 1.8e-03 in 200 — the single-star noise floor |
| `control_pair_pp` ×3 rungs | two positive stars merge at t ≈ 35; their centroid moves ≤ 7.3e-04 all run |
| `control_pair_mm` ×2 rungs + `_frames` | centroid ≤ 3.5e-04; the pair **merges at t = 32.8**, against PP's 33.6 — the driving force is sign-blind. min χ = 1.00000 is the diagnostic's blindness to hills, not flat geometry |
| `massscale_w0804` | phantom ×0.7995 → canonical pull ×0.809 (pred 0.7995) |
| `massratio_w088_r060` | phantom ×0.5974 → canonical pull ×0.597 (pred 0.5974, the tightest point) |
| `massratio_heavyphantom` | mass ordering reversed → the gap **opens** +0.6028; phantom pull ×0.739 (pred 0.747) |
| `amrcheck_lev1` | refinement on: drift Δ 0.0008%, a Δ 0.0020%, level 1 never created |
| `wavezone_L128` | box doubled: drift −4.2% vs L = 64; ψ₄(l=2) ∝ r^−4.0 — no wave zone, no radiation measured |
| `longrun_t400` | a within 2% of 1.43e-04 in every late window; velocity still growing at t = 400 |
| `deepsolve_pair_d10_L64_N128` | solve tolerance 4.4× tighter changes the drift by 0.0015% (2.88117 → 2.88122) — the base rung is grid-limited, not solve-limited |
| `stability_w{075..090}` | the canonical family is stable where it is used |

### The two interaction channels — the paper's defence, as headlines with the runs behind them

The claim under attack will be "something other than gravity moves the pair".
The defence is that the model has exactly two interaction channels, each pair
type isolates one of them, and the campaign measured both. Every bullet names
the run(s) that carry it.

**Channel 1 — the field's own overlap force (same-sector pairs only).**

- *Two lumps of the same field, oscillating in phase, pull together through
  direct field overlap — no gravity required.* Both same-sector pairs fall
  toward each other and merge: the canonical pair at t = 33.6
  (`control_pair_pp_d10_L64_N256_lev0`, well tracking packed as
  `well_tracking.dat`), the phantom pair at t = 32.8
  (`control_pair_mm_d10_L64_N128_lev0_frames`, the frames re-run).
- *This force is ~35× gravity at d = 10.* Newtonian free-fall from 10 apart at
  these masses takes t ≈ 200; both pairs merged in ~35. Force ratio
  ≈ (200/35)² ≈ 35.
- *It is blind to the sign of the gravitational coupling — the single cleanest
  statement in the campaign.* Newtonian gravity's sign is **opposite** between
  the two same-sector cells: two positive masses attract, two negative masses
  mutually repel (Bondi's own result, and the reason "−− should fly apart" was
  the naive expectation). If gravity set the timescale, PP would coalesce while
  MM separated. Measured with the same detector on the same quantity
  (`analysis/track_wells.py` on each cell's χ slice cache), they instead merge
  **within 2.4% of the same time** — 33.6 against 32.8 — with matching gap
  trajectories throughout (8.78/8.57 at t = 0, 7.75/7.56 at t = 20, 6.00/5.54 at
  t = 30). Read 2.4% as "the same timescale", not as a resolved difference: the
  PP cell is N256 and the MM cell N128, so the profiles are sampled at
  dx = 0.25 and dx = 0.5.
- *Consequence: the same-sector cells are not gravity tests and the paper must
  not present them as the ++/−− entries of a Bondi sign matrix.* Their value is
  (a) they calibrate the field channel — fast, short-range, sign-blind — and
  (b) they are momentum nulls: the pair centroid stays pinned to ≤ 7.3e-04
  through a violent merger and ringdown (`control_pair_pp_*`,
  `control_pair_mm_*`, all rungs).

**Channel 2 — gravity (the mixed pair's only channel).**

- *By construction.* Canonical and phantom are two separate complex fields;
  each self-interaction potential is evaluated on its own field's modulus only
  (`Source/Matter/ComplexScalarField.impl.hpp`,
  `Source/Matter/ComplexExoticScalarField.impl.hpp` — `compute_potential`
  takes one field's components; no term in the Lagrangian contains both
  fields). The metric is the only object both fields touch.
- *By measurement, three ways, none requiring trust in the Lagrangian:*
  1. **Same separation, opposite fate.** At the same d = 10 where same-field
     pairs merge in ~35, the mixed pair holds its gap for the whole run,
     10.000 → 10.003 (`runaway_pair_d10_L64_N128_lev0`). A few-percent leak of
     the 35×-gravity overlap force across sectors would have slammed the gap
     shut.
  2. **Range.** A field-overlap force dies exponentially with distance; gravity
     falls as 1/d². The measured force is a clean power law, a ∝ d^−2.05
     across d = 8/10/12/16 (`runaway_pair_d{08,10,12,16}_L64_N128_lev0`), with
     a·d² returning the star's mass to ~1% at d ≥ 12.
  3. **What sets the strength.** A field force scales with field overlap and
     knows nothing about ADM mass. The measured pull tracks the *partner's ADM
     mass* through a factor 2.2 in mass ratio — including the sign flip of the
     gap at the equal-mass point (`massscale_pair_d10_w0804`,
     `massratio_w088_r060`, `massratio_heavyphantom`, all `_L64_N128_lev0`).
- *The phantom's exotic nature enters through exactly one thing: the sign with
  which its energy curves space.* Its field dynamics, internal glue and
  response to curvature are standard. Visible in one pair of frames: the
  phantom pair's merged remnant is a **hill** in χ (χ − 1 up to +0.015,
  `control_pair_mm_..._frames`, `chi_minus_1` at t = 38.4) where the canonical
  remnant is a **well** (χ − 1 ≈ −0.02, `control_pair_pp_..._N256`). Same
  matter behaviour, opposite curvature — the model's one sign, photographed.

**The sentence the two channels buy:**

> Same-sector pairs share a field, and the field's own overlap force — ~35×
> gravity at this separation — dominates their dynamics identically regardless
> of the sign of their gravitational coupling. The mixed pair's two stars live
> in different fields with no coupling term between them, so curvature is
> their only interaction channel: the runaway is gravitational by
> construction, and by measurement.

**Honest residual, so the referee finds nothing we did not state first:** the
canonical star's field tail does reach the phantom's position, but what it
delivers there is energy, and energy curving space *is* the gravitational
channel. Extended sources in each other's curvature give the small excess in
a·d² — +2.9% at d = 8 falling to 0.1% by d = 12 — the fingerprint of a
finite-size gravitational correction fading with distance. That reading was put
to the test by `runaway_pair_d20_L64_N128_lev0`, run precisely because a *floor*
in the measurement would have looked the same at d = 12–16: at d = 20 the ratio
a·d²/GM is 1.000, so it keeps fading and is not a floor.

### Convergence, honestly: what Richardson can and cannot give here

The ladder triple is **non-monotone** — N128 sits below, N192 slightly above
N256 — so the apparent-order ratio (f₁₂₈−f₁₉₂)/(f₁₉₂−f₂₅₆) is negative (≈ −10)
and **no formal convergence order can be quoted**. The reason is visible in the
constraint norms and is worth a sentence in the paper: the t = 0 violation
*rises* with resolution (L2 Ham 1.06e-04 → 1.46e-04 → 1.80e-04), because the
initial-data interpolation noise scales as 1/dx², while the evolution's own
violation *falls* with resolution (late-time 1.06e-05 → 6.6e-06 → 5.6e-06).
Two error sources with opposite resolution behaviour cross somewhere near the
base rung, and the ladder never enters a single-power asymptotic regime.

What survives is still strong:

- The finest two rungs agree to 0.4% (drift) and 0.9% (a), and the evolution
  error falls grid over grid.
- Richardson on the fine pair alone (ratio 4/3): assuming 2nd order gives
  a∞ = 1.577e-04, assuming 4th gives a∞ = 1.589e-04 — both corrections are
  *smaller than the N192–N256 spread itself*.
- **Quote:** drift 3.00 ± 0.01 at t = 200, a = (1.60 ± 0.02)e-04, error bar =
  the fine-pair spread, which dominates any extrapolation correction.

One nuance the coarse grid hides: at N128 the separation returns to 10.003
("rigid pair"), but both fine rungs agree it actually *closes* slightly, to
9.930/9.915. Rigidity is a 1%-level statement and the fine-grid value is the
honest one.

### The N128 offset is not solve residual — the deep-solve twin

The obvious referee question about a non-monotone ladder is "your base rung is
off because its initial data was not converged far enough — solve it harder."
That is now answered with a run rather than an argument.

`deepsolve_pair_d10_L64_N128_lev0` is the headline cell with **one** parameter
changed (verified by diffing against the packed `launch_config.sh`): the
elliptic solver's tolerance is tightened from 8.6e-04% to 2e-04%, and its
stall tolerance with it. The solver duly went deeper — it exited at pass 16 on
Ham 1.94e-04%, a 4.4× tighter residual, and its convergence ladder reproduces
the original digit-for-digit through pass 12, so the two solves are the same
solve followed further.

The evolution does not notice. Midpoint drift, both cells run to t = 200 and
compared on the same time grid:

| t | original N128 | deep-solve N128 | difference |
|---|---|---|---|
| 50 | 0.17283 | 0.17283 | +0.000002 |
| 100 | 0.71272 | 0.71273 | +0.000011 |
| 150 | 1.61519 | 1.61522 | +0.000025 |
| **200** | **2.88117** | **2.88122** | **+0.000044** |

The acceleration over the last two thirds agrees to four digits as well:
1.4634e-04 in both cells (0.001%). The *matter* agrees too, which is the check
worth making before believing any geometry diagnostic: confined fraction
0.276 → 0.269 and rms radius 6.642 → 7.237 over the run, identical in both
cells to three decimals. (That confined fraction is not dispersal — it reads
0.276 already at t = 0, because the diagnostic's confinement sphere is smaller
than a pair separated by 10. The scorer's "matter DISPERSED" note fires on the
headline cell in exactly the same way.)

**A 4.4× deeper solve moves the runaway drift by 0.0015% at t = 200. Refining
the grid N128 → N192 moves it by 4.6% — three thousand times more.** The
resolution offset is therefore a property of the evolution grid, not of how far
the constraints were driven down on the initial slice.

The constraint norms say the same thing and explain why. The deeper solve
lowered the t = 0 violation on the *evolution* grid by only ~1% (L2 Ham
2.022e-04 → 1.999e-04), and by t = 30 the two are identical (5.00e-06). The
t = 0 violation is dominated by the noise the solve-grid → evolution-grid
transfer injects, which scales as 1/dx² and which a tighter solver tolerance
cannot touch. **Solve depth is not what separates the rungs.** The practical
consequence for the paper: there is no cheap fix at fixed resolution, and the
fine rungs are the only honest ladder — which is exactly how the error bar
above is constructed.

### The drift is real — the six-way confirmation

| perturbation | change in the answer |
|---|---|
| cell size ×(3/4)² | +4% then −0.4% (converged) |
| mesh refinement on | 0.0008% — and the tagger never fired: **max_level changes nothing** |
| box volume ×8 | −4.2% |
| sector swap | ×(−1.000022) |
| run 2× longer | a within 2% of the t ≤ 200 value; velocity still growing linearly |
| initial-data solve 4.4× deeper | 0.0015% on drift, 0.001% on a — the base rung is not solve-limited |

And the null side: lone stars ≤ 1.8e-03, same-sign centroids ≤ 7.3e-04 —
2000–4000× below the signal. Total momentum stays consistent with zero while
the pair displaces by 3: px_total = 3.7e-05 at N128 falling to **3.8e-06 at
N256** — displacement without momentum is the Bondi signature, and it converges
*toward* zero, not away from it.

### The force law — five separations, and the excess goes to zero

`a·d²` should return the partner's mass at every separation. Against
GM = 0.014350, fitting each run over its last two thirds:

| d | 8 | 10 | 12 | 16 | **20** |
|---|---|---|---|---|---|
| a | 2.341e-04 | 1.463e-04 | 1.005e-04 | 5.664e-05 | **3.637e-05** |
| a·d² | 0.01498 | 0.01463 | 0.01448 | 0.01450 | **0.01455** |
| vs GM | +4.4% | +2.0% | +0.9% | +1.0% | **+1.4%** |

and on the `t ≥ 5` convention the pack README uses, which is the cleaner
statement of the same thing:

| d | 8 | 10 | 12 | 16 | **20** |
|---|---|---|---|---|---|
| a·d² | 0.01476 | 0.01448 | 0.01437 | 0.01434 | **0.01436** |
| vs GM | +2.9% | +0.9% | +0.1% | −0.1% | **+0.07%** |

**Power-law slope across all five: d^−2.028**, against −2 exact (the four-point
value was −2.046; adding d = 20 moved it *toward* the exact exponent).

`runaway_pair_d20_L64_N128_lev0` was run specifically to decide whether the
close-pair excess is physics or a systematic floor — if `a·d²` had plateaued
around 1% above GM, that would have pointed at a floor in the measurement. It
does not plateau. The excess decays monotonically and lands within 0.1% of GM
from d = 12 outward: the deviation is a finite-size correction, largest exactly
where the two stars' envelopes overlap most, and it behaves like one.

Total momentum halves with each step out, tracking the acceleration rather than
sitting at a fixed floor: px_total = 5.8e-05, 3.7e-05, 2.3e-05, 9.8e-06,
**4.9e-06** at d = 8/10/12/16/20 — the pair displaces by 0.7 to 4.6 while its
momentum stays four to five orders of magnitude below that.

One caveat on quoting a single number: the two fit conventions differ by
0.8–1.5% (a late-window fit weights the part of the run where the pair has
drifted furthest). The *trend* is the robust statement; any single `a·d²`
should be quoted with its window named.

### The mass ladder — the sign flip and its systematics

Gap change at t = 200 against phantom/partner mass ratio: +0.6028 (1.333),
+0.0030 (1.000), −0.5919 (0.7995), −1.3816 (0.5974). A straight line crosses
zero at ratio **1.06**; the gap opens, sits still, or closes strictly according
to which star is heavier.

The separation-corrected per-star fit carries a measurable systematic, and the
paper should state it rather than hide it. The three *changed-partner* ratios —
the scaling law itself — land at 0.809 / 0.597 / 0.739 against 0.7995 / 0.5974
/ 0.7471 predicted (1.2%, 0.05%, 1.1%). The *unchanged*-partner ratio in each
cell is the internal control, and it degrades as the pair deforms: 0.973
(gap −0.59), 1.052 (gap +0.60), 0.725 (gap −1.38). The fit is trustworthy while
the pair stays near-rigid and visibly is not once the gap moves by >1;
ratios of the changed partner are the quotable numbers, constants never are
(on the matched cell the same fit returns 0.0161 and 0.0119 for two constants
both known to be 0.01435).

### What t = 400 adds

a stays inside 2% of 1.43e-04 in every window (133–266: 1.451e-04, 300–400:
1.422e-04, 350–400: 1.417e-04) — no runaway-of-the-runaway, and the mild late
decline has the right sign for its cause: the gap opens to 11.74 over the long
run, and a weaker pull follows a wider pair. The pair is rigid on the t = 200
window and measurably non-rigid on twice that; both facts belong in the text.

### Finding: the flat pulse in the shift movies is emitted by the outer boundary at t = 0 — it is not a reflection

(Absorbed from the retired `Findings.md`, 2026-08-23; originally recorded
2026-08-22 on the runaway N=256 cell's preview movies.)

A flat front is visible crossing the box early in every `shift1` movie. Tracked
on the y-averaged planar detector (planes survive the average, spherical arcs
wash out):

| t | 1.6 | 4.8 | 8.0 | 12.8 | 16.8 |
|---|---|---|---|---|---|
| x of front | 63.12 | 58.12 | 53.62 | 46.88 | 41.12 |

A straight-line fit gives speed **1.389** inward and back-extrapolates to
**x = 64.45 at t = 0** — it starts *at the wall*, at the start. Three reasons it
cannot be a reflection of anything physical: nothing has reached the wall (the
stars sit 27 from the nearest face and the fastest signal in the run moves at
≈ 1.41, so the earliest possible touch is t ≈ 19); it is superluminal,
consistent with the 1+log gauge mode √(2/α) ≈ 1.41, and matter cannot source
it; and it is flat, where anything sourced by the stars would be curved.

What it is: the initial shift is exactly zero everywhere, which is
inconsistent with the outer boundary condition, and the boundary's correction
propagates inward at the gauge speed from t = 0. It appears in every field but
is only *visible* in shift because shift's own dynamic range is tiny:

| field | pulse amplitude | field's own scale | pulse/scale |
|---|---|---|---|
| shift1 | 1.9e-06 | 2.0e-05 | 10% |
| chi | 3.1e-06 | 1.0e-02 | 0.03% |
| local_speed | 3.2e-06 | 1.0e+00 | 0.0003% |

It leaves no trace on the result: the N128 Ham norm is monotone through the
crossing (4.03e-06 at t=9 → 4.49e-06 at t=23 → 7.01e-06 at t=69, no kink), and
the fitted acceleration is unchanged (1.428e-04 at t=63 vs 1.440e-04 at
t=157–199). For the movie figures: label it, don't crop it.

Reproduce (against the run tree's slice cache):

```python
.venv/bin/python - <<'PY'
import numpy as np, glob
d = 'runs/bondi/staging/runaway_pair_d10_L64_N256_lev0/bondi_sg_pair_pm_eqm_w075'
xs = np.arange(256) * 0.25 + 0.125
for f in sorted(glob.glob(f'{d}/frames/_slice_cache/shift1_z/*.npz')):
    z = np.load(f, allow_pickle=True)
    a = z['arr'].astype(float)
    prof = np.abs(a - a.mean()).mean(axis=0)   # y-average: planes survive, arcs wash out
    R = xs > 40
    j = np.argmax(prof[R])
    print(f"{float(z['time']):6.2f}  x={xs[R][j]:6.2f}  amp={prof[R][j]:.2e}")
PY
```

### Finding: the same-sign "flooding" was the pair merging — a misdiagnosis and its correction

(Absorbed from the retired `Findings.md`, recorded 2026-08-23.)

The same-sign controls grow their total box scalar activity sevenfold, onset
t ≈ 40, peak t ≈ 95, at every resolution (7.1× / 7.0× / 6.9×). The first
diagnosis blamed the elliptic solve's outer boundary condition emitting a wave
at t = 0 that reaches the centre at the light-crossing time t ≈ 32 and floods
the box; the cells were declared "quotable only for t < 32" and that
restriction propagated into this document and both results READMEs.

**The diagnosis was wrong, and one column already on disk refuted it:**
`boundary_flux.dat`. Net inward flux in the PP N=256 cell is at round-off
(≤ 1e-08) through the entire early run — nothing arrives from the wall, the
sponge (r = 24 → 32, strength 4) absorbs as designed. The only significant flux
is *outward*, ~1e-04 near t = 100: the story leaving the box.

What actually happens (well tracking now packed as
`campaign/control_pair_pp_d10_L64_N256_lev0/well_tracking.dat`): the two wells
close 8.78 → 8.00 (t = 20) → 5.25 (t = 32) and are merged by t = 33.6. The two
stars pull together and coalesce — predominantly via the same-field overlap
force, ~35× gravity at this separation (free-fall alone would need t ≈ 200) —
and the ×7 growth is merger ejecta; the ringing remnant follows. The centroid stays pinned to ≤ 7.3e-04 through all of
it, which makes the same-sign cells a *stronger* null than the t < 32
restriction allowed — valid over their whole run.

Two mistakes worth naming so they are not repeated:

1. **A timescale coincidence was taken as a mechanism.** t ≈ 32 matched the
   light-crossing time and the search stopped there; the merger completing at
   t ≈ 34 fits the same window.
2. **Resolution-flatness was read backwards.** Identical growth at three grids
   was cited as proof of an initial-data artefact; converged *physics* is flat
   in resolution too. Flatness distinguishes "not a grid effect" from "grid
   effect" — it does not distinguish physics from artefact.

The MM pair's per-star fate, measured the same way by
`control_pair_mm_d10_L64_N128_lev0_frames` on 2026-08-24 (`track_wells.py`,
tracking |χ−1| hills rather than wells): **it merges too, and on the same
clock.**

| t | PP gap (N256) | MM gap (N128) |
|---|---|---|
| 0 | 8.78 | 8.57 |
| 10 | 8.25 | 8.57 |
| 20 | 7.75 | 7.56 |
| 25 | 7.25 | 7.43 |
| 30 | 6.00 | 5.54 |
| 32 | 5.25 | 4.54 |
| **merged** | **t = 33.6** | **t = 32.8** |

The two merger times differ by 2.4%. Read that as "the same timescale", not as
a resolved difference: the only PP cell with a slice cache is N256 and the only
MM cell with one is N128, so the two profiles are sampled at dx = 0.25 and
dx = 0.5 — which is already visible in the t = 0 gaps (8.78 vs 8.57, both from
a nominal separation of 10; the |χ−1| extrema do not sit exactly on the star
centres, and a coarser grid smooths them inward).

This is the single most useful control in the campaign for the mechanism
argument. Newtonian gravity's sign is *opposite* between these two cells —
two positive masses attract, two negative masses mutually repel (that is
Bondi's own result, and the sign flip is why "−− should fly apart" was the
naive expectation). If gravity set the timescale, PP would coalesce and MM
would separate. Instead they collapse together within 2.4% of each other, so
whatever drives them is **blind to the sign of the mass** — the same-field
overlap attraction, ~35× gravity at d = 10. Bondi's −− repulsion is real but
far too weak to see underneath it. See "the two interaction channels". What survives of the
original concern, both t = 0 statements: the same-sign elliptic solve genuinely
floors near 5.4e-04 % (hence the flat 0.002 gate), and the shift pulse above is
real but ~1e-06.

### The star-family floors (from the branch scan, `star_family_massratio.csv`)

Both branches bottom out: the lightest canonical star is M+ = 0.00776 at
ω ≈ 0.92, the lightest phantom |M−| = 0.00791 at ω ≈ 0.94 — about 0.54/0.55 of
the matched mass — and neither branch has a bound star at ω = 1.0. This is why
the mass ladder has no 0.40 rung and why the reversed cell lightens the
canonical instead of fattening the phantom. The bound is a property of the
model worth one paragraph on its own.

## 0. Decisions taken up front (so they are not re-litigated per cell)

### Mesh refinement: `max_level = 0` everywhere. Not negotiable, and here is why

1. **The tagger would never fire anyway.** Refinement triggers at a threshold
   of 0.02; the corrected runs peak at |chi − 1| ≈ 0.005. A `max_level=1` cell
   would evolve on level 0 for its entire life, paying AMR bookkeeping for
   nothing. (The old d=8 data *did* tag at t ≈ 47.5 — because the potential bug
   was collapsing the star to chi ≈ 0.43. That physics no longer exists.)
2. **The convergence ladder must mean exactly one thing** — the cell size. A
   tagging criterion re-evaluated per rung refines a different region on each
   rung, and the ladder stops measuring convergence.
3. **The initial-data fix is proven at matched uniform spacing.** The whole
   artefact was the metric arriving displaced when grids disagree; refined
   levels would put interpolation back into the transfer path, and
   `check_gridinit_alignment.py` only certifies level 0.
4. **Memory.** A refined N=128 cell measured ~35–41 GB against 6 GB unigrid;
   refined N=256 would not fit an 80 GB card. Unigrid never comes close:
   6 GB (N=128) → 20 GB (N=192) → 48 GB (N=256), all on one card.
5. **The known RadialRecipe MPI+CUDA segfault is in the AMR fill path** —
   staying unigrid is also what keeps the multi-GPU option alive.

If a referee asks for an AMR cross-check: one `max_level=1` d=10 cell, same
threshold, is cheap and is predicted to be *identical* because the tagger never
fires. Keep it in the back pocket; do not spend it now.

### MPI and the GPUs

- **Elliptic solves: multi-rank MPI, verified in production.** Every archived
  cell already solved on 32 ranks (~20 min at 256³). The node has 128 CPUs:
  up to three 32-rank solves may overlap; do not run four.
- **Evolutions: one GPU per run is the default.** No cell in this matrix needs
  more than 48 of a card's 80 GB, so multi-GPU is a wall-clock option, not an
  OOM necessity. The four cards are best spent running four *cells* at once.
- **Multi-GPU evolution is unverified on RadialRecipe** (the AMR segfault is
  documented; the unigrid path has simply never been tried). Phase 0 smoke-tests
  it. If it passes, the two N=256 rungs may take two cards each
  (`BONDI_EVO_RANKS=2 BONDI_GPU="0,1"`) to halve the longest runs; if it fails,
  nothing in the plan changes.

All four cards (0–3) stay busy throughout; the gates order the *families*,
not the queue:

| wave | GPU 0 | GPU 1 | GPU 2 | GPU 3 |
|---|---|---|---|---|
| smoke | 2-rank smoke test (cards 0,1) | ← | — | — |
| A | PP d10 N128 | MM d10 N128 | MP mirror N128 | mass-scale (after the CPU family scan) |
| B — phase-1 gate green | PM N192 | PP N192 | wave-zone L128 | long-run t400 (optional) |
| C — N=192 gate green | PM N256 | (+0 if 2-rank passed) | PP N256 | (+2 if 2-rank passed) |

### Time window

Every quantitative claim is fitted on **t ≤ 200**, the window all archived
cells share (fit from t ≥ 5; before that the gauge is settling). The stability
survey stays at its archived t = 120. Anything running longer is illustration,
never a fit window.

### Frames and movies — only where a figure needs them

Drawing frames is the dominant per-plotfile CPU cost, and most cells' product
is numbers, not pictures. The matrix default is therefore `GRTECLYN_FRAMES=0`
(metrics-only consumer; plotfiles are still deleted on the fly — that gate
does not depend on frames). Frames are switched on only for the cells a figure
or movie actually comes from, and those launch with
`GRTECLYN_FRAMES_CACHE_SLICES=1` so every kept series can be re-rendered on one
fixed colour scale afterwards.

**The cells that draw, and the cells that do not** (settled 2026-08-22):

| draws frames | why |
|---|---|
| `runaway_pair_d10_L64_N256_lev0` | the headline movie, at the best resolution the campaign has |
| `control_pair_pp_d10_L64_N256_lev0` | the null beside it, same grid and same colour scale |
| `longrun_pair_d10_t400_L64_N128_lev0` | the sustained-acceleration money plot; it exists for the picture |
| P/M singles | archived, already have frames |

| frameless | why |
|---|---|
| `control_mirror_mp_d10_L64_N128_lev0` | **decided while it was already running.** Its product is the reversed acceleration, a number, and it is the mirror of a cell that already has a movie |
| `control_pair_pp/mm_d10_L64_N128_lev0`, `massscale_…` | numbers only |
| both N=192 rungs, `wavezone_…` | ladder and extraction cells; nothing is read off a picture |

**Frames are ON unless a cell switches them off.** `consumer_frames_enabled()`
returns true for anything that is not literally `GRTECLYN_FRAMES=0/off/no/false`
— an unset variable draws frames. So a frameless cell needs the flag written out
explicitly, while a frames cell works either way. Every frames cell in this
campaign nevertheless states `GRTECLYN_FRAMES=1` beside its slice cache, because
"frames on" should be visible in the launcher rather than inferred from a
default: verified across all ten launchers on 2026-08-22 before wave B started.

**This cannot be revisited after the fact** (README rule 6): the plotfiles a
frame is drawn from are deleted as they are consumed, so a frameless cell can
never be re-rendered — it can only be re-run. A movie of the mirror cell would
therefore cost a fresh ~1.1 h cell, which is cheap if the paper turns out to
want it and wasted if it does not.

The re-render, once a frame-bearing cell lands:

```bash
grteclyn-wrapper/scripts/plot/rerender_frames.py <episode>/frames --movies
```

Plotfiles are still consumed and deleted on the fly — the cache keeps only the
2-D slice behind each frame (~1.4 GB per pair cell; delete
`frames/_slice_cache/` once the movies are final).

### Where the GW shells may sit — outside the stars, inside the sponge

Measured at t=0 on the archived d=10 cell: star cores at ±5 from the box
centre with rms radius 4.34, so ~90% of the matter lives inside r ≈ 11.5 and
the tails reach further. The template's inner extraction shell at r = 8 passes
straight through both stars — it is a near-zone monitor, never a GW
measurement. And the in-code Weyl extraction has been OFF in every cell so
far: all psi4 to date is the consumer's coarse per-plotfile shell sampling.

| box | matter ends | sponge (inner→outer) | shells that mean radiation |
|---|---|---|---|
| L=64 | ~11.5 + tails | 24→32 | r = 16 only — quote nothing else |
| L=128 wave zone | ~11.5 + tails | **48→60, must be moved out** | 16, 24, 32, 40 |

The sponge radii are absolute numbers, not box-scaled: left at their 24→32
defaults, the doubled box would run its dissipation band straight through
three of the four planned shells. Phase 4 sets them explicitly, and switches
the dense in-code extraction on with the new `BONDI_EXTRACTION_RADII` knob.

### Naming and staging

Cell names follow the archive's pattern exactly —
`runaway_pair_d08_L64_N128_lev0` reads as: what the run is, star separation,
box side, cells per side, refinement depth. Every new cell keeps that shape:
`<what>_<dN>_L<box>_N<cells>_lev0`, with extra qualifiers (a retuned omega, a
longer stop time) slotted in before the grid part, e.g.
`massscale_pair_d10_w0790_L64_N128_lev0`, `longrun_pair_d10_t400_L64_N128_lev0`.
A name alone must tell you what the folder holds and on which grid. New runs
land in `runs/bondi/staging/<cell>` and are **moved into
`runs/bondi/staging/archive/` only after their pass gate is green**, so the
archive never contains an unchecked run. After a cell's alignment and t=0
gates pass, delete its `initial_data.gridinit` (0.5–4.4 GB; regenerable from
the cell's own `launch.sh`).

### Running order — launched by hand, one solve at a time

**Cells are started manually, one command each. Do not write a queue script
that launches them by itself.** An orchestrator outlives the session that made
it: on 2026-08-22 a leftover queue from an earlier session was found orphaned
to init, still launching cells on its own, two 32-rank solves at once, while a
replacement queue was being started. Nothing in the run directories shows that
is happening. Before any launch, look for leftovers and kill them first:

```bash
python3 grteclyn-wrapper/scripts/ops/sweep_ranks.py            # what is alive
python3 grteclyn-wrapper/scripts/ops/sweep_ranks.py --kill solves
```

`sweep_ranks.py` walks `/proc` directly because `ps`, `top`, `pgrep` and `free`
are all broken on this node. Kill any orchestrator **first**, then the workers,
then re-run the bare command to verify — killing workers alone just makes an
orchestrator advance to the next cell.

Four cards, but only **one 32-rank elliptic solve may run at any moment**
(README rule 10: a second solve starves every evolution in flight down to a
fifth speed). So each launch waits for the previous cell's solve to hand over
to the GPU — that wait is a decision to take, not something to automate.

Wave B, in the intended order with its pinned card:

| # | cell | card | ~GPU-hours |
|---|---|---|---|
| 1 | `runaway_pair_d10_L64_N192_lev0` | 1 | 5.5 |
| 2 | `wavezone_pair_d10_L128_N256_lev0` | 2 | 9 |
| 3 | `control_pair_pp_d10_L64_N192_lev0` | 3 | 5.5 |
| 4 | `runaway_pair_d10_L64_N256_lev0` | 0 | 17, frames |
| 5 | `control_pair_pp_d10_L64_N256_lev0` | 1 | 17, frames |
| 6 | `longrun_pair_d10_t400_L64_N128_lev0` | 3 | 2.2, frames |

A relaunch needs the cell's half-built episode directory deleted first, or the
wrapper refuses to start with "already exists" and exits within a second. Check
for that before assuming a launch took.

Wave B went up 2026-08-22 03:2x–03:32 on all four cards. The last two wait for
a card and are launched by hand when one frees:

| when this finishes | the card it frees | launch this |
|---|---|---|
| `runaway_pair_d10_L64_N192_lev0` | 1 | `control_pair_pp_d10_L64_N256_lev0` |
| `control_pair_pp_d10_L64_N192_lev0` | 3 | `longrun_pair_d10_t400_L64_N128_lev0` |

```bash
# card 1, after runaway N192 lands  (~4 h solve + ~17.6 h evolve, frames on)
bash runs/bondi/staging/control_pair_pp_d10_L64_N256_lev0/launch.sh

# card 3, after PP N192 lands       (~0.4 h solve + ~2.2 h evolve, frames on)
bash runs/bondi/staging/longrun_pair_d10_t400_L64_N128_lev0/launch.sh
```

Run them from the repository root — the launchers take that path as given, and
a stale working directory is the one way they fail instantly with
"No such file or directory". Preflight both first, and check no solve is
already running before starting a second one.

**Projected finish of the whole stack: Sun 23 Aug, 08:00 at best, ~15:00 with
contention margin.** The critical path is `control_pair_pp_d10_L64_N256_lev0`,
which cannot start until a card frees at ~Sat 10:20 and then needs 21.6 h of
its own; everything else is done by Sun 01:10. The only way to pull the stack
in is to give that cell a card sooner.

---

## 1. What already exists and is reused as-is

Paper-ready, in `runs/bondi/staging/archive/` (details in its README):

| cells | feeds |
|---|---|
| `runaway_pair_d{08,10,12,16}_L64_N128_lev0`, t=200 | separation scan, a·d² = GM table, point-mass limit, headline figures |
| `control_lone_canonical…` / `control_lone_phantom…` (off-centre, x=37), t=200 | single-star null rows of the run matrix |
| `stability/canonical_w{075,080,085,090}…`, t=120 | stability survey section |

Measured anchors these provide: lone-star drift ≤ 1.8e-03 over t=200; the
separation scan gives exponent −2.051 on its original four points (−2.028 once
d = 20 was added, exact −2) and a·d² → 0.014350; d=10 acceleration constant
across four disjoint fit windows.

**What is missing** is everything below: the same-sign pair nulls, the mirror,
the resolution ladder, the mass-scaling lever, and the wave-zone box.

---

## 2. The phases

Baseline launch environment (identical to the archived d=10 cell; every phase
edits only the lines it names). One block = one cell = one command:

```bash
GRTECLYN_FRAMES=0 \
BONDI_GPU=<card> BONDI_STOP_TIME=200 \
BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=10 \
BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004 \
BONDI_S0=0 BONDI_S1=1 BONDI_S0_OMEGA=0.75 BONDI_S1_OMEGA=0.7603 \
BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 \
BONDI_GRTRESNA_RANKS=32 BONDI_GRTRESNA_TIMEOUT=21600 \
BONDI_RUNS_DIR="runs/bondi/staging/<cell>" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

The solve grid always obeys the alignment law
`GRTRESNA_N = NFULL × (GRTRESNA_DOMAIN_L / LFULL)` with no solve refinement —
that is the fix that killed the drift artefact, and it is re-verified per cell.

### The universal pass gate (every cell, before its GPU time is trusted)

1. **Alignment** — `research/bondi_dipole/check_gridinit_alignment.py` on the
   fresh `initial_data.gridinit`: metric-minus-matter centroid offset `0.0000`.
2. **Solve exit** — the tolerance is a request, not a guarantee (rule 8): the
   solver can also leave by stalling or by running out of its 50 iterations,
   and it says "converged" either way. Read the verdict, not the label:

   ```bash
   grteclyn-wrapper/.venv/bin/python research/bondi_dipole/check_solve_exit.py \
       runs/bondi/staging/<cell>
   ```

   Exit 0 or the cell does not count. The d=10 archive landed at 0.00086 /
   0.00080 against 0.002 — only ~2x of headroom, so this is not a formality.
   Remember the gate binds the phantom side, not the canonical.
3. **Matter first** — t=0 row of `small_data/sector_barycenters.dat` matches
   the archived d=10 cell's t=0 row (totals and rms per sector). A dissolved or
   half-painted star makes every later geometry number meaningless.
4. **During the run** — sector total weight flat; scratch dir not growing
   (consumer keeping up).
5. **After the run** — write the measured wall time into the ledger in section
   3. Every episode records it, so this is a lookup, not a stopwatch:

   ```bash
   python3 -c "import json,sys;print(json.load(open(sys.argv[1]))['simulation_elapsed_seconds']/3600)" \
       runs/bondi/staging/<cell>/*/metadata.json
   ```

   The estimates in this plan came from the archive; the ledger is what the
   next campaign should cost its runs from.

### Phase 0 — preflight (minutes, no science GPU time)

- [ ] Binary current: rebuild if any source changed since the archive was
      produced; `git status` clean, note the commit.
- [ ] `BONDI_DRYRUN=1` on the Phase 1 PP command: printed grid, omegas, solve
      N, tolerance all as intended.
- [ ] Disk: ≥ 100 GB free for staging.
- [ ] **MPI evolution smoke test** (the only new machinery question):

```bash
BONDI_EVO_RANKS=2 BONDI_GPU="0,1" BONDI_STOP_TIME=1 \
BONDI_NFULL=64 BONDI_GRTRESNA_N=128 BONDI_NL_TOL=0.1 BONDI_NL_STALL_TOL=0.002 \
BONDI_GRTRESNA_RANKS=8 BONDI_SCRUTINY=0 \
BONDI_RUNS_DIR="runs/bondi/staging/smoke_mpi_evo" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

  Pass: exit 0, no segfault at the first RK4 advance. Records whether
  `BONDI_EVO_RANKS=2` is available for the N=256 rungs. Either outcome is fine.
  Delete the smoke cell afterwards.

  **Result (2026-08-22): passed.** The solve converged in 8 nonlinear
  iterations (Ham 0.082%, Mom 0.086%), the evolution ran all 50 steps split
  evenly across the two cards (identical FAB footprint per rank), and the run
  exited 0 with nothing in the log resembling a fault. So the AMR segfault is
  confined to the refined fill path, exactly as suspected, and the two N=256
  rungs may use `BONDI_EVO_RANKS=2 BONDI_GPU="0,1"` to halve their wall time.

### Phase 1 — the sign matrix (three cells, N=128, ~1.1 h each on the GPU)

### Result — the mirror reverses, and the same-sign pairs do not move

Re-measured 2026-08-22 from `sector_barycenters.dat`, **after both cells
finished and their consumers drained**. "Drift" is the mean of the two sector
barycentres relative to birth; the acceleration is the quadratic coefficient of
a fit over the full shared window `t ≤ 200`. An earlier version of this table
quoted a mid-run snapshot read off the core tracker at `t = 147` and is
superseded — the numbers below are the ones that reproduce from the archived
files.

| | archive (canonical left, phantom right) | mirror (sectors swapped) | ratio |
|---|---|---|---|
| drift at t=200 | **+2.88147** | **−2.88153** | **−1.000022** |
| acceleration (fit, t ≤ 200) | **+1.44877e−04** | **−1.44878e−04** | **−1.000010** |
| worst disagreement over the whole run | — | — | **2.2e−05** relative |

Swapping which star carries positive mass and which carries negative mass
flips the direction of travel and changes nothing else — same speed, same
acceleration, same wobble, agreeing to two parts in a hundred thousand. That is the
point of the control: a drift produced by the grid, the boundary or the solver
would have kept pointing the same way when the physics was swapped. The
direction follows the matter, not the machine.

The same-sign nulls, same date, same files, same window:

| | MM (two phantoms) | PP (two canonicals) | mirror (mixed) |
|---|---|---|---|
| drift of the pair, birth → t=200 | **−0.00028** | **+0.00073** | −2.88153 |
| sector barycentre at t=200 | 31.99972 | 32.00073 | moves 2.9 |
| peak field amplitude, birth → t=200 | 0.02443 → 0.02435 | 0.02457 → 0.02419 | 0.02452 → 0.02443 |
| min χ over the run (1 = flat, 0 = horizon; blind to hills) | **1.00000** | 0.97947 | 0.98889 |
| sector field weight, birth → t=200 | 8.06 → 31.53 (3.9×) | 7.83 → 32.37 (4.1×) | 3.91 → 4.73 (1.2×) |

**Four** orders of magnitude separate the mixed pair from either same-sign pair
— 2.88 against 3e−04 and 7e−04. Neither same-sign pair has a mass dipole to
drive, so neither centroid moves. **Do not read MM's min χ = 1.00000 as "flat
geometry"**: min χ only sees wells, and a phantom star makes a *hill* (χ up to
1.011 on the R = 8 shell). The MM geometry is as curved as PP's, opposite in
sign — and far from static: those two stars merge at t = 32.8, see below. Nothing in the phase is anywhere near collapse — the
lowest χ on the board is 0.979, and a horizon needs it near zero.

The last row is the one to keep an eye on. Both same-sign cells **quadruple**
their sector field weight by t = 200 while the mixed pair grows only 1.2×.
That is the same-sign halo spreading, already documented — the cores stay put
and steady, the envelope does not. It is why the same-sign nulls are quoted on
the barycentre and the core, and why nothing from these two cells is quoted
past t ≈ 60 without saying so.

The phase gate is met on the barycentre and the mirror. The gap half of it is
not measurable with the present tracker; see below.

**Measured 2026-08-22, and it is the documented behaviour**
(`MatterDebugg.md`): the two same-sign cells grow a large halo at late times —
tracked activity up ~6x, rms radius reaching the domain, confined fraction
falling to ~0.10 — while their cores survive intact (peak amplitude steady to
within 3%). This is not a fault. In a same-sign cell both lumps live in *one*
field, so the potential's cross terms give them a direct scalar interaction on
top of gravity; the mixed pair's two sectors share no potential coupling and
meet only through gravity, which is exactly why it is the clean case. The
sponge cannot help: it is zero inside its inner radius, where essentially all
of the halo sits, and it is Kreiss-Oliger dissipation tuned to grid-scale
noise rather than an absorbing layer. **Do not quote a same-sign cell past
t ~ 60.**

The barycentre null is unaffected and is the number this phase exists for: at
t=80 both same-sign pairs sit within 3e-04 of their birth position, six times
below the lone-star noise floor and 1500x below the mixed pair's 0.467 over
the same window.

One measurement this phase does *not* deliver: the **gap** between two
same-sign stars. Both live in one sector, so the core tracker locks onto the
interference peak at their midpoint and reports the separation as `nan` on
every row. The attraction/recession half of the gate needs a two-lump-aware
tracker, which does not exist yet — the barycentre half stands on its own.


The falsifiable core of the paper: only the mixed pair may move.

| cell | change vs baseline | must show |
|---|---|---|
| `control_pair_pp_d10_L64_N128_lev0` | `BONDI_S0=0 BONDI_S1=0`, drop `BONDI_S1_OMEGA`, **`BONDI_GRTRESNA_MAXIMAL_SLICING=1`** (see below) | gap **shrinks** (attraction), pair barycentre still (≤ 5e-03 at t=200) |
| `control_pair_mm_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=1`, `BONDI_S0_OMEGA=0.7603`, drop `BONDI_S1_OMEGA` | gap **grows** (two negative masses recede), barycentre still |
| `control_mirror_mp_d10_L64_N128_lev0` | `BONDI_S0=1 BONDI_S1=0`, `BONDI_S0_OMEGA=0.7603 BONDI_S1_OMEGA=0.75` | acceleration equal and **opposite** to the archived d=10 cell within a few % |

**The all-canonical cell needs one extra knob, and it is the easiest thing in
this campaign to get wrong.** The solver turns maximal slicing on *by itself*
whenever any lump carries negative energy — the CTTK ansatz
`K = sign·sqrt(24πGρ)` is imaginary for `ρ < 0`, so it has no choice. Every
phantom-bearing cell therefore starts from `K = 0`. PP is the only cell in the
whole matrix with no phantom star, so left alone it starts on an
already-collapsing slice and is not comparable to the very cell it is the null
for. Forcing the flag also picks up the rest of the matched path (psi
relaxation 0.6, psi floor 0.1, arithmetic coefficient averaging), so the two
cells then differ only in the sign of the mass.

This was caught on 2026-08-22 by comparing the solve records of the four wave-A
cells: PP was the only one showing `maximal_slicing=0`, and its residuals
oscillated over four orders of magnitude between iterations instead of falling
by a clean factor of ~2.5 like every other cell. The run was stopped before it
reached the GPU and relaunched with the flag. **Every same-sign-canonical cell
in this plan carries it — including the PP ladder rungs in phase 2.**

**The archive was audited for the same fault and is clean** — nothing there
needs re-running. The lone-canonical cell, the only other all-canonical cell
ever produced, set the flag explicitly in its own `launch.sh`. The stability
survey reused pre-solved `K = 0` data. And the decisive check is measured
rather than inferred: the trace of the extrinsic curvature at the innermost
shell reads exactly `0.00000` at birth in every one of the ten archived cells,
which is only true of a flat start. A cell born on the CTTK slice would carry
`K` of order `1e-01` there and show the lapse swinging within the first few
time units.

All three can run at once (three cards, three staggered 32-rank solves — start
them ~5 min apart). Gate for the phase: the two same-sign barycentres sit at
the lone-star noise floor while their gaps move in the *predicted opposite
directions* — attraction where both masses are positive, recession where both
are negative — and the mirror reverses the runaway. Any same-sign cell whose
barycentre drifts like the mixed pair kills the claim; stop and debug, do not
proceed.

### Phase 2 — the resolution ladder (the expensive phase, ~45 GPU-hours)

Headline cell d=10 plus the PP null, each at N=192 and N=256. The tolerance
scales as dx⁴ — this is the one-input fix the old paper could only promise:
every rung previously exited the solve at the *same* residual, so the ladder
had no order to measure. Now the initial violation must fall ~16× per
resolution doubling, and the evolution's own convergence becomes visible.

The tolerances below are the **mixed-cell** values. The same-sign rungs cannot
use them — see "The same-sign rungs keep the N=128 gate" below.

| rung | grid changes | solve changes | GPU cost (t=200) |
|---|---|---|---|
| N=192 | `BONDI_NFULL=192 BONDI_PLOT_INTERVAL=120` | `BONDI_GRTRESNA_N=384 BONDI_NL_TOL=0.000395 BONDI_NL_STALL_TOL=0.0000079` | ~5.5 h |
| N=256 | `BONDI_NFULL=256 BONDI_PLOT_INTERVAL=160` | `BONDI_GRTRESNA_N=512 BONDI_NL_TOL=0.000125 BONDI_NL_STALL_TOL=0.0000025 BONDI_GRTRESNA_TIMEOUT=43200` | ~17 h (or 2 cards if Phase 0 passed) |

Cells: `runaway_pair_d10_L64_N{192,256}_lev0` (baseline env otherwise) and
`control_pair_pp_d10_L64_N{192,256}_lev0` (Phase 1 PP env otherwise — which
means they carry `BONDI_GRTRESNA_MAXIMAL_SLICING=1` too; a PP rung without it
is not a null for anything). The two
N=256 cells are the figure cells: they replace `GRTECLYN_FRAMES=0` with
`GRTECLYN_FRAMES_CACHE_SLICES=1`; every other new cell in the matrix stays
frameless.

Order: run the two N=192 cells first (~5.5 h, overnight is generous), check,
then commit to the two N=256 cells. At N=256 watch the consumer lag — if the
scratch dir grows without bound, stop, double `BONDI_PLOT_INTERVAL`, relaunch.

Gate: (a) the **mixed** rungs exit *at* the tighter tolerances (measured
2026-08-22: N=192 converged at step 14, Ham 3.52e-04 %, still improving 58% a
step — the dx⁴ ladder is sound for these cells); (b) fitted acceleration a
agrees across N=128/192/256 — the spread is the paper's error bar and ~5% or
less is the expectation; (c) the PP barycentre residual *shrinks* rung to rung,
but any residual that does *not* shrink must be checked against the same-sign
floor below before it is called physical; (d) t=0 Ham falls ~16× per rung —
**mixed rungs only**; the PP rungs are flat across resolution by construction.
Remember the standing trap: refining the evolution grid alone makes t=0
*worse* — it is the tolerance scaling that pays for the ladder, which is why
the two are locked together in the table above.

#### The same-sign rungs keep the N=128 gate — measured 2026-08-22

The dx⁴ ladder is right for the mixed cells and **unreachable** for the
same-sign ones. `control_pair_pp_d10_L64_N192_lev0` was launched at
`BONDI_NL_TOL=0.000395`, ran 16 nonlinear steps, and stopped improving:

| step | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|
| Ham [%] | 9.55e-04 | 6.27e-04 | 5.57e-04 | 5.46e-04 | 5.44e-04 |
| improvement on the step | 53% | 34% | 11% | 2.1% | 0.4% |

It floors at 5.4e-04 %, *above* the 3.95e-04 % it was asked for. The cell was
killed at step 16; left alone it would have ground through all 50 steps and
exited by `max_NL_iterations` — the door that does not count — because the
stall detector requires **both** residuals to flatten and Mom was still falling
59% a step (7.02e-05 at step 15, 2.86e-05 at step 16). Note the trap: a stall
tolerance tuned to the Hamiltonian will never fire while the momentum residual
is healthy.

**Cause, confirmed in the solver source.** `BoundaryConditions::fill_constraint_box`
fills the psi correction on the outer boundary with a hard constant 0.0 (the
comment reads "zero for psi"), and `CTTKHybrid::initialise_multigrid_vars` sets
`c_psi_reg` to the constant `regularised_part_psi`. Together these pin the
conformal factor at the solve boundary to exactly 1 — flat space — for the whole
solve. That is the correct condition only for a configuration with zero net ADM
mass. A same-sign pair carries 2M, whose true conformal factor at the boundary is
1 + M/(2R); the solve is structurally unable to represent it. The resulting error
is set by the box, not by the grid, so refinement cannot remove it.

The floor tracks net mass, not resolution — this is the discriminating evidence:

| cell | net ADM mass | lowest Ham reached | still improving? |
|---|---|---|---|
| `runaway_pair_*` (mixed) | 0 | 3.52e-04 % | yes, 58%/step |
| `control_mirror_mp_*` (mixed) | 0 | 8.62e-04 % | yes, 59%/step |
| `massscale_pair_*` (light phantom) | +0.0029 | 8.58e-04 % | yes, 59%/step |
| `control_pair_pp_*` N=192 | +0.0287 | 5.44e-04 % | **no — floored** |

Ten times less leftover mass and no floor appears at all; the same leftover mass
at two resolutions gives the same floor. Note also that `control_pair_pp_*` and
`control_pair_mm_*` at N=128 passed only because their 0.002 gate caught them
(1.27e-03 and 1.33e-03) before the floor was revealed — one step earlier both
were still above it. The N=128 pass was luck, not headroom.

**Fixes considered and rejected.** Moving the boundary out shrinks the error as
1/R while the solve cost grows as R³: merely reaching the dx⁴ gate at N=192 needs
R×1.4 (2.6× cost, zero margin), and a factor-3 margin needs R×2.7 (~20× cost, a
~630 GB intermediate). Replacing the hard zero with a 1/r falloff for the
correction is the physically correct fix and costs nothing at runtime, but it is
new code in the shared elliptic solver and would have to be revalidated before it
could produce paper data.

**Decision: accept the floor and say so.** The same-sign rungs run at a flat
`BONDI_NL_TOL=0.002 BONDI_NL_STALL_TOL=0.00004` at *every* resolution — the same
gate the N=128 rung already used — because their achievable accuracy is set by
the solve boundary and does not improve with cell size. All three PP rungs
therefore begin from initial data of the same quality (~1e-03 % Ham). That is the
honest description of what this solver delivers for a net-mass configuration, and
it is what the paper should state.

What this costs, and what it does not:

| claim | affected? |
|---|---|
| PP/MM as null controls (no self-acceleration) | **no** — the boundary error is spherically symmetric and cannot manufacture a dipole |
| the headline runaway ladder at N=128/192/256 | **no** — mixed pairs have zero net mass, the flat boundary is *correct* for them, and they show no floor |
| PP initial violation falling with resolution | **yes** — it does not fall; do not quote a convergence order for the PP rungs |
| PP evolution convergence rung to rung | partly — the rungs differ in evolution grid only, initial data held at a common quality, which is the cleanest reading available here |

### Result — the ladder converges, and the runaway survives it

All three rungs reached t = 200. Fitting the pair's midpoint to a quadratic over
the last two thirds of each run:

| rung | cells | drift at t=200 | separation | acceleration |
|---|---|---|---|---|
| N=128 | 128³ | +2.8815 | 10.000 → 10.003 | 1.463e-04 |
| N=192 | 192³ | +3.0139 | 10.000 → 9.930 | 1.611e-04 |
| N=256 | 256³ | +3.0016 | 10.000 → 9.915 | 1.596e-04 |

The two finest rungs land on top of each other — 0.4% apart on total drift and
0.9% on acceleration — while N=128 sits about 4% low. That is the shape a
converging quantity has, and it is the opposite of the shape a discretisation
artefact has: refining the grid does not reduce the effect toward zero, it
settles it onto a value. The headline number for the paper is the N=256 rung,
with N=192 as its convergence partner.

Two independent checks agree. Doubling the box (phase 4, L=128) moves the drift
by 4%, so the outer boundary is not driving it. Running four times longer
(t = 400, phase 5) leaves the acceleration flat at 1.44e-04 through the final
quarter, so it is a steady acceleration and not a late-time blow-up.

### The same-sign pairs merge — the activity growth is the merger, not a boundary artefact

**Corrected 2026-08-23.** An earlier version of this section claimed the
same-sign cells were flooded by a wave emitted from the outer boundary at t = 0,
arriving at the centre at t ≈ 32, and were therefore "quotable only for t < 32".
That was wrong, and the measurement that kills it is the boundary flux: the net
inward flux in the PP cell is at round-off (1e-08 and below) through the entire
early run. Nothing arrives from the wall — the sponge (r = 24 → 32, strength 4)
absorbs exactly as designed. The t ≈ 32 ≈ light-crossing coincidence was just
that, a coincidence, and the resolution-flatness of the growth (7.1×/7.0×/6.9×)
was evidence of converged *physics*, not of an initial-data artefact.

What actually happens, established by tracking the two wells frame-by-frame
in the PP N=256 cell (mechanism note, added after the mm frames run: the
attraction driving this infall is predominantly the *field-overlap* force
between two in-phase lumps of the same field, ~35× gravity at this separation
— Newtonian free-fall from d = 10 would take t ≈ 200, not 35 — see "the two
interaction channels" above):

| t | 0 | 13 | 20 | 26 | 32 | ~35 |
|---|---|---|---|---|---|---|
| gap between the two wells | 8.78 | 8.25 | 8.00 | 7.25 | 5.25 | **merged** |

The two stars pull together — predominantly by same-field overlap, ~35×
gravity here (see "the two interaction channels") — and **merge at t ≈ 35**. The
merger ejects scalar field — total box activity rises ×7, peaking at t ≈ 97 —
and the ejecta then drains out through the sponge, which is the *outward*
boundary flux (~1e-04) appearing near t = 100. The well deepens from 1.5e-02 to
2.0e-02 at coalescence and the remnant rings down around 1.8e-02 afterwards:
the "collide and bounce" visible in the chi movie is a merger plus remnant
oscillation, i.e. real physics end to end. The score card's *matter dispersed —
only 15–17% confined* is the merger ejecta; the surviving peak is the remnant
core, not two buried stars.

| | runaway (mixed) | null (same-sign) |
|---|---|---|
| total scalar field in box | 7.9 → 9.2 (+16%) | 7.8 → **54.7** → 27.6 (merger ejecta, then sponge drain) |
| **peak** field strength | 0.0247 → 0.0246 | 0.0247 → 0.0233 (remnant core) |
| net inward boundary flux | round-off | round-off — **nothing arrives** |
| constraint error | 3.7e-06 → 6.6e-06 | 2.1e-05 → 2.3e-05, bounded |

**The MM cells are now explained — by the frames re-run.**
`control_pair_mm_d10_L64_N128_lev0_frames` (physics identical to the archived
MM cell, frames on) shows the phantom pair falling together on the PP pair's
trajectory (lump gap tracked frame by frame with `track_wells.py`: 8.57 → 7.56
at t = 20 → 4.54 at t = 32, merged at **t = 32.8**, against PP's t = 33.6 — a
2.4% difference, and the two are measured on different grids). Two negative masses do NOT repel here,
and that is the *correct* outcome: Bondi's −− repulsion is a statement about
gravity, and the same-field overlap attraction outguns gravity ~35× at d = 10.
The identity of the mm and pp trajectories is itself the proof — gravity is the
only force whose sign differs between those two cells. (The tracker still
reports separation = nan for same-sign pairs — both stars share one sector —
which is why the frames were needed; the flux diagnostic likewise integrates
canonical-sector energy and reads zero in a phantom-only box.)

**Consequence for the paper — better than before.** The same-sign nulls are
valid over their whole run, not just to t ≈ 32: the pair centroid moves by only
+0.00073 (PP) and −0.00028 (MM; −0.00022 in the frames re-run) over the full
t = 200 against the runaway's
+2.88 — four orders of magnitude, sustained through a violent merger and
ringdown. A configuration that cannot manufacture net momentum even while
merging is a stronger null than one that merely sits still. What the same-sign
cells cannot give is a per-star trajectory (single-sector tracking), so the
star-level statement rests on the frames cells. Two residual caveats survive
from the old section, both about t = 0 rather than the evolution: the elliptic
solve's outer boundary condition is genuinely wrong for a same-sign pair (net
2M in the box), which is why the solve gate floors at ~5.4e-04 % and is held
flat at 0.002 — and the tiny t = 0 gauge pulse in the shift documented in
the findings section above is real but ~1e-06 and invisible outside the shift.

### Phase 3 — gravity scales with the source (one cell, ~1.5 h)

The cleanest differential test, and one the grid bug never touched: change the
phantom's mass by ~20%, nothing else. Each star is accelerated by the *other's*
mass, so the canonical star's acceleration must drop by exactly the mass
ratio while the phantom's own acceleration must not move.

First compute the exact retuned mass (no GPU): extend the omega list in
`results/bondi-dipole-runaway/analysis/star_family_scan.py` to cover
0.775–0.80 and pick the phantom omega nearest |M−| ≈ 0.0115 (≈ 0.8 × matched).
The branch anchors: omega 0.750 → −0.015131, 0.775 → −0.013226, 0.7603 →
−0.014350 (the matched point).

**Done 2026-08-22 (CPU, no GPU): the frequency is omega = 0.8040**, giving
M = −0.011472, which is 79.95% of the matched phantom mass. Neighbours for
reference: 0.7900 → −0.012267, 0.8000 → −0.011692, 0.8100 → −0.011162.

Cell `massscale_pair_d10_w<omega>_L64_N128_lev0`: baseline with
`BONDI_S1_OMEGA=<omega>`. Gate: canonical a ratio (this cell / archived d=10) =
mass ratio within a few %; phantom a unchanged within a few %; total momentum
now sums to the expected *non*-zero rate, which is itself a check.

### Result — the pull does follow the partner's mass, but only ratios are quotable

**Measured 2026-08-22.** With unequal masses the pair is no longer rigid: the
gap closes 10.000 → 9.408 over the run, and acceleration goes as 1/d², so a
plain quadratic fit to either star is biased. Fit instead

    x_i(t) = x0 + v0·t + C_i · ∫∫ dt/d(t)²

so that C_canon returns the *phantom's* mass and C_phantom the *canonical's*.
Taking the ratio against the archived matched cell:

| | measured | predicted | agreement |
|---|---|---|---|
| canonical star's pull | **0.809** | 0.7995 (the mass ratio) | 1.2% |
| phantom star's pull | **0.973** | 1.000 (unchanged) | 2.7% |

**Quote ratios, never the constants.** Run the same fit on the matched cell,
where both answers are known to be 0.014350, and it returns 0.0161 and 0.0119 —
+12% and −17%, though their mean is right to 2%. The two stars breathe against
each other over the first ~50 of coordinate time and that oscillation
contaminates the per-star fits *identically in both runs*, so it cancels in the
ratio and does not cancel in the absolute value. The uncorrected midpoint
number, for the record, is 0.928.

### Phase 3b — the mass-ratio ladder (two cells, ~1.5 h each) — branch scanned, cells launched 2026-08-23

Phase 3 leaves the scaling claim resting on two points: ratio 1.000 and ratio
0.7995. Two points fit any monotone law through them, so as it stands the
campaign cannot separate "the pull follows the partner's mass" from "the pull
follows some sublinear function of it". This phase closes that, cheaply — every
cell is `L64_N128_lev0`, about 0.4 h of solve plus 1.1 h of GPU.

Run them in this order; the first is worth more than the other one.

| cell | what changes | prediction | why it matters |
|---|---|---|---|
| `massratio_heavyphantom_d10_L64_N128_lev0` | **canonical** lightened to ω = 0.81, M+ = 0.010721, so \|M−\|/M+ = 1.333 | the gap **OPENS** instead of closing | qualitative sign flip |
| `massratio_w088_r060_d10_L64_N128_lev0` | phantom to ω = 0.88, M− = −0.008573 = 0.597 × matched | canonical pull ratio 0.597 | lever arm 2× Phase 3, clears the ~3% fit systematics |
| ~~`massratio_*_r040_*`~~ | phantom ≈ 0.40 × matched | — | **cannot be built — see the branch bound below** |

**Why the reversed cell lightens the canonical rather than fattening the
phantom.** Every uneven cell so far has the gap closing, so a sceptic can argue
the closing is the artefact and the drift is a consequence of it. Reversing the
mass ordering predicts the opposite — the gap opens — and no artefact story
predicts a sign flip that tracks which star is heavier. Getting there by making
the phantom heavier would need \|M−\| ≈ 0.0187, off the end of the scanned
branch and possibly past its maximum; making the canonical lighter reaches the
same ordering inside frequencies already surveyed (the stability run covers
canonical omega 0.75–0.90). Cheaper and safer.

**Prerequisite, CPU only — done 2026-08-23.** The scan lives in
`results/bondi-dipole-runaway/analysis/star_family_massratio_scan.py` and writes
`stars/star_family_massratio.csv`. It is a separate file from
`star_family_scan.py` on purpose: that table is what the matched pairing was
chosen from, and every published number traces to it, so it is not rewritten to
answer a later question.

**The phantom branch has a floor on how light its stars get, and that bound
belongs in the paper.** Walking ω upward, |M−| falls, flattens, turns around,
and then the branch ends:

| ω | 0.81 | 0.84 | 0.88 | 0.90 | 0.92 | **0.94** | 0.96 | 0.98 | 1.00 |
|---|---|---|---|---|---|---|---|---|---|
| \|M−\| | 0.01116 | 0.00983 | 0.00857 | 0.00816 | 0.00791 | **0.00791** | 0.00835 | 0.01005 | none |
| ÷ matched | 0.778 | 0.685 | 0.597 | 0.569 | 0.551 | **0.551** | 0.582 | 0.700 | — |

The lightest phantom star that exists on this branch is |M−| ≈ 0.00791, about
**0.55 of the matched mass**, near ω = 0.94; past that the star gets heavier
again and by ω = 1.0 there is no bound star at all. A 0.40 rung would need
|M−| ≈ 0.00574, which is below the floor and does not exist. That is why this
phase has two cells and not three, and the bound is a result in its own right —
it is the reason the reversed cell is built by lightening the canonical rather
than by fattening the phantom.

**Both branches have a floor**, at almost the same place. The canonical bottoms
out at M+ ≈ 0.00776 near ω = 0.92 (0.541 of matched) and the phantom at
|M−| ≈ 0.00791 near ω = 0.94 (0.551); both then rise, and neither has a bound
star at ω = 1.0. For reference the canonical branch over the window used here:
0.81 → 0.010721, 0.84 → 0.009495, 0.88 → 0.008336, 0.90 → 0.007963 — still well
above its floor at 0.81, which is what makes ω = 0.81 a comfortable choice for
the reversed cell.

**If a third point is wanted later**, the lever must come from the canonical
side, and it must go *down* in ω, not up: hold the phantom at its lightest
(ω ≈ 0.94, |M−| = 0.00791) and make the canonical heavier than 0.014350 by
dropping its ω below 0.75, away from its own floor. That reaches ratios under
0.55, at the cost of changing both stars rather than one. Raising the canonical
ω instead would run into its floor at 0.92 and buy nothing.

Gate for every cell in this phase: the per-star pull ratio from the
separation-corrected fit above lands within a few % of the mass ratio, and the
three ratios together are linear in the mass ratio rather than merely monotone.
For the reversed cell the gate is the **sign** of d(separation)/dt, not its
size.

### Result — the gap's sign follows the mass ordering, and that closes the artefact argument

Both cells reached t = 200. Placing them beside the two earlier rungs gives four
points spanning mass ratios 0.597 to 1.333:

| cell | \|M−\|/M+ | separation, 0 → 200 | change |
|---|---|---|---|
| `massratio_heavyphantom` | **1.333** | 10.000 → **10.603** | **+0.603 — OPENS** |
| `runaway_pair_d10` (matched) | 1.000 | 10.000 → 10.003 | +0.003 — flat |
| `massscale_pair_w0804` | 0.799 | 10.000 → 9.408 | −0.592 — closes |
| `massratio_w088_r060` | 0.597 | 10.000 → 8.618 | −1.382 — closes more |

**The gate for the reversed cell was the sign of d(separation)/dt, and it
passes.** The gap opens, sits flat, or closes strictly according to which star is
heavier, and it crosses zero at the equal-mass point. This is the answer to the
sceptic's reading of phase 3 — that the gap closing is the artefact and the drift
follows from it. No artefact story predicts a sign reversal that tracks the mass
ordering, and none predicts the crossing landing exactly where the masses match.

**Read the transient before reading the trend.** Every cell dips around t ≈ 50
and the dip is not the signal: the matched cell falls to −0.255 at t = 50 and
recovers to +0.003 by t = 200. At t = 54 the reversed cell sits at −0.207 — still
negative — and only crosses over later. Any verdict taken before t ≈ 150 will be
wrong, in this cell and in the others.

### The quantitative scaling holds in the middle and degrades at the extremes

Fitting each star separately to `x_i(t) = x0 + v0·t + C_i·∫∫dt/d(t)²` and
quoting `C_i` as a ratio to the matched cell:

| cell | canonical pull | predicted | err | phantom pull | predicted | err |
|---|---|---|---|---|---|---|
| `massratio_w088_r060` | 0.597 | 0.597 | **−0.1%** | 0.725 | 1.000 | **−27.5%** |
| `massscale_pair_w0804` | 0.809 | 0.799 | +1.2% | 0.973 | 1.000 | −2.7% |
| `massratio_heavyphantom` | 1.052 | 1.333 | **−21.1%** | 0.739 | 0.747 | **−1.0%** |

Four of the six land within 3%, which is the gate. The two that miss are not
scattered — they are both the star whose partner's mass was changed *most*, both
miss in the same direction (under-reporting the pull), and both belong to the
cells whose separation swings hardest (10 → 8.62 and 10 → 10.60). The fit assumes
the pair is near-rigid; once the separation moves by more than about 10% that
assumption is spent, and the double integral no longer absorbs the geometry.

**So quote the extremes for their sign and the middle rungs for their
magnitude.** The linear-in-mass-ratio claim is supported over
0.8 ≲ \|M−\|/M+ ≲ 1.0 and should not be extended to the ends of this table
without a fit that models the changing separation properly. That is a limitation
of the analysis, not evidence against the scaling — the same fit returns 0.0161
and 0.0119 on the matched cell for two constants both known to be 0.014350.

### Phase 4 — the wave zone (one cell, ~9 h)

Doubled box, same cell size, so the light-crossing distance grows without
touching resolution — the wave-zone question is distance, not dx.

Cell `wavezone_pair_d10_L128_N256_lev0`:

```
BONDI_NFULL=256 BONDI_LFULL=128 BONDI_RADII="16 24 32 40" \
BONDI_EXTRACTION_RADII="16 24 32 40" BONDI_PSI4_HIGHER_L=1 \
BONDI_SPONGE_INNER=48 BONDI_SPONGE_OUTER=60 \
BONDI_GRTRESNA_DOMAIN_L=256 BONDI_GRTRESNA_N=512 BONDI_GRTRESNA_TIMEOUT=43200
```

(the solve box doubles with the evolution box, keeping both the spacing match
and the far outer boundary). 48 GB, one card. Gate: r·psi4 (l=2) amplitude
constant across the four shells; features arrive at retarded time between
shells; the l=2 signal at R=16 consistent with the L=64 runs' extraction.
The dipole-radiation subsection stands or falls here — for a pair with zero
total momentum the mass dipole is static, so the prediction is *no* l=1-type
growth; whatever is measured is reported.

### Result — no wave zone is reached, and that is the answer

The cell ran to t = 200 in the doubled box. Its drift, +2.7606, is within 4% of
the L = 64 N = 128 cell, so doubling the light-crossing distance does not change
the runaway — the boundary is not driving the effect.

The r·ψ₄ gate fails, and the way it fails is the result. Amplitudes of the l = 2
mode at the four shells, at t = 199:

| R | 16 | 24 | 32 | 40 |
|---|---|---|---|---|
| ψ₄ | 3.23e-05 | 6.52e-06 | 2.17e-06 | 7.99e-07 |
| r·ψ₄ | 5.16e-04 | 1.57e-04 | 6.93e-05 | 3.20e-05 |

r·ψ₄ is meant to be flat in the wave zone. It falls by a factor of 16 across the
four shells. Fitting a power law gives ψ₄ ∝ r^−4.0 (it steepens from r^−3.2 at
t = 50 as the near field organises). A radiative tail would give r^−1.

So every shell available in this box is in the near zone, and no gravitational
radiation is measurable out to R = 40. This is consistent with what phase 4 was
set up to predict: for a pair with zero total momentum the mass dipole is
static, so there is no l = 1 growth to find, and the l = 2 content that exists
is quasi-static near field rather than radiation.

**How to report it.** As a null with a number, not as a missing measurement. The
runaway is not powered by, and is not accompanied by, gravitational radiation
detectable at these radii; the energy bookkeeping is between the two stars and
the field that binds them. Reaching a true wave zone would need extraction radii
well beyond this box, which is a different cell and not one this paper needs.

### Phase 5 — optional garnish (only if the paper wants the figures)

- `longrun_pair_d10_t400_L64_N128_lev0` — baseline with `BONDI_STOP_TIME=400`,
  ~2.2 h: the sustained-acceleration money plot (velocity still growing
  linearly at t=400). It exists for the picture, so it launches with frames
  and the slice cache on. Illustration only; fits stay on t ≤ 200.
- `control_pair_mm_d10_L64_N{192,256}_lev0` — completes the null-ladder figure
  with MM alongside PP. Same cost as the PP rungs; skip unless the convergence
  figure looks thin without it.

### Phase 6 — analysis, movies, packing (no GPU)

1. Update the hardcoded cell lists in `results/bondi-dipole-runaway/analysis/`
   (`separation_scaling.py`, `convergence_check.py`, `momentum_balance.py`,
   `newtonian_reference.py`, `make_tables.py`, …) from the old `convA_*` names
   to the archive cell names, then regenerate every table.
2. Re-render every movie on its fixed colour scale (`rerender_frames.py
   <episode>/frames --movies`), then delete each `frames/_slice_cache/`.
3. `research/bondi_dipole/pack_runaway.sh` to refresh the git-tracked extract
   pack in `results/bondi-dipole-runaway/campaign/`; update both READMEs with
   the new cells.
4. Rewrite the results sections of `bondi_dipole.tex` from the new tables —
   the old text's gap-closing/contact-contamination narrative does not survive
   the fix and must not be patched around.

---

### Result — mesh refinement changes nothing, measured rather than argued

`amrcheck_pair_d10_L64_N128_lev1` ran the headline cell with `max_level = 1`
against the campaign's uniform-grid rule, to answer the referee question once.

**Level 1 was never created.** No step ran on it; the 213 log lines reading
`Now regridding at level lbase = 0` are the periodic check, not a refinement. The
tagger threshold is |χ − 1| = 0.02 and these runs peak near 0.005, so there was
nothing to tag — exactly the prediction in section 0.

The cell is therefore a bit-level re-run of the uniform cell, and reports as one:

| | drift at t=200 | separation | acceleration |
|---|---|---|---|
| uniform (`lev0`) | +2.8815 | 10.0030 | 1.4634e-04 |
| AMR enabled (`lev1`) | +2.8815 | 10.0030 | 1.4634e-04 |
| difference | +0.001% | — | +0.002% |

The `max_level = 0` decision costs the campaign nothing and can be defended with
a number instead of an argument.

## 3. Budget

| phase | cells | GPU-hours | wall (4 cards) |
|---|---|---|---|
| 1 sign matrix | 3 | ~3.5 | ~1.5 h |
| 2 ladder | 4 | ~45 | ~1.5 days |
| 3 mass scale | 1 | ~1.5 | with phase 2 |
| 4 wave zone | 1 | ~9 | with phase 2 |
| 5 optional | 0–3 | 0–25 | — |
| **total (required)** | **9** | **~59** | **~2–3 days** |

### Measured wall time — the ledger, and the number to quote in the paper

**The campaign cost about 73 GPU-hours on one H100 per cell.** That is the
figure to quote: it is the card-time the 27 cells actually consume, computed
from each cell's own end time and the fastest clean rate measured for its grid.

Two numbers exist and they are not the same, so it is worth being explicit about
which is which:

| | hours | what it means |
|---|---|---|
| elapsed wall time, summed over cells | 81.6 | what each cell's own clock recorded |
| **exclusive-equivalent card-time** | **73.0** | **the compute actually consumed** |

They differ because several cells shared a card with another cell, and a shared
cell's clock keeps running while the other one holds the GPU. Summing those
elapsed times counts the same card-hour twice. The 73 h figure removes that
double count by pricing every cell at the clean single-occupancy rate for its
grid. Quote 73; mention 81.6 only if wall-clock duration is what is being
discussed.

| cell | N | L | t | elapsed h | card-h | h / 1000 t |
|---|---|---|---|---|---|---|
| `runaway_pair_d08_L64_N128_lev0` | 128 | 64 | 200 | 1.10 | 1.09 | 5.44 |
| `runaway_pair_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.09 | 1.09 | 5.44 |
| `runaway_pair_d12_L64_N128_lev0` | 128 | 64 | 200 | 1.09 | 1.09 | 5.44 |
| `runaway_pair_d16_L64_N128_lev0` | 128 | 64 | 200 | 1.10 | 1.09 | 5.44 |
| `runaway_pair_d20_L64_N128_lev0` | 128 | 64 | 200 | 1.69 † | 1.09 | 5.44 |
| `runaway_pair_d10_L64_N192_lev0` | 192 | 64 | 200 | 4.75 | 4.75 | 23.75 |
| `runaway_pair_d10_L64_N256_lev0` | 256 | 64 | 200 | 14.60 | 14.60 | 72.99 |
| `control_lone_canonical_L64_N128_lev0` | 128 | 64 | 200 | 1.09 | 1.09 | 5.44 |
| `control_lone_phantom_L64_N128_lev0` | 128 | 64 | 200 | 1.10 | 1.09 | 5.44 |
| `control_pair_pp_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.40 † | 1.09 | 5.44 |
| `control_pair_pp_d10_L64_N192_lev0` | 192 | 64 | 200 | 5.82 † | 4.75 | 23.75 |
| `control_pair_pp_d10_L64_N256_lev0` | 256 | 64 | 200 | 14.63 | 14.60 | 72.99 |
| `control_pair_mm_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.41 † | 1.09 | 5.44 |
| `control_pair_mm_d10_L64_N128_lev0_frames` | 128 | 64 | 200 | 1.70 † | 1.09 | 5.44 |
| `control_pair_mm_d10_L64_N192_lev0` | 192 | 64 | 200 | 6.03 † | 4.75 | 23.75 |
| `control_mirror_mp_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.41 † | 1.09 | 5.44 |
| `massscale_pair_d10_w0804_L64_N128_lev0` | 128 | 64 | 200 | 1.40 † | 1.09 | 5.44 |
| `massratio_w088_r060_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.09 | 1.09 | 5.44 |
| `massratio_heavyphantom_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.22 † | 1.09 | 5.44 |
| `amrcheck_pair_d10_L64_N128_lev1` | 128 | 64 | 200 | 2.99 † | 1.09 | 5.44 |
| `deepsolve_pair_d10_L64_N128_lev0` | 128 | 64 | 200 | 1.11 | 1.09 | 5.44 |
| `wavezone_pair_d10_L128_N256_lev0` | 256 | 128 | 200 | 7.37 | 7.35 | 36.76 |
| `longrun_pair_d10_t400_L64_N128_lev0` | 128 | 64 | 400 | 3.82 † | 2.18 | 5.44 |
| `stability_canonical_w{075,080,085,090}` ×4 | 128 | 64 | 120 | 2.63 | 2.61 | 5.44 |
| **total** | | | | **81.6** | **73.0** | |

† shared a card with another cell for part of its run, so its elapsed time
overstates the compute it consumed.

**Cost per unit of evolution, measured.** At N = 128, L = 64 the cost is
**5.44 GPU-hours per 1000 units of t**, and it does not vary by more than 2%
across sixteen cells — separation, sector signs, star frequency, mesh
refinement and frame rendering all cost the same. Run length is the only thing
that sets the bill at this grid.

**The resolution scaling is shallower than N⁴.** The naive expectation is three
spatial dimensions plus the shorter timestep, so `(192/128)⁴ = 5.1×` and
`(256/128)⁴ = 16×`. Measured against the N = 128 rung: **4.4×** and **13.4×**.
The GPU is better occupied at the larger grids, which is where the difference
comes from.

The wave-zone cell is a useful cross-check on the whole model: it has the same
256³ cell count as the finest rung but twice the cell size, so its timestep is
twice as long and it should cost exactly half per unit of t. Measured, 36.76
against 72.99 — half, to within a part in a thousand.

**Refinement was free here.** `amrcheck` carries the same 5.44 h/1000 t as the
uniform cells because its tagger never fired and level 1 was never created; the
2.99 h on its own clock is card sharing, not AMR overhead.

**The elliptic solves are CPU, not GPU, and are not in the 73 h.** Each cell is
preceded by a GRTresna solve on 32 MPI ranks: roughly 20 minutes at 256³ up to
about 4 hours at 512³. Those overlap the GPU work of *other* cells, which is
why cells were launched while a card was still busy — the evolution only claims
the GPU once its solve has landed. Counted as node time rather than card time,
the campaign's solves add on the order of 30 CPU-core-days.

**For the paper's methods section**, the defensible sentence is: *the campaign
comprises 27 numerical-relativity evolutions totalling ≈73 GPU-hours on NVIDIA
H100 hardware, one card per evolution, preceded by constraint solves on 32 CPU
ranks each.*

---

## 4. Beyond t=200 — the 0.3c chase (moved out of this document)

The outlook analysis (how far the dipole can be pushed toward the speed of
light) and the recentring-box ("treadmill") realisation plan that used to be
sections 4 and 5 here have moved to
`research/bondi_dipole/docs/CHASE_TO_03C.md`, and the treadmill is no longer a
plan: it is implemented (`Source/Grids/GridTreadmill.hpp`,
`Examples/RadialRecipe/RecentringBox.hpp` + `SectorCoreTracker.hpp` +
`TreadmillLog.hpp`, tests in `Tests/GridTreadmillTest/`) and its four
implementation-test legs all passed — see
`runs/bondi/staging/treadmill_pair_d10_L64_N128_lev0/README.md`. This document
is the record of the t = 200 campaign; the chase is the follow-up paper's.
