# Matter debugging — getting to a star that does not collapse

> **Status — 2026-08-20: the `ω = 0.55` matter model is rejected.**  The
> follow-up campaign completed all six cells, and its own data shows the
> positive-mass star sits on the unstable branch of its family and collapses on
> the timescale the measurement needs.  The next step is a stable star, not
> another cell on this one.
>
> **Status — 2026-08-21: that rejection is now in doubt, and the reason is
> numerical.**  The replacement stars at `ω = 0.80` and `0.75` failed the same
> way, and the failure has an under-resolution signature rather than a
> stability one: the cores never moved, while the initial data's own constraint
> violation got *worse* when the grid was refined.  See
> [The initial data is the error floor](#the-initial-data-is-the-error-floor--2026-08-21)
> — read it before trusting any resolution result in this document, including
> queue item 1.
>
> **Status — 2026-08-21 (later the same day): the cause is found, and it is a
> bug in how the initial data is built, not the matter model.**  The canonical
> star was being handed a slice that is *already collapsing* at birth, while
> its phantom counterpart was handed one at rest — an asymmetry with no
> physical content, arising from a single line that decides which
> constraint-solving method to use.  See
> [Critical bug: the canonical star is born mid-collapse](#critical-bug--the-canonical-star-is-born-mid-collapse--2026-08-21).
> **Both rejections above are now suspect**, because every canonical cell ever
> run — `ω = 0.55`, `0.75`, `0.80`, at every resolution — was built this way.
>
> This is the single working document for that loop: **what the campaign found**
> → **why the matter is wrong** → **what a stable replacement looks like** →
> **how to run and check it**.  The original GPU queue is kept at the end,
> flagged as historical.  Results are in
> [What the campaign returned](#what-the-campaign-returned--2026-08-20) directly
> beneath this section, together with the one finding that changes the article's
> interpretation.  Everything from *Before launching anything* onward is kept
> exactly as it was written *before* the campaign, so the reasoning that
> motivated each cell can be read against what actually came back.

*Written before the campaign:*  The published campaign is complete: the article
stands on the runs already in `results/bondi-dipole-runaway/campaign/`, and
nothing below is needed to reproduce a number in it.  This is the **queue for
the next campaign** — the open items that cannot be answered by re-analysing
what is on disk, why each one needs the device, and what it unblocks.

*That last clause did not survive.*  Every published number still reproduces,
but the campaign showed the canonical star to be on the unstable branch of its
own family, which changes what those numbers mean.  See below.

Read [`../../../results/bondi-dipole-runaway/LAUNCH.md`](../../../results/bondi-dipole-runaway/LAUNCH.md)
first for the launch mechanics: the `/usr/bin/env` shadowing trap, how to
verify a detached launch actually started, and `stop_campaign.sh`.  Costs quoted
here are the measured ones from Appendix B of the article.

---

## Critical bug: the canonical star is born mid-collapse — 2026-08-21

**Every canonical star this project has ever evolved was given a large inward
velocity of the geometry at `t = 0` that nobody asked for.  Every phantom star
was not.  The canonical stars are the only ones that collapse.**

### The mechanism

Constructing initial data means choosing, at the first instant, not only where
the matter is but also whether space is expanding, contracting or static.
GRTresna's default CTTK ansatz ties that choice to the local energy density
(`Source/Methods/CTTKHybrid.impl.hpp`):

```cpp
Real K_0_squared = 24.0 * M_PI * G_Newton * emtensor.rho;
multigrid_vars_box(iv, c_K_0) = m_method_params.sign_of_K * sqrt(K_0_squared);
```

`K` is the trace of the extrinsic curvature — the rate at which volume
elements shrink.  The ansatz is a legitimate solvability trick: it dumps the
matter energy into `K` so the Hamiltonian constraint becomes algebraic instead
of elliptic.  But it takes **the square root of the energy density**, so for
exotic matter (`ρ < 0`) it is imaginary and cannot be used at all.  A second
path exists for that case — `maximal_slicing`, which fixes `K = 0` and pushes
the whole matter energy into the elliptic `ψ`-solve (the standard
York/Lichnerowicz CMC construction, valid for either sign of `ρ`).

**The bug is the rule that chooses between them**
(`grteclyn-wrapper/src/grteclyn_wrapper/grtresna/fit/motif.py`):

```python
maximal_slicing=has_exotic or motif.exotic_needed,
```

The careful `K = 0` path is switched on **if and only if exotic matter is
present**.  It is a solvability guard, written to stop `sqrt()` seeing a
negative number — but its side effect is that the two sectors of a Bondi pair
are built by *different methods*.  Nothing in the code, the launchers or the
parameter files ever says so.

### The evidence

Measured `max|K|` at `t = 0.01`, same rung, same box, same matched solve grid,
`data/collapse_diagnostics.dat` column 4:

| cell | slicing path | `max\|K\|` at birth |
|---|---|---|
| canonical `ω = 0.80` | CTTK ansatz | `1.005e-01` |
| canonical `ω = 0.75` | CTTK ansatz | `9.918e-02` |
| phantom `ω = 0.80` | maximal, `K = 0` | `1.034e-05` |
| phantom `ω = 0.75` | maximal, `K = 0` | `3.867e-05` |

**Four orders of magnitude, and it is exactly the ansatz formula.**  Evaluating
`K = sqrt(24 π G ρ_core)` from the catalogue's own central amplitude, with
`ρ = ½ω²A² + V(A)` on this rung (`m = 1`, `λ = 10240`, `μ = 21845333`):

| `ω` | `φ_c` | `ρ_core` | predicted `K` | measured `max\|K\|` | agreement |
|---|---|---|---|---|---|
| 0.80 | `0.019126` | `1.356e-04` | `0.1011` | `0.10055` | 0.5 % |
| 0.75 | `0.019720` | `1.308e-04` | `0.0993` | `0.09918` | 0.1 % |

The number in the data is the formula, not physics.  `K ≈ 0.10` against a field
frequency of `ω = 0.80` is a coherent contraction at ~12 % of the star's own
internal rate, applied uniformly across the core at the instant of creation.

### Why the star cannot survive it

The star has no way to push back.  Its central amplitude is pinned to the
potential's preferred value `A_* = sqrt(3λ/4μ) = 0.01875` at *every* frequency
in the family, so it cannot stiffen in response to being squeezed.  What is
measured instead, in every canonical cell:

1. lapse crashes to ~0.93 within `t ≈ 5` — the initial squeeze, about two field
   periods long
2. rebound and ring-down to `t ≈ 20`, shedding the outer envelope
3. the shed field is sub-luminal against the scalar mass and cannot escape, so
   it fills the box: `total_activity` **doubles** (`3.43 → 7.34` by `t ≈ 42`)
   while the confined fraction falls `0.714 → 0.316`
4. the core does the reverse — thins by 37 % while the central height climbs
   15 % — narrower, taller, denser, compounding with a ~20–25 `t` doubling time
5. horizon.  `min_lapse → 0.21` by `t = 115` at `ω = 0.80`; `→ 0.12` by
   `t = 110` at `ω = 0.75`

Heavier stars die sooner (`ω = 0.55` at `t ≈ 60`, `0.75` at `t ≈ 103`, `0.80`
at `t ≈ 108`).  **That ordering is why this looked like physics** — it is the
signature of a genuine instability.  It is equally the signature of this bug:
heavier means larger `ρ`, which means a larger `K` kick.  The two hypotheses
predict the same ordering, which is precisely why three campaigns failed to
separate them.

### What this rules out — everything else was eliminated first

| suspect | verdict | evidence |
|---|---|---|
| under-resolved / interpolated initial data | **ruled out** | rebuilt on the matched solve grid, 100–195× cleaner interior Hamiltonian error at birth; the decline changed by < 0.5 % |
| evolution resolution | **ruled out** | `N = 192` tracks `N = 128` to ~0.1 % |
| outer boundary, reflections | **ruled out** | net boundary flux ≤ `1.1e-06` all run; the phantom shares the box and is flat |
| the matter pump / trajectory well | **ruled out** | `rl_pump_stop_time = 0` gates the pump off for all `t ≥ 0` (`num_sites = 0`, gains zeroed, `RadialRecipeMatterDispatch.hpp`); the evolution is unforced Einstein–Klein–Gordon |
| measurement artefact | **ruled out** | `chi` falls with the lapse — space is really compressing — and horizons actually form |
| past the family's turning point | **not applicable** | the family's mass never turns; it falls monotonically `0.0794 → 0.0080` across `ω = 0.53 → 0.90` |

The A/B that settled the initial-data question, `min_lapse`, clean-born vs
original, same grid and star:

| `t` | canon 0.80 clean / dirty | canon 0.75 clean / dirty | phantom clean / dirty |
|---|---|---|---|
| 10 | 0.9686 / 0.9682 | 0.9647 / 0.9643 | 0.9998 / 0.9997 |
| 20 | 0.9764 / 0.9783 | 0.9730 / 0.9755 | 0.9998 / 0.9997 |
| 30 | 0.9623 / 0.9653 | 0.9563 / 0.9600 | 0.9998 / 0.9997 |
| 40 | 0.9467 / 0.9505 | 0.9366 / 0.9415 | 0.9998 / 0.9998 |
| 50 | 0.9276 / 0.9321 | 0.9119 / 0.9180 | 0.9998 / 0.9998 |

A hundredfold cleaner birth changed nothing — because the birth `K` was never
an accuracy problem.  It is the *intended* output of the method that was
selected.

### The fix — slicing is now a choice, not a side effect of being exotic

`maximal_slicing` is exposed as a knob so canonical matter can be built the
phantom's way.  Nothing about the star changes — same profile, mass and
frequency — only the state of the slice at `t = 0`.

```bash
BONDI_GRTRESNA_MAXIMAL_SLICING=1   # -> replay_eval.py --grtresna-maximal-slicing
```

The flag forces `maximal_slicing = True` and runs `apply_exotic_safe_solver`
on the base config before the matter fit, which only ever turns the flag *on*,
so the choice survives.  It also carries the same under-relaxation the phantom
uses (`psi_relaxation = 0.6`, `psi_floor = 0.1`) so the two constructions match
in every respect except the sign of `ρ`.

**The `K = 0` path converges better for canonical matter than the ansatz does**
— there was never a technical reason to avoid it, only the `sqrt()` guard that
made it mandatory for exotic matter.  Measured on the `ω = 0.80` star at the
matched grid:

| NL iteration | 0 | 1 | 3 | 6 | 9 | 12 | 13 |
|---|---|---|---|---|---|---|---|
| Ham error [%] | `100` | `24.7` | `3.25` | `0.200` | `0.0127` | `1.03e-3` | `7.21e-4` |

Monotonic, no stall, exits under the `0.001 %` gate.

**Confirmed at birth.**  The same `ω = 0.80` canonical star, matched grid, with
and without the flag — and the phantom it is meant to match:

| quantity at `t = 0.01` | canonical, CTTK | canonical, `K = 0` | phantom, `K = 0` |
|---|---|---|---|
| `max\|K\|` | `1.005e-01` | **`1.023e-05`** | `1.034e-05` |
| `L2_Ham_amr` | `2.622e-06` | `1.538e-06` | `1.558e-06` |
| `Linf_Ham_amr` | `1.992e-04` | `1.095e-04` | `1.058e-04` |
| `min_lapse` | `0.99009` | `0.99209` | `0.99981` |

The birth kick is gone — down by a factor of 9800, landing on the phantom's
value — and the constraint violation improves at the same time, so nothing was
traded away to get it.  The canonical and phantom sectors are now built
identically in every respect except the sign of `ρ`.

**And the initial squeeze — step 1 of the death sequence above — does not
happen.**  `min_lapse` over the first two time units, same star, same grid:

| `t` | canonical, CTTK | canonical, `K = 0` |
|---|---|---|
| 0.5 | `0.89889` | `0.99225` |
| 1.0 | `0.81946` | `0.99282` |
| 1.5 | `0.77172` | `0.99408` |
| 2.0 | `0.77992` | `0.99551` |

The CTTK star loses 22 % of its lapse in a violent squeeze inside `t < 2`,
before its field has completed two periods; the `K = 0` star stays flat and is
drifting gently *upward*.  Whether that holds to `t = 120` is what the campaign
below measures.

### The test

`grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_slicing_test.sh` — four
canonical cells, all built with `K = 0`, `t = 120` on the published grid.
`ω = 0.80` and `0.75` are one-flag A/B partners of the completed CTTK cells;
`0.85` and `0.90` walk up the frequency ladder so a partial cure is still
interpretable.

- **flat** → the collapse was an artefact of the initial-data method, and this
  rung is usable for the pair campaign after all
- **still sinks** → the birth kick is not the cause, the rung is genuinely
  unstable to compression, and the mass dial or a new rung is the next move

**Until this returns, treat every canonical-sector conclusion in this document
as provisional** — including the `ω = 0.55` rejection, the turning-point
argument, and the stable-replacement recommendation.

---

## The initial data is the error floor — 2026-08-21

**Refining the evolution grid makes the starting constraint violation worse,
not better.**  Measured on the canonical `ω = 0.80` star, same box, same
physics, at matched time `t = 0.02`, in `L2_Ham` — a *volume-weighted RMS*
(`RadialRecipeLevel.cpp`: `sqrt(sum_ham2 / sum_vol)`), so it is
resolution-independent by construction and the two numbers are directly
comparable:

| quantity | `N = 128` | `N = 192` | ratio |
|---|---|---|---|
| `L2_Ham` (level 0, whole domain) | `5.628e-4` | `1.020e-3` | `1.81` |
| `L2_Ham_amr` (AMR composite) | `2.515e-4` | `4.710e-4` | `1.87` |
| `Linf_Ham_amr` (worst single point) | `2.367e-3` | `5.180e-3` | `2.19` |

Fourth-order convergence predicts `1/1.5⁴ = 0.20`.  Grid-scale noise
differentiated twice predicts `1.5² = 2.25`.  The data lands on `1.8–2.2`.
**The starting data carries noise at the cell scale, and the Hamiltonian
constraint — which takes second derivatives — amplifies it as `1/Δx²`.**

### Where the noise comes from

The elliptic solve ran on a box of `L = 128` while the evolution box is
`L = 64`, and the solve's cell count was hardwired to the evolution's `N`:

```
run_single_selfgrav.sh   --grtresna-domain-l 128      (hard-coded)
replay_eval.py           grtresna_nx = grtresna_ny = grtresna_nz = n_full
```

So the solve cell was **always exactly twice the evolution cell, at every
resolution**, and refining both together never closed the gap.  The evolution
then resolved the coarse solve's interpolation seams more sharply and
amplified them harder.  This is a fixed structural handicap, not a tuning
mistake — no choice of `N` escaped it.

### Why this matters more than any single run

**A resolution ladder on the old wiring is uninterpretable.**  Each finer rung
starts with a *larger* constraint violation than the coarser one it is being
compared against, so "the finer grid failed too" does **not** acquit
resolution, and a convergence order fitted through those rungs measures the
initial-data pipeline rather than the evolution.  Always read the `t = 0` row
of `constraint_norms.dat` before reading anything else in a convergence study.

### The fix — the solve grid is now configuration, not a constant

Three environment knobs on `run_single_selfgrav.sh`, backed by a new
`--grtresna-n` in `replay_eval.py`.  All default to the previously hard-coded
values, so an un-set environment reproduces every earlier cell bit-for-bit:

| knob | default | what it is |
|---|---|---|
| `BONDI_GRTRESNA_DOMAIN_L` | `128` | solve box width — wider than the evolution box on purpose, so the outer boundary condition sits far from the star |
| `BONDI_GRTRESNA_N` | `BONDI_NFULL` | solve cells per axis |
| `BONDI_GRTRESNA_MAXLEVEL` | `3` | refinement depth inside the solve |

Set `BONDI_GRTRESNA_N = BONDI_NFULL × (DOMAIN_L / LFULL)` to make the solve
cell equal the evolution cell.  The cell centres then coincide and the handoff
becomes a **straight copy instead of an interpolation**, which removes the
noise source rather than shrinking it.  `BONDI_GRTRESNA_MAXLEVEL=0`
additionally removes the solve's own refinement seams, the other candidate
source.  The verdict cell now running:

```bash
# solve dx = 128/384 = 0.333 == evolution dx = 64/192.  Uniform, seam-free.
BONDI_NFULL=192 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
  BONDI_GRTRESNA_DOMAIN_L=128 BONDI_GRTRESNA_N=384 BONDI_GRTRESNA_MAXLEVEL=0 \
  BONDI_GRTRESNA_RANKS=32 BONDI_NL_TOL=0.001 BONDI_GRTRESNA_TIMEOUT=21600 \
  BONDI_EXOTIC=0 BONDI_OMEGA=0.80 BONDI_STOP_TIME=120 \
  BONDI_SCRUTINY=1 BONDI_SPONGE=1 \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_single_selfgrav.sh
```

The number that decides it is the `t = 0` violation: the old cells started at
`5.6e-4` (`N = 128`) and `1.0e-3` (`N = 192`).  Well below `5.6e-4` confirms
the diagnosis; unchanged means the noise is coming from somewhere else.

### Correction: the solve tolerance is not the floor it was thought to be

Queue item 1 below states that every rung "exits at the same `0.0832%`
Hamiltonian violation" and concludes the tolerance is the floor.  That is true
of the **phantom** sector only.  The two sectors converge by different paths
and their percentages are normalised differently, so they are **not comparable
to each other**:

| sector | how it converges | exits at | gate `0.1%` binding? |
|---|---|---|---|
| canonical | alternating Ham/Mom, ping-pongs (`4.6e-15 → 12.5 → 1.4e-3 → 0.108 → 1.4e-3 %`) | `≈0.0014%` | **no** — beats it `70×` |
| phantom | geometric, `≈0.4` per iteration | `≈0.0844%` | **yes** — just barely |

Despite the `60×` gap in reported percent, the **absolute** violation entering
the evolution is nearly identical for both (`L2_Ham` `5.7e-4` … `8.2e-4`) —
which is the clue that pointed at the initial-data pipeline rather than the
stopping rule.  Tightening `BONDI_NL_TOL` is still worth doing for phantom
cells (it costs `≈5` CPU iterations of a permitted `50`, and phantom cores were
the ones that drifted off station), but it will not move the canonical star.

---
## What the campaign returned — 2026-08-20

All six cells ran, on four cards, `10:06` to `18:51`: `8.7` hours of wall clock,
`32.6` GPU-hours summed, no failures, `1.4` GB of output, every initial-data
file pruned on success.  Results are packed under
`results/bondi-dipole-runaway/campaign/next_*`.

Four of the five questions came back answered.  The fifth — the constant-gap
test — came back with something else instead, and it is more important than
anything the queue was designed to ask.

### The result that matters: the canonical star was never stable

**The positive-mass star collapses into a black hole under its own gravity, in
every run in the campaign, on the same timescale over which the acceleration is
measured.**  It is not being pushed into collapse by its phantom partner and it
is not a resolution artefact.  It was picked from the unstable half of the star
family and was going to do this regardless.

The long `d₀ = 16` cell is where it becomes unmissable, because that run goes to
`t = 120` instead of stopping at `60`:

| `t` | 0 | 40 | 55 | 60 | 70 | 75 | 85 |
|---|---|---|---|---|---|---|---|
| lapse at the centre | `0.98` | `0.82` | `0.57` | `0.40` | `0.11` | `0.047` | `0.08` |
| conformal factor `χ` | `0.98` | `0.81` | `0.49` | `0.27` | `0.027` | `0.036` | `0.030` |
| canonical peak amplitude | `0.0226` | `0.0233` | `0.0259` | `0.0290` | `0.0183` | `0.0000` | `0.0000` |
| phantom peak amplitude | `0.0223` | `0.0221` | `0.0221` | `0.0221` | `0.0221` | `0.0221` | `0.0221` |

The canonical core brightens as it squeezes down, then the field goes to zero —
swallowed.  The phantom does not move by `1%` across the entire run.  On
`movie_scalar_activity_z.mp4` this reads as the right-hand blob vanishing while
the left-hand one is still barely moving.

**Three independent lines of evidence say the collapse is intrinsic.**

*One — a star on its own, with no partner anywhere in the box:*

| alone, no partner | lapse `t = 0 → 40` | rms radius `t = 0 → 40` |
|---|---|---|
| canonical | `0.976 → 0.867`, falling monotonically | `5.05 → 13.40` |
| phantom | `1.001 → 0.998`, flat | `5.43 → 6.32` |

The canonical star destabilises with nothing near it.  Its lapse at `t = 40`
alone (`0.867`) is the same as in the `d₀ = 16` pair (`0.824`) and the `d₀ = 8`
pair (`0.853`) — the partner changes almost nothing.  Cells:
`published/single_p`, `published/single_m`.

*Two — stacking the deck one way or the other:*

| pair | worst lapse | worst `χ` | outcome |
|---|---|---|---|
| canonical + canonical | `0.014` | `0.003` | horizon almost at once (`lapse < 0.5` by `t = 4`) |
| canonical + phantom | `0.047` | `0.024` | horizon at `t ≈ 57–75` |
| phantom + phantom | `0.974` | `0.989` | nothing happens, ever |

More canonical matter collapses faster; more phantom matter is more stable; two
phantoms never collapse at all.  The ordering is also visible *within* the mixed
runs — `d₀ = 16` crosses `lapse < 0.5` at `t = 57` but `d₀ = 8` not until
`t = 62`, i.e. bringing the negative mass *closer* delays the collapse, exactly
as cancelling part of the attractive field should.

*Three — the star family table already said so.*  Along the canonical sequence,
`ADM` mass falls as the star gets denser at the frequency every run uses:

| `ω` | `φ_c` | `M_ADM` | `dM/dφ_c` | branch |
|---|---|---|---|---|
| `0.900` | `0.0159491` | `0.007963` | — | stable |
| `0.800` | `0.0191264` | `0.011209` | `+2.4` | stable |
| `0.775` | `0.0194775` | `0.012612` | `+4.0` | stable |
| **`0.550`** | **`0.0196947`** | **`0.063951`** | **`−182`** | **unstable ← used in every run** |
| `0.570` | `0.0197682` | `0.052411` | `−143` | unstable |
| `0.615` | `0.0199119` | `0.035175` | `−110` | unstable |

A negative slope here is the standard turning-point signature of the unstable
branch.  The stable branch is up at `ω ≈ 0.775–0.90`, where the stars are three
to five times lighter and roughly `50%` larger.  Source:
`results/bondi-dipole-runaway/stars/star_radius.csv`, already in the pack — no
new computation was needed to see this, only to look.

The phantom sits on the equivalent negative-slope branch (`−225` at `ω = 0.55`)
and is nonetheless perfectly stable, which is the physical point: a
negative-energy object's own gravity pushes *outward*, so it has no collapse
channel to lose to.  The instability is one-sided by construction.

### Two more things the campaign exposed

**The stars overlap at every separation the article uses.**  Each star's `99%`
radius is `8.72` (canonical) and `9.16` (phantom) at `ω = 0.55`:

| `d₀` | centres apart | sum of `r₉₉` | sum of `r₉₀` | verdict |
|---|---|---|---|---|
| `8` | `8` | `17.9` | `15.8` | centres closer than *one* radius — a single blob |
| `12` | `12` | `17.9` | `15.8` | heavily overlapping |
| `16` | `16` | `17.9` | `15.8` | envelopes touch; cores just clear |

`d₀ = 8` — the configuration behind the headline number — is not two objects.
And it does not hold a gap: it **merges**, `8.00 → 1.81` by `t = 60` and through
zero by `t = 70`.  `d₀ = 16` is the first separation where the bulk of the two
stars are apart, and there the drift at `t = 60` is `0.44` against `7.11` at
`d₀ = 8` — a factor `16` smaller, where an inverse-square force between
separated masses would predict a factor `4`.

**The constant gap at `d₀ = 16` is a coordinate statement only.**  Coordinate
separation holds to `2.8%` while the measured distance between the two grows by
half:

| `t` | 0 | 20 | 40 | 60 | 80 |
|---|---|---|---|---|---|
| coordinate separation | `16.00` | `15.94` | `15.81` | `15.68` | `15.56` |
| proper separation | `14.99` | `16.96` | `16.82` | `20.74` | `23.08` |

**The phantom's late acceleration is real but the wrong shape.**  Its speed
doubles roughly every `10` time units — `0.0018`, `0.0055`, `0.0157`, `0.0292`,
`0.0504`, `0.0906` over successive windows — taking off precisely as the
canonical star dies.  This is not coordinate drag: the shift at the phantom's
location is `0.003–0.005` against a coordinate speed reaching `0.12`, forty
times too small to explain it.  But a genuine runaway at fixed separation
requires *steady* acceleration, i.e. speed rising linearly.  Compounding growth
is not that.

### d₀ = 12 does show the Bondi signature — with the right diagnostic

The subsection above is about `d₀ = 16`, where the effect is genuinely tiny.  At
`d₀ = 12` it is not, and the campaign's own peak tracker sees it.  Same run, two
ways of locating each star:

| `d₀ = 12`, `N = 256` | canonical moves | phantom moves | ratio | Bondi predicts |
|---|---|---|---|---|
| core tracker, `t = 40` | `+0.404` | `+0.457` | **`1.13`** | `1.205` |
| core tracker, `t = 45` | `+0.637` | `+0.740` | **`1.16`** | `1.205` |
| core tracker, `t = 50` | `+0.912` | `+1.142` | **`1.25`** | `1.205` |
| barycentre, `t = 40` | `+0.094` | `+0.667` | `7.11` | `1.205` |
| barycentre, `t = 45` | `+0.218` | `+0.974` | `4.46` | `1.205` |

Bondi's kinematic prediction for an unequal pair is that both masses accelerate
the *same* way, with their accelerations in the ratio `M₊/|M₋| = 1/0.83 = 1.205`
— the lighter phantom moving more.  **The core tracker reproduces that to within
6 %**, while the separation holds to `1 %` (`12.00 → 11.90` at `t = 45`), the
canonical star is still within `3.7 %` of its initial peak at `t = 40` and
`5.5 %` at `t = 45`, and the displacements are `1.6–1.8` grid cells at
`Δx = 0.25` — comfortably resolved.

**The barycentre does not reproduce it at all** — ratios of `14`, `−5`, `18`,
`7`, `4.5`, converging on nothing.  It integrates the whole domain and drags in
the radiation halo.  Every Bondi signature in this campaign lives in the core
tracker and is invisible to the diagnostic the article currently quotes.  That
is a diagnostic change, not a physics problem, and it makes the case stronger.

**What still does not match is the time dependence.**

| `t` | 20 | 30 | 40 | 45 | 50 |
|---|---|---|---|---|---|
| drift | `0.031` | `0.133` | `0.432` | `0.690` | `1.029` |
| implied acceleration | `1.6e-4` | `3.0e-4` | `5.4e-4` | `6.8e-4` | `8.2e-4` |

The drift grows as `t^3.9`, not `t^2`.  A Bondi runaway is *constant*
acceleration at constant separation, which is `t^2`.  Direction right, mass
ratio right, growth rate wrong — and the separation closing cannot explain it
(`1 %` over the window is worth `2 %` on a `1/d²` force).  The most economical
reading is that the acceleration is still turning on while the star is already
destabilising, so the two effects overlap and the run ends before they separate.
**A stable star would settle this in one cell**, because the window would no
longer be capped at `t ≈ 45`.

### Item by item

**Items 1 + 2 — the `Δx⁴` tolerance ladder, and the peak tracker.  Answered:
the tolerance was not the problem.**  The tightened exit worked exactly as
designed — the solve now floors `15×` lower at `N = 256` instead of stopping at
the same place every published rung did:

| cell | Newton steps | achieved `L2` Ham |
|---|---|---|
| `next_tolB_n128` | 7 | `8.319968e-02` |
| `next_tolB_n192` | 9 | `1.336698e-02` |
| `next_tolB_n256` | 10 | `5.372809e-03` |

And it changed nothing.  Barycentre Richardson ratios came back `1.42 / 1.33 /
1.92` against the published `1.47 / 1.37 / 1.94`, and the spread at `t = 60`
stayed at `6.3%`.  **The campaign's failure to show a convergence order is not
a solve-accuracy problem** — it is that the quantity being converged is the
centroid of a star that is dispersing and collapsing, which has no continuum
limit to converge to.

The peak tracker is the better diagnostic and confirms this: its spread is
`3.3%` at `t = 20` growing to `5.0%` at `t = 60`, against the barycentre's
`12.5%` shrinking to `6.3%` — the barycentre's apparent "improvement" with time
is halo contamination cancelling, not convergence.

**Item 3a — the sponge.  Answered: safe, and unnecessary.**  Against its
sponge-off twin, integrated negative energy differs by `−0.03%` at `t = 20` and
`−0.13%` at `t = 90`; the Hamiltonian norm is unchanged to five decimals at
every output time; the initial-data residual matches the published value to
seven digits (`0.08318538`).  The sponge neither corrupts the interior nor
rescues anything.  (An earlier reading of "no difference at all" was my own
rounding — the difference is real, and it is this small.)

**Item 3b — the constant-gap test.  Not answered; superseded.**  See above.  The
run does what it was built to do — the coordinate gap holds — but the object
whose acceleration it was meant to measure turns into a black hole at `t ≈ 57`.

**Item 5 — the wider lever arm.  Answered, and the answer is uninterpretable as
a force law.**  Drift at `t = 60` against mass ratio is non-monotonic:

| phantom-to-canonical mass ratio | `0.62` | `0.83` | `1.00` |
|---|---|---|---|
| drift at `t = 60` | `5.66` | `6.41` | `6.23` |
| gap `8.00 →` | `0.84` | `1.86` | `1.42` |

All three are measured on pairs that have already merged.  The non-monotonicity
is merger dynamics, not force scaling.

### One correction to the note above

"Before launching anything" justifies `Δt = 0.02 Δx` by saying the star
disperses at `Δt = 0.2 Δx` (rms `5.05 → 19.2` by `t = 40`).  True, but the
comparison is weaker than it reads: at the *chosen* timestep the same lone star
still goes `5.05 → 13.40` over the same interval.  Most of that expansion is the
branch instability, not the timestep.  The Courant choice is still defensible —
lapse holds at `0.867` instead of crashing to `0.191` — but it should not be
quoted as the reason the star stays together, because it does not.

### Verdict — this matter model is a no-go

The setup is sound.  The solver converges, the sponge is clean, the tracker
works, and at `d₀ = 12` the pair reproduces Bondi's kinematic prediction to 6 %.
What fails is the star itself: it is drawn from the unstable half of its own
family and collapses on the same timescale the measurement needs.  Every
question the campaign could not close — no convergence order, no constant
acceleration, no clean force scaling, no run past `t ≈ 45` — traces back to that
single choice.

**Every frequency this project has ever used is on the unstable branch.**
Walking the canonical family from dilute to dense, the central field value rises
with mass until it turns over, and that turn is the stability boundary:

| `ω` | `φ_c` | `M_ADM` | `r₉₉` | branch |
|---|---|---|---|---|
| `0.900` | `0.0159491` | `0.007963` | `13.40` | stable |
| `0.800` | `0.0191264` | `0.011209` | `9.81` | stable |
| `0.775` | `0.0194775` | `0.012612` | `9.34` | stable |
| `0.750` | `0.0197205` | `0.014350` | `8.99` | stable |
| `0.690` | `0.0199808` | `0.020421` | `8.44` | stable |
| **`0.670`** | **`0.0199941`** | `0.023321` | `8.37` | **turning point — `φ_c` peaks here** |
| `0.650` | `0.0199812` | `0.026859` | `8.33` | unstable |
| `0.615` | `0.0199119` | `0.035175` | `8.30` | unstable ← the lever cell |
| `0.566` | `0.0197539` | `0.054459` | `8.53` | unstable ← `stars/phantom_omega0.566.dat` |
| **`0.550`** | `0.0196947` | `0.063951` | `8.72` | **unstable ← every run, published and follow-up** |
| `0.530` | `0.0196187` | `0.079391` | `8.96` | unstable |

Below `ω ≈ 0.67` the central density turns back while the mass keeps climbing.
That sign change is the turning-point criterion, and it is the whole story.
**The stable branch, `ω ≳ 0.67`, has never been run.**

### What a stable replacement looks like — sizes, separation, run length

The natural target is **`ω = 0.80`**, and it lands on three good things at once:
it is comfortably above the turning point, it is where the star's *bulk* is at
its smallest anywhere in the family, and it is where the two sectors happen to
come out nearly equal in mass.

Read the sizes carefully, because the two radii move in *opposite* directions.
`r₉₀` is the bulk — where the star's mass actually is, and what the core tracker
cares about.  `r₉₉` is the bulk plus the diffuse tail.

| `ω` | `M₊` | `\|M₋\|` | ratio | `r₉₀` c/p | `Σr₉₀` | `d₀/Σr₉₀` | `r₉₉` c/p | span | rel. accel | `t` for the same drift |
|---|---|---|---|---|---|---|---|---|---|---|
| **`0.55`** | `0.06395` | `0.07696` | `1.203` | `7.57 / 8.20` | `15.77` | **`0.76`** | `8.72 / 9.16` | `±15.2` | `1.000` | `45` ← now |
| `0.75` | `0.01435` | `0.01513` | `1.054` | `5.24 / 5.34` | `10.58` | `1.13` | `8.99 / 9.04` | `±15.0` | `0.224` | `95` |
| `0.775` | `0.01261` | `0.01323` | `1.049` | `5.14 / 5.24` | `10.38` | `1.16` | `9.34 / 9.39` | `±15.4` | `0.197` | `101` |
| **`0.80`** | `0.01121` | `0.01169` | **`1.043`** | **`5.11 / 5.18`** | `10.29` | **`1.17`** | `9.81 / 9.84` | `±15.8` | `0.175` | `107` ← pick |
| `0.85` | `0.00916` | `0.00946` | `1.033` | `5.71 / 5.70` | `11.41` | `1.05` | `11.07 / 11.10` | `±17.1` | `0.143` | `119` |

*(all rows computed at `d₀ = 12`, `L = 64`, `N = 256`, `Δx = 0.25`)*

**Size.**  The bulk gets *smaller*, not bigger — `r₉₀` falls from `7.57 / 8.20`
to `5.11 / 5.18`, a third off.  `ω = 0.80` is the exact minimum of `r₉₀` for
**both** sectors across the whole scan; this family's radius is famously
non-monotonic in `ω` (it inflates in both the thin-wall and thick-wall limits),
and `0.80` is the bottom of that curve.  Only the diffuse tail grows,
`8.72 → 9.81`.

**Separation: keep `d₀ = 12`.**  Measured against the bulk, `d₀ = 12` at
`ω = 0.55` gives a clearance of `0.76` — the two stars are *inside* each other,
which is the whole reason the published `d₀ = 8` cell is a single blob.  At
`ω = 0.80` the same `d₀ = 12` gives `1.17`, so the bulk is genuinely clear.
**`d₀ = 12` on the stable branch is better separated than `d₀ = 16` is today**
(`16 / 15.77 = 1.01`, only just touching), while keeping the much stronger
signal of the shorter separation.  There is no reason to widen it.

**The tails will always overlap, and that is not fixable here.**  `Σr₉₉ = 19.65`
against `d₀ = 12`.  Clearing the `99 %` envelopes needs `d₀ ≳ 20`, which needs
`L = 128`, which at `Δx = 0.25` is `N = 512` and roughly `370` GB — four times
an `80` GB card.  Accept the tail overlap and measure with the core tracker,
which is halo-free by construction; the barycentre is exactly the diagnostic
that cannot survive this and it is the one being retired.

**Box and resolution are unchanged.**  The pair spans `±15.8` against a `±32`
boundary — comfortable, with room for the sponge band.  `L = 64`, `N = 256`,
`Δx = 0.25`, `max_level = 0`, `Δt = 0.02 Δx`, exactly as now.

**Run length is the only real cost.**  Stars `5.7×` lighter accelerate `5.7×`
more weakly at the same separation, so the same drift needs `t ≈ 107` instead of
`45`.  `spongeB_sep16_long` already ran `t = 120` at this grid size, so the
budget is known.  And this time running long is *allowed*, because the star does
not collapse at `t ≈ 57`.

**If `t ≈ 107` is too long, `ω = 0.75` is the trade:** `t ≈ 95`, bulk clearance
`1.13`, mass ratio `1.054`, all still on the stable side — just with less margin
above the turning point.  Do not go below `ω = 0.70`.

**The mass ratio is the quiet prize.**  Bondi's constant-gap runaway is strictly
an *equal-mass* result.  At `1.203` the present pair is `20 %` away from it and
the gap closes for that reason alone, with no numerical excuse needed.  On the
stable branch the two sectors come out within `4–5 %` of equal mass with no
detuning at all — the configuration the theorem is actually about, for free.

### Is the phantom stable too?

**Yes, and it is the best-established fact in the whole dataset.**  Over the
longest run in the campaign — `t = 0 … 120`, sitting next to a star that turned
into a black hole — the phantom's core peak never moves:

| `t` | 0 | 20 | 40 | 60 | 80 | 100 | 120 |
|---|---|---|---|---|---|---|---|
| phantom peak | `0.0223` | `0.0222` | `0.0221` | `0.0221` | `0.0221` | `0.0222` | `0.0223` |
| deviation | `0.00 %` | `−0.40 %` | `−0.86 %` | `−0.86 %` | `−0.56 %` | `−0.43 %` | `+0.32 %` |

Alone in an empty box its lapse never leaves `0.998 – 1.001`.  It has no
collapse channel: its energy is negative, so its own gravity pushes outward and
there is nothing for field pressure to lose to.

Note that at `ω = 0.55` the phantom is *formally* past its own turning point
(which sits at `ω ≈ 0.67`, essentially the same place as the canonical one) and
is nonetheless perfectly stable.  That is the criterion failing to apply, not
the star getting lucky — the turning-point argument is derived for
positive-energy configurations.  Moving to `ω = 0.80` puts the phantom on the
formally stable side as well, so there is no trade-off in the choice.

**Two caveats, both about the halo rather than the star.**  The phantom sector's
rms radius does grow — `5.43 → 17.6` by `t = 100` — while its core sits
motionless.  That is halo and debris being counted by a domain-integrated
quantity, and it is one more reason every measurement on a long run has to come
from the core tracker.  And the lone-phantom control has only ever been run to
`t = 40`; **nobody has evolved a phantom on its own past that**.  If the plan is
`t ≈ 110`, run the single-star control for *both* sectors, not just the
canonical one.

### Do this first, before booking any GPU time

**One lone stable-branch star, evolved on its own.**  This is the control the
project has never had, it is a single cheap cell, and it is decisive: solve at
`ω = 0.80`, evolve to `t = 120`, and require the lapse to stay flat and the rms
radius to stay put — the way `single_m` behaves (`0.998`, `5.43 → 6.32`) and the
way `single_p` does not (`0.867`, `5.05 → 13.40` by `t = 40`, with no partner
anywhere in the box).

If it holds, the pair campaign is worth running and the whole question reopens
on a sound footing.  If it does not, the potential itself is the problem rather
than the branch, and no separation, box or resolution will rescue it.

Everything downstream — the tolerance ladder, the force-scaling lever, the
energy budget — should wait for that one cell.  Re-running any of them on
`ω = 0.55` matter would only measure the collapse again.

> **2026-08-21 — this campaign now exists.**  `run_stable_campaign.sh` (same
> folder as the other launchers) writes four cells into the new tree
> `runs/bondi_correct` — lone canonical and lone phantom at `ω = 0.80` (the
> design point) and at `ω = 0.75` (the higher-signal backup) — each `N = 128`
> unigrid at `L = 64` (`Δx = 0.5`, the published rung: a yes/no screen at the
> grid every published cell used, ~8× cheaper than the design grid), `t = 120`,
> sponge on at 24/32, core tracker on, published solve tolerance.  Launch with
> `bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_stable_campaign.sh`
> (four cards, one cell each; subset/dry-run knobs in the header).  Pass gate
> per cell, over the full run: core peak amplitude flat to ±2%, lapse steady
> near 1, core parked at `x = +4` to within one grid cell.  A canonical fail
> at **both** frequencies condemns the rung, not the frequency — move rung
> before spending more GPU time.  `runs/bondi_rerun` is closed as the
> `ω = 0.55` archive; nothing new goes there.


---
## How to run a campaign

Everything below is the operating manual for the loop *find a stable star →
evolve it → check it stayed put*.  Three scripts do all of it.

| script | what it is |
|---|---|
| `grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh` | one cell: CPU solve, then GPU evolve |
| `grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_next_campaign.sh` | many cells: writes job files, starts the queue detached |
| `grteclyn-wrapper/scripts/campaigns/stop_campaign.sh` | stops a detached campaign properly |

### One cell, in the foreground

This is the right shape for a debugging cell — you watch it, it dies with your
shell, nothing detaches:

```bash
BONDI_S0=0 BONDI_S1=1 BONDI_GPU=0 \
  BONDI_SEP=12 BONDI_NFULL=256 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
  BONDI_STOP_TIME=120 BONDI_SCRUTINY=1 \
  BONDI_RUNS_DIR="$PWD/runs/bondi_rerun/next/mycell" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

Add `BONDI_DRYRUN=1` to print the assembled command and run nothing.  A cell
does its own CPU solve first (`8` MPI ranks, `2–20` min depending on `N`) and
only then takes the card.

### The knobs

| variable | default | what it does |
|---|---|---|
| `BONDI_S0` / `BONDI_S1` | `0` / `1` | sector per lump: `0` canonical, `1` phantom.  `0 1` is the mixed pair, `0 0` two canonical, `1 1` two phantom |
| `BONDI_SEP` | `8` | initial centre-to-centre separation `d₀` |
| `BONDI_NFULL` | `128` | cells per side |
| `BONDI_LFULL` | `64` | box side; `Δx = L/N` |
| `BONDI_MAXLEVEL` | `1` | **set to `0`** — refinement blows the memory budget, see below |
| `BONDI_STOP_TIME` | `60` | final time |
| `BONDI_DT_MULT` | `0.02` | Courant factor |
| `BONDI_SCRUTINY` | `0` | **set to `1`** — switches on the halo-free core tracker (`sector_dynamics.dat`).  Costs `~1.6 s` per plotfile and it is the only diagnostic that shows the Bondi signature |
| `BONDI_S1_OMEGA` | *(unset)* | per-lump frequency for lump1 only, used for equal-mass cells |
| `BONDI_RADII` | `8 16` | radiation extraction shells |
| `BONDI_PSI4_HIGHER_L` | `0` | also record `ℓ = 3, 4` |
| `BONDI_SPONGE` | `0` | boundary damping; band via `BONDI_SPONGE_INNER` / `_OUTER` (`24`/`32` default, use `40`/`60` at `L = 128`) |
| `BONDI_NL_TOL` | `0.1` | elliptic solve exit tolerance, in percent |
| `BONDI_NL_STALL_TOL` | `0.002` | give-up threshold; keep it at `1/50` of the exit tolerance |
| `BONDI_GRTRESNA_RANKS` | `8` | MPI ranks for the solve.  `8` is `7×` faster than `1` and gives the same answer to `4` digits |
| `BONDI_GRTRESNA_TIMEOUT` | `7200` | seconds; raise to `21600` when tightening the tolerance |
| `BONDI_GPU` | `3` | CUDA device |
| `BONDI_RUNS_DIR` | derived | where the cell lands |
| `BONDI_DRYRUN` | `0` | print, do not run |

### Several cells, queued across the cards

`run_next_campaign.sh` writes one job file per cell into
`runs/bondi_rerun/next/_queue/pending/` and starts the queue runner detached.
Each job runs in the foreground *inside* its slot; the runner is the only thing
detached, once.  To add a cell, add a `case` arm to its `cell_env()` and a name
to `ALL_CELLS`.

```bash
# see what it would enqueue, start nothing
BONDI_NEXT_DRYRUN=1 bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_next_campaign.sh

# a subset, on chosen cards
BONDI_NEXT_CELLS="mycell" BONDI_NEXT_GPUS="0 1" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_next_campaign.sh
```

A cell whose run directory already exists is skipped, so re-running the launcher
is safe and resumes rather than duplicates.  Jobs dispatch in filename order —
the launcher prefixes `010_`, `020_`, … and lists cells longest-first so the
short ones backfill.

### Watching it

```bash
pgrep -af gpu_queue.sh                              # the runner must exist
tail -f runs/bondi_rerun/next/_queue/queue.log      # dispatch events
tail -f runs/bondi_rerun/next/_queue/logs/*.log     # per-cell output
nvidia-smi                                          # one cell per card
ls runs/bondi_rerun/next/_queue/{pending,running,done,failed}
```

**A silent launch looks exactly like success.**  Always confirm with `pgrep`.

### Stopping it

```bash
bash grteclyn-wrapper/scripts/campaigns/stop_campaign.sh --dry-run runs/bondi_rerun/next
bash grteclyn-wrapper/scripts/campaigns/stop_campaign.sh runs/bondi_rerun/next
pgrep -af 'gpu_queue|main3d|GRTresna'               # verify, then escalate
```

Order matters and the script enforces it: freeze the queue first, then sweep the
workers.  **Killing a running cell on its own does not stop the campaign** — it
looks like a finished step and the runner immediately starts the next one.  By
hand, `touch` the queue's `STOP` file first, then kill the runner from
`runner.pid`, then the workers.

### Traps that have actually bitten

- **The wrapper's Python venv can silently lose its interpreter.**  The venv
  is built by `uv` and its `python` is a symlink into a `uv`-managed
  interpreter under the user's home; a cleanup or node change can delete that
  interpreter while the venv's `lib/` survives.  Symptom:
  `.venv/bin/python: No such file or directory` from every launcher, on a tree
  that ran fine the day before (bit on 2026-08-21 — the campaign of 2026-08-20
  had used the same path successfully).  Fix from the wrapper directory:
  `uv venv --python <a surviving interpreter> && uv sync --frozen`, which
  rebuilds from the lockfile in under a minute.

* **`env` is shadowed on this node.**  Other users' `bin` directories precede
  `/usr/bin` on `PATH` and hold an `env` that is a `PATH`-setup snippet meant to
  be sourced.  Run as `env VAR=x cmd` it prepends its own directory and exits `0`
  **without executing `cmd`** — a launch that silently does nothing and reports
  success.  Spell out `/usr/bin/env`, or drop it.
* **`scripts/lib/env.sh` exports its own `SCRIPT_DIR` and `cd`s to the repo
  root.**  Any script that sources it must re-derive `SCRIPT_DIR` from
  `BASH_SOURCE[0]` straight afterwards, or every sibling path resolves against
  `scripts/lib/`.  It has to be sourced, though — it is the only thing that puts
  `mpirun` on `PATH`, so a multi-rank solve fails without it.
* **One cell per card.**  Two `192³` cells sharing an H100 each ran `2.2×`
  slower, so the pair finished later than running them back to back.  Never
  co-schedule anything with an `N = 256` cell.
* **Memory tracks cell count, not box size.**  `5.8` / `19.5` / `46.2` GB of FAB
  at `N = 128` / `192` / `256`, peaking near `57` GB — identical whether
  `L = 64` or `L = 128`.  A refined `N = 128` cell costs `41` GB against `5.8`
  uniform, which is why `BONDI_MAXLEVEL=0`.
* **Prune the initial data.**  `initial_data.gridinit` is `4.1` GB per `N = 256`
  cell and is regenerable from the cell's own `params.txt`.  Nothing downstream
  reads it.  The campaign launcher deletes it on success; do the same by hand.

---
## The matter configuration

This is the thing being debugged, so all of it is here.

### What the star is

A self-gravitating **sextic (solitonic) Q-ball** — a complex scalar field with a
sixth-order potential, solved at fixed frequency `ω` by amplitude shooting, then
dressed with its own gravity.  The same solver produces both sectors: the sign
of the kinetic term flips, giving the **canonical** star (positive energy,
attractive self-gravity, `lapse < 1` in the core) or the **phantom** star
(negative energy, repulsive self-gravity, `lapse > 1` in the core, negative ADM
mass).

### The exact parameters every run has used

Set inside `run_pair_selfgrav.sh` as `--extra-override` flags:

| parameter | value | meaning |
|---|---|---|
| `grtresna_scalar_lambda` | `10240` | sextic potential coupling |
| `grtresna_scalar_mu` | `21845333` | scalar mass parameter |
| `grtresna_bs_omega` | **`0.55`** | star frequency — **hard-coded, see below** |
| `grtresna_bs_selfgrav` | `1` | dress the star with its own gravity |
| `trajectory_num_lumps` | `2` | two stars |
| `trajectory_lump0_R0` | `d₀/2` | lump 0 at `+x` |
| `trajectory_lump0_phase0` | `0` | …on the `+x` axis |
| `trajectory_lump0_exotic` | `BONDI_S0` | `0` canonical, `1` phantom |
| `trajectory_lump1_R0` | `d₀/2` | lump 1 at `−x` |
| `trajectory_lump1_phase0` | `π` | …on the `−x` axis |
| `trajectory_lump1_exotic` | `BONDI_S1` | `0` canonical, `1` phantom |
| `trajectory_lump1_bs_omega` | `BONDI_S1_OMEGA` | per-lump frequency override, lump 1 only |
| `trajectory_lump{0,1}_well_depth` | `0.15` | confinement well |
| `trajectory_well_width` | `1.2` | confinement well |
| every velocity / rotation / breathing knob | `0` | the pair starts **at rest** |

So in the default mixed pair, **lump 0 is the canonical star at `+x` (the
right-hand blob) and lump 1 is the phantom at `−x` (the left-hand blob).**

### The one change the stable branch needs

**Done, 2026-08-21 — nothing about the matter is hard-coded any more.**  Both
launchers take the full matter configuration from the environment, with the
published values as defaults so an un-set environment reproduces every old
cell bit-for-bit:

| knob | launcher | default | meaning |
|---|---|---|---|
| `BONDI_S0_OMEGA` | pair | 0.55 | base star frequency, both lumps |
| `BONDI_OMEGA` | single | 0.55 | the lone star's frequency |
| `BONDI_SCALAR_LAMBDA` | both | 10240 | potential rung |
| `BONDI_SCALAR_MU` | both | 21845333 | potential rung |

A non-default frequency appends `_wNNN` to the run name (`bondi_sg_single_m_w080`),
so stable-branch cells can never clobber the published ones.  The single-star
launcher was also brought up to the pair launcher's level while it was open:
multi-rank elliptic solve (default 8, it was still hard-wired to 1 with the old
node's "mpirun segfaults" note), sponge knobs, solve-tolerance knobs, and
`BONDI_DRYRUN=1`.

Commit before launching anything, so `metadata.json`'s `git_commit` still
identifies the tree that produced the cell.

### The family, and which half of it is safe

`results/bondi-dipole-runaway/stars/star_radius.csv` holds the whole scan and is
the reference for every number in the *Verdict* section above.  Regenerate or
extend it on the CPU — no GPU time:

```bash
PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
  results/bondi-dipole-runaway/analysis/star_radius_scan.py   # radii + compactness
PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
  results/bondi-dipole-runaway/analysis/star_family_scan.py   # M(omega), both sectors
```

The stability rule, applied to that file: walk the family from dilute to dense
and watch the central field value `φ_c`.  While it **rises** with mass the
branch is stable; once it **turns back** while mass keeps climbing, everything
beyond is unstable.  The turn is at **`ω ≈ 0.67`**.  Use `ω ≳ 0.70`, prefer
`0.775–0.80`.  Never go below `0.67` again.

Note the radius is *not* monotonic in `ω` — sextic Q-balls inflate in both the
thin-wall and thick-wall limits, so "higher `ω` means more compact" is only true
on one side of the minimum.  `star_radius_scan.py` exists precisely to check
that before anyone relies on it.

### How to tell a star is stable, in one cell

Solve one star, evolve it **alone**, and compare against the two references
already on disk:

| | lapse `t = 0 → 40` | rms radius `t = 0 → 40` |
|---|---|---|
| `published/single_m` — phantom, stable | `1.001 → 0.998` | `5.43 → 6.32` |
| `published/single_p` — canonical, unstable | `0.976 → 0.867` | `5.05 → 13.40` |

A stable canonical star must look like the first row, and must still look like
it at `t = 120`.  Watch:

| file | column | what a stable star does |
|---|---|---|
| `data/collapse_diagnostics.dat` | `min_lapse` | stays within a few `%` of `1` |
| `data/collapse_diagnostics.dat` | `min_chi` | stays near `1`; heading for `0` is a horizon |
| `small_data/sector_barycenters.dat` | `rms_radius_*` | flat.  A `2×` growth by `t = 40` is the instability |
| `small_data/sector_dynamics.dat` | `peak_canon` / `peak_phantom` | flat.  A **rise** means the core is squeezing down — collapse begins there, before the lapse shows it |
| `data/constraint_norms.dat` | `L2_Ham` | does not run away |

**The peak amplitude is the earliest warning.**  In the `d₀ = 16` cell it was
already `+28 %` at `t = 60`, while the lapse still read `0.40` and the run still
looked alive.

### Pairing two stars of equal mass

Bondi's constant-gap runaway is strictly an equal-mass result, and the two
sectors do not weigh the same at the same frequency.  Two ways to match them:

* **Detune one lump** — `BONDI_S1_OMEGA`.  At `ω = 0.55` the phantom must run at
  `ω = 0.56598` to weigh the canonical star's `0.0640`.  This is what the
  `*_eqm` cells do.
* **Move up the family** — on the stable branch the sectors come out naturally
  within `4–5 %` of equal mass (`|M₋|/M₊ = 1.049` at `ω = 0.775`, `1.043` at
  `0.80`), against `1.203` at `ω = 0.55`.  No detuning needed.

### Sizing a pair so the stars are actually separate

Each star's `99 %` radius is `≈ 8.7–9.9` across the useful range, so `d₀` has to
be read against `r₉₉(canonical) + r₉₉(phantom)`, not chosen for roundness.  At
`ω = 0.55` that sum is `17.9`, which is why `d₀ = 8` is a single blob and
`d₀ = 12` is the smallest separation with a defensible core-to-core measurement.
Keep `d₀ / Σr₉₉ ≳ 0.67`, and keep the pair clear of the boundary: stars span
`±(d₀/2 + r₉₉)`, which must stay well inside `±L/2`.

---
## Before launching anything

> *Everything from here to the end of the document is the pre-campaign queue,
> kept verbatim.  It was written against `ω = 0.55` matter, so treat the cell
> definitions as historical: the reasoning is still worth reading, the
> frequencies are not.*

**Code reaches the node by `git push` then `git pull` — never `scp`.**  Every
item below that needs a code or launcher change must be committed on
`feature/interstellar`, pushed, and pulled on the node.  A cell launched from a
tree that does not match a commit cannot be reproduced, and the campaign's
provenance (`metadata.json` records `git_commit` per cell) becomes a lie.

**One run per card.**  Two `192³` cells sharing an H100 each ran `2.2×` slower,
so the pair finished `9%` later than running them back to back — at these grid
sizes a single evolution already saturates memory bandwidth.  Budget `5.8`,
`19.5` and `46.2` GB at `N = 128`, `192`, `256` (peak `≈57` GB at `N = 256`),
and do not co-schedule anything with an `N = 256` cell.

**Keep `max_level = 0` and `Δt = 0.02 Δx`.**  Both are measured constraints, not
caution: a refined `N = 128` cell takes `41` GB against `5.8` GB uniform, which
at `N = 256` overruns an `80` GB card; and at `Δt = 0.2 Δx` the star disperses
outright (rms `5.05 → 19.2` by `t = 40`).

---

## The queue, in priority order

### 1. A resolution ladder with the solve tolerance scaled as `Δx⁴`

> **Superseded in part — 2026-08-21.**  The premise below (that the shared
> tolerance is the error floor) holds for the phantom sector but not the
> canonical one, and a more basic problem sits underneath it: on the old
> wiring the solve cell was always twice the evolution cell, so each rung of
> this ladder starts *dirtier* than the one below it.  Fix the solve grid
> first — see
> [The initial data is the error floor](#the-initial-data-is-the-error-floor--2026-08-21)
> — or the ladder measures the initial-data pipeline instead of the evolution.

**Why.**  This is the weakest technical link in the paper and the one a referee
goes for first.  The campaign quotes no convergence *order*, and the reason is
now measured: each rung of a ladder gets its own elliptic solve at its own
resolution, but all three stop on the same criterion —
`NL_exit_tolerance = 0.1` — and therefore all three exit at the same `0.0832%`
Hamiltonian violation, identical to four digits at `N = 128`, `192` and `256`.
The floor is the **tolerance**, not a shared grid, so Richardson extrapolation
against it measures the solver's stopping rule rather than the evolution.

**What it unblocks.**  A quotable fourth-order convergence order — turning
Sec. V C's honest caveat into a standard NR validation, and the article from a
proof-of-concept into a reference.

**Blocker — a launcher change first.**  The tolerance is hard-coded at
`run_pair_selfgrav.sh:121` (`--grtresna-nl-exit-tolerance 0.1`).  It needs a
`BONDI_NL_TOL` knob alongside the other `BONDI_*` overrides.  Watch
`--grtresna-nl-stall-tolerance 0.002` as well: the stall guard may halt the
iteration before a tight exit tolerance is reached, in which case it has to come
down with it.  The solves currently converge in `7` of a permitted `50`
iterations, so there is headroom.

| rung | `Δx` | tolerance | evolution | solve (CPU, 1 rank) |
|---|---|---|---|---|
| `N = 128` | 0.50 | `0.1` (unchanged) | 20 min | 14 min, longer once tightened |
| `N = 192` | 0.33 | `0.019` | 86 min | 32 min, ditto |
| `N = 256` | 0.25 | `0.00625` | 266 min | 70 min, ditto |

**Cost.**  Three cells, `≈6.2` GPU-hours for the evolutions, whatever the
tightened solves add on the CPU side.  Cheap for what it buys.

**Run it together with item 2.**  Both are the same configuration --- `d₀ = 12`,
`L = 64`, uniform, the `128/192/256` ladder --- differing only in which flag is
set, and `BONDI_SCRUTINY` adds a diagnostic stream without touching the
evolution.  One set of three cells therefore delivers both results.  The only
thing lost is that the halo bias then comes from cells solved at the new
tolerance rather than the published one; the tolerance change perturbs the
initial data slightly, and core-against-barycentre is a within-run comparison
anyway, so this is a footnote rather than an obstacle.

```bash
# after BONDI_NL_TOL exists, is committed, pushed and pulled;
# one cell per card, tolerance per the table above
BONDI_SEP=12 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_NL_TOL=0.00625 \
  BONDI_S0=0 BONDI_S1=1 BONDI_GPU=0 \
  BONDI_RUNS_DIR="$PWD/runs/bondi_tol/convB_pm_sep12_n256" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

### 2. The reference cell with the peak tracker enabled

**Why.**  Sec. VII now bounds the halo bias with a halo-free tracker — the
sub-grid-resolved field peak in `sector_dynamics.dat` — which puts the pair
drift `7–10%` *above* the barycentric one, reading the canonical star `40%` low
at `t = 20`.  But that stream was only ever enabled on four cells, all at
`d₀ = 8`, so the bound is measured at a separation the article does not quote
and carried over to `d₀ = 12` as an inference.

**What it unblocks.**  The halo bias measured at the reference configuration
itself, and headline drifts quoted core-to-core rather than as lower bounds.
This is the cheapest credibility gain in the queue.

**No blocker.**  `BONDI_SCRUTINY=1` already exists; it costs `≈1.6` s per
plotfile.

**Cost.**  `4.4` GPU-hours for the `n256` cell alone, `6.2` for the full
ladder --- or nothing extra, if it rides on item 1's three cells.

```bash
BONDI_SCRUTINY=1 BONDI_SEP=12 BONDI_NFULL=256 BONDI_MAXLEVEL=0 \
  BONDI_S0=0 BONDI_S1=1 BONDI_GPU=1 \
  BONDI_RUNS_DIR="$PWD/runs/bondi_core/convA_pm_sep12_n256_scrutiny" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

### 3. Switch the sponge on, then the constant-gap test

**Why.**  The one headline the paper cannot claim is the textbook Bondi
runaway — a pair accelerating forever at *fixed* separation.  It is not reached
at any separation run: at `d₀ = 8` the stars overlap, at `d₀ = 12` the gap
closes mildly, and at `d₀ = 16` the gap is nearly constant but the drift is not
grid-stable.  The blocker is time, not resolution.  Massive-field radiation
cannot leave through massless-wave (Sommerfeld) boundaries, washes back over the
stars and caps every mixed run at `t = 60`.

**What it unblocks.**  `t ≫ 60` at `d₀ = 16`, which is the only configuration
where the gap stays open long enough for the constant-gap claim to be tested.
Also the energy budget (item 4), whose extraction window closes at `u ≈ 25` for
the same reason.

**Not blocked — the sponge already exists, and this campaign never switched it
on.**  `Source/Grids/SpongeZone.hpp` has been in the tree since `1b567cf1`
(2026-07-10), five weeks before these runs, and it is wired into
`Examples/RadialRecipe` — the very executable the Bondi cells use
(`RadialRecipeMatterDispatch.hpp:348`).  It is a radially-ramped band of extra
Kreiss--Oliger dissipation between `sponge_inner_radius` and
`sponge_outer_radius`, off by default, exposed as five parameters, and already
in production elsewhere: the `gw_beam` campaign runs it at `24/32`,
`strength = 4.0`, quartic ramp, precisely so it can reach `stop_time = 40`
without boundary reflections, and the wormhole campaign has a `--sponge` flag.
It reaches the Bondi launcher through the same `--extra-override` channel that
already carries `dt_multiplier`.

**Two things to establish on the first cell rather than assume.**  The sponge is
extra dissipation, not a true outgoing boundary, and it was validated against
*massless* wave reflections in a canonical-matter campaign — whether it clears
this campaign's *massive* scalar bath is the empirical question the first run
answers.  And the default `24/32` band is sized for `L = 64`: the Weyl shells at
`R = 8` and `16` sit safely inside it, but the canonical envelope's rms radius
crosses `r = 16` by `t ≈ 51`, so at late times a sponge placed there starts
eating the star's own halo.  That moves the barycentric diagnostic, which
integrates over the domain — so run item 2's peak tracker alongside any sponge
cell, or the halo the sponge removes cannot be told apart from the halo the
diagnostic was mis-counting.

**Cost.**  `≈3.4` GPU-hours per `L = 128` cell to `t = 90` at `Δx = 0.5`,
`≈4.5` to `t = 120`.  Launchable today, and with four cards it can run alongside
items 1 and 2.  Two cells are the minimum: one repeat of `boxC_pm_L128_n256`
with the sponge on, which is the validation against its existing sponge-off
twin, and one long `d₀ = 16` cell, which is the test itself.

**A decision to take before counting further cells.**  `d₀ = 16` is quoted as a
bound today because its drift spreads `53%` across the resolution ladder, and
making the constant-gap result a *measurement* rather than a bound means a
grid-stable `d₀ = 16` --- which does not simply cost two more cells.  A ladder
at `L = 128` runs out of card: `N = 384` needs `≈156` GB against the `80` GB
available.  Either the ladder stays at `L = 64` with the sponge moved inward,
and accepts that it will absorb some of the stars' own halo, or the `L = 128`
ladder is built downward in resolution from `Δx = 0.5`, which is coarser than
the reference cell.  Neither is free; pick one before booking the cells.

```bash
BONDI_SCRUTINY=1 BONDI_SEP=16 BONDI_LFULL=128 BONDI_NFULL=256 \
  BONDI_MAXLEVEL=0 BONDI_RADII="16 24 32 40" BONDI_STOP_TIME=120 \
  BONDI_S0=0 BONDI_S1=1 BONDI_GPU=2 \
  BONDI_RUNS_DIR="$PWD/runs/bondi_sponge/pm_sep16_sponge" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
# with the launcher passing, alongside its existing overrides:
#   --extra-override sponge_enabled=1 \
#   --extra-override sponge_inner_radius=40 \
#   --extra-override sponge_outer_radius=60 \
#   --extra-override sponge_strength=4.0 \
#   --extra-override sponge_ramp_power=4
```

### 4. ADM surface integrals for the energy budget

**Why.**  The article establishes that a null signal leaves the system —
`1/r` falloff across a factor `2.5` in radius, propagation speed `0.87–1.03`,
box-independent to `0.18%` — but not the *sign* of the energy it carries.  That
matters more than its size: the pair starts at `M₊ + M₋ = −0.013`, the
news-squared flux is positive-definite, and the theorem that would bound the
descent is the positive-mass theorem, whose hypothesis this matter violates by
construction.

**Blocker — mostly not GPU work.**  The odd-`ℓ` half is *already on disk*:
`boxC_pm_L128_n256/psi4_mode_higher_l.dat` carries `ℓ = 3` and `ℓ = 4`, all
`m`, on all four shells, `227` rows to `t = 90`, and nothing reads it yet.
**Analyse that before launching anything.**  The longer extraction window is
now item 3's parameter rather than a code change, so the only code left here is
the ADM surface-integral diagnostic itself.

**Cost.**  Zero until the diagnostic exists; then it rides along with item 3.

### 5. A wider mass lever arm

**Why.**  The force-scaling test in Sec. V D compares `PM` against `PM-eq`,
which differ by `17%` in the phantom's mass — enough to establish the scaling
but a short lever arm, and the obvious referee question.  The limit is what was
run, not what is solvable: the dressed phantom family spans
`|M₋|/M₊ = 0.13` to `1.55`, so a factor-of-two contrast is available downward
(the heavy branch ends just below `ω = 0.53`).

**Cost.**  One cell, `≈4.4` GPU-hours at `N = 256` --- the partner it is
measured against, `convA_pm_n256`, is already in the pack --- plus a CPU-side star solve
at the new frequency — `ω = 0.615` gives `|M₋|/M₊ = 0.62`.  The trade to state
up front: lighter phantoms are more diffuse, so a longer lever arm buys a worse
point-mass comparison.

```bash
BONDI_SEP=8 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_S1_OMEGA=0.615 \
  BONDI_S0=0 BONDI_S1=1 BONDI_GPU=2 \
  BONDI_RUNS_DIR="$PWD/runs/bondi_lever/pm_half" \
  bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
```

### 6. The contact end state

**Why.**  Every mixed run was stopped before contact, so what a
canonical–phantom collision produces — including whether a horizon forms — is
untouched.  Explicitly out of scope for this article; listed so it is not
mistaken for an oversight.

---

## Do **not** spend GPU time on these

Everything here runs on the CPU, against data already in the pack:

* **the `ℓ = 3`/`ℓ = 4` waveforms** — see item 4, already recorded and unread;
* **every figure and quoted number** — `SCRATCH=/tmp python3 make_article_figures.py`
  regenerates all of it, including the propagation speeds and the Sec. VII
  caveat numbers (`caveat_stats()`);
* **the derived tables** — `analysis/{make_tables,momentum_balance,separation_scaling,newtonian_reference,constraint_check,convergence_check,wave_check}.py`;
* **the star family and radius scans** — `analysis/star_{family,radius}_scan.py`,
  which is where a new frequency for item 5 comes from;
* **initial-data solves** — single-rank CPU work throughout (`10.4` CPU-hours
  for the whole published campaign).  They gate the GPU cells but never occupy
  a card.

## Total if the queue is run in full

| item | cells | GPU-hours | blocked on |
|---|---|---|---|
| 1 + 2. tolerance ladder with the tracker on | 3 | 6.2 | `BONDI_NL_TOL` knob |
| 3a. sponge validation against the existing twin | 1 | 3.4 | nothing — `sponge_enabled=1` |
| 3b. long `d₀ = 16` constant-gap cell | 1 | 4.5 | nothing |
| 4. energy budget | 0 | rides on 3a/3b | the surface-integral diagnostic; read the `ℓ = 3` data first |
| 5. wider lever arm | 1 | 4.4 | a CPU star solve at `ω = 0.615` |
| **total** | **6** | **≈18.5** | |

**Six cells, `≈18.5` GPU-hours** — a third of what the published campaign cost,
and about five hours of wall clock with one cell per card.  All six are
launchable today or after a one-line launcher change; only item 4 needs code
written first, and even there the analysis that has to come before the code is
already paid for.

That count is the *minimum that answers each question once*.  It does not
include making `d₀ = 16` grid-stable, which is what would turn the constant-gap
result from a bound into a measurement — see item 3 for why that is a planning
decision rather than two more cells.
