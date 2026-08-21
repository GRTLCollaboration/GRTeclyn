# What still has to be run on a GPU, and why

> **Status — 2026-08-20: the queue below has been run.**  All six cells
> completed; results are in
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

### What this means for the article

The trustworthy window is roughly **`t ≲ 20–30`**, while both stars are still
intact.  Inside that window, at a separation where they are actually separate
objects, the drift is `≈ 0.01` — consistent with nothing measurable.  Everything
past `t ≈ 45` is a collapsing star.

Three options, in order of cost:

1. **Re-scope the claim.**  Present the result as what it is — an unstable
   positive-mass soliton paired with a stable negative-mass one, where the
   asymmetry in stability is itself the finding, and the one-sidedness of
   gravitational collapse for negative energy is a clean physical statement
   worth making.  No new GPU time.
2. **Re-run on the stable branch.**  `ω ≈ 0.80–0.85` puts the canonical star
   where `dM/dφ_c > 0`.  The stars are `3–5×` lighter and `~15%` larger, so the
   separations all have to grow with them — `d₀ ≈ 25–30` to get the same
   clearance — which means a larger box at the same resolution.  Budget the full
   campaign again, not a cell or two.
3. **Both.**  Keep the present runs as the unstable-branch case and add a
   stable-branch pair as the control the article currently lacks.

Whichever is chosen, the pre-campaign framing of these cells — "the published
campaign is complete and nothing below is needed to reproduce a number in it" —
no longer holds.  The numbers reproduce; their interpretation does not survive.

---
## Before launching anything

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
