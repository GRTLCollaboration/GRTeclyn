# What still has to be run on a GPU, and why

The published campaign is complete: the article stands on the runs already in
`results/bondi-dipole-runaway/campaign/`, and nothing below is needed to
reproduce a number in it.  This is the **queue for the next campaign** — the
open items that cannot be answered by re-analysing what is on disk, why each
one needs the device, and what it unblocks.

Read [`../../../results/bondi-dipole-runaway/LAUNCH.md`](../../../results/bondi-dipole-runaway/LAUNCH.md)
first for the launch mechanics: the `/usr/bin/env` shadowing trap, how to
verify a detached launch actually started, and `stop_campaign.sh`.  Costs quoted
here are the measured ones from Appendix B of the article.

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
