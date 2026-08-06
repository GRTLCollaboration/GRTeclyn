# Bondi Dipole Runaway — Debug Log & Road Forward

*Status as of 2026-08-05. Covers campaigns `bondi_dipole_v1` (strong-field),
`bondi_dipole_weakfield_v1`, `bondi_dipole_weakfield_v2`, and
`bondi_dipole_midfield_v1` — all stopped. Goal unchanged: a clean,
non-dispersing (+,−) runaway with dense frames for smooth movies (N1,
Appendix B of `research/nextsteps.md`).*

> **UPDATE 2026-08-06 — root cause found; §4 diagnosis and the §5 Step-1 scan
> are superseded by §7.** The blob was never wrong-shaped: it is the exact
> Q-ball eigenstate at 95.45% of its amplitude, because the painter clamps the
> amplitude to the thin-wall estimate. One-line opt-in fix landed
> (`grtresna_qball_exact_amplitude=1`). The proposed well-depth scan would
> have been three identical runs (the knob saturates); do not run it.
>
> **UPDATE 2026-08-06 (later) — §7.3 step 1 ran and FAILED; superseded by §8.**
> The exact-amplitude flat ball still blows off, because at the weak rung the
> gravitationally dressed equilibrium near the target ball DOES NOT EXIST.
> Fixed-frequency dressed-star initial data landed; first rung where the
> dressed ω=0.55 star exists is λ=10240 ("ultraweak"). Dressed single-lump
> calibration running as `bondi_dipole_selfgrav_v1`.

## TL;DR

- **Four campaigns, four failures — and the failure modes triangulate the
  root cause.** Heavy lumps collapse, light lumps evaporate, and the
  geometric-mean rung does *both at once* (core crunch + envelope blowoff).
  No amplitude rung can fix that: the constraint-solved initial blob is far
  from any self-holding soliton shape at these couplings, so it violently
  re-arranges no matter its weight.
- **The two "shape" knobs we reached for first are dead ends**:
  `grtresna_bs_phi_c` and `grtresna_bs_profile_width` are only the solver's
  initial *guesses*. Changing them (v1 → v2) left the t=0 sector integrals
  identical to 9 significant digits. The solver relaxes to its own state.
- **The overlap hypothesis is refuted.** A lone lump with nothing to overlap
  dispersed on exactly the pair's curve (v2 `single_p`).
- **Next move: a cheap calibration scan over the dials that actually change
  the solver output** — binding frequency `grtresna_bs_omega` and the seed
  amplitude `trajectory_well_depth` — before burning any more pair cells.
  Details in "Road forward" below.
- Dispersal is **not** an exotic-sector disease: every failed calibration
  lump was canonical, and in the one side-by-side test (v1 `pair_pm`) the
  canonical lump dispersed *faster* than the phantom.

## 1. The experiment

Bondi (1957): a positive/negative-active-mass pair self-accelerates — the
phantom falls toward the canonical lump's well while the canonical lump
rolls off the phantom's hill, so both accelerate the same way (phantom
chasing canonical) with P_ADM ≈ 0. Never evolved in full 3+1 NR with
dynamical, constraint-solved matter. The bicomplex model supplies exactly
the required sign structure: both sectors obey the same Klein-Gordon
equation (positive inertial/passive mass); the sign flip lives only in the
Einstein source (negative active mass).

The control matrix — each cell a falsifiable prediction:

| cell | prediction |
|---|---|
| `pair_pm` (+,−) | runs away along +x, phantom chasing canonical |
| `pair_pp` (+,+) | attracts — merges or orbits, no net drift |
| `pair_mm` (−,−) | mutually repels |
| `single_p` / `single_m` | drift nowhere; calibrate per-sector dispersal |

Setup: two lumps at rest on the x axis, evolution pump off from t=0
(`rl_pump_stop_time=0`), no boosts/breathing/tilt/rotation — pure
self-gravity. Live per-sector barycentres
(`GRTECLYN_SECTOR_BARYCENTERS=1` → `small_data/sector_barycenters.dat`),
Ψ4 on. Plotfiles are consumed and purged, so the barycentre stream is the
only recoverable trajectory record.

**Precondition for any pair claim:** a lone lump must first hold its size.
That gate is what every campaign so far has failed.

## 2. Campaign history

| campaign | couplings (λ, μ) | ω | core amp √(3λ/4μ) | outcome |
|---|---|---|---|---|
| `bondi_dipole_v1` (strong) | 640, 85333 | 0.8 | 0.075 | **collapse** — min_chi → 0.009 by t≈25–30 |
| `bondi_dipole_weakfield_v1` | 2560, 1365333 | 0.55 | 0.0375 | **dispersal** — rms 5 → ~30 by t≈50 (pairs `pm`, `pp` ran to t=100; singles never ran) |
| `bondi_dipole_weakfield_v2` (corrected seeds) | 2560, 1365333 | 0.55 | 0.0375 | **dispersal anyway** — lone `single_p` on v1's exact curve; stopped after calibration cell |
| `bondi_dipole_midfield_v1` | 1280, 341333 | 0.55 | 0.053 | **split**: core collapse *and* envelope blowoff; stopped at t≈34 of `single_p` |

All rungs keep λ²/μ = 4.8 so ω_min = 0.316 is fixed; ω = 0.55 means 45%
binding throughout. Separation 8 (R0 = 4) from the weak-field ladder on.
All runs sequential on one GPU, single-rank (mpirun segfaults on this
node — `RANKS=1` always).

### 2.1 Key data fingerprints

**v2 `single_p` (weak rung, corrected seeds)** — breath in, then evaporate:

```
rms: 5.22@t0 → 4.98@t4 → 4.65@t8 (min) → 5.44@t12 → 6.32@t16 → 7.46@t20
min_chi: settles ~0.88, 0.838@t20            (grip stays feeble throughout)
total canonical activity climbs 32.7 → 41.3  (field spreading into the box)
```
Matches v1's dispersal curve (8.7@t19). The corrected seed knobs changed
nothing measurable.

**mid-field `single_p`** — deeper breath, then both failures at once:

```
rms:  5.221@t0 → 4.467@t8.4 (min) → 5.54@t15.6 → 6.41@t20 → 8.19@t24
      → 9.58@t26 → 11.1@t28 → 16.4@t34          (accelerating blowoff)
min_chi: 0.96@t6.8 → 0.715@t9.2 (AMR level-1 regrid resolves the true
      minimum) → plateau ~0.73 to t≈20 → 0.59@t26 → 0.51@t28 → 0.36@t30
      → 0.22@t32                                 (through the 0.3 line)
```
Read: a small central knot self-compresses toward collapse while the bulk
of the field is flung off. The heavier rung gave gravity a real grip
(0.73 vs v2's 0.88) — and the grip went to the wrong place.

**v1 `pair_pm` sector comparison** (the only completed mixed pair):

```
t:          0     11.2   22.4   33.6   44.8   100
rms canon:  5.22  5.88   10.5   20.0   31.6   30.7
rms phant:  5.22  5.94   9.6    15.7   23.7   32.1
```
Both sectors disperse, but the *canonical* lump spreads faster — consistent
with the Bondi sign structure even in dispersal (phantom is pulled toward
its partner and gathered up; canonical is pushed off the phantom's hill and
smeared out).

## 3. Issues found (chronological), and what each taught us

1. **Stopping detached campaigns by killing workers doesn't work** (strong
   field, 2026-08-05). Killing the worker just advances the orchestrator to
   the next cell; `$!` after `setsid` is the dead parent; GRTresna solvers
   sit in their own session. Fix (permanent infrastructure):
   `scripts/campaigns/stop_campaign.sh` kills orchestrators/ancestors
   first, then sweeps workers, then pgrep-verifies with SIGKILL
   escalation. Every launcher now registers via
   `campaign_register_launcher` (`launcher.pid`). One tool, one way —
   deliberately no per-campaign stop scripts.

2. **v1 seed bug (real, but not the root cause).** v1 never overrode the
   GRTresna bound-state seed, so `bs_profile_width` fell back to 8.0 (as
   wide as the whole pair separation) and `bs_phi_c` to 0.08 (2.1× the
   thin-wall plateau amplitude). Also only 64 frames/run
   (plot-interval 160).

3. **The seed knobs are cosmetic — the falsification that mattered.**
   v2 set `bs_phi_c=0.0375`, `bs_profile_width=2.0` (the "correct" Q-ball
   values) and produced t=0 canonical-sector integrals identical to v1's
   to 9 significant digits (total 34.6488659648 vs …602; rms 5.2208928675
   vs …591). The constraint solver relaxes to its own preferred state;
   those two knobs only seed its iteration. **Any future shape control
   must go through dials that survive the solve** (couplings, ω,
   `trajectory_well_depth` — the well amplitude *is* the matter source the
   solver sees; `trajectory_well_width` shapes only the pump wells, which
   are off).

4. **The overlap hypothesis is refuted.** Hypothesis (from watching v1's
   `pair_pm` movie): the two lumps'/fields' overlap causes the failure.
   Test: v2 ran `single_p` first — one lump, nothing to overlap — and it
   dispersed on v1's exact curve. Overlap and seed shape were never the
   cause.

5. **The mass ladder is falsified as a one-parameter fix.** Strong rung
   collapses, weak rung evaporates, geometric-mean rung *splits* (core
   crunch + envelope blowoff simultaneously). Conclusion: the solver's
   initial blob is not near any stable-soliton shape at these couplings —
   it sheds/re-arranges violently at every weight. Amplitude is the wrong
   lever; the *profile* the solver converges to is the problem.

6. **Environment constraints (permanent):** mpirun segfaults cluster-wide →
   GRTresna and the evolution both run single-rank, always. Campaigns
   launch detached (`setsid nohup … </dev/null &`); never edit a script
   bash is currently executing; the depth-search campaign owns one slot on
   each of GPUs 0–3, so Bondi runs sequentially on GPU 1's spare capacity.

7. **Diagnostics that earned their keep:** the live per-sector barycentre
   stream (rms_radius per sector = col 6/11 of `sector_barycenters.dat`)
   and min_chi (col 18 of `confinement.dat`) called every verdict within
   ~20 min of launch; dense frames (plot-interval 40 → 0.4 t-units/frame,
   ~250 frames/run) make mid-run movie snapshots practical via symlinked
   frame subsets (movies built this way for v2 `single_p` t=0–13 without
   touching run output).

## 4. Diagnosis

The initial data is a constraint-solved blob whose profile is fixed by the
solver, not by us, and at ω=0.55 / λ²/μ=4.8 that profile is far from the
true Q-ball ground state. The evolution then does what mismatched initial
data always does: it re-arranges on a dynamical timescale — ejecting the
mismatch as an outgoing shell and, if the remainder is heavy enough,
over-compressing the core. Self-gravity only tips which failure dominates:

- amp 0.0375 → ejection wins entirely (evaporation),
- amp 0.075 → compression wins entirely (collapse),
- amp 0.053 → both, split by radius.

Gravity is a small correction on top of the scalar self-interaction; the
holding-together job belongs to the Q-ball binding, which the solver's
blob does not realize. **Fix the starting shape, not the weight.**

## 5. Road forward

### Step 1 — calibration scan (recommended, ~an afternoon of GPU-1 time)

Small scan of short `single_p` runs at fixed mid-rung couplings
(λ=1280, μ=341333), varying only the dials that demonstrably change the
solver output:

- `grtresna_bs_omega` ∈ {0.45, 0.55, 0.65, 0.75} — how much internal
  binding the bound state carries (ω_min = 0.316, ω_max = 1);
- `trajectory_well_depth` ∈ {0.08, 0.15, 0.25} — the seed-well amplitude
  the solver actually sources from.

Cheap settings: stop-time 40, plot-interval 80, same grid. Pass gate per
cell: rms_radius flat (±10%) to t=40 **and** min_chi stays above 0.3 with
a stable plateau. The verdict is readable from `sector_barycenters.dat` /
`confinement.dat` within ~20–30 min per cell; a failing cell can be cut
early. Winner becomes the recipe for the full matrix.

### Step 2 — the real matrix at the winning recipe

Five cells in v2/mid-field order (`single_p` first as live calibration,
then `pair_pm`, `pair_pp`, `pair_mm`, `single_m`), stop-time 100,
plot-interval 40 (~250 frames/run) for smooth movies. Deliverables:

- (+,−) drift curve (both sector barycentres vs t) + chase movie;
- controls: (+,+) attracts, (−,−) repels, singles static — the falsifiable
  quartet that seals the runaway claim;
- Ψ4 record: what an accelerating dipole radiates is an open question;
- lone exotic lump (`single_m`) — *never yet run in any campaign*; also
  closes the "do exotic lumps evaporate differently?" question properly.

### Step 3 — if the scan finds no holding cell

Options, in order of preference:

1. **True Q-ball profile initial data**: solve the 1D Q-ball ODE for
   (λ, μ, ω) offline and feed the actual profile to the solver instead of
   its generic ansatz — principled, but needs a solver-side hook for a
   custom radial profile.
2. Keep the pump on at low amplitude as an external "holding hand"
   (accepting that the runaway must then be demonstrated *against* a
   documented, sector-symmetric external potential).
3. Shorter-baseline claim: run the pair during the transient window
   (t ≲ 8, while lumps are still compact) with high frame density and
   measure differential drift before dispersal — weaker, but honest.

### Housekeeping

- Sync `campaign_register_launcher` into `qball_trajectory/run.sh` (the
  one unwired launcher) when the depth campaign is idle.
- Proposed, not yet requested: "surfer test" — a geodesic tracer at the
  pair gap midpoint riding the runaway.

## 6. File inventory

| what | where |
|---|---|
| matrix scripts (strong → mid) | `grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_matrix.sh`, `run_matrix_weakfield.sh`, `run_matrix_weakfield2.sh`, `run_matrix_midfield.sh` |
| exact-amplitude fix (2026-08-06) | `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/config.py` (`grtresna_qball_exact_amplitude`), tests in `tests/grtresna/test_qball_bicomplex_campaign.py` |
| stop tool (the only sanctioned way) | `grteclyn-wrapper/scripts/campaigns/stop_campaign.sh` |
| run dirs | `runs/bondi_dipole_v1`, `runs/bondi_dipole_weakfield_v1`, `runs/bondi_dipole_weakfield_v2`, `runs/bondi_dipole_midfield_v1` |
| launch logs | `runs/bondi_wf2_launch.log`, `runs/bondi_midfield_launch.log` |
| trajectory stream | `<run>/small_data/sector_barycenters.dat` (cols: t, total_canon, bary_xyz_canon, rms_canon, total_phant, bary_xyz_phant, rms_phant) |
| grip / collapse monitor | `<run>/small_data/confinement.dat`, min_chi = col 18 |
| v2 intermediate movies | `runs/bondi_dipole_weakfield_v2/bondi_wf2_single_p/movies/` (19 views, t=0–13) |
| v1 pair movie (overlap-hypothesis source) | `runs/bondi_dipole_weakfield_v1/bondi_wf_pair_pm/movies/movie_scalar_activity_z_t0-32.mp4` |

## 7. Root cause (2026-08-06) — the amplitude clamp, not the profile

Three parallel code traces + an independent re-solve of the 1D Q-ball
equation + digit-level reconstruction of the recorded t=0 diagnostics.
Every claim below was verified against the actual run artifacts.

### 7.1 What actually happens at t=0

1. **The painted shape was always the TRUE Q-ball.** The campaign path
   (`grtresna_qball_ode_profile=1` → lump profile 3) shoots the exact radial
   eigenstate and paints it from `grtresna/qball_profile.dat`. An independent
   re-solve reproduces the table's central amplitude to 6 digits, and the
   painted table reproduces the recorded t=0 stream *exactly*
   (predicted total 49.0006 / rms 5.2209 vs measured 49.0006 / 5.2209,
   mid-field `single_p`). The blob was never "far from any soliton shape" —
   §4 is wrong on that point.
2. **The elliptic solver never touches the matter.** Matter is painted once
   before the nonlinear loop; the solve adjusts only the metric. "The solver
   relaxes to its own state" is not physically possible in this code.
3. **Why v1 = v2 to 9 digits:** `grtresna_bs_phi_c` and
   `grtresna_bs_profile_width` are *dead knobs* on the lump path — written to
   params.txt, read by nobody (Python config takes the amp from
   `cap_well_depth`, the C++ legacy branch that would read them is gated on
   having zero lumps). The two campaigns painted byte-identical profile
   tables (same md5). Cosmetic knobs, not "solver relaxation".
4. **The one real defect:** the amplitude is clamped to the thin-wall
   estimate √(3λ/4μ), but the true eigenstate centre is 4.55% higher
   (ratio 0.95452 at every rung, since all rungs share λ²/μ). Both painters
   rescale the exact soliton by amp/φ_c, so every campaign started a
   **0.9545 × eigenstate** — an off-shell state that must breathe.
5. **Why the three rungs failed three ways:** the breathing ball's fate is
   set by self-gravity, which the flat-space profile ignores entirely.
   Compactness 2E/R of the ω=0.55 ball: weak 0.14 (bounce unbinds →
   evaporation), mid 0.27 (core crunch + envelope blowoff), strong 0.55
   (deep infall → collapse). Matches §2 outcomes cell by cell, including the
   first observed breath *inward* (rms 5.22 → 4.65 by t≈8: gravity switching
   on over a zero-gravity profile).

Also falsified en route: `trajectory_well_depth` saturates at the same
√(3λ/4μ) cap — the §5 Step-1 scan over {0.08, 0.15, 0.25} would have been
**three bit-identical runs**. And ω=0.45 is a bad scan cell regardless: the
true ball there has rms ≈ 11, too big for separation 8 in this box.

### 7.2 The fix (landed)

`grtresna_qball_exact_amplitude=1` (opt-in override, default off so every
existing campaign keeps bit-identical seeds): lump amp := the ODE table's own
φ_c, making the painter's rescale exactly 1 — the seed *is* the stationary
eigenstate. Regression tests cover both flag states. For the next matrix,
add to `common_overrides`:

```
--extra-override grtresna_qball_exact_amplitude=1
```

(the `grtresna_bs_phi_c` / `grtresna_bs_profile_width` lines can be deleted —
dead knobs). **Launch-time verification** (first row of
`sector_barycenters.dat`): total_canon must read ≈ 36.30 (weak rung) or
≈ 51.34 (mid rung) with rms still 5.221. If it reads 34.65 / 49.00 the flag
did not reach the solver.

### 7.3 Corrected road forward

1. **Re-run `single_p` at the WEAK rung with the exact amplitude** (ω=0.55,
   stop 40, same pass gate). Weak rung first *by design*: lowest compactness
   (0.14) = smallest error from the remaining flat-space-profile
   approximation. The strong rung can likely never hold a flat-space ball
   (compactness 0.55); treat it as out of scope for this matrix.
2. If weak `single_p` holds → full five-cell matrix at the weak rung
   (§5 Step 2 unchanged otherwise).
3. If it *still* breathes beyond gate: the remaining defect is the gravity
   dressing. Escalation: the self-gravitating ODE solver
   (`grtresna/profiles/boson_star_ode.py`) already handles the sextic
   potential and emits a 3-column table (φ₀, lapse) that profile 3 accepts —
   but its mode currently (a) forces exotic lumps to canonical (hard veto in
   config), (b) shares one table + one ω across lumps, so a Bondi pair needs
   small wiring work first. Phantom side would also need a sign-flipped
   variant of the metric ODEs (a phantom lump's self-gravity is *repulsive*;
   equilibrium should still exist at low compactness as a perturbed Q-ball,
   but it is not the gravity-bound star branch).
4. **Do NOT switch to plain boson stars for this experiment**: a
   gravity-bound phantom star cannot exist (its own gravity pushes it
   apart) — the code vetoes it for that reason. Q-balls bind through the
   scalar interaction, which is sector-symmetric: the right matter for a
   (+,−) pair. (Unchanged conclusion: dispersal was never an exotic-sector
   disease — the lone canonical lump failed alone, §3.4.)

Dead-knob inventory and full plumbing trace (10 silent-drop points, incl.
the sech fallback when a profile table is missing, and the Python-vs-C++
profile-enum collision 4↔4) — kept out of this journal; ask for the
2026-08-06 trace if needed.

## 8. Exact amplitude is not enough: no dressed equilibrium at the weak rung (2026-08-06, later)

### 8.1 The §7.3 gate run — verified seed, same blow-off

`bondi_ex_single_p` (runs/bondi_dipole_exact_v1, weak rung, ω=0.55,
`grtresna_qball_exact_amplitude=1`, stop 40, GPU 1):

- **Seed verified end-to-end**: `lump0_amp = 0.03928682593468194` (the ODE
  φ_c, not the 0.0375 clamp); t=0 stream read total 36.2998 / rms 5.2209 vs
  the §7.2 prediction 36.300 / 5.221. The fix works exactly as designed.
- **Physics failed anyway**:

```
rms: 5.221@t0 → 4.714@t7.2 (min, −9.7%) → 5.304@t11.2 → 5.604@t12.8
     → 6.167@t16.8 → 6.610@t20        (gate ±10% broken at t≈12.5; stopped)
```

  Differences from the clamped runs: the first breath is a *gravitational
  contraction* (total ∫|φ| fell 36.3 → 29.3 while the core brightened ~2×,
  i.e. the field concentrated instead of leaking), the bottom is shallower,
  and the blow-off is ~10% slower — but it is the same terminal blow-off.

### 8.2 Why: the dressed heavy star does not exist at these couplings

Mapping the self-gravitating (dressed) star family at the weak rung with the
ODE solver (`cached_selfgrav_profile`, ω-shooting at fixed φ_c):

```
phi_c    omega    ADM     alpha(0)
0.0385   0.786    0.045   0.967     <- feather-light stable branch
0.0393   0.753    0.053   0.963
0.0399   0.311    2.50    0.50      <- jump: ultra-compact (unstable) branch
0.0430   0.364    2.22    0.375
0.0500   0.369    2.30    0.207
```

There is **nothing near the flat ω=0.55 ball (E=0.281, compactness
2E/R≈0.14)**. The exact-amplitude seed contracted under gravity, found no
nearby equilibrium, and unbound — consistent with every observation in §8.1.
This also retro-explains the whole §2 ladder: compactness ∝ 1/λ, so heavier
rungs were further from any dressed equilibrium, not closer.

Second finding: the ω-shooting parameterization *cannot* reach heavy-branch
sextic stars even where they exist — the heavy ball sits with φ_c a fraction
of a percent below the effective-potential top, making its eigenvalue an
exponentially thin needle in ω (requesting the ω=0.55 star at φ_c=0.0393
silently returns the 5× lighter ω=0.753 star; verified by scanning the
over/under classification: a single fat transition at ω_int≈0.78).

### 8.3 Code landed (all tested; 22+2 pass)

1. **Evolution repaint for dressed lumps** — `profiles/envelope.py`
   `phi0_at_radius` and `fields/boson_star.py _raw_lump_phi_grid` had no
   branch for `PROFILE_SELFGRAV_BOUND`: a dressed lump fell through to the
   **Gaussian** envelope, so the evolution's t=0 matter would have
   contradicted the constraint solve's tabulated star (likely the residual
   "slow leak" in SELFGRAV_HANDOFF's own result). Both now resolve the same
   cached star table.
2. **Fixed-frequency dressed solve** — `boson_star_ode.py
   solve_selfgrav_at_omega` / `cached_selfgrav_at_omega`: bisect φ_c at fixed
   integration-frame ω (mirroring the proven flat-space parameterization),
   outer-iterate ω_int so ω_phys = ω_int/α_inf hits the request. Note the
   dressed twist: α grows outward so a start exactly at the flat-space
   barrier top rolls off the *outward* side — the "over" bracket edge sits
   below f_top and is found by scanning down.
3. **Params writer** — `solver/params.py _write_selfgrav_profile`: when a
   selfgrav lump carries `qball_omega>0` and sextic couplings, solve at that
   frequency; φ_c becomes an OUTPUT and is back-written onto `lump.amp` (the
   painters rescale the table by amp/φ_c — only amp==φ_c keeps the seed on
   the eigenstate; same bug class as §7's clamp).
4. Regression tests: Gaussian-fallback repaint, at-omega back-write +
   repaint consistency, exact-amplitude flag (both states).

### 8.4 The ultraweak rung — first rung where the dressed star exists

Keeping λ²/μ = 4.8 (ω_min = 0.316 unchanged), ω = 0.55:

| λ | μ | dressed ω=0.55 star? | φ_c | ADM | α(0) |
|---|---|---|---|---|---|
| 2560 (weak) | 1365333 | **no** | — | — | — |
| 5120 | 5461333 | **no** | — | — | — |
| **10240 (ultraweak)** | **21845333** | **yes** | 0.019695 | 0.0640 | 0.977 |

Cost of the rung: the lump is ~4.4× lighter than the weak-rung flat ball, so
the pair drift is gentler — a ≈ M/sep² ≈ 1e-3, displacement ≈ 5 length units
by t=100. Still far above the barycentre stream's resolution.

### 8.5 Running now + revised road forward

- `bondi_sg_single_p` (runs/bondi_dipole_selfgrav_v1, GPU 1, stop 40,
  `run_single_selfgrav.sh`): ultraweak rung, `grtresna_bs_selfgrav=1`,
  ω=0.55. Seed verified at launch: solved eigenvalue 0.5500068,
  lump0_amp = 0.0196947, 3-column (lapse) table. Expected t=0 fingerprint
  total ≈ 15.92 / rms ≈ 5.05; pass gate unchanged (rms ±10% to t=40,
  min_chi > 0.3 plateau). Expect ~1% ripple, not ±10–18% breathing.
- If it holds → pair matrix needs the §7.3-step-3 wiring first:
  per-lump profile tables + per-lump ω (Python emitter + C++
  `ComplexScalarField.cpp:169` global-ω for non-winding lumps), and the
  phantom-dressed profile (sign-flipped gravity in `_ode_rhs` — the phantom
  star's own gravity is repulsive, equilibrium slightly puffed; lift the
  exotic→canonical veto in `config.py:472-482` for the dressed path).
- The `grtresna_qball_exact_amplitude` flag stays: it is the correct (and
  necessary) fix for any FLAT-profile path; it is just not sufficient where
  no dressed equilibrium exists.

### 8.6 Dressed-star calibration v1 verdict (2026-08-06, evening)

`bondi_sg_single_p` (runs/bondi_dipole_selfgrav_v1) ran to t=40. Verdict:
**bound breather, dirty launch — the first lump that did NOT disperse.**

- Seed verified digit-perfect: t=0 total 15.924 / rms 5.045 vs predicted
  15.92 / 5.05; solved eigenvalue 0.5500068; amp = the star's own φ_c.
- **The core survived and breathed**: peak amplitude oscillated with period
  ≈14 (crests 0.0254 / 0.0247 / 0.0247 at t≈6/20/36 — cycles 2-3 undamped),
  confined fraction swung around ~0.5 in rhythm. No dispersal: the flat
  seeds' cores unravelled by this point; this core is intact at t=40.
- **The rms stream is NOT a star diagnostic here**: it climbed 5.0 → 13.2
  monotonically while the core breathed — it tracks the shed-radiation
  bath (amplitude-linear weight × r² leverage). Gate needs rethinking for
  dressed runs: use peak + confined_frac + min_chi, not rms.
- **Launch was dirty**: GRTresna exited at Mom residual 0.64% (exit tol
  1.0%); the residual radiated as a visible χ>1 metric ring (frames at
  t≈18) that sloshed the star's envelope. Fix landed: replay_eval now
  exposes NL exit/stall tolerances; both selfgrav scripts use 0.1 / 0.002.
- **Watch item**: min_chi deepened monotonically 0.99 → 0.86 over t=40,
  accelerating late — the massive-field radiation bath cannot exit through
  massless-wave boundaries and keeps washing over the star (plus possible
  gauge drift). If it persists in the clean run, options: sponge layer, or
  accept and keep pair runs ≤ t≈60.
- v2 (tight solve, runs/bondi_dipole_selfgrav_v2) launched; if the breath
  shrinks and chi stabilises → five-cell matrix
  (`run_matrix_selfgrav.sh`, runs/bondi_dipole_selfgrav_matrix_v1).
