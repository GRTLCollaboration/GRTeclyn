# Bondi Dipole Runaway — Debug Log & Road Forward

*Status as of 2026-08-05. Covers campaigns `bondi_dipole_v1` (strong-field),
`bondi_dipole_weakfield_v1`, `bondi_dipole_weakfield_v2`, and
`bondi_dipole_midfield_v1` — all stopped. Goal unchanged: a clean,
non-dispersing (+,−) runaway with dense frames for smooth movies (N1,
Appendix B of `research/nextsteps.md`).*

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
| stop tool (the only sanctioned way) | `grteclyn-wrapper/scripts/campaigns/stop_campaign.sh` |
| run dirs | `runs/bondi_dipole_v1`, `runs/bondi_dipole_weakfield_v1`, `runs/bondi_dipole_weakfield_v2`, `runs/bondi_dipole_midfield_v1` |
| launch logs | `runs/bondi_wf2_launch.log`, `runs/bondi_midfield_launch.log` |
| trajectory stream | `<run>/small_data/sector_barycenters.dat` (cols: t, total_canon, bary_xyz_canon, rms_canon, total_phant, bary_xyz_phant, rms_phant) |
| grip / collapse monitor | `<run>/small_data/confinement.dat`, min_chi = col 18 |
| v2 intermediate movies | `runs/bondi_dipole_weakfield_v2/bondi_wf2_single_p/movies/` (19 views, t=0–13) |
| v1 pair movie (overlap-hypothesis source) | `runs/bondi_dipole_weakfield_v1/bondi_wf_pair_pm/movies/movie_scalar_activity_z_t0-32.mp4` |
