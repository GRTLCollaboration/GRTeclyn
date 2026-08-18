# Bondi dipole runaway — publishable results pack

A positive-active-mass and a negative-active-mass soliton, released at rest in
full 3+1 numerical relativity with dynamical, constraint-solved matter,
**self-accelerate in the same direction** while both same-sector controls stay
put. This pack holds every light artefact behind that claim: per-cell time
series, the dressed-star tables, the solve/evolution parameters, the code
patches, curated frames and movies, and the derived article tables.

Heavy artefacts (~120 MB of PNG frames per cell, plotfiles, gridinit, HDF5
metric stacks) stay in the gitignored `runs/` tree on the machine that produced
them. Everything here is small enough to live in git.

## Headline result

Six cells, all solved to Ham ≤ 0.09 % / Mom ≤ 0.08 %, all born on their
predicted t = 0 fingerprints:

| cell | sectors | t_end | drift canon | drift phantom | sep 0 → end | verdict |
|---|---|---|---|---|---|---|
| `single_p` | canonical | 40 | +0.06 | — | — | stable star, drift = noise floor |
| `single_m` | phantom | 40 | — | +0.07 | — | stable star, drift = noise floor |
| **`pair_pm`** | **canonical + phantom** | **60** | **+3.35** | **+9.56** | **8.00 → 1.79** | **runaway** |
| `pair_pp` | canonical ×2 | 60 | +0.16 | — | — | null (attracting merger) |
| `pair_mm` | phantom ×2 | 60 | — | +0.04 | — | null (mutual repulsion) |
| `pair_pm_eqm` | canonical + phantom, \|ADM\| matched | 60 | +2.94 | +9.58 | 8.00 → 1.36 | runaway persists |

The co-drift is exclusive to the mixed cell — 20–80× above either control's
residual wobble and 50–150× above the single-star noise floor. Full numbers,
the quantitative gravity tests, and the honest caveats are in
[`FINDINGS.md`](FINDINGS.md).

## What is where

| path | contents |
|---|---|
| [`FINDINGS.md`](FINDINGS.md) | the physics: results, quantitative tests against a point-mass reference, caveats, open questions |
| [`DEBUGGING.md`](DEBUGGING.md) | how the matter was made to stop dispersing — two root causes, the falsifications that found them |
| [`MATTER_MODEL.md`](MATTER_MODEL.md) | the bicomplex model, where the sign flip lives, dressed-star initial data, code map |
| [`LAUNCH.md`](LAUNCH.md) | exact launch command and full configuration for every cell |
| `data/<cell>/` | per-cell time series + provenance (see the data dictionary below) |
| `campaign/<cell>/` | the error-bar campaign: same streams at three resolutions + a double-size box (see the campaign section below) |
| `stars/` | dressed-star profile tables `r φ₀(r) α(r)` + the M(ω) family scan |
| `analysis/` | derived tables and the scripts that regenerate them |
| `figures/<cell>/` | curated frames: matter activity and geometry at t = 0 → end |
| `movies/<cell>/` | the views that carry the result (matter motion, geometry sign, signed energy density) |
| `patches/` | the matter-model modifications this campaign required |
| `debug_log/` | the full campaign journal, verbatim |

### Data dictionary — `data/<cell>/`

Every `.dat` stream carries its own `#` header line naming each column.

| file | what it is | key columns |
|---|---|---|
| `sector_barycenters.dat` | **the trajectory record** — per-sector integrals | 1 `t`, 2 `total_canon`, 3 `bary_x_canon`, 6 `rms_canon`, 7 `total_phantom`, 8 `bary_x_phantom`, 11 `rms_phantom` |
| `confinement.dat` | core health / grip | 3 `peak_activity`, 5 `confined_frac`, 18 `min_chi` |
| `psi4_mode_l2m0.dat`, `psi4_mode_l2_all.dat`, `psi4_directional.dat` | radiation extraction at R = 8, 16 | — |
| `shell_profiles.dat` | metric on extraction shells (χ, lapse, K) | — |
| `constraint_norms.dat` | **constraint violation during the evolution** (downsampled to Δt = 0.5) | 1 `t`, 2 `L2_Ham`, 3 `L2_Mom`, 7 `L2_Ham_rel`, 8 `L2_Mom_rel` |
| `energy_conditions.dat`, `curvature_invariants.dat` | NEC/WEC monitors and invariants (downsampled) | — |
| `collapse_diagnostics.dat` | **horizon watch** — how deep the gravity well gets and whether a trapped surface appears (downsampled; column names are supplied by the pack script, the raw stream has none) | 1 `t`, 2 `min_lapse`, 3 `min_chi`, 4 `max_abs_K`, 8 `max_ah_r`, 9 `min_theta_plus` |
| `boundary_flux.dat`, `areal_radius.dat`, `ftl_timeseries.dat` | outflow, minimal areal radius, geodesic diagnostics | — |
| `grtresna_params.txt` | the constraint-solve input (couplings, per-lump seeds, tolerances) | — |
| `evolution_params.txt` | the GRTeclyn evolution input (grid, AMR, boundaries, dt) | — |
| `Ham_and_Mom_errors.txt` | per-iteration constraint residuals (last row = final) | — |
| `metadata.json` | provenance: git commit, overrides, solve convergence | — |
| `initial_data.matter.json` | the matter configuration actually seeded | — |

Grid geometry for every cell: `L = 64`, `N = 128`, centre at `32 32 32`, so a
star seeded at `x = ±4` reads as `bary_x = 36 / 28`. Barycentre coordinates in
the streams are absolute; drifts quoted anywhere in this pack are
`x(t) − x(0)`.

### `analysis/`

| file | contents |
|---|---|
| `summary.csv`, `summary.md` | one row per cell: birth checks → final state |
| `trajectories.csv` | drift / separation / core series for every cell, sampled every Δt = 4 |
| `newtonian_reference.csv` | point-mass Bondi integration with the measured ADM masses, aligned to the NR output |
| `convergence_check.csv`, `convergence_check.md` | drift / gap / control artefact across grid resolutions, with the spread between grids |
| `wave_check.csv`, `wave_check.md` | gravitational-wave amplitude on each extraction shell in retarded time, and whether the shells agree |
| `make_tables.py` | regenerates `summary.*` and `trajectories.csv` from `data/` |
| `newtonian_reference.py` | regenerates the point-mass reference (pure stdlib RK4) |
| `convergence_check.py` | regenerates `convergence_check.*` from `campaign/` |
| `wave_check.py` | regenerates `wave_check.*` from `campaign/` |
| `star_family_scan.py` | regenerates `stars/star_family.csv` (needs the wrapper venv) |

## campaign/ — the error-bar campaign (2026-08-17, complete)

The headline numbers above come from adaptive-mesh runs at a single resolution.
To attach error bars, the same physics was rerun on **uniform** grids (no mesh
refinement — the convergence math needs a single cell size everywhere) at three
sharpness levels, plus a double-size box for wave extraction. Every
`campaign/<cell>/` folder carries the same streams and provenance files as
`data/<cell>/` (same data dictionary), plus reference frames under `frames/` —
matter (`scalar_activity_z`) and geometry (`chi_minus_1_z`) on the z slice,
one every Δt = 10 plus the final state, named by simulation time
(`scalar_activity_z_t0030.png` is that field at t = 30). That is 7 per field
for the t = 60 cells and 10 for the t = 90 double-box cells.

### The eleven runs

`pm` = one positive + one phantom star (the runaway pair). `pp` = two positive
stars (control — it must not run away). `eqm` = phantom retuned so both stars
have the same |ADM| mass. All were released from rest, 8 apart.

| cell | box | grid | cell size | matter | what it is for |
|---|---|---|---|---|---|
| `convA_pm_n128` | 64 | 128³ | 0.50 | runaway pair | coarse leg of the resolution ladder |
| `convA_pm_n192` | 64 | 192³ | 0.33 | runaway pair | middle leg |
| `convA_pm_n256` | 64 | 256³ | 0.25 | runaway pair | fine leg — **quoted values come from here** |
| `convA_pp_n128` | 64 | 128³ | 0.50 | control, two positive | coarse leg of the artefact ladder |
| `convA_pp_n192` | 64 | 192³ | 0.33 | control, two positive | middle leg |
| `convA_pp_n256` | 64 | 256³ | 0.25 | control, two positive | fine leg — the residual drift here is the error floor |
| `convA_pm_eqm_n128` | 64 | 128³ | 0.50 | equal-mass pair | does the runaway need a mass difference? coarse leg |
| `convA_pm_eqm_n192` | 64 | 192³ | 0.33 | equal-mass pair | middle leg |
| `convA_pm_eqm_n256` | 64 | 256³ | 0.25 | equal-mass pair | fine leg |
| `boxC_pm_L128_n256` | 128 | 256³ | 0.50 | runaway pair | waves on four shells (r = 16, 24, 32, 40) and the box-size systematic; run to t = 90 |
| `boxC_pp_L128_n256` | 128 | 256³ | 0.50 | control, two positive | the same for the control |

The two `boxC` cells trade resolution for reach: 256³ spread over twice the
width gives cell size 0.50, so compare them against the `n128` column, not
`n256`. Their drifts are measured from the box centre `64 64 64`, not `32 32 32`.

A second phase is running (launched 2026-08-17 evening) and will appear here as
`convA_mm_n128/192/256` — the two-phantom control on the same ladder, which so
far rests on a single adaptive-mesh run — and `convA_pm_sep12_n128/192/256`,
the pair started 12 apart instead of 8, to test whether the runaway survives at
a constant gap. Re-running the pack script picks them up automatically.

### What the campaign shows

Full tables: [`analysis/convergence_check.md`](analysis/convergence_check.md).

The runaway drift is resolution-independent — the three grids agree to 2.2 % at
t = 15 and to 0.5 % at t = 60 — while the control's spurious drift falls by
roughly half at every resolution step. The effect therefore pulls away from the
artefact as the run proceeds: on the finest grid the runaway beats the control
by 8.4× at t = 15, 9.8× at t = 20, 56× at t = 40 and 179× at t = 60. Real
physics survives grid refinement; the artefact does not. (That is also why the
control table shows a huge "spread": for an artefact, disagreement between
grids is the expected, healthy outcome.)

The equal-mass series reproduces the runaway with the mass asymmetry removed:
its drift is 3 % below the unequal-mass pair (+6.227 vs +6.407 at t = 60),
well outside the 0.3–0.5 % grid spread — a resolved physical difference, not
noise. The runaway does not depend on one star outweighing the other.

**No successive-difference order fit is quotable.** The three-grid ratios sit
at 0.2–1.6, below the 1.41 that any positive convergence order requires, so the
error is dominated by initial-data offsets (~1e-3) rather than by a clean
truncation term. Quote the grid-to-grid spread as the error bar, not an order.

**Trust window: t ≲ 50.** Both lumps stay localized (spread ≈ 10, separation
≈ 5) up to there. Past t ≈ 55 their spreads exceed their separation, and in the
double box they cross at t ≈ 65. The spreading is *converged* (the three grids
agree to 2.5 %), so it is a property of this matter model, not grid error — but
past it the per-sector barycentre no longer means "where the star is", so the
apparent slow-down of the drift after t ≈ 65 in `boxC_pm` is a diagnostic
artefact of merged matter, not a physical deceleration.

**No horizon forms during the measurement window — but the pair does collapse
after it.** Every cell that stops at t = 60 stays horizon-free: the expansion
θ₊ never goes negative, no apparent horizon is found, and the lapse only dips
to ≈ 0.55. The runaway is therefore not a collapse artefact. The two double-box
cells, which run to t = 90, keep going and both form a black hole — a trapped
surface appears at **t = 68.5** in the mixed pair (horizon radius ≈ 2.5) and at
**t = 69.5** in the two-positive control (≈ 4.5), with the lapse collapsing to
0.06 and 0.01 by t = 80. So a canonical–phantom contact ends in collapse, not
in mutual annihilation. That end state comes from a single cell size (0.50) and
has not been convergence-tested; treat it as observed, not measured.

### Gravitational waves

Full tables: [`analysis/wave_check.md`](analysis/wave_check.md).

A wave amplitude is only quotable if every extraction shell sees the same
outgoing wave once each is read at its own retarded time `u = t − r`, so that
`r × Ψ₄` is radius-independent. For the runaway pair it is:

| check | result |
|---|---|
| shells agree (4 shells, r = 16 → 40) | within 1.1–1.4× over `u = 10 … 25` |
| across resolutions (n128 / n192 / n256) | **0.3 %** |
| across box sizes (L = 64 vs L = 128) | **0.1 %** |
| equal-mass pair | emits ≈ 6 % less than the unequal pair |

So the l = 2, m = 0 amplitude over `u = 10 … 25` is a converged, box-independent
measurement. Outside that window it is not: before `u ≈ 10` the shells still
carry start-up junk from the initial data, and after `u ≈ 25` the expanding
matter reaches the shells (it crosses r = 16 at t ≈ 51), which is what ends the
window — not the box size.

**The two-positive control has no usable wave window at all**: its shells never
agree on one wave for two consecutive samples, because its matter reaches the
innermost shell at t ≈ 31. Do not quote a control wave amplitude.

**A bigger box is not worth running.** The double box already reproduces the
small box's wave amplitude to 0.1 %, so reach is not the limitation — the
matter expanding onto the shells is, and that happens at the same physical time
in any box. A quadruple-size box at this cell size is also not possible on this
machine: 512³ needs roughly 460 GB against 81 GB per GPU, and multi-rank runs
are broken here. The only thing more box would buy is the post-break-up phase,
which is not a claim this pack makes.

### Reading notes

- These are uniform-grid runs; the AMR twin of the mixed pair is `data/pair_pm`.
- `boxC` drifts are `bary_x − 64`, not `− 32` as everywhere else in this pack.
- Re-running the pack script rebuilds every cell folder from scratch; nothing
  outside `campaign/` and the `analysis/*_check.*` tables is touched.

## Reproducing

Regenerate this whole pack from the run tree:

```bash
bash research/bondi_dipole/pack_results.sh    # published cells -> data/ stars/ figures/ movies/
bash research/bondi_dipole/pack_campaign.sh   # campaign cells  -> campaign/ + convergence tables
```

Both copy only small artefacts and scrub absolute machine paths at runtime
(`research/bondi_dipole/scrub_paths.py`) — no host, user, or site identity
enters git. To re-run the physics itself, see [`LAUNCH.md`](LAUNCH.md).

## Provenance

- Initial data: **GRTresna** (CTTKHybrid, `BosonStarBH` example, complex scalar
  matter with per-lump signs), single rank.
- Evolution: **GRTeclyn** (`RadialRecipe`, CCZ4 + bicomplex scalar matter),
  single rank, one GPU per cell.
- Campaign journal: [`debug_log/bondi_dipole_debug.md`](debug_log/bondi_dipole_debug.md)
  (working copy: `research/bondi_dipole_debug.md`).
- Code state: wrapper commit recorded per cell in `data/<cell>/metadata.json`;
  the GRTresna matter modifications are in [`patches/`](patches/).
