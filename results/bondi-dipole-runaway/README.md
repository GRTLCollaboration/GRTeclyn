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
[`FINDINGS.md`](../../research/bondi_dipole/docs/FINDINGS.md).

> **Read the drift magnitudes with the separation in mind.** Every cell in the
> table above starts the stars 8 apart, and each star is about 5 across, so
> their envelopes overlap from t = 0. That inflates the drift: against a
> point-mass calculation of the same configuration it is 2.14× too large at
> separation 8, falling to 1.36× at 12 and 1.17× at 16. The runaway is real
> and gravitational — the excess shrinks exactly as the stars are pulled
> apart — but quote separation 12–16 for magnitudes. See the campaign section
> and [`FINDINGS.md`](../../research/bondi_dipole/docs/FINDINGS.md) §5.5.

## What is where

| path | contents |
|---|---|
| [`MATTER_MODEL.md`](MATTER_MODEL.md) | the bicomplex model, where the sign flip lives, dressed-star initial data, code map |
| [`LAUNCH.md`](LAUNCH.md) | exact launch command and full configuration for every cell |
| `campaign/<cell>/` | every run: the `convA_*`/`boxC_*` convergence campaign at three resolutions plus a double-size box, and the four unprefixed `pair_*`/`single_*` production cells at Δx = 0.5 (see the data dictionary below) |
| `stars/` | dressed-star profile tables `r φ₀(r) α(r)` + the M(ω) family scan |
| `analysis/` | derived tables and the scripts that regenerate them |
| `movies/<cell>/` | the views that carry the result (matter motion, geometry sign, signed energy density). Folders named `pair_*` / `single_*` are the original adaptive-mesh cells; folders named `convA_*` are the campaign cells, including the separation series |
| `patches/` | the matter-model modifications this campaign required |

### Data dictionary — `campaign/<cell>/`

Every `.dat` stream carries its own `#` header line naming each column.

| file | what it is | key columns |
|---|---|---|
| `sector_barycenters.dat` | **the trajectory record** — per-sector integrals | 1 `t`, 2 `total_canon`, 3 `bary_x_canon`, 6 `rms_canon`, 7 `total_phantom`, 8 `bary_x_phantom`, 11 `rms_phantom` |
| `confinement.dat` | core health / grip | 3 `peak_activity`, 5 `confined_frac`, 18 `min_chi` |
| `psi4_mode_l2m0.dat`, `psi4_mode_l2_all.dat`, `psi4_directional.dat` | radiation extraction at R = 8, 16 | — |
| `shell_profiles.dat` | metric on extraction shells (χ, lapse, K) | — |
| `constraint_norms.dat` | **constraint violation during the evolution** (downsampled to Δt = 0.5; column names are supplied by the pack script, the raw stream has none) | 1 `t`, 2 `L2_Ham`, 3 `L2_Mom`, 16 `Linf_Ham_amr` (worst single point). Columns 7–8 are relative measures — see the caution in the constraint section |
| `energy_conditions.dat`, `curvature_invariants.dat` | NEC/WEC monitors and invariants (downsampled; column names supplied by the pack script) | — |
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
| `separation_scaling.csv`, `separation_scaling.py` | the separation series: measured drift against a point-mass integration of the *same* configuration, at separations 8 / 12 / 16 — the test that isolates the overlap |
| `convergence_check.csv`, `convergence_check.md` | drift / gap / control artefact across grid resolutions, with the spread between grids |
| `wave_check.csv`, `wave_check.md` | gravitational-wave amplitude on each extraction shell in retarded time, and whether the shells agree |
| `constraint_check.csv`, `constraint_check.md` | constraint violation at both stages — the initial-data solve and the evolution — plus whether refining the grid reduces it |
| `make_tables.py` | regenerates `summary.*` and `trajectories.csv` for the Δx = 0.5 production cells |
| `newtonian_reference.py` | regenerates the point-mass reference (pure stdlib RK4) |
| `convergence_check.py` | regenerates `convergence_check.*` from `campaign/` |
| `wave_check.py` | regenerates `wave_check.*` from `campaign/` |
| `constraint_check.py` | regenerates `constraint_check.*` from `campaign/` |
| `star_family_scan.py` | regenerates `stars/star_family.csv` (needs the wrapper venv) |

## campaign/ — the error-bar campaign (2026-08-17 → 18, complete)

The headline numbers above come from adaptive-mesh runs at a single resolution.
To attach error bars, the same physics was rerun on **uniform** grids (no mesh
refinement — the convergence math needs a single cell size everywhere) at three
sharpness levels, plus a double-size box for wave extraction. Every
`convA_*`/`boxC_*` folder carries the same streams and provenance files as
the production cells (same data dictionary), plus reference frames under `frames/` —
matter (`scalar_activity_z`) and geometry (`chi_minus_1_z`) on the z slice,
one every Δt = 10 plus the final state, named by simulation time
(`scalar_activity_z_t0030.png` is that field at t = 30). That is 7 per field
for the t = 60 cells and 10 for the t = 90 double-box cells.

### The twenty runs

`pm` = one positive + one phantom star (the runaway pair). `pp` = two positive
stars (control — it must not run away). `mm` = two phantom stars (the other
control). `eqm` = phantom retuned so both stars have the same |ADM| mass.
`sep12` / `sep16` start the pair 12 and 16 apart instead of 8. All were
released from rest.

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
| `convA_mm_n128/192/256` | 64 | 128/192/256³ | 0.50 / 0.33 / 0.25 | control, two phantom | the other null test, on the same ladder |
| `convA_pm_sep12_n128/192/256` | 64 | 128/192/256³ | 0.50 / 0.33 / 0.25 | runaway pair, 12 apart | envelopes no longer interpenetrating |
| `convA_pm_sep16_n128/192/256` | 64 | 128/192/256³ | 0.50 / 0.33 / 0.25 | runaway pair, 16 apart | cleanly separated — the point-mass limit |
| `boxC_pm_L128_n256` | 128 | 256³ | 0.50 | runaway pair | waves on four shells (r = 16, 24, 32, 40) and the box-size systematic; run to t = 90 |
| `boxC_pp_L128_n256` | 128 | 256³ | 0.50 | control, two positive | the same for the control |

The two `boxC` cells trade resolution for reach: 256³ spread over twice the
width gives cell size 0.50, so compare them against the `n128` column, not
`n256`. Their drifts are measured from the box centre `64 64 64`, not `32 32 32`.

The second phase (2026-08-17 evening → 2026-08-18 morning) added nine cells:
the two-phantom control on the full ladder, and the separation series at 12 and
16 — which turned out to be the most consequential result of the whole campaign
(see **Separation** below). Re-running the pack script picks them all up.

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

**Separation: the stars were touching, and that mattered.** In the baseline
cell the two stars are 8 apart but each is about 5 across, so their envelopes
overlap from the first moment — you can see it in the very first frame of the
matter movies. Repeating the run with the pair 12 and then 16 apart shows what
that overlap was doing. Against a point-mass calculation of the same
configuration, the measured drift is 2.14× too large at separation 8, 1.36× at
12 and 1.17× at 16 (all at t = 30). Pull the stars apart and the run converges
on the textbook answer, which is exactly what an overlap effect should do and
what a numerical artefact should not.

The trade-off is that the clean geometry is also the faint one: separation 16
drifts about 40× less than separation 8, so the grid-to-grid spread grows from
1.9 % to 40 % at t = 30. Separation 12 is the usable compromise. **Quote
separation 12–16 for magnitudes, separation 8 for illustration.** Details and
the full ladder are in [`FINDINGS.md`](../../research/bondi_dipole/docs/FINDINGS.md) §5.5.

**The two lumps do not stay the same size.** Through every mixed run the
phantom lump stays compact while the canonical lump fades and spreads. Two
separate things are doing that, and it is worth keeping them apart:

- The canonical star's **core is collapsing**. Its peak field strength rises
  36 %, three-quarters of its material ends up within a radius of 4.8, and the
  lapse and the conformal factor fall from 0.98 to 0.35 and 0.20 by t = 60.
  This is the same star that goes on to form a black hole at t ≈ 68 in the
  double-size box. That is physics.
- The canonical star's **halo is leaving the box**. It is the lump being pushed
  toward +x, so its shed matter reaches the edge first. Losing the leading edge
  also drags its measured centre backwards — which is why the widest-separation
  cell reports a small *negative* canonical drift.

This does not conflict with the earlier result that the phantom sector is the
better-behaved one. That claim is about shape and geometry, and it still holds
everywhere: the canonical star spreads 2.7× against the phantom's 1.2× when
each is alone in the box, and the two-canonical control drives the gravity well
down to 0.44 while the two-phantom control holds it at 1.00. The total amount
of field being tracked is *not* a conserved mass — it rises for almost every
star here, and it rises faster for the canonical one. Practical consequence:
the phantom side carries the quantitative tests, and late-time canonical
numbers are upper bounds. See [`FINDINGS.md`](../../research/bondi_dipole/docs/FINDINGS.md) §5.6 and §6.

**No successive-difference order fit is quotable.** The three-grid ratios sit
at 0.2–1.6, below the 1.41 that any positive convergence order requires, so the
error is dominated by initial-data offsets (~1e-3) rather than by a clean
truncation term. Quote the grid-to-grid spread as the error bar, not an order.
The constraint section below measures that directly: the violation is the same
on all three grids because they share one elliptic solve.

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

### Constraint violation

Full tables: [`analysis/constraint_check.md`](analysis/constraint_check.md).

Einstein's equations carry four constraints that a correct solution must
satisfy at all times. They are never exactly zero numerically, and their size —
plus whether refining the grid shrinks them — is what says how far the reported
physics can be trusted. Two stages fail independently and are reported
separately.

**Initial data.** Every cell's elliptic solve converged below the 0.1 % exit
tolerance before the evolution began:

| family | iterations | Hamiltonian | momentum |
|---|---|---|---|
| runaway pair (all separations) | 7 | 0.083 % | 0.075 % |
| equal-mass pair | 7 | 0.082 % | 0.079 % |
| two-positive control | 4 | 0.004 % | 0.045 % |
| two-phantom control | 7 | 0.094 % | 0.035 % |

**Evolution.** Box-averaged violation for the runaway pair, identical on all
three grids to three or four digits:

| t | Hamiltonian (L2) | momentum (L2) | worst single point |
|---|---|---|---|
| 15 | 1.07e−3 | 1.20e−4 | 1.9e−2 |
| 30 | 2.48e−3 | 1.51e−4 | 5.5e−2 |
| 45 | 4.07e−3 | 2.43e−4 | 1.8e−1 |
| 60 | 6.91e−3 | 7.46e−4 | 8.4e−1 |

Momentum violation runs about ten times below the Hamiltonian throughout. The
sharp rise in the worst-point column past t ≈ 45 tracks the canonical star
compactifying toward the collapse documented above; it is not a code failure,
and it sits outside the t ≲ 50 trust window in any case.

**Refining the grid does not reduce the violation, and that is the reason no
convergence order is quotable.** Dividing the coarsest grid by the finest gives
1.00× at every sampled time for both mixed-pair ladders, where a second-order
scheme over a 2× refinement would give 4×. The cause is structural rather than
a code defect: the initial data is solved once on its own fixed elliptic grid
(box 128, 3 levels) and handed unchanged to all three evolution grids, so the
solve's residual is a floor that evolution refinement cannot get under. The
first-step numbers confirm it — the violation *rises* with resolution there
(1.0e−2 → 2.3e−2 → 4.8e−2), exactly what happens when a fixed error field is
drawn more sharply. Getting a genuine order out of this campaign would require
refining the elliptic solve alongside the evolution grid, which is a fresh set
of runs.

**The two-phantom control is the quietest cell in the campaign**, as it is for
drift: at t = 60 its box-averaged violation is half the mixed pair's (3.3e−3
against 7.0e−3) and its worst single point is 27× smaller (3.0e−2 against
8.3e−1). Positive masses attract and pile into a sharp lump; negative masses
repel and stay smooth.

Two cautions when reading the raw stream:

- `L2_Ham_rel`/`L2_Mom_rel` (columns 7 and 8) sit near 1 for these cells and are
  **not** a useful measure here. Both numerator and denominator are averaged
  over a box that is almost entirely empty space, where the constraint terms
  are individually at noise level and there is nothing to cancel. Use the
  absolute norms.
- The box-averaged norms are not comparable across box sizes: `boxC` looks 3×
  quieter than the small box at t = 60 purely because it averages over eight
  times the empty volume.

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

- These are uniform-grid runs; the AMR twin of the mixed pair is `campaign/pair_pm`.
- `boxC` drifts are `bary_x − 64`, not `− 32` as everywhere else in this pack.
- Re-running the pack script rebuilds every cell folder from scratch; nothing
  outside `campaign/` and the `analysis/*_check.*` tables is touched.

### Movies

Up to four views per cell, all on the z slice through the pair (the older
`pair_*` folders carry two or three of them; the `convA_*` folders carry all four):

| file | what it shows |
|---|---|
| `scalar_activity_z.mp4` | the matter itself, as a slice — the clearest view of the two lumps moving |
| `scalar_activity_proj_z.mp4` | the same matter, but with the whole depth stacked into one image (brighter, and it makes the overlap obvious) |
| `chi_minus_1_z.mp4` | the gravity well — positive and negative curvature have opposite colours, so the two sectors are visually distinguishable |
| `rho_req_z.mp4` | the signed energy density, i.e. where the negative-mass matter actually is |

Three campaign cells were added for the separation result: `convA_pm_n256`
(8 apart, envelopes overlapping), `convA_pm_sep12_n256` (12 apart) and
`convA_pm_sep16_n128` (16 apart, cleanly separated). Playing the three matter
movies side by side is the quickest way to see the overlap argument of
[`FINDINGS.md`](../../research/bondi_dipole/docs/FINDINGS.md) §5.5.

Two cautions when reading them. **The colour scale on `chi_minus_1_z` is
re-chosen for every frame**, because the gravity well starts almost flat and
becomes very deep — so a colour means a different number at different times,
and that movie shows *where* structure is, not *how big* it is. The other three
hold a fixed scale and can be compared across time. And in the mixed pairs the
canonical lump (the one at +x, on the right) genuinely fades and spreads as the
run proceeds; that is the boundary-loss effect of §5.6, not a rendering
artefact.

The complete set — 19 fields per cell rather than four — is left outside this
pack under `runs/bondi_rerun/<cell>/<run>/movies/`, together with the PNG frames
they were built from.

## Reproducing

Regenerate this whole pack from the run tree:

```bash
bash research/bondi_dipole/pack_results.sh    # dx=0.5 cells   -> campaign/ stars/ movies/
bash research/bondi_dipole/pack_campaign.sh   # campaign cells  -> campaign/ + convergence tables
```

Both copy only small artefacts and scrub absolute machine paths at runtime
(`grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py`) — no host, user, or site identity
enters git. To re-run the physics itself, see [`LAUNCH.md`](LAUNCH.md).

## Provenance

- Initial data: **GRTresna** (CTTKHybrid, `BosonStarBH` example, complex scalar
  matter with per-lump signs), single rank.
- Evolution: **GRTeclyn** (`RadialRecipe`, CCZ4 + bicomplex scalar matter),
  single rank, one GPU per cell.
- Working notes (narrative, not part of this pack): [`research/bondi_dipole/docs/`](../../research/bondi_dipole/docs/)
  (working copy: `research/bondi_dipole_debug.md`).
- Code state: wrapper commit recorded per cell in `data/<cell>/metadata.json`;
  the GRTresna matter modifications are in [`patches/`](patches/).
