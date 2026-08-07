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
| `make_tables.py` | regenerates `summary.*` and `trajectories.csv` from `data/` |
| `newtonian_reference.py` | regenerates the point-mass reference (pure stdlib RK4) |
| `star_family_scan.py` | regenerates `stars/star_family.csv` (needs the wrapper venv) |

## Reproducing

Regenerate this whole pack from the run tree:

```bash
bash research/bondi_dipole/pack_results.sh
```

It copies only small artefacts and scrubs absolute machine paths at runtime
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
