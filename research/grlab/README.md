# Gravitational Splash — MAP-Elites Implementation

Quality-diversity search for **critical collapse**: configurations where matter at the grid origin shows inward wave focusing, rising central energy density, and (eventually) horizon formation. Implemented in `grteclyn-wrapper` as a first-class campaign alongside the FTL / boson-star studies.

## Scientific goal

Find GRTresna initial data that evolves into **non-trivial splash dynamics** at the origin:

- Central `rho_req` grows toward collapse-scale values (~`10⁻²`)
- Lapse at the origin falls toward a collapse band
- Optional corroborated trapped-surface proxy (horizon bonus)
- Diverse wave profiles across the MAP-Elites archive (`wave_focusing` descriptors)

FTL / warp-motor terms are **not** used in this objective.

## Campaign launcher

```bash
cd grteclyn-wrapper
QD_NAME=scalar_splash_gpu_v1 QD_TARGET_EVALS=12 \
  GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=1 \
  bash scripts/campaigns/splash/run.sh
```

**Current defaults** (`scripts/campaigns/splash/run.sh`):

| Variable | Value | Role |
|----------|-------|------|
| `GRTRESNA_MATTER_SECTOR` | `scalar` | Validated real-scalar lumps (`ScalarFieldBH`) |
| `GRTRESNA_ANSATZ` | `shell` | Spherical shell of momentum-carrying lumps |
| `GRTRESNA_SHELL_PROFILE` | `compact` | Small lump width (1.8–3.0) so frames are not filled by a single blob |
| `GRTRESNA_MATTER_COUPLING` | `canonical` | Positive-energy matter only |
| `OBJECTIVE_MODE` | `critical_collapse` | Splash scorer (see below) |
| `DESCRIPTOR_MODE` | `wave_focusing` | MAP-Elites behavior space |
| `SPLASH_MODE` | `discovery` | Broad collapse search (vs `threshold` lapse band) |
| `GRTECLYN_FRAMES` | `1` | Slice / projection frames per eval |
| `GRTECLYN_FRAMES_AUTO_ZLIM` | `1` | Auto colorbar scale (scalar lumps are much fainter than boson defaults) |
| `GRTECLYN_FRAMES_ZOOM` | `28` | Zoom into origin on L=64 domain |

`search_common.sh` auto-enables `GRTECLYN_CENTRAL_TIMESERIES=1`, disables 4D geodesic FTL, and skips FTL-specific preflight when `OBJECTIVE_MODE=critical_collapse`.

## What was added

### 1. Central field metrics (origin time series)

Plotfile consumer extracts samples at the grid origin each dump:

- **File:** `small_data/central_timeseries.dat`  
  Columns: `time`, `rho_req`, `lapse`, `scalar_activity`, `phi_re`, `phi_im`
- **Module:** `grteclyn_wrapper/visualisation/process_wave/consume_plotfiles/extraction/central.py`
- **Types / reader:** `metrics/types/central.py`, `metrics/diagnostics/central.py`
- **Wiring:** collector, episode metrics, metrics catalog (`central` group)
- **CLI:** `--central-timeseries` on the plotfile consumer; enabled via `GRTECLYN_CENTRAL_TIMESERIES=1`

Derived quantities include `peak_rho_req_at_origin`, `focusing_efficiency` (peak/initial ρ at origin), and `wave_chromaticity` (FFT Q-factor of central scalar activity).

### 2. Splash score components

**Module:** `metrics/score/splash.py` — `compute_splash_components()`

| Component | Meaning |
|-----------|---------|
| `central_energy_peak` | Peak origin ρ, normalized by target `10⁻²` |
| `focusing_efficiency` | Relative ρ growth at origin (capped at 5×) |
| `wave_focusing_quality` | Chromaticity, down-weighted if activity fades (dispersion) |
| `collapse_lapse_progress` | Graded reward as min lapse → 0.05–0.2 (discovery) |
| `central_lapse_collapse` | Band reward/penalty around lapse 0.01–0.05 (`threshold` mode) |
| `dispersion_penalty` | Penalty for fading blobs with weak central peaks |
| `pre_collapsed_penalty` | Penalty if already dense/collapsed at t = 0 |
| `horizon_formation_time` | Bonus when corroborated trapped surface appears |

### 3. Objective: `critical_collapse`

**Module:** `metrics/score/objectives.py` — `_critical_collapse_total()`

Discovery-mode scalarization (survival gates everything):

```
score = 1000 × peak × survival
      + 300 × min(focus, 5) × peak × survival   # gated: relative growth needs absolute ρ
      + 200 × wave_focusing_quality × survival
      + 200 × collapse_lapse_progress × survival
      + 500 × horizon_formation_time × survival   # if horizon forms
      + 100 × pre_collapsed_penalty
      + 100 × dispersion_penalty
```

`SPLASH_MODE=threshold` swaps lapse progress for the sharper `central_lapse_collapse` band term (+500×).

Registered in CLI: `--objective-mode critical_collapse`.

### 4. MAP-Elites descriptors: `wave_focusing`

**Module:** `search/qd_search/descriptors.py`

- **X:** `wave_chromaticity` from central metrics, or fallback from search params (`grtresna_bs_omega` / `grtresna_shell_omega`)
- **Y:** Normalized lump width (`grtresna_bs_profile_width` or `grtresna_shell_width`)

### 5. Search spaces

- **`grtresna_boson_splash_search_space()`** — 7-D boson star with unpinned ω (compact Gaussian defaults). Kept for future boson-splash experiments; not used by the current campaign launcher.
- **Active campaign** uses the existing **`grtresna_shell_search_space(profile=compact)`** scalar shell space (18-D).

Ansatz `splash` in CLI still maps to boson-star matter for legacy compatibility; the campaign script uses `shell` + `scalar` explicitly.

### 6. Tests

Splash-related unit tests live under `grteclyn-wrapper/tests/`:

- `metrics/diagnostics/test_central_timeseries_parser.py`
- `metrics/aggregation/test_central_collector.py`
- `visualisation/test_central_extraction.py`
- `metrics/score/test_splash_components.py`
- `metrics/score/test_critical_collapse_objective.py`
- `metrics/score/test_scorer_splash_integration.py`
- `search/test_splash_descriptors.py`
- `grtresna/test_boson_splash_search_space.py`
- `scripts/test_splash_campaign_env.py`

Included in QD preflight when `OBJECTIVE_MODE=critical_collapse`.

## Run artifacts

Per eval under `runs/grtresna_qd/<QD_NAME>/eval_NNNNNN/`:

- `small_data/central_timeseries.dat` — origin diagnostics
- `score.json` — full metrics + `critical_collapse` components
- `frames/` — field slices and projections (when `GRTECLYN_FRAMES=1`)
- `data/collapse_diagnostics.dat` — global horizon proxies

Campaign-level: `trajectory.jsonl` (score, descriptors, overrides, status per eval).

## Design notes

- **Scalar vs boson:** Early GPU runs used boson-star (`BosonStarBH`) splash ansatz; production campaign switched to **canonical scalar shell** because GRTresna boson convergence was brittle and produced symmetric center blobs with inflated scores. Boson splash search space and tests remain for later work.
- **Score calibration:** Initial runs scored ~800–1000 for modest 3× ρ bumps at the origin; gating `focusing_efficiency` by `peak_norm` and adding dispersion / wave terms fixes that.
- **Compact lumps:** Shell `compact` profile (width 1.8–3.0) keeps matter off the frame edges; boson splash space uses the same idea via `bs_profile_width` default 3.5 when re-enabled.

## Related paths

```
grteclyn-wrapper/
  scripts/campaigns/splash/run.sh
  scripts/campaigns/lib/search_common.sh      # critical_collapse env block
  src/grteclyn_wrapper/metrics/score/splash.py
  src/grteclyn_wrapper/metrics/score/objectives.py
  src/grteclyn_wrapper/search/qd_search/descriptors.py
  src/grteclyn_wrapper/search/optimize/spaces.py   # grtresna_boson_splash_search_space
```
