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
QD_NAME=boson_shell_gpu_v1 QD_TARGET_EVALS=12 \
  GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=1 \
  bash scripts/campaigns/splash/run.sh
```

**Current defaults** (`scripts/campaigns/splash/run.sh`):

| Variable | Value | Role |
|----------|-------|------|
| `GRTRESNA_MATTER_SECTOR` | `boson_star` | Complex U(1) scalar (`BosonStarBH`) |
| `GRTRESNA_ANSATZ` | `shell` | Five lump sites on a shell → one superposed boson field |
| `GRTRESNA_SHELL_PROFILE` | `compact` | Shell geometry profile (radius, layout, etc.) |
| `GRTRESNA_MATTER_COUPLING` | `canonical` | Positive-energy matter only |
| `PIN_DIMS` | `grtresna_scalar_sign=1`, `grtresna_shell_static=1` | Canonical coupling; static shell (zero lump velocities) |
| `ITERATIONS` | `80` | GRTresna NL iterations (boson shell needs more than scalar) |
| `OBJECTIVE_MODE` | `critical_collapse` | Splash scorer (see below) |
| `DESCRIPTOR_MODE` | `wave_focusing` | MAP-Elites behavior space |
| `SPLASH_MODE` | `discovery` | Broad collapse search (vs `threshold` lapse band) |
| `GRTECLYN_FRAMES` | `0` | Off during QD (scoring uses `central_timeseries`); set `1` for viz |
| `PLOT_INTERVAL` | `160` | Plotfile dump cadence (~8 samples over 16 s) |
| `GRTECLYN_SPLASH_EARLY_TERM` | `1` | Stop when matter disperses (~t≈10–12 s); set `0` for full 16 s |
| `GRTECLYN_FRAMES_AUTO_ZLIM` | `1` | Auto colorbar scale (when frames enabled) |
| `GRTECLYN_FRAMES_ZOOM` | `28` | Zoom into origin on L=64 domain |

Production example:

```bash
cd grteclyn-wrapper
QD_NAME=boson_shell_gpu_v1 QD_TARGET_EVALS=100 QD_ITERATIONS=30 \
  GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 BATCH_SIZE=8 \
  bash scripts/campaigns/splash/run.sh
```

Validation smoke (16 evals, Jun 2026): `boson_shell_validate_v8` — 15/16 pass GRTresna gate, 9/16 full pipeline, best score **+8.65** (eval 015).

When `OBJECTIVE_MODE=critical_collapse`, `search_common.sh` also sets:

| Variable | Role |
|----------|------|
| `GRTECLYN_CENTRAL_TIMESERIES` | Origin / central-ball diagnostics |
| `GRTECLYN_CENTRAL_BALL` | Sphere-averaged central sampling (resolution-robust) |
| `GRTECLYN_CENTRAL_RADIAL` | Radial ρ/lapse profiles → `central_radial_profile.dat` |
| `GRTECLYN_INCREMENTAL_SCORE` | Live prefix scores in `score_timeseries.jsonl` |
| `GRTECLYN_SPLASH_EARLY_TERM` | Stop sim when collapse/dispersion predicates fire |

`search_common.sh` disables 4D geodesic FTL and skips FTL-specific preflight when `OBJECTIVE_MODE=critical_collapse`.

## What was added

### 1. Central field metrics (origin time series)

Plotfile consumer extracts samples at the grid origin each dump:

- **File:** `small_data/central_timeseries.dat`  
  Columns: `time`, `rho_req`, `lapse`, `scalar_activity`, `phi_re`, `phi_im`, `noether_charge`, `phase_coherence`, `ham_abs`, `mom_abs` (trailing columns optional)
- **Radial profile:** `small_data/central_radial_profile.dat` — shell-averaged ρ/lapse vs radius per dump; drives `splash_width`, `peak_radius`, `compression_ratio`, `cusp_unresolved`
- **Sampling:** With `GRTECLYN_CENTRAL_BALL=1`, fields are averaged over a small sphere (~2× finest dx) instead of a single AMR cell
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
| `wave_focusing_quality` | Chromaticity (MAP-Elites descriptor only; not in `critical_collapse` total) |
| `geometric_curvature_well` | Initial-relative χ drop at origin (target 0.20) |
| `geometric_crunch` | Peak \|K\| at origin (target 0.35) |
| `geometric_wave_arrival` | Peak GW signal at origin — `Weyl4_Re` or A_ij proxy `\|h\|` (target 1e-2) |
| `collapse_lapse_progress` | Graded reward as min lapse → 0.05–0.2 (discovery) |
| `central_lapse_collapse` | Band reward/penalty around lapse 0.01–0.05 (`threshold` mode) |
| `dispersion_penalty` | Penalty for fading blobs with weak central peaks |
| `pre_collapsed_penalty` | Penalty if already dense/collapsed at t = 0 |
| `horizon_formation_time` | Bonus when corroborated trapped surface appears **and** center χ-well ≥ 0.15 |
| `constraint_quality` | Origin Ham vs ρ at peak (splash-only plot vars) |
| `peak_radius`, `splash_width`, `compression_ratio`, `cusp_unresolved` | Radial-profile diagnostics |

### 3. Objective: `critical_collapse`

**Module:** `metrics/score/objectives.py` — `_critical_collapse_total()`

Discovery-mode scalarization (survival gates everything; geometric terms primary since v12):

```
score = 800 × geometric_curvature_well × survival
      + 600 × geometric_wave_arrival × survival
      + 300 × geometric_crunch × survival
      + 400 × peak × survival
      + 200 × min(focus, 5) × peak × survival
      + 200 × collapse_lapse_progress × survival
      + 500 × horizon_formation_time × survival   # gated on center χ-well
      + 100 × pre_collapsed_penalty + 100 × dispersion_penalty
      − 50 × (1 − constraint_quality) × survival
```

**GW wave proxy (splash path):** RadialRecipe campaigns dump native `A11 A12 A22` in plotfiles (via `splash_wiring.py`). Central extraction computes `|h| ≈ sqrt((A11−A22)² + (2A12)²)` and stores it in the `weyl4` column of `central_timeseries.dat`. SupportedWormholeCollapse uses true `Weyl4` derive instead — different example.

**Offline rescore:** `uv run python scripts/search/rescore_splash_campaign.py runs/grtresna_qd/<campaign>/`

`SPLASH_MODE=threshold` swaps lapse progress for the sharper `central_lapse_collapse` band term (+500×).

Registered in CLI: `--objective-mode critical_collapse`.

### 4. MAP-Elites descriptors: `wave_focusing`

**Module:** `search/qd_search/descriptors.py`

- **X:** `wave_chromaticity` from central metrics, or fallback from search params (`grtresna_bs_omega` / `grtresna_shell_omega`)
- **Y:** Normalized splash timing `peak_rho_time / stop_time`, fallback: normalized `grtresna_shell_radius`

### 5. Search spaces

- **`grtresna_boson_shell_search_space()`** — active campaign space. Scalar shell geometry minus exotic dims, with boson-specific convergence bounds:
  - `grtresna_shell_amp` 0.04–0.12 (not scalar 0.08–0.16)
  - `grtresna_scalar_mass` 0.05–0.35, `grtresna_scalar_lambda` 0–0.05
  - `grtresna_bs_omega` 0.05–0.35
  - No velocity dims (campaign pins `grtresna_shell_static=1`)
- **`grtresna_boson_splash_search_space()`** — 7-D centered boson star (legacy single Gaussian). Kept for experiments; not used by the launcher.
- **`grtresna_shell_search_space(profile=compact)`** — scalar shell (18-D + exotic). Still valid for canonical real-scalar splash without boson work.

Shell expansion for boson: Python `build_grtresna_config()` expands `grtresna_shell_*` → `cfg.lumps`, then `apply_boson_star_overrides()`. GRTresna receives `num_lumps` + `lump{k}_amp/width/center/velocity/omega/mode/profile` in `params.txt`.

Ansatz `splash` in CLI still maps to boson-star matter for legacy compatibility; the campaign uses `shell` + `boson_star`.

### 6. Tests

Splash-related unit tests live under `grteclyn-wrapper/tests/`:

- `metrics/diagnostics/test_central_timeseries_parser.py`
- `metrics/aggregation/test_central_collector.py`
- `visualisation/test_central_extraction.py`
- `metrics/score/test_splash_components.py`
- `metrics/score/test_critical_collapse_objective.py`
- `metrics/score/test_scorer_splash_integration.py`
- `search/test_splash_descriptors.py`
- `grtresna/test_boson_shell_ansatz.py`
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

## Bosonic shell (Jun 2026)

Shell lump geometry is now wired to **bosonic matter**: five Gaussian sites superpose into one U(1) complex scalar solved by `BosonStarBH`, not five independent real scalars.

### Architecture

| Layer | What it does |
|-------|----------------|
| **Python** (`config.py`, `spaces.py`) | Expand `grtresna_shell_*` → `cfg.lumps`; force exotic=0 for boson |
| **GRTresna** (`BosonStarParams.hpp`, `ComplexScalarField.cpp`) | Read `num_lumps` + indexed lump params; paint `phi_re` from superposed Gaussians |
| **Evolution** | Single complex scalar: Re=`phi`/`Pi`, Im=`phi_lump0`/`Pi_lump0` |

One complex field carries all shell sites — not five independent complex fields.

### GRTresna paint (Phase 2)

Initial-data paint at each grid point:

- `phi_re` = Σ lump Gaussians (respects `mode`, `profile` Gaussian vs top-hat)
- `phi_im` = 0
- `Pi_re` = Σ lump boost/rotation terms (`velocity`, `omega` per lump — same formula as `ScalarFieldBH`)
- `Pi_im` = −ω × `phi_re` (global U(1) phase velocity from `bs_omega`)

Params read from `params.txt`: `lump{k}_amp`, `width`, `center`, `velocity`, `omega`, `mode`, `profile`.

### Blockers encountered (v6/v7)

| Issue | Symptom | Fix |
|-------|---------|-----|
| Geometry-only paint | 100% `grtresna_rejected`; Mom 12–30% | Add `lump_pi1()` / `total_pi1()`; read full lump kinematics from params |
| Scalar shell bounds on boson | Ham/Mom never converged | Boson-specific amp/mass/λ caps in `grtresna_boson_shell_search_space()` |
| Moving shell + incomplete Π | Momentum constraint stall | Pin `grtresna_shell_static=1` in campaign (static shell for first production sweeps) |
| Too few NL iterations | Residual plateau | `ITERATIONS=80` in `splash/run.sh` |

### Validation (`boson_shell_validate_v8`, 16 evals)

| Gate | Pass rate |
|------|-----------|
| GRTresna Ham/Mom (5%/5%) | **15/16** (was 0/20 on pre-fix v7) |
| Postload Ham L2 (2e⁻²) | 10/16 |
| Full GPU pipeline | **9/16** (56%) |

Best elite: eval **015**, score **+8.65**. Frames show clumpy shell structure in `rho_req_z` / `phi_z` (not a single centered blob).

Remaining filter: **postload gate** on ~6 configs (Ham L2 ≈ 0.02–0.04 after gridinit load). Optional tuning: `POSTLOAD_MAX_HAM_L2=3e-2`.

### Branches

- **GRTeclyn:** `feature/interstellar` — `grtresna_boson_shell_search_space`, config expansion, campaign pins, tests
- **GRTresna:** `feature/interstellar` — multi-lump `BosonStarParams`, `ComplexScalarField` superposition paint

Run pytest via uv: `cd grteclyn-wrapper && uv run pytest tests/grtresna/test_boson_shell_ansatz.py -q`

Push note: if SSH fails with OpenSSL mismatch, use `env -u LD_LIBRARY_PATH git push myfork feature/interstellar`.

## Design notes

- **Scalar vs boson:** Early runs used centered boson splash (single Gaussian → diffuse center blob). Scalar shell (`ScalarFieldBH`) was the interim production path. **Bosonic shell** (Jun 2026) combines shell geometry with `BosonStarBH` superposition; validated on `boson_shell_validate_v8`. Legacy scalar shell remains available via `GRTRESNA_MATTER_SECTOR=scalar`.
- **Score calibration:** Initial runs scored ~800–1000 for modest 3× ρ bumps at the origin; gating `focusing_efficiency` by `peak_norm` and adding dispersion / wave terms fixes that. `exotic_penalty` added to `_critical_collapse_total()` when exotic matter is present.
- **Compact lumps:** Boson shell caps amp at 0.12 and width 2–4; static pin avoids unconstrained velocity search until moving-shell Π is tuned further.

## Related paths

```
grteclyn-wrapper/
  scripts/campaigns/splash/run.sh
  scripts/campaigns/lib/search_common.sh      # critical_collapse env block
  src/grteclyn_wrapper/metrics/score/splash.py
  src/grteclyn_wrapper/metrics/score/objectives.py
  src/grteclyn_wrapper/search/qd_search/descriptors.py
  src/grteclyn_wrapper/search/optimize/spaces.py   # grtresna_boson_shell_search_space
  src/grteclyn_wrapper/search/optimize/config.py   # boson + shell lump expansion
  tests/grtresna/test_boson_shell_ansatz.py

GRTresna/  (separate repo)
  Source/Matter/BosonStarParams.hpp
  Source/Matter/ComplexScalarField.cpp
```
