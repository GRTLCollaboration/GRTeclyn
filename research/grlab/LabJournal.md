# Spacetime Splash — Lab Journal

Reverse-chronological log of experiments, decisions, and results for the bosonic-shell critical-collapse campaign.  Newest entries first.

---

## Current state (2026-06-23)

### What works

- Bosonic shell search space with canonical coupling, 5-lump superposition
- Geometric splash scoring (χ, K, GW wave proxy) as primary reward, matter (ρ) secondary
- GW strain proxy from native `A_ij` plot vars — no Weyl4 derive on splash path
- Scoring fixes: initial-relative χ-drop, horizon bonus gated on center well, scalar wave floor removed from objective
- **V13 scoring**: survival forgiveness on deep collapse, `GW_WAVE_PEAK_TARGET=0.20`, tighter velocity bounds, frames off by default
- Offline rescore: `grteclyn-wrapper/scripts/search/rescore_splash_campaign.py`
- Inward-converging radial velocity with sphere layout (layout=0)
- NFS-resilient postload gate with retry loop
- AMR regrid threshold guard prevents full-domain refinement segfault
- Early-terminated runs scored correctly (`gpu_ok`)
- **`spacetime_splash_v12` complete** — 31/50 gpu_ok, 0 segfaults (see entry below)

### Outstanding

- **`spacetime_splash_v13` launching** — v13 scoring (geometric-only forgiveness, late GW gate, tighter velocities)
- Confirmed **ingoing GW → center focus** splash (late converging wave, not just early A_ij spike)
- **ham_abs NaN** at origin → −50 constraint penalty (deferred)

### Branches

- **GRTeclyn:** `feature/interstellar`
- **GRTresna:** `feature/interstellar`

---

## 2026-06-23 — V13 readiness: survival paradox & GW recalibration

### Analysis of v12 leaderboard

Analysis of the v12 results revealed that the scoring function actively punishes the best physical candidates.

1. **The survival paradox:** eval 004 and eval 016 triggered rapid collapse (χ dropped, |K| spiked). Physical success caused early termination (survival 0.09 and 0.26). Because all rewards were multiplied by `survival`, scores were crushed. Meanwhile eval 008 (static dud) survived to t≈16 and ranked high on saturated base points.
2. **GW target saturation:** the A_ij proxy amplitude is an order of magnitude above true Ψ₄. `GW_WAVE_PEAK_TARGET = 1e-2` saturated on noise, giving a free 600 points to every run.

### Fixes implemented for v13

- **Survival forgiveness (geometric only):** in `objectives.py`, `geometric_survival = 1.0` when `geometric_curvature_well > 0.4` or a horizon forms; matter terms still use raw `survival`.
- **GW proxy recalibration:** `GW_WAVE_PEAK_TARGET` raised from `1e-2` to `0.20`; peak must occur at **t ≥ 2** (`GW_WAVE_MIN_FOCUS_TIME`).
- **Velocity bounds:** radial/poloidal velocity tightened to `[-0.2, 0.2]` in `spaces.py` to reduce GRTresna Mom rejections.
- **Frame overhead:** `GRTECLYN_FRAMES=0` default in `run.sh` (scoring uses central timeseries only).

### v13 scoring tune (post–eval 16 review)

Eval 16 rescore at 2410 exposed stacked rewards: full geometric saturation + matter peak 0.70 + survival forgiveness on *all* terms. Follow-up tuning before launch:

- Forgiveness scoped to **χ/K/GW/lapse only** — ρ/focus still × raw survival.
- Forgiveness threshold lowered **0.8 → 0.4** so eval 004-class splashes qualify.
- GW credit uses peak **after t = 2** only (early A_ij noise zeroed).

### Next step

Launch **`spacetime_splash_v13`** (50 evals, 8 GPUs) with tuned scoring.

**Offline rescore (3 evals with `small_data/`):** eval 016 **1687** (was 738), eval 044 **917** (was 1020), eval 008 **450** (was 786). Geometric-only forgiveness + late GW gate removes eval 008 false lead; eval 016 stays #1 without the 2410 inflation from matter terms.

---

## 2026-06-23 — spacetime_splash_v13 launch

Launch: `QD_NAME=spacetime_splash_v13 QD_TARGET_EVALS=50`, 8 GPUs, v13 scoring (geometric-only survival forgiveness, `GW_WAVE_PEAK_TARGET=0.20`, late GW gate t≥2, velocity bounds ±0.2, `GRTECLYN_FRAMES=0`).

Run dir: `runs/grtresna_qd/spacetime_splash_v13/`.

---

Campaign complete (`50/50` evals, exit 0). Run dir: `runs/grtresna_qd/spacetime_splash_v12/`.

### What was added (vs v11)

| Change | Effect |
|--------|--------|
| A_ij → GW proxy in central timeseries | Finite `weyl4` column on all gpu_ok; no NaN Ψ₄ |
| Initial-relative χ-drop | eval 008 no longer gets false χ-well credit |
| Horizon bonus gated on χ-well ≥ 0.15 | Zero horizon bonuses across all gpu_ok |
| `wave_focusing_quality` removed from objective | No ~100 pt scalar FFT floor |
| Recalibrated χ/K targets (0.20 / 0.35) | Better separation of modest vs strong geometry |

### Outcome mix (50/50)

| status | count | notes |
|--------|-------|-------|
| `gpu_ok` | 31 | 62% success rate |
| `postload_rejected` | 18 | Ham L2 > 0.03 |
| `grtresna_rejected` | 1 | Mom/Ham convergence |
| `gpu_failed` | 0 | no segfaults |

### Top gpu_ok (raw score — still miscalibrated)

| eval | score | χ-well | \|K\| | survival | physics note |
|------|-------|--------|-------|----------|--------------|
| 044 | 1020 | 0.52 | 0.81 | 0.78 | **new #1** — strong geometry + best survival among elites |
| 008 | 786 | 0.00 | 0.60 | 1.00 | GW floor + lapse; χ **rises** at origin |
| 016 | 738 | 1.00 | 1.00 | 0.26 | **best geometry** — full χ-well + K saturation |
| 039 | 585 | 0.00 | 0.62 | 0.72 | high survival, weak χ-well |
| 034 | 473 | 0.00 | 0.74 | 0.50 | K-driven mid-tier |
| 004 | 66 | 0.46 | 0.88 | 0.09 | **most splash-like** (χ drops, ρ peaks early); survival kills rank |

19/31 gpu_ok have χ-well > 0.15 — genuine curvature dynamics, not pure matter pile-up.

### Still missing

- Recalibrated **GW_WAVE_PEAK_TARGET** and offline v12 rescore (expect 004/016/044 to reorder vs 008/039)
- Confirmed **ingoing GW → center focus** splash (late converging wave, not just early A_ij spike)
- Stable long-run survival for highest-geometry elites (016/015 stall ~t≈11 via early termination)

---

## 2026-06-18 — v11 score inflation fix + spacetime_splash_v12 launch

### v11 validation

`spacetime_splash_v11` (50/50 evals) completed but leaderboard was untrustworthy:

| Issue | Fix |
|-------|-----|
| `#1 eval 008 @ 947` dominated by ungated horizon bonus + false χ-well | Gate horizon on `geometric_curvature_well ≥ 0.15`; χ-drop = initial − min(χ) |
| `geometric_wave_arrival = 0` everywhere | Wire A_ij → GW proxy into `central_timeseries.dat` col 12 (`weyl4`) |
| ~100 pt scalar floor from `wave_focusing_quality` | Removed from `critical_collapse` objective (descriptor only) |
| Targets miscalibrated | `CHI_DROP_TARGET` 0.6→0.20, `ABS_K_TARGET` 1.0→0.35 |

Offline rescore on v11 evals with `small_data/`:

| eval | old | new |
|------|-----|-----|
| 015 | 505 | **799** (best physics; higher well target saturation) |
| 008 | 947 | **187** (horizon + false well removed) |
| 026 | 325 | **19** |

Launch: `QD_NAME=spacetime_splash_v12 QD_TARGET_EVALS=50` with corrected pipeline.

---

## 2026-06-22 — Segfault fix + early-termination classification

### Symptom

`spacetime_splash_v11` crashed on the first GPU timestep with a `MultiFab` segfault in `storeRKCoarseData`.  Root cause was the AMR chi-tagger refining **100% of the domain** (512 grids, 256³ Level 1) on broad boson-shell initial data, exhausting GPU memory.

### Fix

- `search_common.sh`: raised `regrid_threshold` from 0.01 → 0.1 so refinement stays on the central collapse region (~0.1% of domain, 1–2 Level 1 grids).
- `spaces.py`: capped `grtresna_shell_amp` at 0.12 (was 0.15) to keep matter compact.
- `splash/run.sh`: reduced `PLOT_INTERVAL` from 40 → 80 to cut GPU/disk load.
- `evaluation.py`: fixed trajectory classification — early-terminated runs via `.stop_sim` (e.g. `dispersion_complete`) were recorded as `gpu_failed` due to SIGTERM exit code (-15); now treated as `gpu_ok` so they are scored normally.

**Commits:** `2ae95b7`, `3d1ac8e`

Full test suite: 426 passed, 3 skipped.

### Campaign results

`spacetime_splash_v11` relaunched (50 evals, 8×H100, 20 iterations).  First 24 evals:

| status | count | notes |
|--------|-------|-------|
| `gpu_ok` | 5 | eval 8: 947, eval 15: 505, eval 17: 151, eval 13: 27, eval 11: -37 |
| `postload_rejected` | 10 | legitimate Ham L2 > 0.03 |
| `gpu_failed` | 7 | early batch from before SIGTERM fix; future early-terminated runs will be `gpu_ok` |
| `grtresna_rejected` | 2 | GRTresna Mom/Ham convergence > 5% |

No segfaults.  Strongest candidates (eval 8, eval 15) show genuine geometric dynamics.

---

## 2026-06-22 — NFS race condition fix

Fixed NFS read-after-write race causing spurious postload rejections (`constraint_norms.dat missing or empty`).  Added retry loop with cache invalidation in `run_postload_gate()` (5 retries × 0.5s).  Added tests for delayed materialization and missing file.

**Commit:** `ab88333` — `fix: NFS read-after-write race in postload constraint gate`

---

## 2026-06-20 — Geometric splash scoring

### Motivation

Matter density (ρ) alone is an ambiguous collapse indicator — a slowly accreting blob scores well on `central_energy_peak` but is not a spacetime splash.  The genuine signal is a **converging gravitational wave** that: (1) crushes the conformal factor (χ → 0), (2) spikes the extrinsic-curvature trace |K|, and (3) emits a Weyl/Ψ₄ pulse at the center.

### New diagnostics (central timeseries extraction)

Added `chi_min`, `abs_K_max`, and `weyl4_re` columns to the central timeseries extraction (`extraction/central.py`).  Extended `CentralFieldMetrics` with derived geometric properties:

| Property | Definition |
|----------|-----------|
| `chi_drop` | `1 − min(chi_min)` over the evolution |
| `peak_abs_K` | `max(abs_K_max)` over the evolution |
| `peak_abs_weyl4` | `max(abs(weyl4_re))` over the evolution |
| `has_geometric_data` | True when at least one geometric column was extracted |

### Geometric splash scoring components (`splash.py`)

Three new components in `_compute_geometric_components()`:

| Component | Target scale | Meaning |
|-----------|-------------|---------|
| `geometric_curvature_well` | `CHI_DROP_TARGET = 0.6` | χ collapse: strong curvature well at center |
| `geometric_wave_arrival` | `WEYL4_PEAK_TARGET = 1e-2` | Ψ₄ pulse at center from focused GW |
| `geometric_crunch` | `ABS_K_TARGET = 1.0` | Crunch rate from extrinsic curvature spike |

Each is normalized to `[0, 1]` via `min(1, value / target)`.

### Objective reweighting (`_critical_collapse_total`)

Geometric terms are now **primary**; matter terms are **secondary**:

```
total = 800 × curvature_well × survival    # Primary (geometric)
      + 600 × wave_arrival × survival
      + 300 × crunch × survival
      + 400 × peak × survival              # Secondary (matter)
      + 200 × focus_effective × survival
      + 200 × wave × survival
      + 200 × lapse_progress × survival    # (discovery mode)
      + 500 × horizon_bonus × survival     # (if horizon forms)
      + 100 × pre_collapsed_penalty
      + 100 × dispersion_penalty
      −  50 × (1 − constraint_quality) × survival
      +  50 × exotic_capped × survival
```

Maximum theoretical geometric contribution: 1700 / ~2900 total = 59% of budget.

### Inward-converging shells

Radial velocity bounds adjusted to prioritize inward motion: `[-0.5, 0.5]` with default seed `−0.2` (inward).  Toroidal/poloidal bounds tightened to ±0.2 / ±0.15 to keep kinetic energy mainly radial.

### Critical bug: matter layout

Discovered that `matter_layout=1` (channel layout) applies a **uniform** velocity to all lumps — net linear momentum, not inward convergence.  Pinned `matter_layout=0` (sphere layout) so each lump moves inward from its own position.  This is essential for spacetime splash.

**Commit:** `dbbd4fa` — `feat: spacetime splash — geometric (chi/K/Weyl4) diagnostics + scoring`

---

## 2026-06-18 — Baseline: exotic penalty fix + moving shell

### Exotic penalty cap

The `exotic_penalty` in `_critical_collapse_total()` was capping at −1.6 (costing −320 points) because evolved geometry naturally develops ρ < 0 regions near the focusing core — gauge artifacts, not bad initial data.  Capped the penalty at `max(exotic, −0.3)` with weight 50, reducing the maximum cost to −15 points.

### Moving boson shell

Added velocity dimensions to the search space (`grtresna_shell_toroidal_velocity`, `grtresna_shell_poloidal_velocity`, `grtresna_shell_radial_velocity`, `grtresna_shell_omega`).  When `static=False`, the campaign includes velocity search.  Initial bounds were later tightened due to GRTresna convergence issues.

### Amplitude & postload gate

- Upper bound of `grtresna_shell_amp` widened from 0.12 → 0.15 to push closer to collapse.
- `POSTLOAD_MAX_HAM_L2` relaxed from `2e-2` → `3e-2`, recovering ~15% of configs rejected due to grid interpolation noise.

### Temporal sampling

`PLOT_INTERVAL` in `splash/run.sh` reduced from 320 → 40, yielding ~8 central-timeseries samples per evolution instead of ~1.

### Validation campaigns (v9)

| Campaign | GPUs | Evals | Notes |
|----------|------|-------|-------|
| `boson_shell_v9_static_validate` | 0–3 | 50 | Static shell baseline |
| `boson_shell_v9_moving_validate` | 4–7 | 50 | Moving shell — killed early due to Mom convergence |

**Commit:** `d933cf1` — `feat: broader splash campaign — exotic fix, moving shell, higher amp`
