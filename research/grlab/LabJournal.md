# Spacetime Splash — Lab Journal

Reverse-chronological log of experiments, decisions, and results for the bosonic-shell critical-collapse campaign.  Newest entries first.

---

## Current state (2026-06-22)

### What works

- Bosonic shell search space with canonical coupling, 5-lump superposition
- Geometric splash scoring (χ, K, Ψ₄) as primary reward, matter (ρ) secondary
- Inward-converging radial velocity with sphere layout (layout=0)
- NFS-resilient postload gate with retry loop
- AMR regrid threshold guard prevents full-domain refinement segfault
- Early-terminated runs scored correctly (`gpu_ok`)
- 426 tests passing

### Outstanding

- **Geometric target scales** (χ_drop=0.6, |K|=1.0, |Ψ₄|=1e-2) are educated guesses — may need recalibration from v11 data.
- **Moving shell convergence**: some high-velocity configs still fail GRTresna (Mom > 5%).  May benefit from adaptive iteration count or velocity-dependent rejection thresholds.

### Branches

- **GRTeclyn:** `feature/interstellar`
- **GRTresna:** `feature/interstellar`

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
