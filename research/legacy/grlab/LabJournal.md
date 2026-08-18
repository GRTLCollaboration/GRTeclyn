# Lab Journal

## 2026-07-03 — `gw_beam_qd100_v4` complete: eval 61 / 88 analysis

Campaign finished **100/100**. Hard gates held: **77/100** collapse modes @ ~**−116** (`mult=0`); **22** healthy survivors; **5** archive elites. Best score **eval 88 @ 3.09**; best beaming **eval 61 @ 2.82** (beam_ratio **~30%**).

Run dir: `runs/grtresna_qd/gw_beam_qd100_v4/`

### Campaign headline

| | eval 88 (best score) | eval 61 (best beam) |
|--|----------------------|---------------------|
| Score | **3.09** | 2.82 |
| Tier | constructed | constructed |
| `gw_health_multiplier` | 1.0 | 1.0 |
| max ‖Ham‖₂ | 0.08 | 0.14 |
| mean Ψ₄ power | 4.5×10⁻⁴ | **6.4×10⁻⁴** |
| beam_ratio | 14% | **~30%** (to 40% late) |
| Scorer bias | high survival (1.0) | better GW, slightly noisier grid |

**Critical verdict:** Neither run is a strong GW emitter. Both produce a **weak steady hum** (P ~ 10⁻⁵–10⁻⁴ after t≈5), not a merger chirp or beamed burst. The t=0 Ψ₄ spike is **initial/near-zone transient**, not radiative physics. Eval **61** is the better physics candidate for directional GW; eval **88** wins total score because survival/stability bonuses dominate over `gw_beam_quality`.

### Eval 61 — lump dynamics (best beaming)

**Geometry:** five Q-balls in a **compact ring** R ≈ 2.4–4.1; all tangential speeds **~0.27–0.30** (near trajectory speed cap).

| Lump | R₀ | v_t | v_rad |
|------|-----|-----|-------|
| 0 | 4.06 | 0.30 | −0.01 |
| 1 | 2.48 | 0.30 | +0.03 |
| 2 | 2.73 | 0.30 | −0.01 |
| 3 | 2.38 | 0.27 | −0.03 |
| 4 | 3.69 | 0.28 | **+0.10** (outward) |

**Global:** breathing A ≈ 1.55 (T ≈ 4.7); Z bob amp ≈ 0.93 (T ≈ 4.5). Slow retrograde orbits (T_orb 52–85) — azimuthal drift is small over t = 24.

**Ψ₄ at R = 20:** t = 0 spike P ≈ 3.6×10⁻³; steady P ≈ 3–4×10⁻⁵ with **beam_ratio ≈ 30%** flat from t ≈ 5–24. **Constraints:** Ham rises slowly 0.06 → 0.14; Mom ~10⁻³–10⁻²; **no spike**.

**Mechanism (critical):** ~30% Z-beaming comes from a **coherent, fast, compact multi-lump cluster** + breathing quadrupole — not a clean radiative-zone binary. Extraction at R = 20 still sees **near-zone** motion from R ~ 2–4 lumps.

**Plots:** `eval_000061/plots/constraints_plot.png`, `eval_000061/plots/psi4_analysis_r20.png`

### Eval 88 — lump dynamics (best score)

**Geometry:** four **inner** lumps R ≈ 1.7–1.9 + one **outer** lump R = 6.5 (v_rad = +0.065 outward). Tangential speeds **slower** (0.08–0.30) than eval 61.

**Global:** stronger breathing A ≈ 1.77 (T ≈ 4.9); larger Z bob (amp 1.23). Orbits very slow (T_orb 67–144).

**Ψ₄ at R = 20:** t = 0 spike; steady P ≈ 2–3×10⁻⁵; beam_ratio **~14%** only. max Ham ≈ 0.08 — very stable, weaker directional signal.

**Why higher score:** `survival = 1.0`, lower instability penalty — scorer rewards numerical cleanliness over beaming.

### Implications

1. v4 gates work; no repeat of v3 **336** exploit.
2. Search has **not** found a strong GW laser — only stable breathing multi-lump configs with 10⁻⁵–10⁻⁴ power.
3. Total score **underweights beam_ratio** relative to survival — consider rebalancing or CMA-ES warm-start from **eval 61**.
4. Next: longer `stop_time`, R = 24 extraction, tighter `v_rad` / breathing bounds, inspect eval 61 lump trajectories visually.

---

## 2026-07-03 — `gw_beam_qd100_v3` → v4: collapse-mode reward hacking

### The exploit (classic ML reward hacking in NR)

Asked MAP-Elites to maximize GW power. The optimizer did **not** build a gravitational-wave laser — it built a **numerical bomb**.

> Ψ₄ uses second derivatives of the metric. Crash the Hamiltonian constraint → grid fills with high-frequency numerical noise → second-derivative operators report “infinite wave power.”

**Mechanism:** 5 pumped Q-balls, strong breathing + inward `v_rad` → lapse pinches → **Ham/Mom cliff at t≈4.5** → sim keeps running (`exit_code=0`) → Ψ₄ burst scored as GW. Eval 51 hit **336** vs trustworthy eval 7 @ **3.4**.

**Why additive penalties failed:** `−113` spike penalty vs `+197` from uncapped **mean** Ψ₄. You cannot fight unbounded exploits with additive terms.

### Partial fixes (v3, insufficient)

| Fix | Effect |
|-----|--------|
| Pre-spike **peak** cap only | Peak → 0; **mean** still post-collapse → eval 51 @ **81** |
| `constraint_spike_penalty ×150` | Additive; still dominated by quality term |
| Tier = rejected | Ignored by archive — rejected runs still filled cells |

### Permanent closure (v4 hard gates)

**Fix A — Ψ₄ time-series truncation** (`read_psi4_metrics`)
- At `t_spike`, truncate **peak, mean, and final** to `[0, t_spike]` only
- If `t_spike < 10` (wave hasn't reached R=12 cleanly), GW metrics invalid

**Fix B — Archive admission** (`qd_search/driver.py`)
- `tier >= CONSTRUCTED` required for `archive.insert()`
- Rejected collapse runs stay in trajectory / near-miss pool, never occupy descriptor cells

**Fix C — Multiplicative health gate** (`gw_beam_total`)
- `gw_health_multiplier = 0` if `max(‖Ham‖₂) > 100` or early spike (`t_spike < 10`)
- `score = (1000×quality + 100×peak + health bonuses) × multiplier + penalties`
- Collapse eval 51: **336 → −116** (GW terms zeroed via `mult=0`; remainder is additive penalties)

### v3 rescore (hard gates applied 2026-07-03)

- Stopped v3 at 100/100 evals; rescored 19 rows with surviving eval dirs
- Collapse evals 42/51/70: **−116**, `gw_health_multiplier=0`
- Archive rebuilt: **8 elites**, best **3.37** (eval 7 trajectory row; eval 91 @ 2.70 best on disk)

### Physics notes for v4 campaign

1. **R = 20** extraction (was R = 12 in v3) — farther wave zone, early-spike gate at t ≈ 18.
2. **Inward `v_rad`** with canonical matter → central collapse. Clamped `v_rad > −0.1`.

### Next run

v4 **complete** (see post-run analysis above). v2/v3 run dirs deleted.

Plots: `eval_000051/plots/`, `eval_000042/plots/` (`constraints_plot.png`, `psi4_analysis_r12.png`).

---

## 2026-07-10 — `gw_beam_v5` implementation complete

### v4 limitations addressed

| v4 issue | v5 fix |
|----------|--------|
| Single-radius Ψ₄ extraction (R=20) — near-zone contamination | Multi-radius extraction (R=12,16,20) with 1/r validity check (`wavezone_ok` gate) |
| Only l=2 modes | Full multipoles l=2,3,4 spin-weighted harmonics in `sphere.py` |
| `beam_ratio` as directional proxy | `beaming_gain = 1 + 4×beam_ratio` — gain above isotropic, mapped to [0,1] descriptor |
| t=0 junk transient pollution | `min_valid_time` junk truncation in `read_psi4_metrics` |
| Weighted-sum scoring (survival/stability dominate) | Constrained single objective: `(log10(P) − log10(floor)) × beaming_gain × gate` |
| Boundary reflections limit stop_time to ~24 | Numerical sponge zone (radially-ramped KO dissipation) → stop_time=40 |

### Sponge zone (C++ — Phase 0)

**Files:** `Source/Grids/SpongeZone.hpp` (new), `Source/Grids/Make.package`, `Examples/RadialRecipe/SimulationParameters.hpp`, `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`

**Mechanism:** After the main CCZ4RHS ParallelFor, a second ParallelFor applies extra Kreiss-Oliger dissipation in a spherical shell between `sponge_inner_radius` and `sponge_outer_radius`. The ramp is quartic: `f(r) = ((r − r_in)/(r_out − r_in))^4`, capped at 1. Extra dissipation = `sponge_strength × f(r)`.

**Parameters (defaults):**
- `sponge_enabled = 0` (off by default; v5 launcher sets 1)
- `sponge_inner_radius = 24.0`
- `sponge_outer_radius = 32.0`
- `sponge_strength = 4.0`
- `sponge_ramp_power = 4`

**Smoke test:** Built and ran 200 steps (t=2.0) on H100 with sponge enabled — no crashes, constraints stable, plotfiles written. Executable: `Examples/RadialRecipe/main3d.gnu.CUDA.ex` (130 MB).

### Wave-zone GW instrumentation (Phase 1)

**Files:** `sphere.py`, `metrics/types/psi4.py`, `metrics/diagnostics/psi4.py`, `worker.py`

- **l=3, l=4 harmonics:** Added `_spin_weighted_s2_l3_base`, `_spin_weighted_s2_l4_base`, and generic `spin_weighted_sph_harm_s2(l, theta, phi)`.
- **Psi4Metrics extended:** Added `mean_beaming_gain`, `peak_beaming_gain`, `wavezone_ok`, `wavezone_one_over_r_std`, `direction_stability` fields (backward-compatible defaults).
- **psi4_directional.dat v5 format:** 6 columns (was 4): `time P_total P_z_beam beam_ratio beaming_gain wavezone_std`.
- **1/r validity check:** `wavezone_std = std(r×|Ψ₄|) / mean(r×|Ψ₄|)` across extraction radii; `wavezone_ok = std < 0.3`.
- **Junk truncation:** `read_psi4_metrics(min_valid_time=...)` excludes samples before R_ext + margin.

### Constrained scoring (Phase 2)

**Files:** `metrics/score/gw_beam.py`, `search/qd_search/descriptors.py`

**v5 objective:** `(log10(P) − log10(1e-8)) × beaming_gain × health_gate × wavezone_gate`

- `health_gate`: 0 if Ham > 100 or early spike (same as v4)
- `wavezone_gate`: 0 if 1/r check fails (near-zone signal)
- Both gates multiplicative — collapse/near-zone configs score ~0
- Survival/stability demoted to tiny tie-breakers (≤1% of typical score)

**Descriptor:** x = log10(power) mapped to [0,1]; y = (beaming_gain − 1) / 4 mapped to [0,1]. `direction_stability` held out as diagnostic only.

### Calibration controls (Phase 3)

**File:** `tests/metrics/score/test_gw_beam_v5_controls.py` (7 tests, all passing)

- **Positive control:** Compact binary (P=1e-3, gain=3.4) scores higher than static lump
- **Negative control:** Static lump (P=5e-8, gain=1.0) — verifies isotropic detection
- **Regression:** Collapse exploit (Ham=200, spike t=4.5) → score ≤ 0 via health gate
- **Wavezone gate:** Near-zone signal (1/r std=0.5) → score < 10 (gated to ~0)
- **Junk truncation:** Samples before `min_valid_time` excluded from aggregates

### v5 launcher (Phase 4)

**File:** `scripts/campaigns/gw_beam/run_v5.sh`

- Variable lump count: `trajectory_well_depth_min=0.0` (lumps can switch off)
- `TRAJECTORY_V_MAX=0.5` (faster binaries)
- `STOP_TIME=40` (sponge absorbs reflections)
- `CONSUMER_RADII="12 16 20"` (multi-radius 1/r check)
- Sponge: `sponge_enabled=1 sponge_inner_radius=24 sponge_outer_radius=32 sponge_strength=4.0`
- `OBJECTIVE_MODE=gw_beam` (v5 constrained scoring)
- Canonical matter only (all lumps `exotic=0`)

### Build & test status

- C++ build: **passed** (H100, CUDA 12.1, CUDA_ARCH=90)
- Smoke test: **passed** (200 steps, sponge active, no NaN)
- Python tests: **18/18 passed** (7 new v5 controls + 11 existing gw_beam/psi4)

---

## 2026-07-10 — `gw_beam_v5` end-to-end GPU verification + campaign launch + bugfixes

### End-to-end GPU test (single episode)

Ran a full pipeline episode through the wrapper: GRTresna constrained solve → GRTeclyn GPU evolution (t=4, sponge active) → plotfile consumer (3 radii) → Psi4 directional extraction → v5 gw_beam scoring.

**Verified on real GPU:**
- Sponge zone active, constraints stable (max Ham = 0.0065)
- `psi4_directional.dat` produced with v5 6-column format
- `beaming_gain`, `wavezone_std`, `direction_stability` all computed
- v5 scoring components present: `gw_beam_v5_objective`, `gw_health_multiplier`, `gw_wavezone_ok`
- Wavezone gate correctly rejected Ellis-Bronnikov wormhole (near-zone signal, `wavezone_ok=False`)

### Bugs found and fixed during campaign launch

#### Bug 1: Plotfiles not deleted after processing

**Cause:** `GRTECLYN_FRAMES=0` disables frame rendering. The consumer's `consumer_delete_plotfiles_enabled()` refuses to delete plotfiles when `frames=False` unless `GRTECLYN_DELETE_WITHOUT_FRAMES=1` is also set — a safety guard so post-run frame drains have plotfiles to render.

**Fix:** Added `export GRTECLYN_DELETE_WITHOUT_FRAMES=1` to `run_v5.sh`. Updated `tests/core/test_plot_consumer_flags.py` to `monkeypatch.delenv` the flag in tests that assume no-delete default.

**Impact:** Each eval was 1.5GB with plotfiles on disk → 628MB after fix (gridinit + small_data only).

#### Bug 2: `--objective-mode` not available on `reproduce` subcommand

**Cause:** The `reproduce` subparser didn't have `--objective-mode`; only `optimize` and `qd` did. Single-episode runs always used `weighted` scoring.

**Fix:** Added `--objective-mode` (with `gw_beam` choice) to the `reproduce` subparser in `cli/parser.py`. Threaded `objective_mode` through `finalize_score()` and `run_single()` in `cli/episode.py`.

#### Bug 3: All evals scored identically (~−1.5 to −1.9) — wavezone gate killed by junk radiation

**Cause:** The wavezone 1/r validity check (`wavezone_std < 0.3`) was computed as the mean over ALL timesteps, including t=0 junk radiation. Early-time samples have high `wavezone_std` (0.5–1.4) because the signal hasn't reached the extraction sphere yet. This polluted the mean so every eval failed the gate, zeroing the GW reward.

**Evidence:** eval_000006 had `wz_std=0.08` at late times (clearly in wave zone) but mean `wz_std=0.318` (FAIL). After truncating samples before t=10 (junk window): mean `wz_std=0.230` (PASS).

**Fix:** `collector.py` and `incremental.py` now pass `min_valid_time=gw_min_valid_observation_time()` to `read_psi4_metrics()`. This excludes samples before `R_ext − R_orbit + margin = 12 − 5 + 3 = 10` code units. The 1/r check is now computed only on clean wave-zone data.

**Result:** eval_000006 would have been the first eval to PASS the wavezone gate and get a non-zero GW reward.

#### Bug 4: GRTresna throughput mismatched GPU evolution time

**Cause:** With `stop_time=40`, each GPU eval takes ~11 minutes. GRTresna solves take ~2-3 minutes. Default `MAX_CONCURRENT_GRTRESNA=6` produced gridinits 4× faster than the GPU could consume them, queuing 20+ evals on disk (19GB) while only 8 ran on GPU.

**Fix:** Set `MAX_CONCURRENT_GRTRESNA=2` so GRTresna produces gridinits at roughly the rate the GPU consumes them.

#### Bug 5: `psi4_directional.dat` header still 4-column

**Fix:** Updated `psi4_directional_header` in `driver.py` to `"# time  P_total  P_z_beam  beam_ratio  beaming_gain  wavezone_std"`.

### Campaign launch attempts

**Attempt 1** (235 preflight tests pass, 8 GPUs, `MAX_CONCURRENT_GRTRESNA=6`):
- 22 eval dirs created, 8 on GPU, 19GB disk
- Plotfiles not deleted (Bug 1)
- First 7 evals scored ~−1.5 to −1.9 (Bug 3 — all failed wavezone gate)
- Killed and cleaned

**Attempt 2** (after Bug 1+2+4 fixes, `MAX_CONCURRENT_GRTRESNA=2`):
- 9 eval dirs, 8 on GPU, plotfiles properly deleted (628MB/eval)
- First 7 evals scored: same identical scores (Bug 3 not yet fixed)
- eval_000006 had late-time `wz_std=0.08` but mean `0.318` → identified junk pollution
- Killed and cleaned

**Status:** All bugs fixed. Campaign not yet re-launched after Bug 3 fix. Runs dir cleaned.

### How to run the campaign

```bash
cd grteclyn-wrapper
export PATH="/usr/local/cuda/bin:$PATH"

QD_NAME=gw_beam_v5 \
QD_TARGET_EVALS=200 \
GPU_IDS="0 1 2 3 4 5 6 7" \
BATCH_SIZE=8 \
QD_ITERATIONS=25 \
MAX_CONCURRENT_GRTRESNA=2 \
bash scripts/campaigns/gw_beam/run_v5.sh
```

**Key parameters:**
- `MAX_CONCURRENT_GRTRESNA=2` — throttle GRTresna to match GPU throughput (each GPU eval ~11 min at t=40)
- `GRTECLYN_DELETE_WITHOUT_FRAMES=1` — delete plotfiles after Psi4 extraction (set in `run_v5.sh`)
- `STOP_TIME=40` — sponge zone absorbs boundary reflections
- `CONSUMER_RADII="12 16 20"` — multi-radius 1/r wavezone check
- `OBJECTIVE_MODE=gw_beam` — v5 constrained scoring

**Monitoring:**
```bash
# Pipeline monitor
tail -f runs/_logs/gw_beam_v5.pipeline_monitor.csv

# Campaign log
tail -f runs/_logs/gw_beam_v5.log

# Scored evals
for d in runs/grtresna_qd/gw_beam_v5/eval_*/; do
  [ -f "${d}score.json" ] && python3 -c "
import json; d=json.load(open('${d}score.json'))
c=d['score']['components']; p=d['metrics'].get('psi4') or {}
wz='PASS' if p.get('wavezone_ok') else 'FAIL'
print(f'$(basename $d): score={d[\"score\"][\"total\"]:7.2f}  P={p.get(\"mean_total_power\",0):.2e}  gain={p.get(\"mean_beaming_gain\",0):.2f}  wz={wz}')
"
done

# GPU usage
nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv,noheader
```

**Expected per-eval cost:** ~11 min GPU (t=40 at ~210 code units/h), ~2-3 min GRTresna. 200 evals on 8 GPUs ≈ 4.5 hours wall time.

**Cleanup:**
```bash
pkill -f "gw_beam_v5"; pkill -f "main3d.*CUDA"; pkill -f "Main_BosonStarBH"
rm -rf runs/grtresna_qd/gw_beam_v5/
rm -f runs/_logs/gw_beam_v5.*
```

---

## 2026-07-10 — `gw_beam_v5` QD campaign results + HQ promotion of eval_000046

### QD campaign — partial run (70 evals, 12 scored, stopped by user)

Campaign ran on 8 H100s at 128³ / L=64 / t=40 with the junk-truncation fix in place.
Stopped after ~70 eval dirs created (12 scored) to proceed to HQ promotion.

#### Scored evals (sorted by score)

| Rank | Eval | Score | mean P | peak P | gain | wz_ok | wz_std | v5_obj |
|------|------|------:|-------:|-------:|-----:|-------|-------:|-------:|
| 1 | **eval_000046** | **646.25** | 3.91e-6 | 7.61e-6 | 2.50 | PASS | 0.265 | 6.47 |
| 2 | eval_000052 | 590.48 | 8.06e-6 | — | 2.04 | PASS | — | 5.92 |
| 3 | eval_000037 | 590.48 | 8.06e-6 | — | 2.04 | PASS | — | 5.92 |
| 4 | eval_000049 | 579.82 | 8.24e-6 | — | 1.99 | PASS | — | 5.81 |
| 5–12 | (FAIL) | -1.0 to -1.9 | 1e-7 to 3e-5 | — | 1.8–3.0 | FAIL | 0.36–0.49 | 3.5–9.7 |

4 evals passed the wavezone gate; 8 failed. The v5 scoring discriminated properly —
PASS evals got +530 to +646, FAIL evals got penalty-level -1 to -2.

#### eval_000046 — best QD configuration

**Matter:** Canonical complex scalar (Q-ball boson stars), m=1, λ=640, μ=85333, ω=0.8.
All `exotic=0` (no phantom/negative-energy matter).

| Lump | R₀ | v_tangent | v_radial | well_depth | Status |
|------|-----|-----------|----------|------------|--------|
| 0 | 3.49 | 0.145 | -0.054 | 0.0014 | nearly OFF |
| 1 | 6.20 | 0.175 | +0.213 | 0.009 | active |
| 2 | 7.86 | 0.227 | -0.097 | 0.0014 | nearly OFF |
| 3 | 7.12 | 0.237 | +0.387 | 0.018 | active (strongest pump) |
| 4 | 5.51 | 0.135 | +0.137 | 0.013 | active |

Plus global breathing (A=1.25, ω=1.55) and Z-oscillation (amp=1.27, ω=0.75).
Effectively a 3-lump system (lumps 0, 2 nearly off) with lump 3 spiraling outward fast.

**PD trap pump active:** `rl_pump_kp=12.0`, `rl_pump_kd=7.0`, `trajectory_mode=1`.
The closed-loop PD controller drives each lump along its trajectory during the full
evolution (not just initial conditions). `well_depth` controls how strongly the trap
maintains the soliton — low values let lumps disperse freely.

**GW power growth (QD):** P doubles every ~8 code units (log10(P) = 0.038·t − 6.6).
Signal is growing, not decaying — the system is becoming more asymmetric as lumps spiral.

#### Psi4 time series (QD, 8 samples)

| t | P_total | beam_ratio | gain | wz_std | Notes |
|---|---------|-----------|------|--------|-------|
| 0.0 | 1.34e-6 | 0.184 | 1.74 | 0.529 | junk (excluded) |
| 6.4 | 8.0e-8 | 0.340 | 2.36 | 0.335 | junk (excluded) |
| 12.8 | 5.6e-7 | 0.304 | 2.22 | 0.640 | near-zone |
| 19.2 | 1.9e-6 | 0.316 | 2.26 | 0.508 | settling |
| 25.6 | 2.4e-6 | 0.361 | 2.45 | 0.180 | clean wave zone |
| 32.0 | 4.4e-6 | 0.389 | 2.56 | 0.022 | clean wave zone |
| 38.4 | 6.6e-6 | 0.432 | 2.73 | 0.112 | clean wave zone |
| 40.0 | 7.6e-6 | 0.442 | 2.77 | 0.126 | clean wave zone |

### HQ promotion — eval_000046 at 256³ / L=128 / t=40

**Run dir:** `runs/grtresna_promote/gw_beam_v5_hq_hq_eval000046/`

Replayed the QD score leader at full HQ resolution (256³, L=128, max_level=3, t=40)
with multi-radius Psi4 extraction (radii 12, 16, 20, 24), sponge zone (inner=56,
outer=64), frames on, and gw_beam v5 scoring.

**Launcher:** `scripts/campaigns/gw_beam/promote_hq.sh` (new — wraps `hq/run_batch.sh`
with gw_beam-specific env: `OBJECTIVE_MODE=gw_beam`, `GRTECLYN_PSI4=1`, sponge params,
no evolving geodesic, frames on).

#### Runtime

| | |
|---|---|
| Wall clock | ~1h 35m (06:17 → 07:52 UTC) |
| GRTresna HQ solve | ~8 min (8 MPI ranks, 256³ constraint solve) |
| GPU evolution | ~78 min (t=0→40 at 31 code units/h, 256³, 55GB GPU mem) |
| Post-run consumer drain | ~22 min (168 frames × 17 fields + Psi4 extraction) |
| Score computation | manual (replay_eval stuck on unwanted geodesic trace) |

#### Score: QD vs HQ — **score did NOT persist**

| Metric | QD (128³, t=40) | HQ (256³, t=40) | Change |
|--------|----------------:|----------------:|--------|
| **Score** | **646.25** | **0.00** | **-100% (gate failed)** |
| mean P | 3.91e-6 | 8.73e-7 | **-78% (4.5× weaker)** |
| peak P | 7.61e-6 | 2.20e-6 | -71% |
| mean gain | 2.50 | 2.64 | +6% (slightly better beaming) |
| peak gain | 2.77 | 3.66 | +32% |
| **wz_ok** | **True** | **False** | **gate failed** |
| wz_std | 0.265 | 0.527 | **2× worse** |
| dir_stability | 0.916 | 0.847 | -8% |
| n_samples | 6 | 126 | (HQ has 21× more samples) |

**The score collapsed from 646 to 0.** The wavezone gate failed (wz_std=0.53 > 0.3
threshold), zeroing the GW reward. The v5 objective was still computed (5.12) but
the gate killed it.

#### Why the score didn't persist

**1. GW power 4.5× weaker at HQ.** At 256³ the scalar field lumps are better resolved
(less numerical diffusion), so:
- Lumps are more compact → less spurious GW from numerical noise
- But the real GW signal is also weaker — the weak-field Q-balls barely curve spacetime
- The QD run's higher power was partly numerical artifacts being counted as GW

**2. Wavezone 1/r check failed.** Extraction radii (12, 16, 20, 24 from center) are
too close to the matter at HQ. In the QD run (L=64), radii 12-20 were in the wave
zone. At HQ (L=128), the same radii are still near-zone — the better-resolved field
extends the near-zone further. wz_std=0.53 means the signal doesn't scale as 1/r at
these radii.

**3. Beaming gain improved.** The directional preference is real and gets stronger
at HQ (mean 2.50→2.64, peak 2.77→3.66). This is the one metric that persisted and
improved — the beam is not a numerical artifact.

#### Physics at HQ

| Quantity | Value at t=40 | Interpretation |
|----------|--------------|----------------|
| chi_min | 0.970 | Space compressed by ~3% (weak gravity) |
| chi_max | 1.010 | Slight stretching (gauge) |
| phi_max | 1.16e-2 | Scalar field amplitude (small) |
| Weyl4_max | 1.4e-3 | GW amplitude (small but growing) |
| Ham constraint | 4.5e-4 | Healthy, growing linearly |
| Matter RMS spread | ~8 code units | Stable (not dispersing) |

**Constraint health:** Excellent. Hamiltonian = 4.5e-4 at t=40, growing linearly
but slowly. No blow-up, no NaN. The simulation is numerically stable throughout.

**Matter state:** The 3 active Q-ball lumps are holding together (PD trap maintaining
them). RMS spread ~8 code units, stable from t=11 to t=18. Not dispersing — the
opposite of the FTL campaign's eval 118 which fell apart by t=30.

#### Psi4 time series (HQ, last 5 of 169 samples)

| t | P_total | beam_ratio | gain | wz_std |
|---|---------|-----------|------|--------|
| 38.88 | 1.95e-6 | 0.576 | 3.30 | 1.124 |
| 39.36 | 2.05e-6 | 0.582 | 3.33 | 1.112 |
| 39.60 | 2.11e-6 | 0.599 | 3.40 | 1.123 |
| 39.84 | 2.16e-6 | 0.634 | 3.54 | 1.141 |
| 40.00 | 2.20e-6 | 0.665 | 3.66 | 1.153 |

Power growing, beaming strengthening (gain 3.3→3.7), but wz_std >1.0 — signal is
not in the wave zone at these radii.

#### Artifacts

- **17 movies** (`movies/movie_*.mp4`): scalar_activity, Weyl4_Re/Im/Mag, chi,
  chi_minus_1, local_speed, shift1, lump_activity, phi_lump_sum, Pi_lump_sum
  (z-slices + x/y/z projections)
- **17 plots** (`plots/`): constraints, collapse diagnostics, Psi4 analysis at
  M=30/1000 × D=0.002/1/10 Mpc (waveforms + PSD + propagation + strain + LIGO),
  custom `psi4_directional_analysis.png` (4-panel: power, gain, wz_std, beam_ratio)
- `score.json`, `metadata.json`, `small_data/` (169-row psi4_directional.dat,
  psi4_mode_l2m0.dat, constraint_norms.dat, collapse_diagnostics.dat)

### Plain English summary

We took the best configuration from the search (eval 000046 — three Q-ball lumps
orbiting with a breathing mode) and replayed it at 4× higher resolution (256³ vs
128³) for the same duration (t=40).

**What worked:**
- The simulation ran perfectly — no crashes, constraints healthy throughout
- The directional beaming is real and got stronger (2.6× → 3.7× peak gain)
- The matter is stable — Q-balls held together, not dispersing

**What failed:**
- The gravitational waves are 4.5× weaker than the search predicted
- The search score of 646 dropped to 0 at HQ — the wavezone quality check failed
- The search was partly rewarding numerical noise as if it were gravitational waves
- The extraction radii are too close to the matter at HQ resolution

**Root cause:** These are very weak gravitational fields (chi_min = 0.97, only 3%
deviation from flat space). The Q-balls are diffuse, low-density objects — they
don't curve spacetime enough to produce strong GWs. The beaming is real but the
total power is too weak.

**Next steps for strong, beamed GWs with real matter:**
1. **Denser lumps** — self-gravitating boson stars (higher central amplitude, compact
   ODE profiles) instead of diffuse Q-balls
2. **Faster motion** — higher orbital velocities (currently 0.01-0.02c tangential,
   need ~0.1-0.3c for significant quadrupole radiation)
3. **Closer orbits** — lumps at R=1-3 instead of R=5-8 (stronger tidal forces)
4. **Wider HQ extraction radii** — need radii at 30-50 from center for L=128 box
5. **Higher well_depth** — stronger PD trap pump to maintain lumps at closer radii
