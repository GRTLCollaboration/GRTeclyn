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
cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/grteclyn-wrapper
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
tail -f /home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/_logs/gw_beam_v5.pipeline_monitor.csv

# Campaign log
tail -f /home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/_logs/gw_beam_v5.log

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
rm -rf /home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/grtresna_qd/gw_beam_v5/
rm -f /home/jovyan/nachevsky/test/simulation/GRTeclyn/runs/_logs/gw_beam_v5.*
```
