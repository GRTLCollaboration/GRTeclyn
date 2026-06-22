# Spacetime Splash — Lab Journal

Reverse-chronological log of experiments, decisions, and results for the bosonic-shell critical-collapse campaign.  Newest entries first.

---

## Current state (2026-06-23)

### What works

- Bosonic shell search space with canonical coupling, 5-lump superposition
- Geometric splash scoring (χ, K, GW wave proxy) as primary reward, matter (ρ) secondary
- GW strain proxy from native `A_ij` plot vars — no Weyl4 derive on splash path
- Scoring fixes: initial-relative χ-drop, horizon bonus gated on center well, scalar wave floor removed from objective
- **V13 scoring validated** — survival forgiveness (geometric only), `GW_WAVE_PEAK_TARGET=0.20`, late GW gate t≥2, frames off
- Offline rescore: `grteclyn-wrapper/scripts/search/rescore_splash_campaign.py`
- Inward-converging radial velocity with sphere layout (layout=0)
- NFS-resilient postload gate with retry loop
- AMR regrid threshold guard prevents full-domain refinement segfault
- Early-terminated runs scored correctly (`gpu_ok`)
- **`spacetime_splash_v13` complete** — 125/200 gpu_ok, elite eval 020 @ 1684 (see entries below)
- **`spacetime_splash_v14_moving` complete** — 27/50 gpu_ok, best eval 015 @ 632; moving shells underperform v13 static (see entry below)
- **`spacetime_splash_v12` complete** — 31/50 gpu_ok, 0 segfaults (see entry below)

### Outstanding

- **v14 follow-up** — pin `v_rad ∈ [−0.2, 0]` and/or `matter_layout=0` (sphere) to enforce inward convergence; channel layout + uniform velocity still breaks splash physics
- **ham_abs NaN** at origin → −50 constraint penalty on all gpu_ok (deferred)

### Branches

- **GRTeclyn:** `feature/interstellar`
- **GRTresna:** `feature/interstellar`

---

## 2026-06-23 — spacetime_splash_v14_moving final results

Campaign complete (`50/50` evals, exit 0). Run dir: `runs/grtresna_qd/spacetime_splash_v14_moving/`. First MAP-Elites pass with **moving shells** (`shell_static=0`, 4 velocity dims active).

### Outcome mix (50/50)

| status | count | notes |
|--------|-------|-------|
| `gpu_ok` | 27 | 54% success rate (vs 61% v13) |
| `postload_rejected` | 16 | Ham L2 > 0.03 |
| `grtresna_rejected` | 7 | Mom/Ham convergence — moving IC stress |

No segfaults. MAP-Elites: **9 elites**, 14% bin coverage (`8×8`), best score **632** (eval 015).

### Top gpu_ok (v13 scoring unchanged)

| eval | score | χ-well | GW | \|K\| | survival | v_rad | layout | physics note |
|------|-------|--------|-----|-------|----------|-------|--------|--------------|
| 015 | 632 | 0.43 | 0.27 | 0.61 | 0.10 | **+0.16** | channel (~1) | **#1** — modest χ-well only; early stop t≈1.6 |
| 010 | 150 | 0.37 | 0.30 | — | 0.25 | **+0.17** | bipolar (~2.7) | best GW among mid-tier |
| 050 | 50 | 0.21 | 0.14 | — | 0.15 | +0.15 | sphere (~0.25) | weak geometry |
| 018 | 47 | 0.31 | 0.13 | — | 0.11 | +0.19 | channel | outward velocity |
| 046 | 43 | 0.21 | 0.29 | — | 0.10 | +0.15 | channel | GW-only signal |

Only **1/27** gpu_ok has χ-well > 0.4 (forgiveness threshold). GW spread: 0.00–0.30 (median 0.09) — far below v13 elite range (eval 042 GW=0.75). **No lapse collapse** on any scored run (`collapse_lapse_progress=0`).

### Velocity exploration — wrong direction

Among 27 gpu_ok, radial velocity is overwhelmingly **outward**:

| metric | value |
|--------|-------|
| v_rad range | −0.185 … +0.192 |
| v_rad mean | +0.074 |
| inward (v_rad < 0) | **3/27** (evals 012, 034, 043; scores ≤ 6) |

Top elites all have **positive v_rad** (+0.15–+0.19) despite default seed −0.2 inward. MAP-Elites selected configs that disperse matter early (low survival ~0.10) rather than converge — scorer still rewards modest χ-drop + GW proxy without deep collapse.

Channel layout (`matter_layout≈1`) dominates elites; known issue: uniform velocity on all lumps → net linear momentum, not radial convergence (see 2026-06-20 entry).

### Comparison vs v13 static

| | v13 (200 evals) | v14 moving (50 evals) |
|---|-----------------|----------------------|
| gpu_ok rate | 61% (125/204) | 54% (27/50) |
| Best score | **1684** (eval 020) | **632** (eval 015) |
| Best GW | 0.75 (eval 042) | 0.30 (eval 010) |
| χ-well > 0.4 | 26/125 | 1/27 |
| MAP-Elites elites | 23 | 9 |
| Deep collapse + lapse | yes (020 family) | **none** |

**Conclusion:** unpinned moving shells with ±0.2 velocity bounds did **not** improve toward ingoing GW splash. Kinetic stress raised Mom rejections and early termination; QD drifted to outward velocities and shallow geometry. Static v13 IC (gravity-driven accretion) remains the stronger collapse path under current scoring.

### Score composition (eval 015)

Geometric terms dominate but at modest saturation: χ-well 0.43 → ~345 pts (× survival 0.10 without forgiveness), GW 0.27 → ~41 pts, |K| 0.61 → ~55 pts. No lapse progress, no matter ρ credit, heavy instability penalty. Total 632 reflects **partial geometry before early stop**, not v13-class deep well.

### Recommended next steps

1. Pin **`grtresna_matter_layout=0`** (sphere) so radial velocity is per-lump inward.
2. Restrict **`grtresna_shell_radial_velocity` to [−0.2, 0]** — forbid outward radial motion.
3. Optional: raise Mom gate tolerance or reduce `shell_amp` cap for moving IC stability.
4. HQ frames on best **moving** elite once inward convergence is confirmed.

Retained eval dirs: `eval_000010`, `eval_000015`, `eval_000050`.

Launch command (reference):
```bash
SPLASH_MOVING=1 QD_NAME=spacetime_splash_v14_moving QD_TARGET_EVALS=50 \
  GPU_IDS="0 1 2 3 4 5 6 7" BATCH_SIZE=8 \
  bash scripts/campaigns/splash/run.sh
```

---

## 2026-06-23 — eval 020 deep dive (v13 #1)

Champion genome from `spacetime_splash_v13` (eval 020 / 054 identical, score **1684**). HQ promotion to N=256, t=30 was **stopped** after consumer crash (missing `weyl4` frame field — fixed in `fields.py`); promote dir removed.

### How the score is built

| Term | pts (approx) | Component | Notes |
|------|-------------|-----------|-------|
| χ-well | 800 | 1.00 | survival **forgiven** on geometry |
| Lapse progress | 200 | 1.00 | origin α: 1.0 → 0.018 |
| \|K\| crunch | 298 | 0.99 | max domain \|K\| ≈ 6.1 |
| GW arrival | 194 | 0.32 | late A_ij peak 0.065 (target 0.20) |
| ρ peak | 73 | 0.70 × surv 0.26 | matter secondary |
| Penalties | ≈ −80 | constraint, exotic (capped), instability, stationary | |

Without v13 survival forgiveness, geometric terms × 0.26 would bury this run (~400 total). High score is **geometry + lapse collapse**, not long runtime.

### Initial matter (canonical, not exotic)

- **Coupling:** canonical, `scalar_sign=+1`, all lumps `exotic=0` — no phantom/ghost matter in IC.
- **Layout:** channel (`matter_layout≈1`) — 5 lumps along a tilted axis (θ≈2.7, φ≈1.1), radii ~1.1–4.6 code units, **not** a hollow sphere at R=5.7.
- **Static IC:** `grtresna_shell_static=1` pinned in v13 → all lump velocities and ω **zero at t=0**. No pre-seeded ingoing shell.
- **Dynamics:** ρ at origin still rises (5×10⁻⁴ → 7×10⁻³ by t≈11); `peak_radius=1.0` — accretion via field evolution + gravity, not bulk initial momentum.

### Exotic matter: IC vs evolved

| | t=0 | evolved (t≈11) |
|---|-----|----------------|
| Seeded exotic? | **No** | — |
| min ρ_req | small positive | **−0.049** |
| ∫ negative ρ | 0 | **≈ 33** |
| WEC violation fraction | — | **≈ 97%** |
| `exotic_penalty` | — | −1.6 raw → **−0.3 capped** in `critical_collapse` |

Strong collapse **creates** effective negative-energy regions (NEC/WEC violations, negative ρ_req in the constraint solve). Scorer note “exotic matter required” refers to the **evolved geometry**, not phantom IC. Source of strong gravity: **self-gravitating canonical boson scalar** stress–energy → Einstein feedback (χ-drop, lapse crash, K spike).

### What “splash” means here vs target physics

**Scorer splash (eval 020):** domain-wide χ-well + lapse collapse + |K| crunch + modest late GW proxy + late ρ pile-up.

**Target splash (not fully achieved):** ingoing gravitational wave from **converging moving shells** focusing at the center. Eval **042** is closer (GW=0.75, low ρ); eval 020 wins on **geometry saturation + matter corroboration**.

**Implication for v14:** unpinned velocity with radial default **−0.2**, bounds **±0.2**, `shell_static=0` — test whether MAP-Elites finds true wave-focus elites without static late-ρ inflation.

---

## 2026-06-23 — spacetime_splash_v14_moving launch

50-eval MAP-Elites **moving-shell** test after v13 static campaign. **Completed** — see final results entry above.

**Changes vs v13:**
- `SPLASH_MOVING=1` → pin `grtresna_shell_static=0` (velocity dims active)
- Search adds 4D: `radial_velocity` **±0.2** (default −0.2 inward), `toroidal`/`poloidal` **±0.2**, `shell_omega` **±0.15**
- Same v13 scoring, frames off, 8 GPUs, canonical boson star

**HQ eval 020:** stopped; promote dir removed (consumer NFS stub from crashed `weyl4` frames — fixed in `fields.py`).

Run dir: `runs/grtresna_qd/spacetime_splash_v14_moving/`.

---

## 2026-06-23 — spacetime_splash_v13 final results

Campaign complete (`200/200` evals, exit 0). Run dir: `runs/grtresna_qd/spacetime_splash_v13/`. Launched at 50 evals, extended mid-run via `QD_RESUME=1 QD_TARGET_EVALS=200`.

### Outcome mix (204 trajectory rows)

| status | count | notes |
|--------|-------|-------|
| `gpu_ok` | 125 | 61% success rate |
| `postload_rejected` | 62 | Ham L2 > 0.03 |
| `grtresna_rejected` | 7 | Mom/Ham convergence |
| `grtresna_failed` | 3 | solver failures |
| `gpu_failed` | 3 | resume handoff artifacts |
| `pipeline_interrupted` | 4 | resume handoff artifacts |

No segfaults. MAP-Elites: **23 elites**, 36% bin coverage (`8×8`), best score 1684.

### Top gpu_ok (v13 scoring)

| eval | score | χ-well | GW | \|K\| | ρ | survival | physics note |
|------|-------|--------|-----|-------|---|----------|--------------|
| 020 / 054 | 1684 | 1.00 | 0.32 | 0.99 | 0.70 | 0.26 | **#1** — identical genome; deep domain χ-drop + full lapse collapse |
| 102 | 1607 | 1.00 | 0.39 | 0.88 | 0.54 | 0.32 | same archetype; shallower domain min χ, better GW |
| 042 | 1388 | 0.73 | **0.75** | 1.00 | 0.25 | 0.31 | **best GW focus** — pure geometric splash, low matter |
| 011 | 1292 | 0.72 | 0.68 | 0.97 | 0.15 | 0.41 | balanced geometry; minimal matter pile-up |
| 029 | 1200 | 0.70 | 0.65 | 1.00 | 0.21 | 0.10 | violent crunch, survival forgiven |
| 093 | 1176 | 0.73 | 0.57 | 0.86 | 0.13 | 0.28 | geometry-heavy, low ρ |
| 008 | 417 | 0.00 | 0.44 | 0.60 | — | 1.00 | **demoted dud** — full survival, no χ-well |

26/125 gpu_ok have χ-well > 0.4 (survival-forgiveness eligible). GW spread among scored runs: 0.01–0.78 (median 0.21) — recalibration restored dynamic range.

### Analysis of top elites

**Two distinct archetypes emerged:**

1. **Deep-well + matter corroboration (eval 020 family).** Eval 020 and 054 are bit-identical genomes (QD reproduced the elite). Large outer shell (`radius≈5.7`, high radial jitter 0.86, `amp≈0.076`). Domain min χ drops to 0.29 with max |K|≈6.1; central ball-averaged χ actually *rises* early (0.87→1.11) before late collapse — the well metric captures the domain-wide drop, not the center monotonicity. Lapse at origin hits 0.018 (full `collapse_lapse_progress=1.0`). Matter ρ peaks late (t≈11.2) and contributes ~280 score points (ρ=0.70 × survival 0.26 × 400 weight); geometric terms get survival forgiveness. θ⁺≤0 detected (horizon proxy r≈0.65) but **horizon bonus suppressed** — no lapse-collapse corroboration at the same timestep.

2. **GW-focused geometric splash (eval 042).** Best late GW signal in the campaign (0.75 of target 0.20). Full |K| saturation with moderate χ-well (0.73) and deliberately low matter (ρ=0.25). Thicker shell, smaller radius (4.8), higher shell mode (0.87). Score ~300 pts below #1 because it lacks full χ-well saturation and matter corroboration — but **closest to the target physics** (converging wave → center crunch without matter pile-up).

**Score composition (eval 020):** geometric terms ≈1290 pts (forgiven survival) + lapse progress 200 + matter ≈73 (raw survival) − constraint/exotic/instability ≈−80 → 1684. High-ρ matter is a tiebreaker among deep collapses, not the primary signal.

**Contrast with v12:** v12 #1 was eval 044 (score 1020, high survival, modest geometry). Under v13 rescoring, eval 016-class splashes rank at the top; v13 campaign confirms this ordering in live search. Eval 008-style duds (survival 1.0, χ-well 0) now score ≤417.

### What v13 scoring validated

- Geometric-only survival forgiveness lifts collapse elites without inflating matter pile-ups
- Late GW gate (t≥2) prevents early A_ij noise from dominating; eval 042's peak is genuinely late
- Horizon bonus gate: zero horizon bonuses across all 125 gpu_ok despite θ⁺ detections
- Static-shell search sufficient to find deep collapses; moving shells not yet explored

### Still missing

- Frame/visual confirmation of ingoing wave convergence for eval 042 (eval dirs pruned to top 3: 020, 054, 102)
- True apparent-horizon credit with lapse corroboration
- `ham_abs` NaN at origin → universal −50 constraint penalty
- Moving-shell campaign with tightened velocity bounds (±0.2)

Retained eval dirs: `eval_000020`, `eval_000054`, `eval_000102`.

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

**Extended (2026-06-23):** resumed from trajectory with `QD_RESUME=1 QD_TARGET_EVALS=200` after initial 50-eval batch. Completed 200/200 evals — see final results entry above.

Run dir: `runs/grtresna_qd/spacetime_splash_v13/`.

---

## 2026-06-18 — spacetime_splash_v12 final results

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
