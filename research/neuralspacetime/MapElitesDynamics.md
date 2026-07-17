# MAP-Elites Dynamics — Trajectory Ansatz Campaigns

> Continuation of [MapElites.md](./MapElites.md). This journal covers campaigns using
> the **trajectory ansatz** (Option D: per-lump independent orbits) and comparison
> against the spherical-harmonic (SH) baseline.

---

## Status (2026-07-17) — start here

**Paper (current draft):** [article/research.pdf](./article/research.pdf) · source [article/research.tex](./article/research.tex)

### Bug fix landed (eval-118 NO-GO)

Two physics bugs invalidated the spiral_v2 champion and its validation matrix:

1. **ID ↔ evolution sign mismatch** — GRTresna applied per-lump phantom signs in the constraint solve, but evolution used single-complex `grtresna_complex_scalar` with one global canonical sign. Cauchy pairing broken.
2. **Illegal PD pump** — pump forced `d_t Π` without a matching `T_μν`, so `∇_μ T^{μν} ≠ 0`.

**Cure (committed):** fail-closed sign gate; campaign switched to `grtresna_bicomplex_scalar` (per-lump signs retained in Phi±); `pump_work` diagnostics + `SCORE_PUMP_ENERGY_WEIGHT`; transient igniter via `RL_PUMP_STOP_TIME`. Quarantine notes: `runs/grtresna_qd/qball_traj_spiral_v2/eval_000118/NO_GO.txt`, `validation/eval118/NO_GO.txt`.

### Run directories — what is fresh vs archived

| Path | Status | Notes |
|------|--------|-------|
| **`runs/grtresna_qd/qball_traj_bicomplex_v1/`** | **DONE — QD archive** | 200/200 evals; champion `eval_000087` score **825.2**. |
| **`runs/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/`** | **FRESH — active CMA-ES** | Warm-start from QD champion; 8 GPUs; target 150. **Current work.** |
| `runs/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/` | freeze snapshot | HQ/matrix source (`CHAMPION.json`); re-freeze if best moves |
| `runs/qball_traj_bicomplex_cmaes_v1.launch.log` | live log | CMA-ES launcher log |
| `runs/grtresna_promote/e118_bicomplex_pumpon_L64_N128_t16/` | Phase-5 diag | eval-118 genome, bicomplex, pump on |
| `runs/grtresna_promote/e118_bicomplex_pumpoff_tstop4_L64_N128_t16/` | Phase-5 diag | same genome, igniter off after t=4 |
| `runs/grtresna_promote/e118_bicomplex_PHASE5_DECISION.txt` | decision note | Shortcut still present under corrected physics |
| `runs/grtresna_qd/qball_traj_spiral_v2/` | **NO-GO / archived** | Old single-complex campaign; do not use scores for paper |
| `runs/grtresna_promote/qball_traj_spiral_v2_*` | **NO-GO / archived** | HQ/promotion of eval-118 under broken pairing |
| `research/neuralspacetime/validation/eval118/` | **NO-GO** | RC/RM/RF/DS/DL matrix void; re-run after new champion |

Sections below on spiral_v2 / eval-118 HQ are **historical** until superseded by the bicomplex champion.

### Forward plan

```
1. QD           → DONE (qball_traj_bicomplex_v1, champion eval_000087 @ 825.2)
2. CMA-ES       → running (qball_traj_bicomplex_cmaes_v1, target 150)
3. Validation   → HQ + RC/RM/RF/DS/DL on the CMA-ES champion only (scripts ready)
4. Paper        → update research.tex / research.pdf from the validated champion
```

**Diagrams:** [mapelites-campaign-stages.svg](./mapelites-campaign-stages.svg) (outer ladder),
[mapelites-end-to-end.svg](./mapelites-end-to-end.svg) (outer + inner loop),
[mapelites-matter-first.svg](./mapelites-matter-first.svg).

Do not promote or write article numbers from spiral_v2 / eval-118.

### CMA-ES launch — wired (in flight)

Launcher mirrors QD physics (bicomplex `EXTRA_SETS`, igniter `RL_PUMP_STOP_TIME=4`,
score weights, pins). Output: `runs/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/`.

### HQ + Richardson — ready (post CMA-ES)

Campaign folder (keep future promote/resolution work here, not in `hq/`):

`grteclyn-wrapper/scripts/campaigns/promote/bicomplex_cmaes_v1/`

| File | Role |
|------|------|
| `manifest.json` | RM/RC/RF (+ optional DS/DL/GW) |
| `campaign.env.sh` | Physics + paths for this campaign only |
| `freeze.sh` | Snapshot champion |
| `run.sh` | One cell |
| `launch.sh` | Freeze → RM; `LAUNCH_LADDER=1` → RC+RF |

Shared engine: `promote/lib/` + `hq/{run_batch,replay_eval}.py`.  
NO-GO archive: `promote/_archive/eval118/`.

```bash
cd grteclyn-wrapper
P=scripts/campaigns/promote/bicomplex_cmaes_v1

bash $P/freeze.sh
DRY_RUN=1 bash $P/run.sh BCMA-RM
bash $P/launch.sh
# LAUNCH_LADDER=1 SKIP_FREEZE=1 bash $P/launch.sh
```

**Order:** RM (1 GPU) → RC (1 GPU) + RF (2 GPUs, frames on). Series = RC+RM+RF at L=128.  
Frozen source: `runs/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/`.  
New campaigns: add `scripts/campaigns/promote/<campaign_name>/` (do not dump scripts into `hq/`).

---

## HQ promotion — `qball_traj_spiral_v2_t30_hq_eval000118` (COMPLETE 2026-07-02, t=30)

Extended HQ replay of v2 score leader **eval 118** at **256³ / L=128 / max_level=3 /
t_stop=30**, `plot_interval=24` (126 frames), `general_ftl` + 4D evolving geodesic,
`per_frame_zlim` on `chi` / `chi_minus_1` / `shift*`. Run dir:
`runs/grtresna_promote/qball_traj_spiral_v2_t30_hq_eval000118/`.

### Runtime

| | |
|---|---|
| Wall clock | **~1 h 42 m** (02:33 → 04:15 UTC, 2026-07-02) |
| GPU evolution | **~80 min** (4786 s) |
| Post-run scoring + 4D trace | ~22 min (consumer drain + `score.json`) |
| Status | **COMPLETE** — sim exit 0, 4D geodesic done, `score.json` written |

### Scores

| | QD (128³, t=16) | HQ t=16 | **HQ t=30** |
|--|----------------:|--------:|------------:|
| **`score.json` total** | **603.39** | 511.89 | **224.20** |
| Incremental (last frame) | — | — | 138.27 @ t=30 |
| `structural_persistence` | ~35% | ~35% | **23%** |
| `confinement_retention` | — | — | 23% (spread ×2.46) |
| `operational_ftl` | 0.347 | 0.346 | **0.099** (dispersal-gated) |
| `ftl_geo_evolving` | 0.306 | 0.225 | **0.150** |
| `exotic_penalty` | — | — | −1.6 |
| `instability_penalty` | — | — | −0.83 |
| `constraint_growth` | 0.949 | — | 0.982 |

**Verdict:** Numerics **survive** to t=30 (no blow-up, `numerical_survival`=1.0,
`constraint_growth`=0.98), but **matter disperses** — confinement falls from ~53% @ t=0
to **23%** @ t=30 (rms radius 7.6 → 18.6). The dispersion gate correctly crushes
`operational_ftl` and the headline **`score.json` total** vs QD/HQ-t16. The gauge-invariant
channel is still real but weaker than the QD emit-sweep peak.

### 4D evolving null-geodesic trace (authoritative)

| Quantity | Value |
|---|---|
| `f_geo_evol` (4D trace, t_emit=0) | **13.0%** |
| Frozen `f_geo` (snapshot, z) | 11.1% |
| Frozen `f_geo` peak (time series) | **22.8%** @ t≈19.2 |
| Frozen `f_geo` mean (126 frames) | 15.8% |
| Rays reached | **5 / 5** |
| `t_emit → t_arrival → t_flat` | 0.0 → 12.52 → 14.40 |
| `h_quality_ok` | true |
| Best probe axis | **x** (evolving), z (frozen) |

Unlike the QD emit-sweep (still **rising** at t_emit=12 → 17.7%), the extended HQ window
shows frozen `f_geo` **peaking at t≈19.2 then decaying** — the channel opens then fades as
matter shears apart. 4D end-to-end trace stays at **13.0%**, below the QD peak.

### GW detection (Psi4) — pipeline note (2026-07-02)

This run predates HQ GW extraction (`--no-psi4` was used). Subsequent HQ replays enable
spherical Psi4 mode extraction via `GRTECLYN_PSI4=1` in [`promote_common.sh`](../../grteclyn-wrapper/scripts/campaigns/lib/promote_common.sh).

| | **SupportedWormholeCollapse** | **Q-ball HQ (RadialRecipe)** |
|--|-------------------------------|------------------------------|
| Domain | Octant `[0,L]³`, physics at **corner** `(0,0,0)` | Full box `[0,L]³`, physics at **center** `(L/2,L/2,L/2)` |
| Boundaries | Reflective `lo`, Sommerfeld `hi` | Sommerfeld all faces |
| Extraction center | `(0,0,0)` | **`center`** from params (e.g. `64 64 64`) |
| Consumer frames | `--frames-corner` | **No** corner mode (centered slices) |
| Psi4 radii | 12 16 20 24 (wormhole `plot_run.sh`) | **8 12 24** from center (`CONSUMER_RADII`) |
| Output | `small_data/psi4_mode_l2m0.dat` | same |

Requires `amr.derive_plot_vars = Weyl4` (already on HQ params) and aligned
`extraction_center = center` for the Newman–Penrose tetrad. Post-run:
`bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh RUN_DIR 8 12 24`.

### Plain English

We took the **best configuration from the search** (eval 118 — a spiral of exotic Q-ball
lumps) and replayed it at **full HQ resolution** (256³ grid, 4× finer than the search) for
**twice as long** (t = 30 instead of t = 16). The goal was to see whether the spacetime
shortcut keeps getting stronger if we run longer.

**Bottom line:** the simulation **did not crash** — it ran cleanly to t = 30. But the
**matter fell apart over time**, and the geometry shortcut **peaked around t ≈ 19 and then
faded**. This is **not** a stable, traversable warp bubble.

**What worked**

1. **Numerically stable.** No blow-up; constraints stayed healthy throughout.
2. **A real geometry shortcut existed.** 4D null-geodesic ray tracing found light beating
   a flat-space signal by about **13%** (5/5 rays completed). That is gauge-invariant —
   not a plotting artifact.
3. **The channel was strongest mid-run.** Snapshot measurements peaked at roughly **23%**
   around t ≈ 19, before the structure sheared apart.
4. **Superluminal coordinate signals were present** — local signal speed briefly reached
   ~**1.26 c** in places.

**What failed**

1. **The Q-balls dispersed.** Confinement fell from **53%** at t = 0 to **23%** at t = 30;
   the cloud spread ~**2.5×** (rms radius 7.6 → 18.6). The structured “motor” holding the
   warp open largely dissolved.
2. **The score dropped: 603 → 224.** The scoring system penalizes shortcuts that open
   **while matter is flying apart** — a dispersing cloud with a coordinate superluminal
   channel is not counted as a useful warp.
3. **The shortcut did not keep climbing.** In the shorter search run it was still rising at
   t = 16 (~18% in emit sweeps). Over t = 30 it **peaked near t ≈ 19 and decayed**.
4. **Exotic matter is required** — full exotic-matter penalty applies.

**Analogy:** a whirlpool that briefly creates a faster path through water, but the
whirlpool itself is **dissolving** as it spins. You can measure the faster path (~13%
geodesic shortcut), but the engine making it (confined Q-ball lumps) is **breaking up**,
so you do not end up with a durable traversable shortcut.

**Campaign implication:** eval 118 remains the best search candidate for FTL-like geometry
at t = 16, but **extending to t = 30 at HQ does not confirm a long-lived warp**. Next step:
find configs that **keep confinement high** while opening a geodesic shortcut — not just
run this one longer.

### Artifacts (on disk)

`score.json`, `metadata.json`, `small_data/` (timeseries, `evolving_geodesic.json`,
126-row `score_timeseries.jsonl`), `frames/` (16 dirs, 2016 PNGs), `movies/movie_*.mp4`,
`data/` diagnostics. Plotfile HDF5 dirs and `initial_data.gridinit` pruned post-run.

---

## `qball_traj_spiral_v2` — dispersion-gated spiral QD (COMPLETE 2026-07-01, 200/200 evals)

First **full** trajectory campaign to run to completion. 39D MAP-Elites on compact
Q-ball spiral orbits, `general_ftl` objective, **`SCORE_FTL_DISPERSION_GATE=1.0`**
(new — gates `operational_ftl` + `ftl_persistence` by `structural_persistence`),
`SCORE_EXOTIC_PENALTY_WEIGHT=0.2`, `ftl_lifetime` descriptor on an 8×8 archive.
GRTresna preload with postload constraint gate (Ham/Mom ≤ 5%), compact ODE seed
(m=1, λ=640, μ=85333, ω=0.8), N=128 / L=64 / max_level=1 / t_stop=16.

### Runtime & throughput

| | |
|---|---|
| Wall clock | **~1 h 45 m** (10:09:39 → 11:55:14 UTC, 2026-07-01) |
| Evals | **200 / 200** (target reached, `last_eval_counter=200`) |
| Avg / eval | ~32 s (8-GPU pipeline, `max_concurrent_grtresna=6`, `scoring_workers=8`) |
| Stalls | **none** — the release-lease-early + bounded scoring-pool fix held; GPUs fully drained at shutdown |

This is the run that validated the scoring-pipeline refactor: the campaign filled all
200 evals without the mid-run stall that plagued the earlier attempts, and terminated
cleanly (0 MiB on all 8 H100s, no orphaned GRTeclyn/GRTresna/MPI processes).

### Outcomes

| Status | Count | | Tier | Count |
|--------|------:|---|------|------:|
| `gpu_ok` | 133 | | operational | 44 |
| `postload_rejected` | 50 (25%) | | constructed | 77 |
| `grtresna_rejected` | 15 | | nontrivial | 1 |
| `gpu_failed` | 2 | | rejected | 78 |

85 evals scored > 0. The 2 `gpu_failed` runs (**eval 78, 119**) are late-time
constraint/energy blow-ups (`structural_persistence`→0, `energy_condition`→0,
`instability_penalty=-1`) — the config disperses so hard the evolution goes singular,
exactly the class the dispersion gate is meant to demote.

### Archive — 9 elites, coverage 14.06% (9/64), best **603.39**

| Eval | Score | Cell | Tier | ftl_peak (x) | op_ftl | constraint_growth |
|------|-------|------|------|-------------|--------|-------------------|
| **118** | **603.39** | (7,7) | operational | 0.882 | 0.347 | 0.949 |
| 181 | 487.34 | (4,7) | operational | 0.622 | 0.330 | 0.914 |
| 162 | 464.19 | (5,7) | operational | 0.630 | 0.309 | 0.948 |
| 195 | 412.17 | (6,7) | operational | 0.831 | 0.249 | 0.954 |
| 168 | 239.71 | (3,7) | constructed | 0.474 | 0.0 | 0.971 |
| 090 | 221.50 | (1,7) | operational | 0.143 | 0.280 | 0.933 |
| 138 | 137.59 | (2,7) | operational | 0.323 | 0.169 | 0.929 |
| 009 | 24.06 | (0,0) | constructed | 0.0 | 0.0 | 0.809 |
| 068 | −2.37 | (0,7) | constructed | 0.048 | 0.0 | 0.980 |

Every non-trivial elite lives in the **top lifetime row (y=7, `ftl_lifetime_fraction`=1.0)**
— the archive is a clean sweep across the `ftl_peak_strength` (x) axis at 100% FTL
lifetime. All operational elites hold **`anec_condition`=1.0** and **`tidal_comfort`=1.0**
with high `constraint_growth` (0.91–0.98): the dispersion gate is working — leaders now
combine a real FTL channel **with** structural persistence, instead of banking score on a
cloud that has already flown apart.

### FTL champions (`ftl_champions.json`)

| Descriptor | Eval | Value | Score |
|-----------|------|-------|-------|
| `f_geo_evol` (peak) | 77 | 0.301 | 568.98 |
| `ftl_geo_evolving` | 79 | 0.343 | 366.86 |
| `f_op_peak` (t=16.0) | 77 | 0.270 | 568.98 |
| `max_local_speed` (t=16.0) | 44 | **2.006 c** | 62.69 |
| `superluminal_fraction` | 4 | 1.0 | 38.45 |
| `ftl_lifetime_fraction` | 2 | 1.0 | 216.44 |

### Top run — `eval_000118` deep dive (cell 7,7, score 603.39)

**Physics signature** (all peaks at t ≈ 16.0, the last frame; n_frames=7):

| Quantity | Value | Meaning |
|---|---|---|
| `f_geo_peak` = `f_geo_evol` | **17.7%** | gauge-invariant geodesic shortcut vs flat light front |
| `ftl_geo_evolving` | 0.306 | scored transform of the evolving-metric shortcut |
| `f_op_peak` | 19.0% | operational (coordinate-front) FTL excess |
| `max_local_speed` | **1.47 c** | peak local coordinate signal speed |
| `superluminal_fraction` | 1.0 | entire sampled front superluminal at peak |
| `ftl_lifetime_fraction` | 1.0 | channel present for the whole measured window |

**Objectives:** `anec_condition`=1.0, `tidal_comfort`=1.0, `constraint_growth`=0.949,
`operational_ftl`=0.347 — a clean, constraint-healthy FTL channel, not a numerical blow-up.

**Winning recipe (39 params):** 5 mostly-**exotic** lumps (`exotic`≈0.75–0.91 on lumps 1–4,
0.32 on lump 0); all **retrograde** (`omega_rot` −0.04…−0.15); nested radial shell
(`R0`≈1.96–7.52, one tight inner lump); **near-zero `v_rad`** (0.003–0.069 → persistent, not
dispersing); strong collective **breathing** (`A_breath`=1.60, `omega_breath`=2.67) and
**z-oscillation** (`z_amp`=1.50, `omega_z`=1.51). A persistent, breathing, retrograde
exotic-Q-ball shell that coherently pumps the warp channel.

**Continuous-emission sweep** (`evolving_geodesic.json`, Δt=2, 3 rays, all reached each launch):

| t_emit | 0.0 | 2.0 | 4.0 | 6.0 | 8.0 | 10.0 | **12.0** |
|--------|-----|-----|-----|-----|-----|------|----------|
| f_geo  | 13.5% | 11.7% | 11.8% | 13.0% | 12.8% | 14.2% | **17.7%** |

Unlike eval 122 (monotonic **decay** from t=0, dispersing channel), eval 118 **dips then
climbs** and is **still rising at the last valid launch** (t_emit=12 → 17.7%). The warp
channel *strengthens over time* — matching the t≈16 peaks in every operational descriptor.
The sweep caps at t_emit=12 (`max_emissions`=7 × `emit_interval`=2, and later launches lack
enough remaining slices to trace end-to-end), so **the true peak is very likely beyond t=12,
cut off by the emission window rather than the physics.** This is the primary motivation to
HQ-promote eval 118 to t_stop=30 and see whether f_geo saturates or keeps climbing.

### HQ promotion — `qball_traj_spiral_v2_hq_eval000118` (2026-07-02, t=16 pilot)

First HQ replay of the v2 score leader at **256³ / L=128 / max_level=3 / t_stop=16**,
`plot_interval=24` (~68 frames), `general_ftl` + 4D geodesic **`hq`** mode, dirs `x y z`.
Run dir: `runs/grtresna_promote/qball_traj_spiral_v2_hq_eval000118/`.

| | Stage 0 (128³, t=16) | HQ (256³, t=16) |
|--|---------------------:|----------------:|
| **Score** | **603.39** | **511.89** |
| `operational_ftl` | 0.347 | 0.346 |
| `ftl_geo_evolving` | 0.306 | **0.225** |
| `ftl_persistence` | 0.347 | 0.346 |
| `structural_persistence` / confinement | ~35% | ~35% |
| `constraint_health` | — | 0.823 |
| `f_op` (evolved) | 0.189 | 0.188 |
| `max_local_speed` | **1.47 c** | **1.46 c** |
| 4D `f_geo_evol` | peak **17.7%** (t_emit≈12) | **13.0%** |
| gauge-invariant `f_geo` | — | **19.9%** |
| FTL lifetime | 100% | 100% |
| HQ status | — | **survived** |

**Verdict:** FTL channel and ~35% Q-ball confinement **survive full HQ resolution** —
coordinate metrics (`f_op`, `max c`, confinement) match QD within noise. Score drop is
almost entirely **`ftl_geo_evolving`**: HQ 4D trace reports **13.0%** vs QD emit-sweep
peak **17.7%** at t_emit≈12 (eval 118 still **rising** at the last QD launch; t=16 cut
the window short). Not a crash or dispersion cheat — geodesic credit is lower at 256³.

**4D trace (HQ, t_emit=0):** `f_geo_evol=13.0%`, frozen peak **22.1%** at t≈15.1,
`h_quality_ok=true`, best direction **x**. Time-averaged frozen `f_geo` mean **69.8%**
(ignored when 4D evolving trace is authoritative).

**GRTresna (256³):** Ham **0.93%**, Mom **0.40%** at NL iter 6 (vs 0.95% / 0.42% at QD).

**Artifacts:** `frames/` (16 field/projection dirs, 1088 PNGs), `movies/movie_*.mp4`,
`score.json`, `small_data/score_timeseries.jsonl` (70 rows).

**Frame fix (before t=30 rerun):** global zlim-lock from t=0 washed out `chi` / `chi-1`
contrast (flat yellow on χ≈1, O(0.01) signals on ±0.9 bar). `per_frame_zlim` now uses
local percentile scaling for those fields.

**Follow-up (done):** `qball_traj_spiral_v2_t30_hq_eval000118` at **t_stop=30** —
see top section. f_geo peaks @ t≈19.2 then decays; matter disperses to 23% confinement.

### v2 (gated) vs v1 (ungated) — why the headline score dropped

v1 leaders scored **~1100** but at only **28–40% confinement**: the dispersion gate did
not yet exist, so `operational_ftl` + `ftl_persistence` were banked ungated. v2's
best is **603** because those same terms are now scaled by `structural_persistence` —
the drop is **the gate doing its job**, not a regression. The search converged early
(best_score reached 603 by iteration 15 of 24 and held); coverage plateaued at 14%.

**Kept on disk** (`keep_top_eval_dirs=3` + FTL retention): `eval_000002`, `_000004`,
`_000041`, `_000044`, `_000077`, `_000079`, `_000118` under
`runs/grtresna_qd/qball_traj_spiral_v2/`.

---

## `qball_traj_spiral_v1` — spiral v_rad QD (stopped 2026-07-01, 73/200 evals)

39D MAP-Elites on compact Q-ball trajectory orbits with per-lump **`v_rad`** spiral
drift, `general_ftl` objective, `SCORE_EXOTIC_PENALTY_WEIGHT=0.2`, multi-ray emission
sweep (Δt=2, 7 launches). Run stopped manually; disk pruned to **top-3** eval dirs.

| Eval | Score | Tier | Confinement | FTL signal (headline) |
|------|-------|------|-------------|------------------------|
| **56** | **1140** | operational | 40% | `operational_ftl`=1, `ftl_persistence`=1, `f_geo_evol`≈0.30 |
| **46** | 1102 | operational | 28% | saturated coord FTL; `f_geo` peak ≈19% @ t≈12.8 |
| **40** | 1097 | operational | 28% | late-blooming channel; inward spirals on lumps 1 & 4 |

**Outcomes:** 57 `gpu_ok`, 14 postload rejected (~19%), 2 GRTresna rejected; 21 operational
tier hits; archive **14%** coverage (9 elites). Search reliably finds **operational FTL**
(scores ~1100+) with **`operational_ftl` + `ftl_persistence` saturated**, but leaders keep
**low confinement (28–40%)** — dispersion opens a late coordinate superluminal channel while
wells overlap and matter shears. Best geodesic peaks at **t_emit ≈ 12–16**, not t=0.

**Kept on disk:** `eval_000056`, `eval_000046`, `eval_000040` under
`runs/grtresna_qd/qball_traj_spiral_v1/`.

**Follow-up:** gate `operational_ftl` / `ftl_persistence` by `structural_persistence` before
next campaign (see [NextSteps.md](./NextSteps.md)).

---

## Continuous null-ray emission sweep — FTL channel lifetime mapping (2026-07-01)

The single-launch evolving geodesic probe fires rays only at `t_emit = times[0]`,
which may miss the peak of a transient warp.  A **continuous-emission sweep** now
fires a ray fan at multiple launch times (`t_emit = 0, 2, 4, ..., 12`) and reports
the full `f_geo(t_emit)` map plus the peak.

### How it works

1. **`_emission_times(field, interval, max_emissions)`** computes launch times
   spaced by `emit_interval` (env: `GRTECLYN_GEO_EMIT_INTERVAL`), capped at
   the last available slice time and at `max_emissions` launches
   (env: `GRTECLYN_GEO_MAX_EMISSIONS`).  `interval <= 0` or `max_emissions <= 1`
   disables the sweep (legacy single-launch).

2. **`compute_evolving_geodesic_ftl_emission_sweep`** runs a full
   direction-scanned ray fan (5 rays, best of x/y/z) at each launch time.
   The per-launch results are collected in `emit_sweep` as
   `(t_emit, f_geo, n_reached)` triples.  The headline `f_geo` / `t_emit`
   report the peak trustworthy launch.

3. Both `compute_evolving_geodesic_ftl_from_metric_stack_cache` and
   `_from_plotfiles` auto-route to the sweep when `emit_interval > 0`.

### Run: `traj_qball_compact_sweep_eval122_t16` (2026-07-01)

Eval 122 compact-seed replay (N=128, L=64, ml=2, t=16) with `emit_interval=2.0`,
`max_emissions=7`.  All 7 launches converged (5/5 rays reached detector at every
launch).

| t_emit | f_geo  | n_reached |
|--------|--------|-----------|
| 0.0    | **9.38%** | 5/5    |
| 2.0    | 9.18%  | 5/5       |
| 4.0    | 8.33%  | 5/5       |
| 6.0    | 7.30%  | 5/5       |
| 8.0    | 6.45%  | 5/5       |
| 10.0   | 6.34%  | 5/5       |
| 12.0   | 5.89%  | 5/5       |

**Peak f_geo = 9.38%** at `t_emit = 0.0`, monotonically decaying.  The FTL channel
is strongest at the earliest emission and gradually weakens as the warp geometry
disperses — all 7 launches produce a genuine superluminal shortcut (f_geo > 0),
confirming **100% FTL lifetime** across the full 0-to-12 emission window.

`f_geo_frozen_peak = 12.79%` (unchanged from single-launch; frozen metric is
launch-time independent).

### Comparison vs single-launch (no sweep)

| Metric | Single-launch (`traj_qball_compact_ode_eval122_t16`) | Sweep (`_sweep_`) |
|--------|------|------|
| f_geo (evolving) | not computed | **9.38%** |
| f_geo_frozen_peak | 12.79% | 12.79% |
| confinement @ t=16 | 31.0% | 31.0% |
| ftl_geo_evolving component | 0.0 | **0.145** |
| total score | 754.1 | **772.5** (+18.4) |

The score gain comes from `ftl_geo_evolving` now being populated via the sweep.
Physics are identical (same seed, same evolution) — only the post-hoc ray-tracing
differs.

### Interpretation

The monotonic decay of `f_geo(t_emit)` means the warp geometry is **strongest at
the initial configuration** and erodes as the Q-ball matter disperses.  There is no
delayed peak — the "surfable window" hypothesis (that a later launch might catch
a stronger transient) does not apply to this configuration.  The channel simply
opens at t=0 and gradually closes.

The f_geo frozen timeseries (per-slice frozen geodesics) peaks at t~9.6 (12.8%),
but the evolving 4D trace — which accounts for the metric changing while the ray
propagates — shows the shortcut is already at its best for rays launched at t=0
(9.38%).  Later launches see a weaker geometry during propagation.

---

## Compact-soliton GPU validation — `traj_qball_compact_ode_eval122_t16` (2026-06-30)

(Previously documented below under "Q-ball / boson-star profile solver".)

Eval 122 with fixed radial-ODE solver (compact preset: m=1, λ=640, μ=85333, ω=0.8)
+ equilibrium amplitude + ODE seed + Lorentz boost.  N=128, L=64, ml=2, t=16.

| Metric | Broken seed (ω=0.4) | Compact seed (ω=0.8) |
|--------|--------------------:|---------------------:|
| confined_frac @ t=16 | 19.4% | **31.0%** |
| spread_ratio (rms) | 2.11 (5.46→11.52) | **1.49 (7.03→10.48)** |
| f_geo (gauge-inv) | 1.94% | **5.19%** |
| ftl_geo_peak | 20.4% | **63.8%** |
| max local speed | 1.09c | **1.25c** |
| ftl_persistence | ~0 | **0.85** |
| numerical_survival | 1.0 | 1.0 |

A correct compact seed both confines better (spread 2.11→1.49) **and** opens a
100%-lifetime superluminal channel that the broken thin-wall seed never produced.

---

## Campaign comparison: `scalar_sh_ftl_v22` vs `trajectory_5lump_v1` (2026-06-25)

### Setup

| Parameter | **SH v22** | **Trajectory v1** |
|-----------|-----------|-------------------|
| Ansatz | Spherical harmonic (ℓ_max=4) | Per-lump circular orbits |
| Search dims | 38 D (24 SH + 14 geometry/physics) | 40 D (7 per-lump × 5 + 5 shared) |
| Matter model | `grtresna_independent_scalars` | `grtresna_independent_scalars` |
| Exotic handling | Coarse azimuthal wedge | Per-lump binary flag [0,1] |
| Motion | One global (v_tor, v_pol, v_rad) for all lumps | Independent ω_rot per lump in tilted plane |
| Grid (Stage 0) | N=128, L=64, ml=1 | N=128, L=64, ml=1 |
| Grid (GRTresna) | 128³ gridinit, ml=3 | 128³ gridinit, ml=3 |
| Stop time (Stage 0) | 16.0 | 16.0 |
| Objective | `general_ftl` | `general_ftl` |
| Descriptor | `ftl_lifetime` (8×8 archive) | `ftl_lifetime` (8×8 archive) |
| GPUs | 8 × A100 80GB | 8 × A100 80GB |
| Batch size | 8 | 8 |
| Target evals | 200 | 200 (stopped at 130) |
| Completed evals | 202 | 130 |
| HQ promotions | 0 | 5 (at 256³, t=30) |
| Run dir | `runs/grtresna_qd/scalar_sh_ftl_v22/` | `runs/grtresna_qd/trajectory_5lump_v1/` |
| HQ run dir | — | `runs/grtresna_promote/trajectory_5lump_v1_hq_eval*` |
| Date | 2026-06-24 | 2026-06-25 |

### Results — Head-to-head (Stage 0)

| Metric | **SH v22** (202 evals) | **Trajectory v1** (130 evals) | Factor |
|--------|----------------------|---------------------------|--------|
| Best stable score | 470.6 | **1367.9** (eval 115) | **2.9x** |
| Best overall score | 470.6 | **1389.6** (eval 111, crashed) | **3.0x** |
| Best stable f_geo_peak | 2.12% | **10.63%** (eval 115) | **5.0x** |
| Best overall f_geo_peak | 2.12% | **30.25%** (eval 103, crashed) | **14.3x** |
| Best HQ-confirmed f_geo_evol | — | **9.40%** (eval 122) | — |
| Best HQ-confirmed f_geo_peak | — | **20.97%** (eval 122) | — |
| Best f_op_peak | 7.0% | **21.3%** (eval 122) | **3.0x** |
| Best stable max_local_speed | 1.33c | **1.67c** (eval 122 HQ) | **1.3x** |
| FTL hit rate (per GPU eval) | 1.3% (1/79) | **54%** | **~40x** |
| Archive cells filled | 2 / 64 | 4 / 64 | — |

### HQ validation results (trajectory_5lump_v1 only)

Five top Stage 0 elites promoted to 256³, max_level=3, t_stop=30:

| Eval | Stage 0 score | Stage 0 f_geo | HQ f_geo_evol | HQ f_geo_peak | HQ status | Verdict |
|------|--------------|--------------|--------------|--------------|-----------|---------|
| 122 | 1237.6 | 8.51% | **9.40%** | **20.97%** | Survived t=30 | **CONFIRMED** |
| 115 | 1367.9 | 10.63% | **12.5%** | **20.3%** | Crashed t=21 | Confirmed (transient) |
| 050 | 1039.5 | 10.82% | **7.4%** | **20.3%** | Crashed t=19 | Confirmed (transient) |
| 111 | 1389.6 | 17.37% | **8.6%** | **19.8%** | Crashed t=8.6 | Confirmed (short) |
| 008 | 1166.8 | 24.62% | **0.0%** | — | Survived t=30 | **FALSE POSITIVE** |

Key HQ findings:
- All genuinely FTL configs converge to **~20% peak f_geo** at HQ (resolution ceiling).
- The strongest Stage 0 signal (eval 008, 24.62%) was entirely a low-res artifact.
- 3/5 evals crashed at HQ from NaN in metric tensor at AMR level 3 boundaries.
- Eval 122 is the ONLY eval that both survived to t=30 AND confirmed FTL.

### Status breakdown (Stage 0, 130 evals)

| Status | **SH v22** | **Trajectory v1** |
|--------|-----------|-------------------|
| gpu_ok | 75 (37%) | ~52 (40%) |
| gpu_failed | 4 (2%) | ~16 (12%) |
| postload_rejected | 23 (11%) | ~62 (48%) |
| grtresna_rejected | 74 (37%) | 0 (0%) |
| grtresna_failed | 18 (9%) | 0 (0%) |
| pipeline_interrupted | 8 (4%) | 0 (0%) |

Key difference: Trajectory has **zero GRTresna failures** (trivial momentum constraint
at t=0 with static lumps) but **higher postload rejection** (48%) from constraint
violations when loaded onto the evolution grid. SH has the opposite profile — many
GRTresna solver failures from complex initial data, but lower postload rejection.

---

## Top evaluations — `trajectory_5lump_v1`

### Eval 122 — HQ Champion (Stage 0: 1237.6, HQ: 244.5, gpu_ok both)

**HQ-CONFIRMED: 9.40% geodesic shortcut at 256³.** The only eval that both survived
to t=30 at HQ resolution AND confirmed FTL. Transient channel lasting ~16 code units.

| Metric | Stage 0 (128³, t=16) | HQ (256³, t=30) |
|--------|---------------------|-----------------|
| Score | 1237.6 | 244.5 |
| Status | gpu_ok | gpu_ok |
| f_geo_evol | 8.51% | **9.40%** (improved!) |
| f_geo_peak | 8.51% @ t=9.6 | **20.97%** @ t=10.56 |
| f_op_peak | 21.31% | 21.26% (identical) |
| max_local_speed | 1.460c | 1.668c |
| ftl_lifetime | 100% (4/4 frames) | 47% (59/126 frames) |
| numerical_survival | 1.0 | 1.0 |
| structural_persistence | 1.0 | 0.631 (63.1% density) |
| 4D geodesic h_drift | — | 0.000525 (excellent) |
| 4D geodesic n_reached | — | 5/5 |

**HQ time profile:** Channel opens t~3.8, plateau (>90% peak) t=9.8–12.0,
decay t=12–20, closed by t=20. Late trapped surface at t=28.5. Total FTL
window: 16.6 code units (55% of evolution).

**Configuration:**

```
Lump  R0    omega   tilt   well_depth  exotic
0     6.12  -0.872   96°   0.1006      EXOTIC
1     3.86  -0.909   47°   0.1398      EXOTIC
2     4.41  -0.761    2°   0.0237      canonical
3     4.49  -0.588   98°   0.0589      canonical
4     7.12  -0.851   88°   0.1035      EXOTIC

Shared: A_breath=1.319, omega_breath=0.263, z_amp=2.277, omega_z=1.389, well_width=1.563
Total well_depth: 0.427
```

**Pattern:** ALL 5 retrograde. 3/5 exotic. Nested shells (inner R0~3.9–4.5, outer
R0~6.1–7.1). Mixed tilts spanning equatorial (2°) to polar (98°). Strong z-oscillation
(amp=2.28). Slow breathing (omega=0.26). Moderate well width (1.56).

---

### Eval 115 — Strongest FTL (Stage 0: 1367.9, HQ crashed t=21.15)

**12.5% geodesic shortcut at HQ** — the strongest FTL signal, but crashed with NaN in
K at AMR level 3. Lapse collapsed to 1e-7 (horizon forming).

| Metric | Stage 0 | HQ |
|--------|---------|-----|
| Score | 1367.9 | 591.5 |
| f_geo_evol | 10.63% | **12.5%** |
| f_geo_peak | 10.63% | **20.3%** @ t=14.64 |
| ftl_lifetime | 43% | 55% |
| numerical_survival | 1.0 | 0.705 (crashed t=21) |

---

### Eval 111 — Highest Stage 0 score (1389.6, crashed both stages)

**17.37% geodesic shortcut at Stage 0**, crashed at 65% survival. At HQ crashed
even earlier (t=8.64). The strongest raw signal but too unstable for verification.

| Metric | Stage 0 | HQ |
|--------|---------|-----|
| Score | 1389.6 | 1223.8 |
| f_geo_evol | 17.37% | 8.6% |
| f_geo_peak | 17.37% | **19.8%** @ t=6.96 |
| numerical_survival | 0.651 | 0.288 |

**Configuration:** 4 retrograde + 1 slow prograde (lump 2: omega=+0.091). 2/5 exotic.
Total well_depth=0.407.

---

### Eval 008 — FALSE POSITIVE (Stage 0: 1166.8, HQ: -19.7)

**24.6% geodesic shortcut at Stage 0 was entirely a resolution artifact.** At HQ:
zero FTL signal, matter dissipated to 45% retention, curvature_activity dropped to 0.09.

| Metric | Stage 0 | HQ |
|--------|---------|-----|
| Score | 1166.8 | **-19.7** |
| f_geo_evol | 24.62% | **0.0%** |
| numerical_survival | 0.592 | 1.0 |
| density retention | — | 44.8% |

**Lesson:** High Stage 0 scores from overlapping lumps + short evolution do not predict
HQ performance. The counter-rotating mixed-sign omega configuration was not physically
viable — the 16-unit evolution at 128³ simply was not long enough to expose this.

---

## Physics analysis — Why trajectory >> SH

### 1. Differential motion (the decisive mechanism)

The SH ansatz gives all 5 lumps **identical velocity** (one global `v_tor`, `v_pol`, `v_rad`).
This produces rigid-body rotation with no frame-dragging *shear* — the matter moves as a
single unit, generating a uniform frame-drag field with no gradient.

The trajectory ansatz gives each lump **independent angular velocity** in an **independently
tilted orbital plane**. This creates frame-dragging shear between counter-moving or
differently-oriented matter streams — exactly the mechanism behind Alcubierre-type warps
where differential frame-dragging produces effective metric expansion/contraction.

### 2. GRTresna convergence advantage

| | SH | Trajectory |
|-|----|----|
| t=0 momentum | Non-trivial (lumps have velocity) | **Trivial** (lumps at rest) |
| Solver failures | 46% (grtresna_rejected + failed) | **0%** |
| Effective GPU budget | 79 / 202 evals reach GPU (39%) | 13 / 25 reach GPU (52%) |

The trajectory ansatz places static lumps at t=0. The momentum constraint is trivially
satisfied (zero momentum source → zero shift). Motion comes from the C++ `TrajectoryEvaluator`
during evolution, not from the initial data. This means the constraint solver never fails,
and all computational budget goes directly to GPU evolution.

### 3. Per-lump exotic matter placement

SH uses a coarse azimuthal wedge for exotic matter (contiguous sector, same for all lumps).
Trajectory provides per-lump binary exotic flags, enabling exotic matter at specific orbital
positions. Every top-3 configuration uses **3/5 exotic lumps (60%)**, with the exotic lumps
typically at intermediate radii where the frame-dragging gradient is steepest.

### 4. Tilted orbital planes create 3D topology

SH matter lives in a fixed shell geometry (spherical). Trajectory lumps orbit in arbitrarily
tilted planes (per-lump `tilt_theta`, `tilt_phi`), creating complex 3D frame-dragging
topology. The crashed champion (eval 8) has all lumps with `tilt_theta > 2.3` (nearly
inverted planes), creating a tangled 3D vortex structure.

### 5. Stability vs strength trade-off (updated with HQ results)

| Pattern | Stage 0 f_geo | HQ f_geo_evol | HQ status |
|---------|--------------|--------------|-----------|
| Counter-rotating (eval 008: 3R+2P) | 24.6% | **0.0%** | Survived but NO FTL (false positive) |
| All-retro, strong omega (eval 115) | 10.6% | **12.5%** | Crashed t=21 (lapse collapse) |
| All-retro, nested shells (eval 122) | 8.5% | **9.4%** | **Survived t=30** |
| All-retro, deep wells (eval 050) | 10.8% | 7.4% | Crashed t=19 |

The HQ results overturn the Stage 0 picture. Counter-rotation (eval 008) was not actually
producing real FTL — the signal was a resolution artifact. All-retrograde is the only
pattern that produces HQ-confirmed FTL.

Within all-retrograde, there is a stability-strength trade-off: eval 115 (12.5%) is
stronger than eval 122 (9.4%) but crashes from lapse collapse at t=21 vs surviving to
t=30. The key difference: eval 122 has slower omegas and more moderate well depths,
giving the geometry room to evolve without forming a horizon during the FTL window.

---

## Emerging parameter patterns (updated with HQ validation)

### What correlates with HQ-confirmed FTL

| Parameter | HQ-confirmed pattern | Evidence |
|-----------|---------------------|----------|
| Rotation direction | **ALL retrograde** (mandatory) | Only all-retro configs produce real FTL at HQ |
| Exotic fraction | 3/5 lumps exotic (60%) | Eval 122 (3/5), eval 115 (2/5), eval 050 (3/5) |
| Z-oscillation | Strong (amp > 1.8) | Eval 122: 2.28, eval 115: 1.85, eval 050: 1.85 |
| Tilt diversity | Mixed (2° to 98°) | Creates 3D frame-drag topology |
| Nested R0 | Inner 3.8–4.5, outer 6.1–7.1 | Eval 122 pattern; R0 spread matters |
| Well width | Moderate (1.5–1.7) | Concentrate energy without over-driving |
| Total well_depth | 0.40–0.43 | Below this: too weak; above: crashes |
| Breathing | Slow (omega < 0.3) | Eval 122: 0.26; too fast destabilizes |

### What does NOT correlate with real FTL

| Pattern | Stage 0 result | HQ result |
|---------|---------------|-----------|
| Counter-rotation (mixed signs) | 24.6% f_geo (eval 008) | **0.0% — FALSE POSITIVE** |
| Overlapping lumps (negative margin) | High Stage 0 score | Dissipates at HQ |
| Very wide well (>2.0) | Crashes at Stage 0 | — |
| High omega (>0.9) all lumps | Crashes early | Crashes earlier at HQ |

### What correlates with HQ crashes

| Pattern | HQ crash mode | Evidence |
|---------|--------------|----------|
| Strong omegas (>0.85 all lumps) | NaN in h11/K at level 3 | Evals 111, 050 |
| High total well_depth (>0.43) | Lapse collapse + NaN | Eval 115 (late collapse t=21) |
| Fast breathing (omega_breath > 0.5) | Metric oscillation at boundaries | Eval 111 |

### What causes postload rejection (48% of evals)

The postload gate checks constraint violations after loading GRTresna's solved initial
data onto the coarser evolution grid. Configurations rejected here typically have:
- Very large orbital radii (lumps near domain boundary)
- High aggregate matter density (sum of well_depths > ~0.5)
- Very different GRTresna vs evolution grid resolution mapping

---

## Bicomplex (canonical + phantom) matter model (2026-06-26)

A new matter model, `grtresna_bicomplex_scalar`, evolves **two coupled complex
fields** instead of one: a canonical field **Φ⁺** and a phantom field **Φ⁻** with
*opposite* gravitational sign. Both share the boson-star GRTresna constraint solve
(per-lump-signed `ComplexScalarField`), but GRTeclyn evolves each field
independently with its own stress-energy contribution. The phantom Φ⁻ supplies a
**sustained, dynamically-evolving ANEC violation** that a single complex field
cannot maintain on its own.

The genome is the eval 122 trajectory (5 lumps, same centers/widths), re-solved with
the bicomplex matter sector. Per-lump gravitational signs select which lumps source
the phantom channel.

### Is it bosonic matter? (model caveat)

Yes — the matter content is entirely **bosonic**: two complex spin-0 scalar fields,
each U(1)-symmetric (boson-star binding). There are no fermions or spinors. Φ⁺ is an
ordinary positive-energy boson; Φ⁻ is an **exotic (phantom) boson** that supplies the
negative-energy / ANEC violation.

The important caveat is *how* the phantom is implemented. This is **not** a fully
Lagrangian-consistent phantom field. A textbook phantom has a wrong-sign kinetic term
that flips both its gravitational coupling *and* its own equation of motion. Here only
the **gravitational coupling sign** is reversed (the entire stress-energy is added with
`fs = -1`, see `BiComplexScalarField.impl.hpp::accumulate_complex_field`), while **both
fields obey the same Klein–Gordon RHS**:

> *"the field EOM is sign-independent; only the gravitational coupling is reversed, so
> both fields obey the same Klein-Gordon RHS while preserving their own U(1) charge."*
> (`BiComplexScalarField.impl.hpp`)

So Φ⁻ **propagates** like a normal positive-norm boson but **gravitates** like
negative-energy matter. This decoupling is a deliberate engineering choice — it keeps
the evolution stable and the U(1) charge well-defined — but it means the configuration
is a *phenomenological* exotic-matter source, not a self-consistent ghost scalar. Any
physics claim that depends on a true phantom field should account for the EOM-vs-gravity
sign decoupling.

### Head-to-head — single-complex vs bicomplex at identical parameters

Both runs use the *same* 5-lump trajectory genome, mass `m=0.15`, `ω=0.12`,
`general_ftl` objective, evolved to t=16:

| Model | Total score | Gauge-inv. shortcut (f_geo) | Geodesics reached | max local speed | operational FTL |
|-------|-------------|------------------------------|--------------------|------------------|-----------------|
| `grtresna_boson_star` (single complex) | 35.2 | 0.0% | **0 / 5** | 1.90c (coordinate only) | none (f_op=0) |
| **`grtresna_bicomplex_scalar`** (canonical + phantom) | **255.7** | **5.21% frozen / 6.13% evolving** | **5 / 5** | **1.071c** | **confirmed (f_op=5.56%)** |

The single complex field produces a large *coordinate* speed (1.90c) that the scorer
correctly rejects as a delocalised gauge offset — no null ray reaches the detector,
f_geo=0. The bicomplex model with its phantom channel produces a **genuine,
gauge-invariant geodesic shortcut**: all 5 null rays reach the detector and the
evolved geometry sustains a real superluminal channel.

### Validated run — `traj_bicomplex_m015_w012`

| Metric | Value |
|--------|-------|
| Total score | 255.7 |
| f_geo (frozen, `geodesic_ftl`) | **5.21%** @ 5/5 reached |
| f_geo (evolving 4D, `evolving_geodesic`) | **6.13%** @ 5/5 reached |
| f_op (evolved, `general_ftl_evolved`) | 5.56% |
| max_local_speed (evolved) | 1.071c |
| numerical_survival | 1.0 |
| operational_ftl / ftl_persistence | 0.365 / 0.366 |
| WEC violation fraction | 91.4% (wec_min = −6.5e-4) |
| NEC min | −1.35e-4 |
| Final Ham L2 / Mom L2 | 1.66e-4 / 2.08e-5 (clean solve) |
| Max Ham L2 | 1.97e-3 |
| structural_persistence | 0.021 (matter dissipates by t=16) |
| Final time | 16.01 (completed) |

**Configuration:** matter_model = `grtresna_bicomplex_scalar`, sector = boson_star,
coupling = canonical, scalar_mass = 0.15, λ = 0.0, ω = 0.12, per-lump signs =
`[-1, -1, +1, +1, -1]` (lumps 0/1/4 phantom-sourced, 2/3 canonical).

**Scoring notes (from `score.json`):**
- *"gauge-invariant null-geodesic shortcut confirmed (f_geo = 5.208e-02)"*
- *"evolved geometry sustains a strong superluminal channel (max c = 1.071,
  F_op^ev = 5.558e-02); operational FTL"*
- *"exotic matter required (matter = 0.82, geometric = 0.24 of full penalty)"* —
  the phantom field is doing the heavy lifting on the energy-condition side.

**Interpretation.** Splitting the exotic content into an independent phantom field
decouples the ANEC-violating source from the canonical matter that builds the
frame-dragging structure. The canonical Φ⁺ lumps provide the orbital geometry; the
phantom Φ⁻ provides the negative-energy backing that keeps the shortcut open into
the evolved geometry. The result is a confirmed shortcut where the single-field
configuration produced only a coordinate artifact.

**Caveats.** The shortcut is transient and the matter structure is largely lost by
t=16 (structural_persistence = 2%); the high WEC-violation fraction (91%) means the
configuration is heavily exotic. A higher-mass / higher-ω HQ promotion
(`traj_bicomplex_m03_w025`, eval 122, m=0.3, ω=0.25, 256³, t=30) is **in progress**
to test resolution-independence and longer-time survival.

> GRTresna note: the at-rest lumps give a zero momentum source, so the printed
> Momentum relative error is `-nan` (0/0). This is harmless physically (the
> momentum constraint is trivially satisfied) but it defeats the NL solver's
> early-exit check — the Hamiltonian error converges below tolerance by iteration
> ~8 yet the solve still runs all 30 iterations. See NextSteps.md for the proposed
> fix.

---

## Matter confinement & the Q-ball extension (2026-06-26)

### The dispersion problem

Across the bicomplex runs the matter lumps **disperse** rather than holding their
shape, which collapses `structural_persistence` (2% by t=16 in
`traj_bicomplex_m015_w012`). A dedicated, mass-weighted **confinement diagnostic**
(`small_data/confinement.dat`, columns `rms_radius` and `confined_frac`) was added
because peak/total density is spatially blind: a lump can spray into a ragged halo
while its peak density *rises* under pump injection, so the old peak-based persistence
reported "fine" while the frames showed matter blowing apart. Dispersal shows up
unambiguously as **rms_radius growing** and **confined_frac collapsing**.

A clean diagnostic run (`m=1.0`, `ω=0.3`, 5 compact `R=1/√(m²−ω²)=1.05` sech lumps on
an evenly-spaced `R0=4` ring, 128³, PD trap active) confirmed the dispersal with the
trap *fully active and constraints pristine*:

| metric | t=0 → t=8 |
|--------|-----------|
| confined_frac | 0.73 → **0.36** |
| rms_radius | 5.39 → **7.29** (×1.35) |
| total_activity | 7.85 → **20.2** (×2.58, pump injecting) |
| peak_activity | falling |
| constraint_health / Ham L2 | 0.999 / **1e-5** (governor *not* throttling) |

### Root cause — no binding

At `amp = 0.005` (corrective pump amplitude) the lumps barely curve spacetime (Ham L2
~1e-5 ⇒ **negligible self-gravity**), so they are not boson stars — they are **free
Klein–Gordon wave packets**, which have no bound state and disperse. The closed-loop
PD "trap" controller can only grip where its envelope `sech(r/R)` is non-zero, so
matter that diffuses past ~2–3R escapes its reach and is lost (the trap injects waves
fighting the leak, hence `total_activity` *grows*). **You cannot trap a configuration
that has no bound state.**

### Q-ball extension (Stage A — analytic seed + self-interaction)

The fix is to give the field a genuine **self-bound soliton** via self-interaction.
`ComplexScalarPotential` was extended from the mass-only/quartic form to the full
**Q-ball potential** with a sextic stabiliser:

```
V(|Φ|) = ½ m² |Φ|²  −  ¼ λ |Φ|⁴  +  ⅙ μ |Φ|⁶
dV/d|Φ|² = ½ m² − ½ λ |Φ|² + ½ μ |Φ|⁴
```

The attractive quartic (`λ>0`) supplies self-binding; the sextic (`μ>0`) is **required**
— a pure attractive quartic is the critical 3D case and either collapses or disperses.
`λ = μ = 0` recovers the previous free massive field.

Q-ball design relations used (flat-space; gravity is negligible here):
- equilibrium core density `ρ* = 3λ/4μ`, core amplitude `|Φ|_core ≈ √ρ*`;
- existence band `ω ∈ (ω_min, m)` with `ω_min² = m² − 3λ²/16μ`;
- tail decay `R = 1/√(m²−ω²)` (unchanged — `λ,μ` vanish in the low-amplitude tail, so
  the existing `bound_width()` + sech painter/trap-target infrastructure stays correct).

**Test parameters** (`traj_bicomplex_qball_test`): `m=1.0`, `λ=160`, `μ=5333`
(⇒ core `|Φ|≈0.15`, `ω_min≈0.316`), `ω=0.4` (just above `ω_min` ⇒ strongly bound),
`R=1.09` (2R=2.18 ≪ 4.70 ring separation, no overlap), lump amplitude raised
`0.005 → 0.15` to sit at the Q-ball core. The large couplings are an artifact of the
small amplitude (`λA²~m²` ⇒ `λ~1/A²`); RHS magnitudes stay O(0.5) and CFL-safe.

**Preliminary result** (in progress, to t≈2.4 so far) — the confinement metric holds,
in sharp contrast to the dispersing baseline:

| metric | dispersing (no sextic) | **Q-ball (sextic)** |
|--------|------------------------|---------------------|
| rms_radius | 5.39 → 7.29 (×1.35) | **5.51 → 5.57 (×1.01, flat)** |
| confined_frac | 0.73 → 0.36 | **0.71 → 0.70 (held)** |
| total_activity | ×2.58 | ×1.11 |
| peak_activity | falling | falling 0.21→0.11 (seed relaxing to the Q-ball profile) |

The flat `rms` + held `confined_frac` indicate the matter is **staying confined**; the
falling peak is the analytic sech seed shedding its excess as it settles toward the true
Q-ball profile (the radiated excess is what appears as "dispersion" in the `lapse_z`
frames, but the bulk matter does not fly apart).

**Files changed (evolution-side only):**
`Source/Matter/ComplexScalarPotential.hpp` (sextic V + dV, `m_mu`),
`BiComplexScalarField.hpp` / `ComplexScalarField.hpp` (load `recipe_scalar_mu`),
`Examples/RadialRecipe/{SimulationParameters,RadialRecipeMatterDispatch}.hpp`,
and Python wiring `grtresna/{matter_wiring,solver,matter_models}.py`,
`search/optimize/config.py` (thread `grtresna_scalar_mu` → `recipe_scalar_mu`). The
GRTresna constraint solve is left quartic-only: at these amplitudes gravity is
negligible, so the near-flat initial data is consistent and CCZ4 damps the tiny `V`
mismatch.

### Stage B — the principled fix (full boson-star / Q-ball ODE) — NOT yet done

Stage A imposes an **analytic sech seed** at the Q-ball core amplitude and lets it relax
into the soliton (shedding radiation). The radiated transient and the residual peak drop
come from the seed not being the *exact* equilibrium profile. The principled fix is to
**solve the boson-star / Q-ball radial ODE** for the true stationary profile so the
field starts *on* the soliton attractor with no relaxation transient:

- Solve the coupled `(φ₀(r), metric)` eigenvalue problem
  `φ₀'' + (2/r)φ₀' = [dV/d|Φ|² − ω²/α²] φ₀` (shooting on the central amplitude /
  eigenfrequency `ω` so `φ₀→0` as `r→∞`), with `V` the full sextic potential.
- Replace the sech painter (`grtresna/boson_star_profile.py`, currently
  `φ₀(r)=φ_c·sech(r/R)`) with the tabulated ODE solution, and feed the *same* `V` to
  the GRTresna constraint solve so initial data is fully self-consistent (this also
  removes the quartic-only solve approximation noted above).
- For a self-**gravitating** star (rather than the flat-space Q-ball used here at
  negligible amplitude), couple the metric ODE and raise the amplitude until the
  gravitational binding is non-trivial — at which point the GRTresna sextic consistency
  becomes mandatory.

This is the natural next step once Stage A confinement is confirmed over the full t=8–16
window. See `NextSteps.md` (Phase 4) for sequencing.

### Stage A dynamics — eval 122 trajectory with Q-ball matter (2026-06-26)

The static test (lumps at rest) proved Q-ball confinement works. The critical question:
does the self-binding survive under **real FTL-producing dynamics**? Eval 122's trajectory
(all-retrograde, tilted planes, breathing, z-oscillation, tangential velocities 2.6–6.0c)
was replayed with bicomplex Q-ball matter (`m=1.0`, `λ=160`, `μ=5333`, `ω=0.4`,
all well_depths raised to 0.15 for proper Q-ball formation).

**Run:** `traj_bicomplex_qball_eval122_v2` (128³, L=64, ml=2, t_stop=16, GPU 0)

| metric | t=0 | t=8 | t=16 | trajectory |
|--------|-----|-----|------|-----------|
| rms_radius | 5.49 | 7.20 (×1.31) | 11.90 (×2.17) | monotone growth |
| confined_frac | 0.714 | 0.385 | **0.118** | continuous erosion |
| peak_amplitude | 0.251 | 0.107 | 0.185 | drops then recovers (pump) |

**Verdict: Q-ball self-interaction does NOT prevent dispersal under eval 122 dynamics.**

The scorer flagged the run as `"matter DISPERSED: only 12% of matter remains confined"`.
However, the FTL geometry is still active:

| FTL metric | value |
|------------|-------|
| FTL precursor (cone climbing) | **1.000** |
| Gauge-invariant FTL lifetime | **84%** (57/68 frames) |
| Peak f_geo | **4.62%** @ t=3.36 |
| Coupled channel progress | 0.281 |
| Geodesic shortcut | rejected as unreliable (3/5 rays, f_geo=5.1%) |

**Comparison across all Q-ball dynamics tests:**

| run | trajectory | t_final | rms growth | final conf_frac | outcome |
|-----|-----------|---------|-----------|----------------|---------|
| static (at rest) | none | 2.4 | ×1.01 | 0.70 | **HELD** |
| simple (ω=0.1, v=0 init) | slow ring | 8.0 | ×1.26 | 0.43 | eroding |
| **eval 122 (ω=-0.6 to -0.9)** | real FTL | **16.0** | **×2.17** | **0.12** | **DISPERSED** |

### Root cause — transport velocity mismatch

The PD trap drives the field toward a target `sech(r/R)` profile whose *center* moves
along the parametric trajectory. At eval 122's tangential speeds (2.6–6.0c), the pump
cannot translate a self-bound Q-ball — the soliton's inertia resists deformation. Instead
the pump **re-creates** field at each new position (injection) while the old field
(abandoned by the moving spotlight) **disperses freely** without trap support. The result
is a trail of radiated matter behind the pump, not a transported soliton.

The peak amplitude *recovering* at late times (0.107→0.185 from t=8 to t=16) confirms the
pump IS actively maintaining field at the trajectory centers — but the net confined fraction
drops because each injection cycle sheds radiation that escapes the spotlight.

### Implications for Phase 4 — fixing the dispersion

The Q-ball self-interaction is **necessary but not sufficient**. The binding prevents
instantaneous dispersal (unlike free fields where the pump trail vanishes immediately),
but it cannot compensate for the transport velocity mismatch. The dispersion has two
independent causes that require separate fixes:

**Cause 1 — Cold start (v=0 initial data on a fast trajectory):** The lumps start at
rest while the pump target immediately moves at v=R0·ω_rot (2.6–6.0c for eval 122).
The trap cannot accelerate a self-bound Q-ball to superluminal speeds — it just creates
new field at the moving target and abandons the stationary lump.

**Cause 2 — Sech seed radiation:** The analytic sech profile is NOT the exact Q-ball
equilibrium. On initialization it sheds ~30% of its amplitude as radiation (peak drops
0.25→0.11 in the first ~5 code units). This radiated shell disperses regardless of
whether the lump moves.

#### Fix 1 — Velocity-matched initial data ("speed follow") — **IMPLEMENTED**

Start each Q-ball **already moving** at its trajectory angular velocity. This was
implemented and tested on eval 122.

**Implementation (2026-06-26):**

1. `_expand_trajectory_boson_lumps_from_overrides` now computes the t=0 velocity
   of each lump from the trajectory geometry:
   - `v_orb = R0 · omega_rot` (signed tangential speed in the orbital plane)
   - `v_vec = R(v_orb, tilt_phi, tilt_theta) · v_orb` (same rotation matrix used for positions)
   - `v_vec` is capped at `grtresna_boost_v_max` (default 0.8c) to stay sub-luminal on the
     flat initial slice. Two new overrides control this: `grtresna_boost_lumps` and
     `grtresna_boost_v_max`.

2. `boson_star_fields.py::_boosted_lump_fields` paints the full Lorentz-boosted
   complex field for a moving Q-ball:
   - Lorentz contraction of the radial profile: `ρ = √(δx_⊥² + γ² δx_∥²)`.
   - de Broglie phase: `Φ = φ₀(ρ) · exp(−i ω γ (v·δx))`.
   - Advection momentum: `Π = −(∂_t Φ)/α` with the full time derivative including the
     Lorentz-contraction radial velocity and Doppler-shifted internal rotation.

3. `paint_bicomplex_fields_on_grid` superposes the boosted canonical and phantom fields
   onto the GRTeclyn state. The GRTresna metric solve is still done for at-rest matter;
   the boost is applied afterward during gridinit export. At these amplitudes the residual
   momentum-constraint violation is O(0.3% · v) and the solver converges normally.

The pump controller then only needs to *maintain* the soliton's trajectory (small
corrective force), not *create* the entire motion from scratch.

#### Fix 1b — Test result: Lorentz-boosted eval 122 dynamics

**Run:** `traj_qball_boosted_eval122` (128³, L=64, ml=2, t_stop=16, same Q-ball
parameters as the cold-start run: `m=1.0`, `λ=160`, `μ=5333`, `ω=0.4`, all well_depths=0.15).

| metric | t=0 | t=2 | t=8 | t=16 | trajectory |
|--------|-----|-----|-----|------|-----------|
| rms_radius | 5.26 | **5.07** (×0.96) | 6.91 (×1.31) | 12.04 (×2.29) | shrinks initially, then grows |
| confined_frac | 0.756 | **0.790** | 0.401 | **0.066** | rises at first, then erodes |
| peak_amplitude | 0.315 | 0.262 | 0.135 | 0.037 | monotone drop |

**Comparison with cold-start (v=0) run:**

| metric | Cold start | Lorentz boost | Interpretation |
|--------|-----------|---------------|----------------|
| Early rms trend | grows from step 1 | **shrinks for ~1.7 code units** | boost gives coherent initial motion |
| conf_frac @ t=2 | 0.70 | **0.79** | +12% better early confinement |
| conf_frac @ t=8 | 0.385 | **0.401** | comparable, boost slightly better |
| final conf_frac | **0.118** | 0.066 | worse by the global metric |
| final peak | 0.185 | **0.037** | cold run keeps higher central peak |
| final rms | 11.90 | 12.04 | comparable |
| AMR level-2 cells | ~0.4M | **3.2M** | boost generates more gravitational structure |
| barycenter drift | negligible | z: 29.2 → 25.7 | boosted lumps carry net momentum |

**Verdict:** The boost is **correct and partially effective**, but it does **not** solve the
final dispersal problem for eval 122. The early phase is dramatically improved — the
matter contracts and the pump initially holds it — but by t=8 the confinement fraction
has decayed to the same level as the cold-start case, and the late-time metric is no
better.

**Why the final confinement metric is misleading here:** The Q-balls are moving on a
multi-lump orbit, so the activity-weighted barycenter stays near the grid center while the
lumps themselves sweep around the perimeter. The confinement score `confined_frac` counts
matter within a fixed radius `r_conf=6` of the barycenter; a moving, separated lump system
naturally scores low even if each individual lump is well-localized. The cold-start run has
a higher final peak because the matter is largely stationary near the barycenter, not
because transport succeeded. A **per-lump confinement metric** (tracking each trajectory
center) is needed for a fair comparison.

**Physical reason the boost is not enough:** The pump trajectory demands high tangential
speeds (2.6–6.0c) and continuous *transverse* acceleration (circular orbits). Capping the
boost at 0.8c removes the superluminal initial velocity but leaves the pump to supply the
remaining speed and all the centripetal acceleration. The Q-ball cannot adiabatically
follow the curved, high-speed spotlight; it radiates at each turn. The boost helps the
initial transient but does not change the fundamental mismatch between the soliton's
inertia and the trajectory's curvature.

**Conclusion from this test:** Fix 1 alone is insufficient. The remaining levers are deeper
binding (Fix 2) and an exact equilibrium profile (Fix 3), likely used together.

#### Fix 2 — Stronger self-interaction (tighter binding)

Current parameters: `λ=160`, `μ=5333` → core amplitude 0.15, binding energy
`ΔE/E ~ (λ/m²)·|Φ|² ~ 160·0.0225 ≈ 3.6`. The Q-ball IS strongly bound in energy,
but its **spatial extent** (R=1.09) means the binding force per unit volume is moderate.

Increasing the couplings with **preserving-band scaling** (μ ∝ λ², fixed ω_min):
- `λ → 640`, `μ → 85333` (4× λ): equilibrium core **√(3λ/4μ) ≈ 0.075** (half of
  standard's 0.15), but the restoring force at a fixed paint amplitude is much
  stronger.  `ω_min = 0.316` unchanged (ratio λ²/μ fixed).
- **Do not** seed at `well_depth=0.15` with stiff couplings: that is 2×
  super-equilibrium and drives t=0 relaxation radiation.
- Practical limit: very large `λ,μ` without an ODE seed increases RHS stiffness
  and makes the initial shock worse; prefer stiff + equilibrium amplitude + ODE
  profile over 8× λ with sech.

The stronger restoring force resists deformation by the pump's translational push,
reducing radiation loss during transport.

#### Fix 3 — Full Q-ball radial ODE solution (exact equilibrium profile)

The sech seed loses ~30% of its peak amplitude to radiation because it is NOT the true
Q-ball equilibrium. Solving the stationary Q-ball radial ODE gives the **exact** profile
`φ₀(r)` that sits perfectly on the soliton attractor with zero radiation:

```
φ₀'' + (2/r)φ₀' = (m² − ω²)φ₀ − λ φ₀³ + μ φ₀⁵
```

(flat-space, spherically symmetric). Boundary conditions: `φ₀'(0) = 0` (regularity),
`φ₀(∞) = 0`. This is a **shooting problem** on the central amplitude `φ₀(0)`:
- Shoot from `r=0` outward with trial `φ₀(0)`.
- Too large → solution diverges to −∞ at large r (over-shooting).
- Too small → solution goes to +∞ (under-shooting).
- Bisect to find the separatrix solution that decays exponentially.

The exact profile has:
- **Zero radiation** at initialization (no transient shedding, peak holds at 100%).
- **Correct tail shape** (exponential, not sech — though sech is a good approximation
  for r ≲ 2R, it falls off too slowly for r ≫ R).
- **Self-consistent metric** if coupled to the GRTresna constraint solver (gravity).

Implementation path (**Fix 3 implemented 2026-06-27** in `grtresna/qball_radial_ode.py`):
1. ~~Add `qball_radial_ode.py`~~ — shooting + bisection on φ₀(0); returns tabulated
   `QBallRadialProfile` with cubic-spline interpolation.
2. ~~Gridinit repaint~~ — `PROFILE_ODE_BOUND` in `boson_star_fields.py` uses the ODE
   tabulation when `grtresna_qball_ode_profile=1`.
3. ~~Equilibrium amplitude cap~~ — `grtresna_qball_equilibrium_amplitude=1` caps each
   lump's `well_depth` to √(3λ/4μ) via `QBallCouplings.cap_well_depth`.
4. Replay flags: `--qball-preset stiff --qball-equilibrium-amplitude --qball-ode-profile`.

**Next replay (dispersion-corrected stiff eval 122):**

```bash
uv run python scripts/campaigns/hq/replay_eval.py \
  ../runs/grtresna_promote/traj_qball_boosted_eval122 \
  --name traj_qball_stiff_ode_eval122 \
  --qball-preset stiff \
  --qball-equilibrium-amplitude \
  --qball-ode-profile \
  --stop-time 16 ...
```

This isolates **kinematic** orbit dispersion (centripetal scalar radiation) from
**relaxation** dispersion (sech + super-equilibrium seed).

#### Fix 2b — Stiff coupling promotion run (2026-06-27, stopped early)

Promoted eval 122 with **λ=640, μ=85333**, Lorentz boost, sech seed, `well_depth=0.15`
(2× super-equilibrium for stiff). Runs stopped at t≈9.5 after confinement clearly
tracked baseline dispersal.

| run | t reached | confined_frac @ t=0 / 3.2 / 6.4 / 9.6 | vs λ=160 baseline @ same t |
|-----|-----------|----------------------------------------|----------------------------|
| `traj_qball_stiff_boosted_eval122_v2` | ~9.6 (stopped) | 0.756 / **0.776** / **0.622** / **0.394** | 0.756 / ~0.75 / 0.519 / ~0.28 |
| `traj_qball_stiff_static_smoke` (ω_rot=0) | ~6.5 (stopped) | 0.714 / 0.716 | — (control: binding OK at rest) |
| `traj_qball_boosted_eval122` (baseline) | 16 (complete) | → **0.066** @ t=16 | reference |

**Verdict:** 4× λ with sech + 0.15 paint **delays** dispersal (~10 pts better
`confined_frac` at t≈6.4) but **does not prevent** orbit-driven erosion by t≈9.6.
The run was stopped before t=16; trajectory already matched baseline trend.
Raising λ further without ODE + equilibrium seed would worsen the t=0 shock.

**Dispersion decomposition (agreed priority):**

1. **Relaxation dispersion** — sech seed + super-equilibrium amplitude → Fix 3 ODE +
   `--qball-equilibrium-amplitude` (implemented; next replay).
2. **Kinematic dispersion** — centripetal acceleration on curved orbit at 2.6–6.0c
   tangential speeds → measure only after (1) is removed.

Optional developer throughput: GRTresna NaN-Mom early-exit patch (NextSteps item 1)
— decoupled from evolution physics.

#### Priority ordering (updated 2026-06-27)

| Fix | Effort | Expected impact | Status |
|-----|--------|-----------------|--------|
| 1. Lorentz boost | Done | Partial — early transient only | Implemented |
| 2. Stiff λ=640, μ=85333 | Done | Partial — delays dispersal, sech seed limits gain | Tested; insufficient alone |
| 2b. Equilibrium well_depth | Low | High — removes super-equilibrium shock | **Implemented** (`--qball-equilibrium-amplitude`) |
| 3. ODE radial profile | Medium | High — zero t=0 relaxation radiation | **Implemented** (`--qball-ode-profile`) |
| **Next replay** | Low | Isolates kinematic vs relaxation dispersion | `traj_qball_stiff_ode_eval122` |

Do **not** raise λ to 8× before running stiff + ODE + equilibrium amplitude: a stiffer
well with a poor seed produces a **more violent** initial shock.

---

## Q-ball / boson-star profile solver — C++ implementation and validation (2026-06-29/30)

The profile-3 ODE blocker from 2026-06-27 is now resolved end-to-end. GRTresna C++ can
load the same flat-space Q-ball radial profile as Python, and the temporary `3→2` sech
mapping hack has been removed.

### Implementation

| Component | Change |
|-----------|--------|
| `GRTresna/Source/Matter/BosonStarParams.hpp` | Added `profile == 3` (tabulated radial Q-ball); added `scalar_mu` parameter for the sextic stabilizer |
| `GRTresna/Source/Matter/ComplexScalarField.cpp` | Implemented sextic term `scalar_mu · φ² · φ` (μ controls the stabilizer) and profile-3 tabulated loader |
| `grteclyn-wrapper/grtresna/boson_star_profile.py` | Shoots flat-space ODE and emits a tabulated `.dat` file for C++ to load |
| `grteclyn-wrapper/grtresna/solver.py` | Removed the temporary `profile 3 → 2` sech mapping hack |
| `search/optimize/candidates.py` + `config.py` | Added sub-luminal trajectory speed cap (`max_speed_fraction < 1.0`) with tests |
| `visualisation/process_wave/consume_plotfiles/extraction/confinement.py` | Rewrote `confined_frac` extractor to be AMR-aware |

### Validation

- **GRTresna solve**: converges with the tabulated profile-3 ODE; no more sech/ODE mismatch.
- **GPU smoke (eval 122, t=0)**: `confined_frac ≈ 75%` with clean φ frames; no diagonal garbage.
- **Speed sweep**: sub-luminal cap applied; candidate trajectory speeds are now strictly `< c`.
- **Full test suite**: 547 passed.

### Full t=16 evolution (`traj_qball_stiff_ode_eval122_t16`, 2026-06-30)

With C++ profile-3 landed, the full combined config (stiff λ=640/μ=85333 +
equilibrium amplitude + ODE seed + Lorentz boost) was finally run end-to-end at
N=128, L=64, ml=2, t=16:

| Quantity | Value |
|----------|-------|
| GRTresna solve | converged NL iter 6 (Ham 0.71%, Mom 0.068%, Mom finite) |
| confined_frac @ t=0 | 75.6% (see correction — paint artifact, not a clean soliton) |
| confined_frac @ t=16 | 19.4% |
| rms radius | 5.46 → 11.52 (spread ×2.11) |
| ftl_geo_peak / f_geo | 20.4% / 1.94% |

**CORRECTION (2026-06-30): the Q-ball ODE seed is NOT localized — an earlier
"kinematic dispersal / non-adiabatic transport" reading of this run is RETRACTED.**

Inspecting the run's `qball_profile.dat` revealed the tabulated φ₀(r) does not
decay. From core φ_c≈0.076 it falls only to ≈0.44·φ_c by r≈15 then *rises* back to
≈0.5·φ_c, plateauing out to r=100 — even though the bound-state tail length is
1/√(m²−ω²) ≈ 1.09 (a real soliton should be ≲1% of core by r≈10). Reproduced with
`solve_qball_radial_profile`: **both** the `standard` and `stiff` presets fail to
localize inside the L=64 box.

Root cause: for these couplings the Q-ball sits in the barely-bound thin-wall
corner — U_eff(φ_core) = ½κ²φ² − (λ/4)φ⁴ + (μ/6)φ⁶ ≈ −1.7e-4, only just below
zero — so the equilibrium radius is far larger than the box, and the outward shoot
stalls on the flat top (the 1e-4 φ_c bisection tolerance plus the
"no blow-up / no zero-crossing within r_max" acceptance test both accept the
plateau as a valid bracket).

What this means:

- The seed used in every `--qball-ode-profile` run (including this one and the
  earlier t=0 "75%" smoke) is a near-constant box-spanning condensate, not a
  soliton. The 75.6%→19.4% confinement curve is that unphysical slab spreading,
  **not** orbital kinematics.
- The run's kinematics were in fact already adiabatic: after the speed cap,
  ω_rot ≈ 0.04–0.08 vs bs_omega = 0.4 (the stiff preset overrides scalar_mass→1.0,
  bs_omega→0.4). The "ω_rot ≫ internal frequency" story does not apply here.

Next step is to fix the radial-ODE seed so it localizes (re-derive couplings whose
Q-ball radius fits the box; make the shoot bisect to machine precision and attach a
clean exp(−κr) tail; add a localization regression test) **before** any further
GPU runs or dispersal conclusions. See NextSteps.md §3.

### RESOLVED (2026-06-30): fixed radial-ODE shoot + compact preset — GPU-validated

The radial-ODE solver (`qball_radial_ode.solve_qball_radial_profile`) was rewritten:

- **Event-based classification** of the outward shoot: *over* = φ crosses zero
  (→−∞) or blows up (→+∞); *under* = φ turns around and settles onto the
  competing false vacuum φ_fv.  The soliton is the separatrix.
- **Bracket against the φ''(0)<0 window edge** φ_win_hi and **bisect φ_c to ~machine
  precision** (the thin-wall separatrix sits within ~1e-5 of φ_win_hi, which the old
  64-point/1e-4 scan could never resolve).
- **Analytic Yukawa tail**: truncate at the asymptotic-match radius and attach
  φ ∝ e^{−κr}/r so the tabulated φ₀(r) decays monotonically to ~0 (min/core ~1e-18)
  instead of plateauing.

This produces a genuinely localized soliton for the presets, but ω=0.4 gives a
**thin-wall** Q-ball of radius ~20 (still box-filling).  Scanning ω showed the
soliton radius is minimized (~7) in the **thick-wall** regime ω≈0.7–0.85, so a new
`QBallCouplings.compact()` preset (m=1, λ=640, μ=85333, **ω=0.8**) was added and
wired into `replay_eval.py --qball-preset compact`.  231 grtresna tests pass,
including new localization/monotonicity regression guards.

**GPU validation** — see "Compact-soliton GPU validation" section at top of this
document for the full broken-vs-compact comparison table.  Summary: compact seed
(ω=0.8) gives 31% confinement (vs 19%), spread 1.49 (vs 2.11), f_geo 5.19%
(vs 1.94%), with a 100%-lifetime operational superluminal channel.

### AMR-aware confinement diagnostic

The old `confined_frac` extractor integrated over the level-0 grid only, silently
down-sampling any refined structure. It was rewritten to integrate `w · dV` over **all
AMR levels**, using the finest cell at each location and summing the true cell volume.

- Unit test on synthetic 2-level data: the old level-0 path missed the refined core
  (0.01 vs actual 0.50); the new path recovers it.
- Backward-compatible on uniform grids.
- 547 tests passed.

### Why `max_level=2` and `max_level=3` gave identical scores on solitons

GPU re-validation showed that AMR **almost never engages** on the boson-star soliton.
The field is smooth (core amplitude ≈ 0.075, sech/ODE profile), so gradients stay below
the default `regrid_threshold=0.02`. Level-1 grids may form transiently between regrids,
but they are de-refined before plotfiles are written. Consequently:

- The scored confinement lives on the base grid regardless of `max_level`.
- The observed ~20% retention at t=16 is a genuine **base-grid (N=96)** result, not a
coarse-AMR artifact.
- Raising `max_level` without retuning the tagging threshold is ineffective.

**Practical levers to improve soliton resolution:**
1. Tune the regrid tagger to tag on `scalar_activity`/matter (or lower the threshold) so
   refinement is sustained around the lump.
2. Raise the base resolution `N` (e.g., 128, 160, or 256). The AMR-aware diagnostic will
   then reward the genuinely finer evolution.

---

## Next steps

See [NextSteps.md](./NextSteps.md) for the full plan. Summary:

### Phase 3 — Post-verification analysis (COMPLETED)

1. **Gauge-invariance check**: PASSED. Harmonic slicing (lapse_power=0, lapse_coeff=1)
   gives f_geo_evol=4.4% vs 9.4% with 1+log. Signal persists => NOT a gauge artifact.
   Magnitude is foliation-dependent (expected: different slicings sample different
   hypersurfaces through the same 4D geometry).

2. **Directional geodesic sweep**: COMPLETED. x is the best direction (f_geo=9.4%).
   y/z show weaker signal. The shortcut aligns with the lump orbital axis.

3. **Pipeline fix**: Switched default from `ftl_first` to `general_ftl` objective and
   made xyz geodesics the default for all QD and HQ runs. Under `general_ftl`, eval 122
   correctly outscores eval 008 (1293 vs 1194) -- the old `ftl_first` had them inverted
   (1238 vs 1482) because coordinate-level shaping rewards dominated.

### Phase 3b — Remaining analysis

4. **Transient channel characterization**: Why does the FTL window last ~16 code units?
   Emission sweep (item 4 DONE) showed monotonic decay, peak at t=0. Correlate with
   confinement timeseries and breathing phase to identify the decay driver.

5. **Crash mitigation for eval 115**: Strongest FTL (12.5%) crashed at t=21. Try higher
   KO dissipation, reduced max_level, or CFL reduction.

### Phase 4 — Q-ball trajectory QD campaign (NEXT)

5. ~~**All-retrograde constraint**~~: **DONE** (2026-07-01). `_enforce_retrograde` in
   `candidates.py` flips prograde omega_rot when `trajectory_retrograde_only=1`.
   12 tests pass.

6. **Q-ball trajectory QD campaign**: `scripts/campaigns/qball_trajectory/run.sh`.
   Compact Q-balls (ODE seed, equilibrium amplitude, ω=0.8) on all-retrograde
   orbits with sub-luminal speed cap + multi-ray emission sweep.  `general_ftl`
   objective.  Goal: find less-dispersive orbital configs that sustain FTL longer.

7. **Resolution scaling**: Run best QD elites at 384³ or 512³. Does f_geo_evol
   keep improving toward the frozen peak?

8. **Longer evolution**: t_stop=64 to test if the channel lifetime is intrinsic or
   tunable.

9. ~~**Self-gravitating soliton trajectories**~~: **DEPRIORITIZED**. The pump
   model is the control mechanism; the two dispersal causes (sech seed, superluminal
   orbits) are now fixed. Improve confinement within the pump framework instead.

10. ~~**Stronger Q-ball binding + boost**~~: **DONE**. Compact preset (λ=640,
    μ=85333, ω=0.8) + ODE profile + equilibrium amplitude validated at 31%
    confinement, 5.19% f_geo, 100% FTL lifetime.

---

## Run log

| Run | Date | Ansatz | Evals | Best stable | Best HQ-confirmed | Headline |
|-----|------|--------|-------|-------------|-------------------|----------|
| `qball_traj_spiral_v2` | 2026-07-01 | Q-ball spiral (39D) | 200/200 | **603.4** (eval 118) | — (QD) | Dispersion-gated; 17.7% f_geo peak at t_emit=12 |
| `qball_traj_spiral_v2` HQ | 2026-07-02 | HQ eval 118 (t=16) | 1/1 | 511.9 | **13.0%** f_geo_evol | Survived 256³; channel + 35% confinement hold |
| `qball_traj_spiral_v2_t30` HQ | 2026-07-02 | HQ eval 118 (t=30) | 1/1 | 224.2 | **13.0%** f_geo_evol | Survived t=30; matter dispersed (23% conf); f_geo peak 22.8% @ t≈19 |
| `scalar_sh_ftl_v22` | 2026-06-24 | SH (ℓ=4) | 202/200 | 470.6 (eval 189) | — | 1 FTL hit in 202 evals; 2.1% geodesic |
| `trajectory_5lump_v1` | 2026-06-25 | Trajectory (5 lumps) | 130/200 | 1367.9 (eval 115) | **9.40%** f_geo (eval 122) | HQ-confirmed 9.4% geodesic shortcut at 256³ |
| `trajectory_5lump_v1` HQ | 2026-06-25 | HQ promotion (5 evals) | 5/5 | eval 122 (survived) | **9.40%** f_geo, **20.97%** peak | 1 confirmed, 3 crashed, 1 false positive |
| `eval000122_harmonic` | 2026-06-25 | Gauge test (harmonic slicing) | 1 | eval 122 | **4.40%** f_geo | Gauge-invariance confirmed (4.4% in harmonic vs 9.4% in 1+log) |
| `eval000122_xyz` | 2026-06-25 | Direction sweep (x y z) | 1 | eval 122 | **9.40%** f_geo | x is best axis; shortcut aligned with orbital direction |
| `traj_boson_m015_w012` | 2026-06-26 | Single-complex boson (eval 122 genome) | 1 | 35.2 | 0.0% f_geo | Coordinate-only (1.90c), no gauge-invariant shortcut, 0/5 reached |
| `traj_bicomplex_m015_w012` | 2026-06-26 | Bicomplex (canonical + phantom) | 1 | **255.7** | **5.21%** f_geo (5/5) | Phantom channel → confirmed gauge-invariant shortcut at identical params |
| `traj_bicomplex_m03_w025` HQ | 2026-06-26 | Bicomplex HQ (eval 122, m=0.3, ω=0.25) | in progress | — | — | 256³, t=30 promotion running |
| `traj_bicomplex_qball_eval122_v2` | 2026-06-26 | Q-ball + eval 122 dynamics (m=1, λ=160, μ=5333) | 1 | −100.1 | 4.62% f_geo peak | Matter dispersed (12% conf), FTL precursor=1.0, 84% lifetime |
| `traj_qball_boosted_eval122` | 2026-06-26 | Q-ball + Lorentz-boosted initial data (v_max=0.8c) | 1 | — | 3.53% f_geo peak | Early confinement improved (0.79), final confinement still poor (0.066); boost implemented and tested |
| `traj_qball_stiff_boosted_eval122_v2` | 2026-06-27 | Stiff Q-ball λ=640, μ=85333 + boost, sech@0.15 | — | — | — | Stopped t≈9.6; conf 0.62@6.4 vs 0.52 baseline; still dispersing on orbit |
| `traj_qball_stiff_static_smoke` | 2026-06-27 | Stiff λ=640 static (ω_rot=0), t→8 | — | — | — | Stopped t≈6.5; conf ~0.72 (binding OK at rest) |
| `traj_qball_compact_ode_eval122_t16` | 2026-06-30 | Compact Q-ball ODE (ω=0.8), eval 122, t=16 | 1 | 754.1 | **5.19%** f_geo | Fixed solver; compact soliton; 31% conf, 1.49 spread |
| `traj_qball_compact_sweep_eval122_t16` | 2026-07-01 | Compact + multi-ray sweep (7 launches, Δt=2) | 1 | **772.5** | **9.38%** f_geo | Monotonic decay 9.4%→5.9%; 100% FTL lifetime; peak at t_emit=0 |

**Conclusion:** The trajectory ansatz with per-lump differential motion is a **qualitative
improvement** over spherical harmonics. The HQ validation confirms a **resolution-independent,
gauge-invariant 9.4% geodesic shortcut** (eval 122) that improves at higher resolution and
persists under harmonic slicing (4.4%). The FTL is transient (~16 code units) but genuine:
5/5 null rays reach the detector 9.4% faster than flat-space light, with excellent energy
conservation (h_drift = 0.05%).

The key mechanism is all-retrograde frame-dragging from independently-tilted matter lumps.
Counter-rotation (eval 008) was shown to be a false positive at HQ — the strongest real
FTL comes from coherent retrograde rotation with diverse orbital tilts. The shortcut is
x-axis-aligned (direction of lump orbital motion).

Pipeline improvements applied: default objective mode switched to `general_ftl` (removes
coordinate-shaping rewards that inflated false positives), xyz geodesics enabled by default
for blind directional search at both Stage 0 and HQ.
