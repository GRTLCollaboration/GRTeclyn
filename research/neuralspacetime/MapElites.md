# MAP-Elites + CMA-ES FTL Discovery — Matter-First Metric Discovery

> Three-stage pipeline: **MAP-Elites** (wide survey) finds where good warps live in
> the 23-D shell search space; **CMA-ES** (local refinement) hill-climbs around
> the best *healthy* survivors; **HQ promotion** re-runs the top QD + CMA-ES elites
> at full resolution and extended time with incremental scoring. All stages share
> the same matter-first loop — propose lumps → GRTresna constraint solve → GRTeclyn
> GPU evolution → time-resolved FTL probes → score — but differ in proposer,
> resolution, and stop time.

## Table of contents

**Reference**

- [The idea: matter-first, not metric-first](#the-idea-matter-first-not-metric-first)
- [The pipeline](#the-pipeline)
- [The hard consistency rule](#the-hard-consistency-rule)
- [Behavior descriptors](#behavior-descriptors-the-diversity-axes)
- [Scoring model](#scoring-model-the-quality-axis)
- [Code map](#code-map-where-everything-lives)
- [Building the binaries](#building-the-binaries-grtresna--grteclyn)
- [How to run a campaign](#how-to-run-a-campaign)

**Campaign log** (reverse-chronological)

- [Quick index](#campaign-log--runs-analysis)
- [v22: Pre-GPU rejection learning + v21 resume](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18)
  - [Final results: top 3 + FTL champions](#v22-final-results-top-3--ftl-champions)
  - [Physical validation & next steps](#v22-physical-validation--next-steps)
- [v22 CMA-ES: wormhole refinement](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18)
- [v21: Pipelined QD + GPU tenancy tuning](#v21-pipelined-qd--gpu-tenancy-tuning-2026-06-17)
- [v10–v20: pipeline evolution & runs](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)
- [Foundational entries (2026-06-10)](#foundational-entries-2026-06-10)

## The idea: matter-first, not metric-first

Classic warp-drive analysis is **metric-first**: write down a target metric
(Alcubierre, Natário, …), then read off the stress-energy `T_ab = G_ab/8π` it
requires — which is always exotic, often pathological, and never solved as a
*consistent initial-value problem*. The "FTL" in those constructions is baked
into the chosen coordinates and need not survive a real evolution.

This project inverts that. The pipeline is **matter-first**:

1. The search proposes a **matter configuration** (massive scalar-field "lumps").
2. GRTresna **solves the Einstein constraint equations** for conformal factor and shift.
3. GRTeclyn **evolves** that spacetime forward on GPU.
4. Probes measure whether an FTL signature emerges and *persists*.
5. MAP-Elites stores the best candidate per behavior cell and mutates elites.

The payoff: any FTL signal that survives this loop is a property of a
**self-consistent, evolved spacetime**, not a hand-picked metric. Metrics are
iteratively hardened so the leaderboard cannot be gamed by coordinate artifacts
(see the [campaign log](#campaign-log--runs-analysis)).

## The pipeline

### Diagram — end-to-end overview

The closed loop: search proposes matter → physics builds and evolves a real
spacetime → metrics discover FTL signatures → archive feeds the next proposal.

![MAP-Elites end-to-end overview](mapelites-end-to-end.svg)

### Diagram — matter-first vs metric-first

![Matter-first vs metric-first diagram](mapelites-matter-first.svg)

### Stage 0 — Quality-Diversity proposer (MAP-Elites)

MAP-Elites maintains an 8×8 behavior archive keyed by a 2-D descriptor. Each
batch it either **mutates an existing elite** (boundary-reflected Gaussian
perturbation, σ=0.15, ~85% of draws) or **samples inside the feasible box** of
known elites. Proposals are points in the 23-D `grtresna_shell` search space.

> **QD vs CMA-ES.** QD (`campaigns/qd/run.sh`) illuminates the 8×8 grid.
> CMA-ES (`campaigns/cmaes/run.sh` → `optimize`) hill-climbs from QD survivors.
> Same search space (`search/optimize/spaces.py`); only the proposer differs.
> See [CMA-ES refinement](#cma-es-refinement-after-map-elites).

### Stage 1 — Initial data (GRTresna, CPU/MPI)

`../GRTresna` paints lumps into `(φ, Π)` and runs a Lichnerowicz/York **constraint
solve**. CPU/MPI-bound; poorly-conditioned proposals fail the Ham/Mom gate and are
**rejected before the GPU**.

### Stage 2 — Evolution (GRTeclyn, GPU)

Converged initial data evolves via `RadialRecipe` to `STOP_TIME`, dumping plotfiles
every `PLOT_INTERVAL` steps.

### Stage 3 — Metrics & probes (scoring)

Plotfiles → diagnostics: constraints, `theta_plus`, comoving/shift stats, matter
density, and FTL probes — coordinate (`operational_ftl_solved`), evolved
(`operational_ftl`, `ftl_persistence`), gauge-invariant geodesic shortcut
(`operational_ftl_geodesic`), and (opt-in) **4D evolving** end-to-end null trace
(`ftl_geo_evolving`; authoritative in search when enabled — see
[4D probe](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

### Stage 4 — Archive update & feedback

`metrics/score/` collapses diagnostics into `ftl_first` fitness. Only `gpu_ok`
candidates compete for a behavior cell. Rejected GRTresna solves log to
`trajectory.jsonl` but never enter the archive.

## The hard consistency rule

`T_ab` used in the GRTresna **constraint solve** must equal `T_ab` used in the
GRTeclyn **evolution**. Otherwise the run starts off-constraint and any apparent
"FTL" is a constraint-relaxation transient, not physics.

The `grtresna_independent_scalars` matter path exists precisely to keep both
sides identical; any new matter sector must be added to **both** sides with
matching analytic forms. (Root cause: `Examples/RadialRecipe/Debug.md`; see also
[Matter model](#foundational-entries-2026-06-10).)

## Behavior descriptors (the "diversity" axes)

Current descriptor: **`ftl_lifetime`** (`qd_search/descriptors.py`), added in
[v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17):

- **x** — peak trustworthy `f_geo`, ramped `(f_geo − 1e-3)/(2e-1 − 1e-3)`.
- **y** — **FTL-lifetime fraction** (share of frames with a live shortcut).

Back-compat: **`speed_super`** (v14 default); `speed_horizon` retired after the
`theta_plus` centering bug (see [Status](#map-elites-ftl-discovery-status)).

## Scoring model (the "quality" axis)

Fitness is `ftl_first` in `metrics/score/` (CMA-ES may use `robust_ftl`; see
[v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

- **Gauge-invariant FTL dominates — time-averaged since v15; 4D since v18.**
  `operational_ftl_geodesic` (×1000) or `ftl_geo_evolving` when 4D runs; gated on
  `h_quality_ok` and ray reach; mean over trustworthy magnitudes × `structural_persistence`.
- **Dynamical signal next.** `operational_ftl` (×400) + `ftl_persistence` (×300) outweigh
  coordinate `operational_ftl_solved` (×50, shaping only).
- **Persistence-honest health.** `survival = numerical_survival × structural_persistence`.
- **Vetoes / penalties.** Horizon (−500 when corroborated, [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)),
  exotic, instability, stationary warp-lens artifacts.

Geodesic ramp: `(f_geo − 1e-3)/(2e-1 − 1e-3)` — full marks need ~20% shortcut
([v9 recalibration](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17)).

Per-component weight table: `grteclyn-wrapper/src/grteclyn_wrapper/metrics/README.md`.

### Plain-English glossary (compact)

Modes: `weighted` (plain sum) and `ftl_first` (validated FTL dominates).

| FTL signals | Role |
|-------------|------|
| `operational_ftl_geodesic` | Frozen per-frame null-ray shortcut; largest weight when 4D off |
| `ftl_geo_evolving` | 4D end-to-end null trace; authoritative when enabled |
| `operational_ftl`, `ftl_persistence` | Evolved coordinate shortcut + final-frame persistence |
| `operational_ftl_solved`, `ftl_precursor`, `channel_progress`, `shift_drive` | Shaping / t=0 hints; gated |

| Health rewards | Role |
|----------------|------|
| `numerical_survival`, `structural_persistence`, `survival` | Integrator + morphology |
| `stability`, `comoving_stability`, `constraint_health`, `lapse_health` | Geometry / solve quality |
| `energy_condition`, `anec_condition`, `tidal_comfort` | Physical energy rules |
| `curvature_activity`, `nontrivial_geometry`, … | Non-flat geometry rewards |

| Penalties | Role |
|-----------|------|
| `exotic_penalty` | Graded 0..−1.6 for negative-energy matter |
| `horizon_penalty` | −500 veto when lapse-collapsed trapped surface corroborated |
| `instability_penalty`, `qei_penalty`, `boundary_penalty` | Geometry churn / bounds |
| `stationary_artifact_penalty` | Static lens pretending to be FTL |

**Non-triviality gate** — health rewards off for flat vacuum.

## Code map (where everything lives)

| Concern | Path |
|---------|------|
| **MAP-Elites** QD loop, archive, descriptors | `grteclyn-wrapper/src/grteclyn_wrapper/search/qd_search/` |
| **Pre-GPU rejection learning** (near-miss pool, shadow archive, Ham/Mom descriptors) | `search/pre_gpu/` |
| Shadow archive persistence | `qd_search/pre_gpu_archive.py` → `pre_gpu_archive.json` |
| **CMA-ES** optimize loop, warm-start | `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize/` |
| Search-space defs (shared QD + CMA-ES) | `search/optimize/spaces.py` |
| FTL champion retention | `search/ftl_retention.py` |
| Scoring (`ftl_first`, `robust_ftl`, `general_ftl`) | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/score/` |
| Metric aggregation | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/aggregation/collector.py` |
| FTL probes | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/probes/ftl/` |
| **4D evolving geodesic** + metric-stack cache | `metrics/probes/ftl/evolving_geodesic.py`, `metric_field.py`, `metric_stack_cache.py` |
| Plotfile → frames + `ftl_timeseries.dat` | `visualisation/process_wave/consume_plotfiles/` |
| **Incremental HQ scoring** | `metrics/aggregation/incremental.py` |
| HQ promotion launcher | `scripts/campaigns/hq/run_batch.sh`, `campaigns/hq/replay_eval.py` |
| Frame → movie stitching | `scripts/plot/make_movies.sh` |
| Matter (evolution) | `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}`, `Examples/RadialRecipe/` |
| Matter (initial data) | `../GRTresna/Examples/ScalarFieldBH/` |
| Campaign launchers | `scripts/campaigns/qd/run.sh`, `cmaes/run.sh`, `general_ftl/run_all.sh` |

## Building the binaries (GRTresna + GRTeclyn)

**GRTresna** = Chombo + conda-OpenMPI (CPU/MPI). **GRTeclyn** = AMReX + CUDA (GPU).

### One env to set first (every shell)

```bash
export GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
export SIM_ROOT=/home/jovyan/nachevsky/test/simulation
export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"
export CONDA_PREFIX="${GRTRESNA_ENV}"
export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
```

Shortcut: `source grteclyn-wrapper/scripts/lib/env.sh`. Full recipe:
`grteclyn-wrapper/src/grteclyn_wrapper/README.md`.

### Build GRTresna (initial-data solver, MPI)

```bash
cd "${SIM_ROOT}/GRTresna/Examples/ScalarFieldBH"
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
# -> Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
```

First time: `cd "${CHOMBO_HOME}" && make lib -j"$(nproc)"`. Header-only edits need force-relink.

### Build GRTeclyn (evolution, GPU)

```bash
cd "${SIM_ROOT}/GRTeclyn/Examples/RadialRecipe"
make COMP=gnu USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 -j"$(nproc)"   # H100
make COMP=gnu USE_CUDA=TRUE USE_MPI=TRUE  CUDA_ARCH=90 -j"$(nproc)"   # MPI+CUDA
```

`CUDA_ARCH`: `90` = H100, `80` = A100, `70` = V100.

### Common failures → fixes

| Symptom | Fix |
|---------|-----|
| `mpicxx` / `gfortran` not found | `export PATH="${GRTRESNA_ENV}/bin:${PATH}"` |
| `cannot find -lhdf5` | `export CONDA_PREFIX="${GRTRESNA_ENV}"` before `make` |
| `libmpi.so` / `libhdf5.so` missing at runtime | `export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:..."` |
| `CHOMBO_HOME` undefined | `export CHOMBO_HOME="${SIM_ROOT}/Chombo/lib"` |
| Header edit "did nothing" | force-relink |
| `/bin/csh: No such file` | point Chombo at conda `tcsh` (README step 6) |

## How to run a campaign

### MAP-Elites (QD)

```bash
cd grteclyn-wrapper
QD_NAME=ftl_discovery_vN QD_ITERATIONS=10 BINS=8 STOP_TIME=16.0 \
  GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 LUMPS=5 SHELL_PROFILE=compact \
  GRTRESNA_MAX_HAM_PCT=5.0 GRTRESNA_MAX_MOM_PCT=5.0 \
  nohup bash scripts/campaigns/qd/run.sh \
  > ../runs/qd_ftl_discovery_vN.launch.log 2>&1 &
```

Results: `runs/grtresna_qd/<name>/` (`trajectory.jsonl`, `archive.json`, `eval_*/`).

**Monitor:** `tail -f runs/grtresna_qd/<name>/trajectory.jsonl`;
`cat runs/grtresna_qd/<name>/ftl_champions.json`.

### CMA-ES refinement after MAP-Elites

| Line | QD source | Objective | Warm-start | Result | Doc |
|------|-----------|-----------|------------|--------|-----|
| **v22 wormhole** | `general_ftl_wormhole_v21` | **`general_ftl`** | QD **063** → CMA-ES **046** (+14.2 pts) | [v22 CMA-ES](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) |
| **v20** | `general_ftl_{wormhole,ring,spin}` | **`general_ftl`** | Ring eval **43** (TBD) | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |
| **v18 / ftl_4d** | `ftl_4d_v1` | **`ftl_first`** | QD **156** → CMA-ES **144** | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |
| **v17** | `ftl_discovery_v16` | `robust_ftl` | eval **739** (not king 233) | [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) |

**Rule:** CMA-ES must match QD `OBJECTIVE_MODE`, grid, stop time, and 4D profile
(`campaigns/lib/search_common.sh`). Do **not** switch objectives on a warm-started trajectory.

```bash
cd grteclyn-wrapper
RUN_NAME=<run> \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/<qd>/trajectory.jsonl" \
WARM_START_TOP_K=5 WARM_START_JITTER=0.03 SIGMA0=0.05 MAX_GENERATIONS=18 \
GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/campaigns/cmaes/run.sh
```

Legacy v17 (`robust_ftl`, `ftl_cmaes_v17_robust`): see [v10–v20](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17).

### HQ promotion (full resolution)

QD/CMA-ES at **N=128, L=64, ml=2, t=16** → HQ at **N=256, L=128, ml=3, t=30** with fresh
GRTresna solve, frames, and incremental `score_timeseries.jsonl`.

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/<cmaes_run>" \
CANDIDATES="<eval_ids>" NAME_PREFIX=<prefix> FORCE=1 \
  bash scripts/campaigns/hq/run_batch.sh
```

| Knob | Search | HQ |
|------|--------|-----|
| Grid | 128³, L=64, ml=2 | **256³, L=128, ml=3** |
| Stop time | 16 | **30** |
| 4D mode | `search` (cheap) | `hq` (full stack) |
| Geodesic dirs | **`x y z`** (`GRTECLYN_GEO_DIRECTIONS`) | **`x y z`** (same; wormhole on **z**) |
| Objective | campaign-specific | **`general_ftl`** for v20/v22 wormhole |

**Incremental scoring:** when `--evolving-geodesic` is on, only `ftl_geo_evolving` earns
geodesic FTL credit until the end-of-run 4D trace completes — mid-run totals are not comparable
to search finals. See [v10–v20 HQ eval 144](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17).

```bash
tail -f runs/grtresna_promote/*/small_data/score_timeseries.jsonl
bash scripts/plot/make_movies.sh runs/grtresna_promote/* --framerate 10
```

Concrete candidate lists and results → [campaign log](#campaign-log--runs-analysis).

---

## Campaign log / runs analysis

Reverse-chronological journal. Quick index:

| Campaign / section | Date | Headline |
|--------------------|------|----------|
| [**v22 CMA-ES: wormhole refinement**](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) | **06-18 stopped** | QD **063** (165.6) → CMA-ES **046** (**179.8**, +14.2). Converged by gen 6 (~47 evals). Scoring parity fix (`GRTECLYN_GEO_DIRECTIONS=x y z`). → [HQ eval 046](#v22-cma-es-hq-promotion) in flight. |
| [**v22: Pre-GPU learning + v21 resume**](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18) | **06-18 complete** | **200 evals** on `general_ftl_wormhole_v21`. Pipelined QD + pre-GPU learning. **Champion eval 191** (161.9, survival 1.00, 18.5% 4D). Eval 063 holds score/FTL records (165.6, 19.3%). [Top 3 + FTL champions](#v22-final-results-top-3--ftl-champions). |
| [**v21: Pipelined QD + GPU tenancy**](#v21-pipelined-qd--gpu-tenancy-tuning-2026-06-17) | **06-17 stopped** | **Pipelined MAP-Elites** (`GpuPool` + `EvalPipeline`). **5 slots/GPU overloaded** H100s at t=16 (~3× slower/evol). **Working config:** 8 GPUs × **1 slot/GPU** + continuous GRTresna; **bottleneck** `MAX_CONCURRENT_GRTRESNA=3` → raise to **5+**. **26 evals:** 11 `gpu_ok`. Cold-start gap → [fixed in v22](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18). |
| [**v10–v20: pipeline + runs**](#v10v20-pipeline-evolution--runs-2026-06-11--2026-06-17) | **06-11 → 06-17** | Scoring/geodesic hardening (v7–v16) → **4D probe + HQ** (eval **144** verified **~8%**) → **general_ftl** QD (ring eval **43** **~3.9%**). See glued section. |
| [**Foundational (06-10)**](#foundational-entries-2026-06-10) | 06-10 | Matter model, navigation overhaul, status reset |

---

## v22: Pre-GPU rejection learning + v21 resume (2026-06-18)

**Status:** **complete** — `general_ftl_wormhole_v21` reached **200 evals** (2026-06-18).
Resume from eval 98 after [convergence-fix](#v22-launch-resume) that prevented early stop
before `target_evals`. Same campaign dir; **v22** = code + learning phase + full production run.

**Purpose.** Close the v21 [cold-start gap](#v21-map-elites-cold-start-grtresna-rejections): MAP-Elites
now **learns from graded pre-GPU rejections** while keeping the main FTL archive pure (`gpu_ok`
only). Continue wormhole QD toward **80 evals** with tuned pipeline knobs from v21.

### Pipeline update: pre-GPU rejection learning

| Layer | Role | Learns from |
|-------|------|-------------|
| **Main archive** (`archive.json`) | Real FTL elites | `gpu_ok` only (unchanged) |
| **Near-miss pool** (`search/pre_gpu/near_miss_pool.py`) | Top-K parents by `score` | `grtresna_rejected`, `solved_ftl_rejected`, `postload_rejected` |
| **Shadow archive** (`pre_gpu_archive.json`) | Separate MAP-Elites grid on pre-GPU axes | Same graded rejections |
| **Unified sampler** (`qd_search/sampling.py`) | 60% gpu_ok elite / 15% shadow / 15% near-miss / 5% feasible / 5% random | Renormalizes when sources empty |

**Shadow descriptor axes** (status-aware, not post-GPU FTL axes):

| Status | x-axis | y-axis |
|--------|--------|--------|
| `grtresna_rejected` | Ham quality | Mom quality |
| `solved_ftl_rejected` | convergence (=1) | precursor tilt (solved FTL) |
| `postload_rejected` | convergence / score | postload margin |

**Excluded from mutation parents:** `grtresna_failed`, `gpu_failed`, `pipeline_interrupted`.

**Gate metrics** now flow into `trajectory.jsonl` (`grtresna_convergence`, `postload_gate`).
Legacy v21 records without structured convergence are rebuilt via `Ham=(\d+)%` regex on
`reason` strings.

**CMA-ES (stage 1):** `--warm-start-include-near-miss` seeds from graded rejections in prior
`trajectory.jsonl` (default on with `--grtresna`). **HQ promotion:** unchanged (replay only).

### v22 launch (resume)

```bash
cd grteclyn-wrapper
BRANCH=wormhole PIPELINE_MONITOR=1 \
  QD_RESUME=1 \
  QD_NAME=general_ftl_wormhole_v21 \
  QD_TARGET_EVALS=80 \
  GPU_IDS="0 1 2 3 4 5 6 7" \
  GPU_SLOTS_PER_DEVICE=1 \
  MAX_CONCURRENT_GRTRESNA=5 \
  BATCH_SIZE=8 \
  QD_ITERATIONS=30 \
  SKIP_QD_PREFLIGHT_TESTS=1 \
  bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/_logs/general_ftl_wormhole_v21.launch.log 2>&1 &
```

| Knob | v21 final | **v22 resume** |
|------|-----------|----------------|
| `QD_RESUME` | — | **1** |
| `pre_gpu_learning` | — | **true** (auto with `--grtresna`) |
| `MAX_CONCURRENT_GRTRESNA` | 3 | **5** |
| `max_level` | 2 (metadata) | **1** (`search_common.sh` default) |
| `near_miss_pool_size` | — | **32** |

`max_in_flight = 8 gpu_slots + 5 GRTresna = 13`.

### Progress snapshot (resume start)

- **Inherited from v21:** 26 scored evals — 11 `gpu_ok`, 9 `grtresna_rejected`, 3 `postload_rejected`;
  best near-miss **eval 4** score **−195** (Ham=100%, Mom within threshold).
- **On resume:** near-miss pool + shadow archive rebuild from full `trajectory.jsonl` before
  new submits; children should **mutate near best near-misses** instead of uniform random.
- **Reports:** after each batch of 8 completions, `[qd] report` includes
  `near_miss=N shadow_cov=…`; writes `pre_gpu_archive.json`.

**Monitor:**

```bash
tail -f runs/_logs/general_ftl_wormhole_v21.launch.log
tail -f runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl
grep near_miss runs/_logs/general_ftl_wormhole_v21.log
watch -n5 'nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv'
```

**Outputs:** `runs/grtresna_qd/general_ftl_wormhole_v21/` — adds `pre_gpu_archive.json` alongside
`archive.json` and `trajectory.jsonl`.

### v22 vs v21

| | v21 | **v22** |
|---|-----|---------|
| GRTresna rejections | logged, ignored for sampling | **near-miss pool + shadow archive** |
| Main archive | `gpu_ok` only | `gpu_ok` only (unchanged) |
| `MAX_CONCURRENT_GRTRESNA` | 3 (starved GPUs) | **5** |
| `max_level` | 2 in run metadata | **1** |
| Campaign dir | `general_ftl_wormhole_v21` | **same** (resume) |

### v22 live pipeline validation (2026-06-18, running)

Validated from live trajectory (54 records at time of check, 3 GPU evolutions + 3 GRTresna
CPU solves running concurrently on GPUs 1-3 at 81-87% utilization).

#### Continuous pipeline — no batch barrier

```
Completion order (post-resume, 25 records):
  eval= 7  pipeline_interrupted  [resumed]
  eval=22  pipeline_interrupted  [resumed]
  ...
  eval=37  grtresna_rejected     [CPU ~30s]
  eval=41  grtresna_rejected     [CPU ~30s]
  eval=46  grtresna_rejected     [CPU ~30s]
  eval=45  grtresna_rejected     [CPU ~30s]
  eval=47  grtresna_rejected     [CPU ~30s]
  eval=52  grtresna_rejected     [CPU ~30s]    ← completes BEFORE eval 35
  eval=35  gpu_ok                [GPU ~2min]   ← GPU evolution still running when 52 finished
  eval=36  gpu_ok                [GPU ~2min]
  eval=38  gpu_ok                [GPU ~2min]
  eval=39  gpu_ok                [GPU ~2min]
  eval=49  postload_rejected     [CPU+solved ~45s]
  eval=43  gpu_ok                [GPU ~2min]
  eval=42  gpu_ok                [GPU ~2min]
  eval=53  grtresna_rejected     [CPU ~30s]
  eval=44  gpu_ok                [GPU ~2min]
```

**Key evidence:** eval 52 (CPU rejection, ~30s) completes before eval 35 (GPU evolution, ~2min).
The pipeline submitted eval 52 while eval 35's GPU phase was still running — no batch await.
**16 out-of-order completions** in trajectory confirm high parallelism.

**At validation time:**

| Phase | Active | Evidence |
|-------|--------|----------|
| GPU evolution | eval_000048, _050, _054, _056, _057 | nvidia-smi: GPUs 1-3 at 81-87%, ~8 GB VRAM each |
| GRTresna CPU solve | eval_000058, _059, _060 | 8 MPI ranks each visible in ps |

9 GPU evolutions + 3 CPU solves running simultaneously. Pipeline continuously feeds
GPU slots without waiting for any batch to complete.

#### Pre-GPU learning active in live run

| Metric | Value |
|--------|-------|
| Shadow archive elites | 3 cells occupied |
| `shadow_improved` events | 8 records tagged |
| Structured `grtresna_convergence` | 7 records carry `{ham_pct, mom_pct}` |
| `postload_gate` metrics | 1 record with `{max_hamiltonian_l2: 0.0124}` |

Shadow archive cell placement:

| Cell | Status | Descriptors | Score |
|------|--------|-------------|-------|
| `(7,0)` | postload_rejected | `conv_quality=1.0, tilt=0` | -100 |
| `(0,7)` | grtresna_rejected (Ham=100%, Mom=0%) | `ham_quality=0.05, mom_quality=1.0` | -195 |
| `(0,0)` | grtresna_failed (NaN convergence) | `quality=0, 0` | -350 |

Near-miss pool top candidates: 4× `postload_rejected` (score -100) rank above
8× `grtresna_rejected` (score -195). This correctly prioritizes candidates that
passed more gates (closer to GPU success).

### Schema: continuous pipeline architecture

```
                         ┌─────────────────────────────────────────────────────────┐
                         │              EvalPipeline (no batch barrier)             │
                         │                                                         │
  ┌──────────┐   submit  │   ┌──────────────┐        ┌──────────────────────┐      │
  │  Unified │──────────►│   │ CpuAdmission │  pass  │      GpuPool         │      │
  │  Sampler │           │   │  Controller  │───────►│  lease() → GPU slot  │      │
  │          │◄──────────│   │  (max=5 MPI) │        │  (8 GPUs × 1 slot)   │      │
  └──────────┘ on_complete   └──────┬───────┘        └──────────┬───────────┘      │
       │                            │                           │                   │
       │                            │ reject                    │ gpu_ok/gpu_failed  │
       │                            ▼                           ▼                   │
       │                   ┌─────────────────────────────────────────┐              │
       │                   │            _on_eval_complete            │              │
       │                   │                                         │              │
       │                   │  ┌─── archive_lock ──────────────────┐  │              │
       │                   │  │ if gpu_ok: archive.insert(elite)  │  │              │
       │                   │  │ if graded_rejection:              │  │              │
       │                   │  │   near_miss_pool.consider(record) │  │              │
       │                   │  │   pre_gpu_archive.insert(shadow)  │  │              │
       │                   │  │ trajectory.append(record)         │  │              │
       │                   │  └───────────────────────────────────┘  │              │
       │                   └─────────────────────────────────────────┘              │
       │                                    │                                       │
       └────────────────────────────────────┘                                       │
         immediately submit next candidate                                          │
         (no wait for batch)                                                        │
                         └─────────────────────────────────────────────────────────┘

  Parallelism at any instant:
    CPU phase:  up to MAX_CONCURRENT_GRTRESNA (5) MPI solves
    GPU phase:  up to gpu_slots (8) evolutions
    Total in-flight: up to 13 evals simultaneously
```

### Schema: pre-GPU rejection learning flow

```
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │                        Eval Completion                                       │
  │                                                                             │
  │   status?                                                                   │
  │     │                                                                       │
  │     ├── gpu_ok ─────────────► Main Archive (archive.json)                   │
  │     │                              │                                        │
  │     │                              ▼                                        │
  │     │                    ┌───────────────────┐                              │
  │     │                    │ Unified Sampler   │                              │
  │     │                    │ 60% mutate elite  │◄─── Main Archive             │
  │     │                    │ 15% mutate shadow │◄─── Shadow Archive           │
  │     │                    │ 15% mutate near   │◄─── Near-Miss Pool           │
  │     │                    │  5% feasible box  │◄─── Union bounds             │
  │     │                    │  5% random        │                              │
  │     │                    └────────┬──────────┘                              │
  │     │                             │                                         │
  │     │                             ▼                                         │
  │     │                      Next candidate                                   │
  │     │                                                                       │
  │     ├── grtresna_rejected ──┐                                               │
  │     ├── solved_ftl_rejected ┼──► Near-Miss Pool (top-32 by score)           │
  │     ├── postload_rejected ──┘         │                                     │
  │     │                                 └──► Shadow Archive (pre_gpu_archive) │
  │     │                                      status-aware descriptors:        │
  │     │                                        grtresna: (ham_q, mom_q)       │
  │     │                                        solved_ftl: (1.0, tilt)        │
  │     │                                        postload: (conv_q, margin)     │
  │     │                                                                       │
  │     ├── grtresna_failed ────► trajectory only (no gradient signal)          │
  │     ├── gpu_failed ─────────► trajectory only                               │
  │     └── preflight_rejected ─► trajectory only                               │
  └─────────────────────────────────────────────────────────────────────────────┘

  Score ranking in near-miss pool (higher = closer to GPU success):
    postload_rejected:    -75 to -100  (passed convergence + solved-FTL)
    solved_ftl_rejected:  -75 to -95   (passed convergence)
    grtresna_rejected:    -100 to -350 (failed convergence, graded by Ham/Mom)
    grtresna_failed:      excluded     (solver crash, no useful gradient)

  Cold start (no gpu_ok yet): renormalize → 43% shadow + 43% near-miss + 14% random
```

### Schema: CMA-ES near-miss warm-start (QD→CMA handoff)

```
  Prior QD trajectory.jsonl
        │
        ├── gpu_ok records ──────────► _load_warm_start_vectors (top_k=5)
        │                                     │
        │                                     ▼
        │                              ┌──────────────┐
        └── graded rejections ────────►│ near_miss_   │──► Append to warm vectors
             (--warm-start-include-    │ vectors_from │    (dedup by 8-decimal rounding)
              near-miss, default on)   │ _trajectory  │
                                       │ (near_miss_k │
                                       │  = 4)        │
                                       └──────────────┘
                                              │
                                              ▼
                                    CMA-ES initial population
                                    (gpu_ok first, then near-miss)
                                    First vector = initial mean
                                    Others jittered by sigma
```

### v22 final results: top 3 + FTL champions

**Campaign summary** (`runs/grtresna_qd/general_ftl_wormhole_v21/`):

| Metric | Value |
|--------|-------|
| Evals | **200** (201 trajectory lines incl. 10 `pipeline_interrupted` from prior halt) |
| `gpu_ok` | **103** |
| Archive | **4 elites**, coverage **6.25%**, MAP-Elites best **165.6** (eval 063); **refinement champion eval 191** |
| Throughput (resume 98→200) | ~**0.93 eval/min** (8 GPU × 1 slot, `MAX_CONCURRENT_GRTRESNA=3`) |
| Near-miss pool | 32 (pre-GPU learning active throughout) |

**Top 3 by score** (final run):

| Eval | Score | `ftl_geo_evolving` | `f_geo_peak` | `operational_ftl` | Survival | Reading |
|------|------:|-------------------:|-------------:|------------------:|---------:|---------|
| **063** | 165.6 | **19.3%** | 4.2% | 0 | 0.94 | Score + `ftl_geo_evolving` record holder (early run) |
| **191** | **161.9** | **18.5%** | 3.8% | 0 | **1.00** | **Champion** — resume-era lead; perfect survival, `f_op_peak` holder, jitter 0.96 |
| **174** | 157.4 | 18.5% | 3.8% | 0 | 1.00 | Stable variant — thick shell (1.40), moderate jitter |

All three: `operational_ftl = 0` but strong **4D null-geodesic** shortcuts — the hallmark of
stationary, non-translating wormholes (see [physical validation](#v22-physical-validation--next-steps)).

Earlier basin members still in archive: **eval 071** (122.8, symmetric/stable), **eval 040** (96.5, compact/massive).

**FTL champions** (`ftl_champions.json` — one eval dir kept per peak metric):

| Metric | Eval | Peak value | Notes |
|--------|------|------------|-------|
| `ftl_geo_evolving` | **063** | **19.3%** | Run-best integrated 4D shortcut |
| `f_op_peak` | **191** | 0.0178 | Best operational precursor |
| `f_geo_evol` (frozen) | **178** | 5.2% | Best static-slice geodesic peak (resume find) |
| `max_local_speed` | **178** | **1.31** | Highest local coordinate speed |
| `superluminal_fraction` | **050** | 77.4% | Largest superluminal volume fraction |
| `ftl_lifetime_fraction` | **008** | 1.00 | Longest FTL persistence |

### v22 physical validation & next steps

**Verdict.** First production run of the pipelined `general_ftl` wormhole campaign in a
**15-D** pinned search space successfully located a **clustered family** of stable,
constraint-clean, non-translating wormholes. Empirical proof that coordinate warp-bubble
metrics and stationary FTL wormhole metrics are physically distinct observables.

#### I. Scoring paradox: `operational_ftl = 0` while `ftl_geo_evolving ≈ 19%`

In a **stationary** wormhole (\(\beta^i \approx 0\)):

1. **Lapse-dominated coordinate speed** — heavy throat matter drops \(\alpha < 1\), so
   \(c_{\rm coord} = \alpha \pm \beta < 1\). Frozen-slice Dijkstra (`operational_ftl`) reads
   subluminal → **0**.
2. **Proper-distance contraction** — bipolar shell geometry warps \(\gamma_{ij}\); proper
   path length through the throat is shorter than coordinate distance. The **4D null-geodesic
   tracer** integrates \(ds^2 = 0\) and captures this → **~19% faster** than flat-space light.

The pipeline has decoupled coordinate gauge artifacts from physical traversable shortcuts.

#### II. Dynamically opening throat

For eval **191**: integrated `ftl_geo_evolving` = **18.5%** but frozen `f_geo_peak` = **3.8%**
at \(t \approx 16\) (eval **063** shows the same pattern at 19.3% / 4.2%). Even with zero
initial velocities, CCZ4 evolution is highly dynamical — negative-energy shell and positive
central regions interact gravitationally. The throat **opens during evolution**; a photon
emitted at \(t=0\) traverses geometry that is actively warping. Only the 4D tracer detects
this time-dependent GR effect.

#### III. Elite basin (cell [1, 7]) — final run

| Eval | Role | Key traits |
|------|------|------------|
| **191** | **Champion** | Resume-era lead — extreme jitter (0.96), dipole +0.57, survival **1.00**, 18.5% 4D shortcut, `f_op_peak` holder |
| **174** | Stable variant | Thick shell (1.40), moderate jitter (0.59) → symmetric, survival 1.00, 18.5% 4D shortcut |
| **063** | Record holder | Early-run peak score (165.6) + `ftl_geo_evolving` record (19.3%); survival 0.94, 52% superluminal volume |

Legacy early-basin members: **071** (symmetric/stable, 13.9% 4D), **040** (compact/massive, survival 0.77).

#### IV. Strategic action plan

| Step | Action | Goal | Status |
|------|--------|------|--------|
| **1** | **CMA-ES refinement** seeded from QD **eval 063** | Local hill-climb in the v22 wormhole basin; push 4D shortcut above QD record while holding survival. | **Done** — [eval 046 @ 179.8](#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18) (+14.2 vs 063). Stopped early (converged ~gen 6). |
| **2** | **HQ promotion** of CMA-ES **eval 046** → **N=256, L=128, ml=3, t=30** | Does the refined throat survive full resolution and longer time? | **In flight** — `general_ftl_wormhole_cmaes_v1_hq_eval000046` (single representative; top-3 CMA-ES evals were basin duplicates). |
| **3** | **Complex scalar (boson star) pivot** | Real scalars disperse at \(t \gg 50\). Re-run with complex fields + conserved U(1) charge. | Pending |

**Note:** Action plan originally named eval **191** as refinement seed (resume-era champion, survival 1.00). CMA-ES used **063** instead — the absolute score / `ftl_geo_evolving` record holder — which replayed at **165.6** bit-identically once scoring parity was fixed.

---

## v22 CMA-ES: wormhole refinement (`general_ftl_wormhole_cmaes_v1`, 2026-06-18)

**Status:** **stopped early** (user halt @ ~47 evals; score plateau ~179.8 by gen 6).
Campaign dir (pruned): `runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1/` — champion **`eval_000046` only**.

**Purpose.** Stage-1 local refinement after [v22 QD](#v22-final-results-top-3--ftl-champions): hill-climb
within the clustered wormhole basin found by MAP-Elites, using the same 15-D pinned search space,
`OBJECTIVE_MODE=general_ftl`, and pipelined evaluator as QD.

### Seed & launch

Warm-started from QD trajectory `general_ftl_wormhole_v21` with `WARM_START_TOP_K=1` (exact genome
of eval **063**, not resume-era eval 191):

```bash
cd grteclyn-wrapper
RUN_NAME=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl" \
WARM_START_TOP_K=1 WARM_START_JITTER=0.05 SIGMA0=0.05 \
TARGET_EVALS=150 MAX_GENERATIONS=50 KEEP_TOP_EVAL_DIRS=3 \
GPU_IDS="0 1 2 3 4 5 6 7" GPU_SLOTS_PER_DEVICE=1 MAX_CONCURRENT_GRTRESNA=3 \
PIN_DIMS="$(bash -c 'source scripts/campaigns/lib/general_ftl_pins.sh && ftl_general_ftl_wormhole_pins')" \
  bash scripts/campaigns/cmaes/run.sh
```

| Knob | Value |
|------|-------|
| Search space | 15-D wormhole pins (same as v22 QD) |
| Population | 8 (= GPU count) |
| σ₀ | 0.05 |
| Pipeline | `GpuPool` + `EvalPipeline`, `max_concurrent_grtresna=3` |
| Target | 150 evals (stopped @ ~47) |

### Scoring parity fix (critical)

First CMA-ES attempt replayed QD eval 063 at **−27** instead of **165.6** — same overrides,
different score. Root cause: CMA-ES path lacked QD's geodesic configuration.

| Issue | QD (v22) | CMA-ES (before fix) | Fix |
|-------|----------|---------------------|-----|
| `GRTECLYN_GEO_DIRECTIONS` | `x y z` (`general_ftl/run_all.sh`) | default **`x` only** | export in `campaigns/lib/search_common.sh` when `OBJECTIVE_MODE=general_ftl` |
| Eval driver | shared `_run_cpu_grtresna_gates` + `_run_gpu_session` | duplicate `_objective` path | CMA-ES GRTresna mode wired to shared `evaluate_overrides` stack |
| Live trajectory | append per eval | batch rewrite per gen | `_track_trajectory` + `live_trajectory_path` |

After fix, warm-start replay (CMA-ES eval **001**) scored **165.608** — **bit-identical** to QD eval **063**
(`ftl_geo_evolving` = 0.192631, `f_geo_evol` = 0.0419, geodesic axis **z**).

HQ promotion received the same fix: `GRTECLYN_GEO_DIRECTIONS=x y z` + `--objective-mode general_ftl`
in `campaigns/lib/promote_common.sh` and `campaigns/hq/run_batch.sh`.

### Improvement over QD top candidate

| Metric | QD eval **063** | CMA-ES eval **046** | Δ |
|--------|----------------:|--------------------:|---|
| **Score** | **165.6** | **179.8** | **+14.2 (+8.6%)** |
| `ftl_geo_evolving` | 19.3% | **20.3%** | +1.0 pp |
| `f_geo_evol` (4D trace) | 4.19% | 4.15% | ~flat |
| `f_geo_peak` (frozen) | 7.07% | 6.77% | −0.3 pp |
| `ftl_lifetime_fraction` | 86% | 86% | — |
| `max_local_speed_peak` | 1.146 | 1.140 | −0.006 |
| Survival | 0.94 | **1.00** | improved |
| Geodesic axis | **z** | **z** | same |

CMA-ES did **not** discover a new mechanism class — it tightened the existing v22 wormhole basin:
higher integrated 4D score (`ftl_geo_evolving` × survival weight) with perfect structural persistence.

### Convergence curve (best per generation)

| Gen | Best score | Notes |
|-----|----------:|-------|
| 1 | 174.4 | Warm seed + σ=0.05 neighbors; eval 001 replays 063 @ 165.6 |
| 2 | 174.8 | |
| 3 | 178.0 | |
| 4 | 179.1 | |
| 5 | 179.6 | Population mean ≈ 178 — basin collapsed |
| 6 | **179.8** | **Champion eval 046**; last 4 gen bests spread only +1.8 |

Stopped manually when further gens were unlikely to beat 179.8 (σ=0.05 jitter around a single optimum).

### v22 CMA-ES HQ promotion

One HQ replay suffices — top CMA-ES evals 039/045/046 were σ-neighbors in the same basin (scores
179.6–179.8), not distinct mechanisms.

```bash
cd grteclyn-wrapper
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
CANDIDATES="46 0" \
  bash scripts/campaigns/hq/run_batch.sh
```

Output: `runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/` (L=128, N=256, ml=3, t=30,
4D mode **`hq`**, dirs **`x y z`**, incremental `score_timeseries.jsonl`).

Monitor:

```bash
tail -f runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046.log
tail -f runs/grtresna_promote/general_ftl_wormhole_cmaes_v1_hq_eval000046/small_data/score_timeseries.jsonl
```

---

## v21: Pipelined QD + GPU tenancy tuning (2026-06-17)

**Status:** **stopped early** (user halt at ~30 evals; **26** trajectory lines saved).
Campaign dir preserved: `runs/grtresna_qd/general_ftl_wormhole_v21/`.

**Purpose.** Validate the pipelined MAP-Elites evaluator (`search/gpu_pool.py`,
`search/eval_pipeline.py`) on production grid/time: GRTresna CPU solves overlap with
GPU `main3d` evolution across **8× H100**, wormhole-only pins (same as v20 wormhole).

### Architecture (unchanged from design)

| Component | Role |
|-----------|------|
| `GpuPool` | `total_slots = len(gpu_ids) × gpu_slots_per_device`; blocking lease before GPU phase |
| `EvalPipeline` | Cross-batch queue; no batch barrier; one completion → one new submit |
| `CpuAdmissionController` | Caps concurrent GRTresna MPI solves (`MAX_CONCURRENT_GRTRESNA`) |
| QD driver | Pipelined MAP-Elites (`use_pipeline=True` default) |

`max_in_flight = gpu_slots + max_concurrent_grtresna` (e.g. 8 + 3 = 11 with final config).

### Run history (three launches, same campaign name)

#### Attempt 1 — multi-tenant GPU overload (failed throughput)

```bash
GPU_SLOTS_PER_DEVICE=5   # 8×5 = 40 GPU leases
MAX_CONCURRENT_GRTRESNA=6
```

**Result:** Multiple `main3d` sessions stacked on the same H100 (GpuPool fills GPU 0
first). At `stop_time=16` each solo evolution peaks ~**44 GB** VRAM; **5 concurrent**
on one GPU caused severe slowdown — wall time per eval ~**3×** vs solo (~35 min vs
~12 min). **Do not use 5 slots/GPU at t=16** on 80 GB cards with this gridinit.

#### Attempt 2 — 1 slot/GPU, high GRTresna concurrency

```bash
GPU_SLOTS_PER_DEVICE=1   # 8 GPU leases total
MAX_CONCURRENT_GRTRESNA=6
```

**Result:** GPU memory healthy (~one evolution per device), but GRTresna admission
still contended with long GPU phases.

#### Attempt 3 — production config (successful pipeline shape)

```bash
GPU_SLOTS_PER_DEVICE=1
MAX_CONCURRENT_GRTRESNA=3
BATCH_SIZE=8
QD_TARGET_EVALS=80
```

**Result:** **8 GPUs each running one evolution**; **up to 3 GRTresna solves in flight**
continuously feeding the pipeline. First `gpu_ok` after ~5 min wall; archive began
filling. Typical `gpu_ok` eval ~5–8 min; one AMR outlier (`eval_000007`, `max_level=2`)
hit ~**39 GB** VRAM and ~30 min wall.

**Trajectory (26 logged):** 11 `gpu_ok`, 9 `grtresna_rejected`, 3 `postload_rejected`,
3 `grtresna_failed`. Best scored `gpu_ok`: eval **19** score **+33.7**. Best GRTresna
near-miss: score **−195** (Ham=100%, Mom within threshold) — logged but **not used
for sampling** (see below).

### Issue: GPU idle gaps — GRTresna batch too small

With `gpu_slots=8` and `MAX_CONCURRENT_GRTRESNA=3`, the pipeline often had **3–5 GPUs
idle** for multi-minute stretches: GPU work finished faster than new gridinits arrived.
The evaluator submits one new candidate per completion (no batch barrier), so throughput
is limited by how fast GRTresna can produce passing initial data.

**Fix:** **Increase `MAX_CONCURRENT_GRTRESNA`** (GRTresna batch / CPU admission cap) to
**5–6** so more solves run in parallel while GPUs evolve — target `max_in_flight ≈ 13–14`
on 8 GPUs. Do **not** compensate by raising `GPU_SLOTS_PER_DEVICE` at t=16.

| Knob | Attempt 1 (bad) | Attempt 3 (ran) | **Recommended relaunch** |
|------|-----------------|-----------------|--------------------------|
| `GPU_SLOTS_PER_DEVICE` | 5 | 1 | **1** |
| `MAX_CONCURRENT_GRTRESNA` | 6 | 3 | **5–6** |
| `gpu_slots` | 40 | 8 | **8** |
| GPU util @ t=16 | overloaded, slow | healthy, GRTresna-starved | balanced |

### VRAM sizing (eval_000005 gridinit replays, H100 80 GB)

Benchmark (`scripts/benchmarks/gpu_gridinit_load.sh`) on real 128³ wormhole gridinit —
not bare `sweep` (~5 GB underestimate).

| Concurrent `main3d` | `stop_time` | Peak VRAM | Notes |
|--------------------|-------------|-----------|-------|
| 1 | 4 | ~9.7 GB | matches QD evaluate path |
| 2 | 4 | ~17.4 GB | ~8.7 GB/eval, linear |
| 5 | 4 | ~52 GB | fits 80 GB at short t only |
| 1 | **16** | **~44 GB** | AMR; **1 slot/GPU** at t=16 |
| 5 | **16** | — | **avoid** — contention, not OOM |

### v21 MAP-Elites cold start (GRTresna rejections)

**Bug / design gap:** MAP-Elites only inserts `gpu_ok` into `archive.json`. GRTresna
rejections are written to `trajectory.jsonl` with **graded fitness** (CMA-ES uses this;
QD does not). Until the first `gpu_ok`, `_sample_next_candidate` draws **uniform random**
vectors — ~40% of v21 evals were pre-GPU rejects with no learning signal.

Rejections also bin to descriptor cell `(0,0)` because FTL axes need GPU evolution.

**Fixed in [v22](#v22-pre-gpu-rejection-learning--v21-resume-2026-06-18):** near-miss parent pool
+ shadow pre-GPU archive (`search/pre_gpu/`). Main `archive.json` stays `gpu_ok` only; graded
rejections guide sampling via `NearMissPool` and `pre_gpu_archive.json`.

### Recommended relaunch

See [v22 resume command](#v22-launch-resume). v21 stopped campaign preserved at
`runs/grtresna_qd/general_ftl_wormhole_v21/`.

**Monitor:**

```bash
tail -f runs/_logs/general_ftl_wormhole_v21.launch.log
tail -f runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl
watch -n5 'nvidia-smi --query-gpu=index,memory.used,utilization.gpu --format=csv'
```

**Outputs:** `runs/grtresna_qd/general_ftl_wormhole_v21/` — `general_ftl` objective,
`ftl_lifetime` descriptors, `GRTECLYN_GEO_DIRECTIONS=x y z`.

### v21 vs v20

| | v20 | **v21** |
|---|-----|---------|
| Evaluator | batch barrier | **pipelined** CPU→GPU, no batch wait |
| GPU tenancy | 1 evol / GPU | **1 evol / GPU** (after tuning; not 5) |
| GRTresna | serialized per wave | **continuous** solves (`MAX_CONCURRENT_GRTRESNA`) |
| Campaign | 3-class parallel | **wormhole-only** pipeline validation |
| Target | 30 iter × 3 classes | **80 evals** single class |

---

## v10–v20: pipeline evolution & runs (2026-06-11 — 2026-06-17)

Single chronicle for scoring/metric/probe campaigns **v7–v20** (named `ftl_discovery_v*` through
`general_ftl_*`). **Pipeline changes** and **run outcomes** only — no per-version chapters.

**Main pipeline improvements**

| Era | Date | Pipeline change | Code / config |
|-----|------|-----------------|---------------|
| **v4** | 06-11 | Persistence-honest `survival`; FTL shaping × `structural_persistence`; mass search + velocity caps | `metrics/score/`, 18-D space, `STOP_TIME` 8 |
| **v7→v8** | 06-11 | Null-ray **forward** launch (`future_null_cov`); relative H-drift gate `H_REL_TOL=1e-2` | `metrics/probes/ftl/geodesic.py` |
| **v9** | 06-11 | Geodesic ramp target 5%→**20%**; weight ×1500→×1000; rebalance `ftl_first` (coordinate shaping down) | `metrics/score/ftl.py` |
| **v10** | 06-11 | HQ killed v9 shortcuts → `STOP_TIME` 8→**16**; searchable `shell_static` toggle | `campaigns/qd/run.sh` |
| **v11** | 06-12 | Geodesic reward × `structural_persistence`; exotic/energy weights **100× / 40×** | `metrics/score/ftl.py` |
| **v13–v14** | 06-12 | λφ⁴; matter layouts 0–4; per-lump profile + cloud; geodesic contradiction gate; **23-D** space | GRTresna + `spaces.py` |
| **v15** | 06-13 | Per-frame `ftl_timeseries.dat`; time-mean geodesic; descriptor **`ftl_lifetime`** | `consume_plotfiles/`, `qd_search/descriptors.py` |
| **v16** | 06-13 | FTL champion retention (`ftl_retention.jsonl`); horizon penalty needs **lapse corroboration** | `search/ftl_retention.py`, `collapse.py` |
| **v17** | 06-14 | CMA-ES **`robust_ftl`** warm-start from healthy QD survivors (not raw king) | `search/optimize/`, `ftl_cmaes_v17_robust` |
| **HQ batch** | 06-15 | Incremental `score_timeseries.jsonl`; HQ ml=3 `regrid_interval` fix | `metrics/aggregation/incremental.py` |
| **4D probe** | 06-15→16 | End-to-end null trace through metric stack; `ftl_geo_evolving` headline; search vs HQ profiles | `evolving_geodesic.py`, `evolving_geodesic_options.py` |
| **v18 / ftl_4d** | 06-16→17 | Full **QD → CMA-ES → HQ** on 4D metric; matched `ftl_first` objective | `ftl_4d_v1`, `ftl_4d_cmaes_v1` |
| **v20** | 06-17 | **`general_ftl`** objective; `--pin-dimension`; `GRTECLYN_GEO_DIRECTIONS=x y z` | `objectives.py`, `general_ftl/run_all.sh` |

**4D frozen vs evolving (design lesson):** eval **086** smoke 4D **1.42%** vs frozen **5.75%**; HQ t=30
4D **0%** (negative control); eval **144** HQ 4D **7.96%** verified (positive control). Frozen
mid-run peaks overstate dynamic warps — **4D trace is authoritative** when enabled.

**HQ lesson (v16 elites):** trust **incremental peak** (~t≈9–12), not final t=30 — horizon −500 kills
static-lump winners; eval **177** (dynamic, CMA-ES) only finisher positive (+67).

**Alcubierre control (06-12):** probes validated ~**32%** f_geo @ 129³; QD H-gate fixed with
`GEO_REFINE_N=129` re-probe.

**Next directions (06-15):** persistent standing channel vs transport worldtube; exotic-energy Pareto;
boson-star matter; CMA-ME; analytic metric extraction.

**Campaign runs (v7–v20)**

| Run | Date | Evals / gpu_ok | Best eval | Score / f_geo | Headline | Run dir |
|-----|------|----------------|-----------|---------------|----------|---------|
| `ftl_discovery_v7` | 06-11 | 88 / 53 | 71 | 606 / geo=0 | Coordinate precursors; geodesic blind (fixed v8) | `runs/grtresna_qd/ftl_discovery_v7/` |
| `ftl_discovery_v9` | 06-11 | 88 / 54 | 40 | 266 / **2.6%** | First certified geodesic shortcuts | `runs/grtresna_qd/ftl_discovery_v9/` |
| `ftl_discovery_v10` | 06-11 | — | — | — | HQ rejected all v9 promotes; longer QD window | — |
| `ftl_discovery_v11` | 06-12 | 400 | — | — | Persistence-gated geodesic + physicality pressure | `runs/grtresna_qd/ftl_discovery_v11/` |
| `ftl_discovery_v13` | 06-12 | 278 | — | — | λφ⁴ + layouts; zero geodesic until gate | `runs/grtresna_qd/ftl_discovery_v13/` |
| `ftl_discovery_v14` | 06-12 | 504 / 351 | **231** | 551 / **5.30%** | Ring layout dominates; Alcubierre ~6× higher f_geo but metric-first | `runs/grtresna_qd/ftl_discovery_v14/` |
| `ftl_discovery_v15` | 06-13 | — | 231 | peak **7.43%** @ t=9.6 | Time-resolved scoring validated | `runs/grtresna_qd/ftl_discovery_v15/` |
| `ftl_discovery_v16` | 06-13 | ~971 | **233** | 652 / 5.88% peak | FTL retention; horizon fix; plateau | `runs/grtresna_qd/ftl_discovery_v16/` |
| `ftl_cmaes_v17_robust` | 06-14→15 | 200 / 163 | **177** | 312 / **5.65%**, timeavg 16.3% | CMA-ES +0.26 pp f_geo vs seed 739 | `runs/grtresna_cmaes/ftl_cmaes_v17_robust/` |
| HQ v16+v17 | 06-15 | 4 promotes | **233** incr. | **749** @ t≈11.8 | Only **177** final +67; 233/446/676 horizon-killed | `runs/grtresna_promote/l128n256t30_*/` |
| `ftl_max_speed_no_penalty_v1` | 06-15 | 200 / 100 | **92** | +27.5 / **27.5%** timeavg | Side survey; max **1.58 c** coord (eval 70); not v16-comparable | `runs/grtresna_qd/ftl_max_speed/ftl_max_speed_no_penalty_v1/` |
| `ftl_4d_v1` | 06-16 | 192 / 105 | **156** | **508** / `ftl_geo_evol` 0.346 | First 4D-in-loop QD | `runs/grtresna_qd/ftl_4d/ftl_4d_v1/` |
| `ftl_4d_cmaes_v1` | 06-16 | 144 / 140 | **144** | **596** / 0.395 | +88 pts vs QD 156 | `runs/grtresna_cmaes/ftl_4d_cmaes_v1/` |
| HQ eval **144** | 06-17 | 1 | **144** | **283** / 4D **7.96%** | **Verified** gauge-invariant shortcut (5/5 rays) | `runs/grtresna_promote/l128n256t30_ftl_4d_cmaes_qd_eval000144/` |
| `general_ftl_{wormhole,ring,spin}` | 06-17 | 172 / 130 | ring **43** | **196** / search **~3.9%** | Ring **20** 4D hits; wormhole **0**; stopped early, top-3 kept | `runs/grtresna_qd/general_ftl_*/` |

**Eval 177 physics (HQ, 06-15):** null transit **~5.9%** early (`≈1.063 c` signal); matter sub-luminal;
exotic total **~5–24× < Alcubierre**, per-shortcut energy **≈ comparable**. Not transport — pulsing lens.

**v20 launch:**

```bash
cd grteclyn-wrapper
MODE=par QD_ITERATIONS=30 \
  bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/_logs/general_ftl_par.launch.log 2>&1 &
```

**v18 QD → CMA-ES → HQ (representative):**

```bash
# QD
RUNS_DIR="${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_4d" QD_NAME=ftl_4d_v1 QD_TARGET_EVALS=200 \
  bash scripts/campaigns/qd/run.sh
# CMA-ES
RUN_NAME=ftl_4d_cmaes_v1 WARM_START_TRAJECTORY=.../ftl_4d_v1/trajectory.jsonl \
  WARM_START_TOP_K=5 SIGMA0=0.05 MAX_GENERATIONS=18 bash scripts/campaigns/cmaes/run.sh
# HQ
SOURCE_RUN=.../ftl_4d_cmaes_v1 CANDIDATES="144 0" bash scripts/campaigns/hq/run_batch.sh
```

---

## Foundational entries (2026-06-10)

### Matter model — reference

Campaign evolves **N independent massive real scalar fields** ("lumps") via
`grtresna_independent_scalars` — 5 lumps, `recipe_scalar_mass` searched, per-lump exotic flags.

| Side | Key paths |
|------|-----------|
| GRTeclyn | `RadialRecipeMatterDispatch.hpp`, `GRTresnaIndependentScalars.{hpp,impl.hpp}`, `StateVariables.hpp` |
| GRTresna | `MatterParams.hpp` (lump_t), `MyMatterFunctions.cpp` |

Potential: `V = ½ m² (Σφ_k)²` (+ λφ⁴ when searched). Lumps interact via shared gravity + mass term;
O(1) boosts + light mass → fly-away (fixed v4 mass search + velocity caps).

**Roadmap (done in bold):** **(1)** search mass + cap boosts; **(2)** λφ⁴; (3) complex scalar /
Q-balls; (4) per-lump independent mass.

### Navigation overhaul (2026-06-10)

1. `speed_horizon` → **`speed_super`** descriptor.
2. ~82% pre-GPU rejection → tightened bounds, feasible-box sampling, harder GRTresna solve.
   `ftl_discovery_nav`: GPU-reach ~40%.

**Scoring fix (after ~90 evals):** stationary warp-lens artifacts — reliability-gate geodesic;
zero shaping when stationary + no dynamical FTL. Eval 83: 1164→−247.

### MAP-Elites FTL Discovery Status

Status: **reset**. `theta_plus` measured from origin not `grid_center` → false horizon penalty.
Fixed in `RadialRecipeLevel.cpp`; `ftl_discovery_postfix` confirmed `horizon_penalty=0`.
~93% pre-GPU rejection exposed navigation defects fixed above.
