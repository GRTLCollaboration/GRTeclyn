# grteclyn-wrapper

Python orchestration layer that runs isolated **GRTeclyn** GPU episodes from the
repo root (`GRTeclyn/`). For each candidate spacetime it asks the sibling
**GRTresna** elliptic solver for constraint-satisfying initial data, hands that
data to GRTeclyn for time evolution, streams plotfiles into `small_data/` +
`frames/` during the GPU run, deletes the heavy HDF5 dirs afterward, and scores
the result with a pluggable objective. On top of this episode loop sit
quality-diversity (MAP-Elites) and CMA-ES search campaigns that explore
matter/geometry configurations for FTL shortcuts, gravitational-wave beaming,
spacetime shear, and critical collapse.

All commands run from **`grteclyn-wrapper/`**. Binaries must be built first —
see [Operations](#operations). Site layout (sibling repos, OpenMPI, GRTresna
env) is configured with a gitignored [`.env`](#site-paths-env) — see below.

> **Roadmap:** critical review of the results' validity and the prioritized
> implementation plan (probe calibration, gauge-honest baseline, ANEC/QI,
> convergence, pump accounting, self-grav, wave-zone GW) live in
> [`NextSteps.md`](NextSteps.md).

---

## What is implemented

| Capability | Where | Notes |
|------------|-------|-------|
| **GRTresna bridge** — `params.txt` → MPI solve → Chombo HDF5 → `.gridinit` | `src/grteclyn_wrapper/grtresna/` | Ham + Mom constraint solve; exotic (`rho<0`) auto-switches to K=0 maximal slicing. Deep docs: [`src/grteclyn_wrapper/grtresna/README.md`](src/grteclyn_wrapper/grtresna/README.md) |
| **Matter sectors** — real scalar, complex scalar / U(1) boson star, Q-ball, shell, trajectory, rotating Q-torus eigenstate | `grtresna/matter/`, `grtresna/fields/`, `grtresna/profiles/` | Per-lump signed scalars (`GRTresnaIndependentScalars` C++); exotic wedge; self-gravitating boson-star ODE solver (isotropic coords); 2D spinning Q-ball (`qball_torus.py`, `profile==4`) |
| **Matter-profile contract (rail)** — single source of truth + t=0 cross-code consistency gate | `grtresna/matter/profile_contract.py`, `tests/grtresna/test_profile_contract.py` | Catches stale-binary empty fields / sign flips / wrong-ω / center bugs before evolution; see [Matter-profile contract](#matter-profile-contract--the-rail-for-adding-new-profiles-safely) |
| **Plotfile consumer** — streaming `small_data/` + PNG `frames/` + HDF5 deletion | `scripts/lib/`, `src/.../visualisation/` | `consume_plotfiles` sidecar; **required** for every production run |
| **Ψ₄ / GW extraction** — in-code C++ `WeylExtraction` (spherical-harmonic modes) | `Examples/RotatingWormholeCollapse/`, `src/.../visualisation/process_wave/` | **Primary: in-code GRTeclyn `SphericalExtraction`** → `data/Weyl4_mode_2{0,1,2}.dat`, dense (every coarse step), multi-radius, decoupled from plotfiles. Python `process_wave` sidecar still extracts a coarse cross-check + drives frames |
| **Search algorithms** — MAP-Elites (QD) archive, CMA-ES hill-climb | `src/.../search/qd_search/`, `src/.../search/optimize/` | Shared pre-evolution gates; warm-start from any trajectory |
| **Objectives** — `ftl_first`, `robust_ftl`, `general_ftl`, `f_geo_max`, `f_geo_depth`, `critical_collapse`, `gw_beam`, `spacetime_shear` | `src/.../metrics/score/objectives.py` | See [Campaigns](#campaigns) for which objective each campaign uses |
| **Descriptors** — `ftl_lifetime`, `speed_horizon`, `wave_focusing`, `spacetime_shear`, `gw_beam` | `src/.../search/qd_search/descriptors.py` | Behavior axes for the MAP-Elites archive |
| **4D null-geodesic probe** — gauge-invariant FTL shortcut measurement | `src/.../metrics/probes/ftl/` | `search` (cheap) and `hq` (full verify) profiles; continuous emission sweep |
| **Falsification tiers** — T0 constructed → T6 analytic | `scripts/search/validate_tiers.py` | Offline ladder; no rerun needed |
| **Geometry-first projection** — motif scout → GRTresna solve | `src/.../initial_data/motif.py`, `grtresna/fit/motif.py`, `projection/` | Additive second stage; never push fitted matter directly into GRTeclyn |
| **Pure-geometry MAP-Elites atlas** — Stage-1 stationary metric scout | `src/.../search/geometry_atlas/`, `scripts/campaigns/geometry_atlas/run.sh` | Searches broad asymptotically flat 4-metrics (no matter); scores frozen `f_geo` + stationary `f_ff` vs exotic-energy cost; see [`docs/GeometryFirst.md`](docs/GeometryFirst.md) |
| **Iterative matter adjustment** — CMA-ES loop over lump params (GRTresna-only) | `projection/iterate.py`, `projection/mismatch.py` | `--iterate N` on `project_geometry_motif.py`; L2 geometry-mismatch fitness; closes the fit→solve→compare loop |
| **Post-load constraint gate** — short GPU load check of `.gridinit` | `projection/postload_gate.py` | Rejects bad loads before the expensive main evolution |
| **Solved-FTL gate** — cheap t=0 filter on `.gridinit` | `search/solved_ftl_gate.py` | ~1 s/candidate; rejects flat/degenerate slices |
| **Visualization** — constraint plots, GW panels, frame movies | `src/.../visualisation/`, `scripts/plot/` | Article-style 6-panel Ψ₄ figures; `make_movies.sh` |
| **LIGO matched-filter search** — methodology + reference script | `src/.../gw_search/` | Intermediate-mass collapse templates vs GWOSC; see [`src/.../gw_search/README.md`](src/grteclyn_wrapper/gw_search/README.md) |
| **RL chassis** — pump/transport training launcher | `scripts/campaigns/rl/` | Opt-in stage 1.5; handoff in [`research/RL/LabJournal.md`](../research/RL/LabJournal.md) |
| **Self-gravitating boson star seed** — stationary single-star confinement | `grtresna/profiles/boson_star_ode.py`, GRTresna C++ | Four-bug fix; see [Self-gravitating boson star](#self-gravitating-boson-star) and [`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md) |

---

## Campaigns

The wrapper runs **three-stage** search campaigns (QD → CMA-ES → HQ) plus several
specialized single-stage campaigns. Each campaign has its own launcher under
`scripts/campaigns/`.

| Stage | What it does | Launcher | Output directory |
|-------|--------------|----------|------------------|
| **0 — QD** | MAP-Elites survey (8×8 archive) | `scripts/campaigns/qd/run.sh` | `runs/grtresna_qd/<name>/` |
| **1 — CMA-ES** | Local hill-climb from QD elites | `scripts/campaigns/cmaes/run.sh` | `runs/grtresna_cmaes/<name>/` |
| **2 — HQ** | Full-res replay + frames | `scripts/campaigns/hq/run_batch.sh` | `runs/grtresna_promote/<prefix>_hq_eval*/` |

**Geometry-first Stage-1 scout** (independent of the matter-first ladder above):

| Stage | What it does | Launcher | Output directory |
|-------|--------------|----------|------------------|
| **G0 — geometry atlas** | MAP-Elites over stationary metrics; frozen `f_geo` / `f_ff` | `scripts/campaigns/geometry_atlas/run.sh` | `runs/geometry_atlas/<name>/` |

This is a pure-geometry scout. It does **not** replace matter-first QD/CMA-ES; elites later hand off to `project_geometry_motif.py` for matter synthesis. Stationary atlas scores are screening metrics only — dynamical shortcuts still require GRTeclyn evolution.

Shared search defaults (grid, gates, 4D probe, objective) live in
`scripts/campaigns/lib/search_common.sh`. HQ defaults in
`scripts/campaigns/lib/promote_common.sh`. QD and CMA-ES **must** stay aligned
(same `search_common.sh`) so warm-started scores are comparable. HQ is
intentionally different: higher resolution and longer time stress-test whether
shortcuts survive refinement.

### Campaign catalog

| Campaign | Launcher | Objective | Descriptor | Matter | What it searches for |
|----------|----------|-----------|------------|--------|----------------------|
| **`geometry_atlas`** | `scripts/campaigns/geometry_atlas/run.sh` | stationary `f_geo`/`f_ff` | `f_geo` × log exotic energy | none (pure geometry) | Broad stationary metric atlas; Stage-1 inverse-design scout |
| **`general_ftl` (wormhole/ring/spin)** | `scripts/campaigns/general_ftl/run_all.sh` | `general_ftl` | `ftl_lifetime` | real scalar (pinned 15-D subspace) | FTL shortcut on a wormhole/ring/spin geometry; current production path |
| **`ftl_4d` (generic QD)** | `scripts/campaigns/qd/run.sh` | `ftl_first` | `ftl_lifetime` | real scalar shell/ring/free | Generic FTL discovery |
| **`qball_trajectory` (spiral)** | `scripts/campaigns/qball_trajectory/run.sh` | `general_ftl` | `ftl_lifetime` | complex scalar Q-ball, 5 per-lump orbits (39-D) | FTL from compact solitons on retrograde spiral orbits |
| **`qball_trajectory` (Lentz)** | `scripts/campaigns/qball_trajectory/run_lentz.sh` | `general_ftl` | `ftl_lifetime` | canonical Q-ball only, v_max=0.5c | Pure positive-energy FTL (no phantom matter) |
| **`qball_trajectory` (shear)** | `scripts/campaigns/qball_trajectory/run_shear.sh` | `spacetime_shear` | `spacetime_shear` | canonical Q-ball | Extreme non-collapsing frame-dragging shear |
| **`qball_trajectory` (f_geo)** | `scripts/campaigns/qball_trajectory/run_fgeo.sh` | `f_geo_max` | `ftl_lifetime` | complex scalar Q-ball, 39-D, phantom free | Evolving-geodesic shortcut × matter retention (`qball_traj_fgeo_v1`: 400 evals, best depth 38.3%) |
| **`qball_trajectory` (depth)** | `scripts/campaigns/qball_trajectory/run_fgeo_depth.sh` | `f_geo_depth` | `ftl_lifetime` | complex scalar Q-ball, 39-D, phantom free | Pure DEPTH hunt: raw uncapped evolving f_geo, no survival/stability shaping; stop_time 32, emission sweep to t=18 |
| **`gw_beam`** | `scripts/campaigns/gw_beam/run.sh` | `gw_beam` | `gw_beam` | canonical Q-ball trajectory | Directional gravitational-wave emission (Z-axis beaming) |
| **`splash` (critical collapse)** | `scripts/campaigns/splash/run.sh` | `critical_collapse` | `wave_focusing` | canonical bosonic shell | Gravitational-wave focusing / critical collapse |
| **`boson_star`** | `scripts/campaigns/boson_star/run.sh` | `ftl_first` | `ftl_lifetime` | complex scalar / U(1), 7-D | Single centered Gaussian boson star |
| **`boson_star` (FTL shell, RL chassis)** | `scripts/campaigns/boson_star/ftl_shell_run.sh` | `general_ftl` | `ftl_lifetime` | `BosonStarBH` superposed shell (~18-D, exotic wedge ON) | Bosonic shell FTL with RL pump handoff |

### Objective modes

| Mode | Goal | Dominant terms |
|------|------|----------------|
| `ftl_first` | FTL shortcut, health secondary | geodesic + operational FTL dominate; shaping gradients cut to ~40%; trapped-surface veto |
| `robust_ftl` | Persistent, low-exotic FTL | persistence boosted (300→500), exotic hardened (40→70), coordinate signals trimmed |
| `general_ftl` | Gauge-invariant shortcut + curvature | geodesic + curvature_activity; graded horizon penalty; warp-motor shaping disabled |
| `f_geo_max` | Evolving shortcut × matter retention | 10000×`ftl_geo_evolving` (= geodesic depth × structural persistence); no exotic penalty; pump tax + graded horizon stay on |
| `f_geo_depth` | **Raw shortcut depth, go-deep** | 10000×`ftl_geo_depth` (uncapped raw evolving f_geo, 1% = 100 pts, NOT survival-multiplied); no survival/stability/confinement/exotic terms; pump tax + graded horizon (−200 max) only |
| `critical_collapse` | GW focusing / collapse | geometric splash (χ-well + Ψ₄ wave + K-crunch) primary; ρ/focus/lapse secondary; FTL ignored |
| `gw_beam` | Directional GW power + Z-beaming | `(1000×quality + 100×peak + health) × gw_health_multiplier + penalties`; collapse → multiplier 0 |
| `spacetime_shear` | Max curvature, avoid horizon | 1000×curvature_activity + confinement; horizon veto (−500); FTL-agnostic |

**Depth scaling is uncapped (2026-08-05).** `_geo_magnitude` (1.0 at 20%
path saving) and the operational arrival-time components (`operational_ftl`,
`ftl_persistence`) used to saturate via `min(..., 1.0)`; the caps are
removed, so every mode now pays linearly for depth beyond target. (The old
cap silently turned `qball_traj_fgeo_v1` into a matter-retention contest
above 20% depth.) Artifact-guard compressions -- precursor / shift-drive log
scales, solved-FTL locality gates -- are deliberately still bounded.

---

## How to run

### Prerequisites

From the GRTeclyn repo root:

```bash
cd /path/to/GRTeclyn
uv sync   # Python deps including yt, h5py>=3.10 for the Chombo→gridinit bridge
```

Then configure site paths (required for GRTresna builds / solves) — see
[Site paths (`.env`)](#site-paths-env).

First build (single GPU, no MPI):

```bash
BUILD=1 bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

Binaries (GRTresna MPI solver + GRTeclyn GPU binary) — see [Operations](#operations).

### Site paths (`.env`)

Machine-specific paths stay out of git. Copy the template, edit values, and
either source `scripts/lib/env.sh` (shell campaigns) or rely on automatic load
from Python (`grteclyn_wrapper.core.site_paths`).

```bash
cd grteclyn-wrapper
cp .env.example .env
# edit .env — set at least SIM_ROOT and GRTRESNA_ENV
source scripts/lib/env.sh   # exports vars for this shell + campaign scripts
```

| Variable | Meaning | Example shape |
|----------|---------|---------------|
| `SIM_ROOT` | Parent of GRTeclyn / GRTresna / Chombo | `/path/to/simulation` |
| `GRTECLYN_ROOT` | This checkout | `${SIM_ROOT}/GRTeclyn` |
| `GRTRESNA_ROOT` | Sibling elliptic solver | `${SIM_ROOT}/GRTresna` |
| `CHOMBO_HOME` | Chombo `lib` dir | `${SIM_ROOT}/Chombo/lib` |
| `OPENMPI_ROOT` | Local OpenMPI prefix (multi-GPU) | `${SIM_ROOT}/local/openmpi` |
| `GRTRESNA_ENV` | Conda/venv with `mpirun` for GRTresna | `/path/to/envs/grtresna` |

Rules:

- **`.env` is gitignored** — never commit it. Commit only [`.env.example`](.env.example).
- Already-exported shell variables win over `.env` (safe to override per run).
- `${VAR}` expansion is supported inside `.env` (e.g. `GRTECLYN_ROOT=${SIM_ROOT}/GRTeclyn`).
- Shell scripts that `source scripts/lib/env.sh` pick up the same keys.
- Python resolves the same layout via `site_paths` (loads `.env` on first use).
- If `.env` is missing, `GRTECLYN_ROOT` is auto-detected from the wrapper layout;
  `GRTRESNA_ENV` is **not** guessed — set it in `.env` for GRTresna MPI work.

Quick check: `source scripts/lib/env.sh && echo $GRTECLYN_ROOT`, or from
Python `from grteclyn_wrapper.core.site_paths import grteclyn_root`.

### Stage 0 — MAP-Elites (QD)

**Generic launch** (background, 8 GPUs, pipelined GRTresna + GPU):

```bash
QD_NAME=my_campaign_v1 \
QD_TARGET_EVALS=200 \
OBJECTIVE_MODE=ftl_first \
GPU_IDS="0 1 2 3 4 5 6 7" \
GPU_SLOTS_PER_DEVICE=1 \
MAX_CONCURRENT_GRTRESNA=3 \
  nohup bash scripts/campaigns/qd/run.sh \
  > ../runs/my_campaign_v1.launch.log 2>&1 &
```

**`general_ftl` wormhole** (current production path — pins 15-D subspace, dirs `x y z`):

```bash
BRANCH=wormhole \
QD_NAME=general_ftl_wormhole_v22 \
QD_TARGET_EVALS=200 \
GPU_IDS="0 1 2 3 4 5 6 7" \
  nohup bash scripts/campaigns/general_ftl/run_all.sh \
  > ../runs/general_ftl_wormhole_v22.launch.log 2>&1 &
```

**Resume** an existing QD run (same `QD_NAME`, same dir):

```bash
QD_NAME=general_ftl_wormhole_v21 QD_RESUME=1 \
  bash scripts/campaigns/qd/run.sh
```

**Campaign-family one-liners** (defaults live inside each launcher; same
override env vars as `run.sh` -- `QD_NAME`, `QD_TARGET_EVALS`, `GPU_IDS`,
`QD_RESUME=1`, ...):

```bash
bash scripts/campaigns/qball_trajectory/run.sh            # spiral: 39-D per-lump orbits, general_ftl
bash scripts/campaigns/qball_trajectory/run_lentz.sh      # Lentz: canonical matter only, v_max=0.5c
bash scripts/campaigns/qball_trajectory/run_shear.sh      # spacetime shear
bash scripts/campaigns/qball_trajectory/run_fgeo.sh       # f_geo_max: shortcut x retention (qball_traj_fgeo_v1)
bash scripts/campaigns/qball_trajectory/run_fgeo_depth.sh # f_geo_depth: raw uncapped depth, seeds from fgeo_v1 elites
bash scripts/campaigns/gw_beam/run.sh                     # directional GW emission (GPU_IDS="0 1 2 3")
bash scripts/campaigns/splash/run.sh                      # critical collapse / GW focusing
GRTRESNA_MATTER_SECTOR=boson_star GRTRESNA_MATTER_COUPLING=canonical \
  bash scripts/campaigns/boson_star/run.sh                # boson star, 7-D

# Bosonic shell + FTL (RL chassis), ~18-D, exotic wedge ON:
QD_NAME=boson_shell_ftl_rl_v1 QD_TARGET_EVALS=200 QD_ITERATIONS=30 \
STOP_TIME=16.0 PLOT_INTERVAL=320 GRTECLYN_FRAMES=1 \
GPU_IDS="0 1 2 3 4 5 6 7" MAX_CONCURRENT_GRTRESNA=5 BATCH_SIZE=8 \
  nohup bash scripts/campaigns/boson_star/ftl_shell_run.sh \
  > ../runs/boson_shell_ftl_rl_v1.launch.log 2>&1 &
```

| Knob | Default (search) | Notes |
|------|------------------|-------|
| Grid | N=128, L=64, ml=1–2 | GRTresna solve on 128³ domain |
| Stop time | t=16 | `STOP_TIME=16.0` |
| Archive | 8×8 bins | `BINS=8` |
| Frames | **off** (search) / **on** (GW beam, boson shell) | `GRTECLYN_FRAMES=0` for speed |
| 4D geodesic | `search` profile | cheap stack for leaderboard |
| Prune eval dirs | top 3 + FTL peaks | `QD_KEEP_TOP_EVAL_DIRS=3` |

**Outputs:** `trajectory.jsonl`, `archive.json`, `eval_*/`, `ftl_champions.json`.

**Monitor:**

```bash
tail -f runs/grtresna_qd/<name>/trajectory.jsonl
cat runs/grtresna_qd/<name>/ftl_champions.json
```

### Stage 1 — CMA-ES optimization

Warm-start from a QD (or prior CMA-ES) trajectory. **Must match** the source
campaign's `OBJECTIVE_MODE`, grid, `STOP_TIME`, pins, and 4D profile.

```bash
RUN_NAME=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
WARM_START_TRAJECTORY="${GRTECLYN_ROOT}/runs/grtresna_qd/general_ftl_wormhole_v21/trajectory.jsonl" \
WARM_START_TOP_K=1 WARM_START_JITTER=0.05 SIGMA0=0.05 \
TARGET_EVALS=150 MAX_GENERATIONS=50 KEEP_TOP_EVAL_DIRS=3 \
GPU_IDS="0 1 2 3 4 5 6 7" MAX_CONCURRENT_GRTRESNA=3 \
PIN_DIMS="$(bash -c 'source scripts/campaigns/lib/general_ftl_pins.sh && ftl_general_ftl_wormhole_pins')" \
  nohup bash scripts/campaigns/cmaes/run.sh \
  > ../runs/general_ftl_wormhole_cmaes_v1.launch.log 2>&1 &
```

| Knob | Typical | Notes |
|------|---------|-------|
| Population | **≥ 4 × GPU slots** (e.g. 16 on 4 GPUs) | NEVER `= #GPU_IDS` — see warning below |
| σ₀ | 0.05–0.08 | local basin width |
| Warm-start | top-K elites | `WARM_START_TOP_K`, `WARM_START_JITTER` |
| Target | eval budget | `TARGET_EVALS` or `MAX_GENERATIONS × pop` |

> **⚠ Population must be several × the GPU slot count — never run CMA-ES
> "generationally starved".** CMA-ES has a hard barrier at every `tell()`:
> no next-generation candidate exists until ALL of the current generation
> finish. With `POPULATION = #GPUs` (the old default) any fast-failing
> candidate (post-load gate reject) or straggler leaves GPUs idle for most
> of each generation — measured ~50 % idle on the 2026-08-06
> `qball_traj_fgeo_depth_cmaes_v1` first launch. With `POPULATION ≥ 4 ×`
> slots the within-generation pipeline streams candidates back-to-back
> (like MAP-Elites) and the barrier tail is amortized away. Bonus: the
> CMA-ES textbook population for an n-dim search is `4+⌊3·ln n⌋` (≈ 15 at
> n = 39) — `pop = 4` was statistically undersized as well as slow.

**Outputs:** `runs/grtresna_cmaes/<RUN_NAME>/` — same layout as QD. **Monitor:**
`tail -f runs/grtresna_cmaes/<RUN_NAME>/trajectory.jsonl`

### Stage 2 — HQ promotion

Replays elite genomes at **N=256, L=128, ml=3, t=30** with fresh GRTresna solve,
**frames on**, 4D geodesic in **`hq`** verify mode, incremental
`score_timeseries.jsonl`.

**Single candidate** (`CANDIDATES` = eval/gpu **pairs**):

```bash
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
CANDIDATES="46 0" \
  bash scripts/campaigns/hq/run_batch.sh
```

**Top-K auto-pick** from trajectory (score-sorted `gpu_ok` rows):

```bash
SOURCE_RUN="${GRTECLYN_ROOT}/runs/grtresna_cmaes/general_ftl_wormhole_cmaes_v1" \
NAME_PREFIX=general_ftl_wormhole_cmaes_v1 \
OBJECTIVE_MODE=general_ftl \
TOP_K=1 \
  bash scripts/campaigns/hq/run_batch.sh
```

| Knob | Search (QD/CMA-ES) | HQ |
|------|-------------------|-----|
| Grid | 128³, L=64, ml=1–2 | **256³, L=128, ml=3** |
| Stop time | 16 | **30** |
| 4D mode | `search` | **`hq`** (full stack) |
| Frames | off | **on** (`GRTECLYN_FRAMES=1`) |
| GW / Ψ₄ | off (`--no-psi4`) | **on** (`GRTECLYN_PSI4=1`, HQ default) |
| Extraction radii | `4 8` | **`8 12 24`** from physics `center` |

**Monitor:**

```bash
tail -f runs/grtresna_promote/<name>.log
tail -f runs/grtresna_promote/<name>/small_data/score_timeseries.jsonl
ls runs/grtresna_promote/<name>/frames/
bash scripts/plot/make_movies.sh runs/grtresna_promote/<name> --framerate 10
```

### End-to-end orchestrator (QD → CMA-ES → HQ)

Single script under one campaign root (`runs/campaigns/<CAMPAIGN_NAME>/`):

```bash
CAMPAIGN_NAME=general_ftl_wormhole_v22 \
  bash scripts/campaigns/run_full_campaign.sh
```

Resume / partial: `RESUME=1 STAGE=cmaes CAMPAIGN_NAME=... bash scripts/campaigns/run_full_campaign.sh`

See [`scripts/campaigns/README.md`](scripts/campaigns/README.md) for stage
comparison and full env-var reference.

### Rules (do not skip)

1. **CMA-ES must mirror QD** — same `OBJECTIVE_MODE`, pins, grid, `STOP_TIME`, geodesic config.
   And `POPULATION ≥ 4 × GPU slots` — never `= #GPUs` (generation barrier starves the pipeline; see the Stage 1 warning).
2. **`general_ftl` needs `GRTECLYN_GEO_DIRECTIONS=x y z`** — wormhole shortcuts live on z; x-only scoring replays elites at the wrong fitness.
3. **HQ `CANDIDATES` is eval/gpu pairs** — e.g. `"46 0 39 1"` not a bare eval list.
4. **Search turns frames off, HQ turns them on** — by design (`search_common.sh` vs `promote_common.sh`).
5. **Promotion must use `N > L`** (or same `L` with larger `N`) to refine the grid. `L=N` only enlarges the domain at `dx=1` — no fidelity gain.
6. **Plotfiles go to node-local scratch, never NFS** — automatic since `core/scratch.py`; `output_path` stays on NFS. Only override it (`GRTECLYN_SCRATCH`) if `/tmp` is not node-local on your machine. See [Plotfile scratch MUST be node-local](#plotfile-scratch-must-be-node-local-required).
7. **The pump runs for the ENTIRE simulation** — never stopped mid-run. Enforced at launch, not by checklist; see [Pump convention](#pump-convention-enforced-at-launch) below.

### Pump convention (enforced at launch)

The PD pump stays on for the whole simulation. In the evolution config this
means the `rl_pump_stop_time` key is **absent** — the binary's default (`-1`)
is "never stop" — so a config with no pump key is already correct, and any
non-negative value someone bakes in silently changes the physics of every run
that inherits it. That failure mode is now closed by the launchers themselves:

**Search launchers** (`qball_trajectory/cmaes_run.sh`): `RL_PUMP_STOP_TIME`
must be set explicitly — there is no silent default any more; an unset value
refuses to launch (exit 2).

```bash
RL_PUMP_STOP_TIME=-1     # the convention: pump on for the whole run
GEODESIC_EMIT_MIN_TIME=4 # required whenever the pump value is negative
```

The floor pairing is also enforced: a negative pump value erases the scorer's
fallback emission floor (`metrics/score/ftl.py` skips negatives), so `-1`
without an explicit `GEODESIC_EMIT_MIN_TIME` refuses to launch — otherwise
`f_geo` silently changes meaning between runs.

**Promotion launchers** (`promote/lib/run_matrix.sh`): every launch and every
`--list` runs `promote/lib/validate_pump_convention.py` over the manifest and
the environment. It refuses (exit 3) when:

- the manifest bakes a non-negative `rl_pump_stop_time` anywhere
  (`physics_frozen`, `extra_overrides`, any nesting, `key=value` string form —
  unparseable values are refused, not waved through), or
- the environment carries `RL_PUMP_STOP_TIME ≥ 0`, or
- a pump-on manifest is launched without `GEODESIC_EMIT_MIN_TIME` in the env.

**Deliberate pump-off controls** (e.g. the pump-free twin that proves `f_geo`
persists without ignition) are the only exception, and they must say so in
the manifest, top-level:

```json
"pump_off_control": true
```

Only with that key may a manifest carry a non-negative `rl_pump_stop_time`.
The key is not inherited from templates — each control manifest declares it
itself, so a copy-paste of an old manifest with a baked stop time is refused
loudly instead of quietly stopping the pump.

Covered by `tests/scripts/test_pump_convention.py` (11 tests, including
end-to-end refusal of the retired template manifest and acceptance of the
pump-free twin).

### Stopping detached campaigns — kill the orchestrator first

Campaign launchers run detached (`setsid nohup bash launcher.sh ... &`) so they
survive shell teardown. Stopping one naively **does not work**, for three
reasons found the hard way (2026-08-05, `bondi_dipole_v1` post-mortem):

1. **The `$!` you captured at launch is the wrong PID.** It is the short-lived
   `setsid` parent; the real launcher is its forked child in a *new session*.
   Killing the recorded pid (or its process group) hits nothing.
2. **Killing workers by path pattern silently *advances* the campaign.** The
   orchestrator's argv is just `bash launcher.sh` and per-run drivers use
   `--runs-dir X --name Y` (never the `X/Y` path you grep for). So a pattern
   kill takes out the simulation + consumer, the orchestrator sees a
   "finished" step, and launches the next run — which looks exactly like the
   campaign refusing to die.
3. **GRTresna solvers detach into their own session/pgid** — even a correct
   group-kill of the launcher leaves a solver running to its timeout.

The fix is a **global tool** — one implementation for every campaign type
(QD/MAP-Elites, CMA-ES, HQ replays, Bondi matrices, one-off ladders):

```bash
# Preview what would be killed (touches nothing):
bash scripts/campaigns/stop_campaign.sh --dry-run <runs_dir | campaign_name>
# Stop for real, with verification and escalation:
bash scripts/campaigns/stop_campaign.sh <runs_dir | campaign_name> [...]
```

It works in the only order that does: (1) freeze the queue — kill the
orchestrators and their shell ancestors first, found via
`<runs_dir>/launcher.pid`, driver argv (`--name <campaign>`,
`--runs-dir <dir>`), and a parent-walk; (2) sweep every worker class by runs
dir and scratch path; (3) **verify with pgrep and escalate to SIGKILL** until
nothing survives. A stop without the verification step is a guess.

Launcher-side contract: every campaign launcher registers its true PID by
calling the shared helper (one line, after the runs dir is known):

```bash
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "${RUNS_DIR}"
```

The stop tool's process discovery covers unregistered launchers too, but the
pid file is the fast, unambiguous path — add the call to any new launcher.
There are deliberately no per-campaign stop scripts: one tool, one way to
stop, for every campaign.

### ALWAYS extract frames on the fly (required)

Every GPU evolution run — QD, CMA-ES, HQ promotion, replay — **MUST** stream
plotfiles through `consume_plotfiles` during the simulation. Do not let heavy
`data/plt*` HDF5 directories accumulate; extract PNG frames + `small_data/`
metrics in flight and delete processed plotfiles immediately.

| Requirement | How |
|-------------|-----|
| Sidecar consumer ON | `CONSUME_PLOTFILES=1` (shell) or `consume_plotfiles=True` (Python) |
| Delete heavy HDF5 | on by default with `--keep-last 3` (even if `GRTECLYN_FRAMES=0`); opt out with `GRTECLYN_KEEP_PLOTFILES=1` |
| PNG frames written | Set `GRTECLYN_FRAMES_FIELDS` → outputs in `eval_*/frames/` |
| Verify consumer alive | `ps aux \| grep consume_plotfiles` while GPU is busy |

#### Plotfile scratch MUST be node-local (required)

Plotfiles are write-once, read-once, delete-immediately transients -- never
write them to NFS. `amr.plot_file` / `amr.check_file` are independent of
`output_path`: route the heavy data to node-local NVMe and keep only `.dat`
diagnostics and `small_data/` on the shared filesystem. A 256-cubed ml=3
plotfile is ~3.2 GB every ~288 s per run; on NFS that capped concurrency at
2 runs (consumers stalled in `D` state). Node-local scratch drops NFS traffic
from ~130 MB/s to KB/s and lets every GPU on the node run concurrently
(extraction 15.7 s vs 288 s cadence -- 18x headroom; measured 2026-07-28).

```
output_path    = "<NFS>/runs/<campaign>/<run>"       # .dat + small_data
amr.plot_file  = "/tmp/<scratch>/<run>/RadialRecipePlt"
amr.check_file = "/tmp/<scratch>/<run>/RadialRecipeChk"
```

Consumer: `--data /tmp/<scratch>/<run>` (local) `--out <NFS>/.../small_data`.

Enforced in code for every stage: `core/scratch.py` maps episode dir to
`/tmp/grteclyn_scratch/<campaign>_<eval>_<hash>/`, applied in
`core/params.py::episode_path_overrides`; the consumer's `--data`, the
scoring plotfile lookup and every cleanup path use the same mapping. An
unwritable root falls back to the episode dir with a warning; launchers that
set `amr.plot_file` themselves keep their explicit value.

| Variable | Effect |
|----------|--------|
| *(unset)* | scratch at `/tmp/grteclyn_scratch` -- the default |
| `GRTECLYN_SCRATCH=/path` | move the scratch root |
| `GRTECLYN_SCRATCH=0` | keep plotfiles in the episode directory (old behaviour) |
| `GRTECLYN_SCRATCH_FORCE_PURGE=1` | purge scratch even when the ledger is incomplete |

Rules that keep this safe:

* **Deletion is ledger-gated.** Only plotfiles recorded in
  `small_data/consume_state.json` are GC'd; a file that fails extraction is
  retried, never collected. `[ok] ... kept` means *protected right now*, not
  retained permanently. Consequence: **any extraction not enabled during the
  run is unrecoverable** -- decide `--areal-radius`, `--shell-fields`, frames
  and scoring passes before launch.
* **NEVER run external plotfile-deletion loops** alongside the consumer. An
  ad-hoc `while true; do ... rm -rf ...; done` loop once raced the consumer
  and destroyed complete, unprocessed plotfiles (2026-07). The consumer's own
  `--delete --keep-last N` is the only safe deletion path; if disk pressure
  is a concern, raise `plot_interval` or reduce extraction cost (`-j 1`,
  disable frames).
* **Purge scratch only when every resident plotfile appears in
  `consume_state.json`**, and allow a **600 s drain window** after the sim
  exits -- shorter windows truncated confinement data.
* **Cloning a baseline `params.txt`:** `amr.plot_file` / `amr.check_file` are
  absolute paths into the source run -- strip and re-emit them, or the clone
  writes into (and a `--delete` consumer prunes) the baseline. Verify each of
  the three path keys occurs exactly once.
* **Sizing.** Local steady state = `n_runs x (keep_last + jobs) x
  plotfile_size` (six 256-cubed runs = ~70 GB); check `df -h /tmp` first;
  checkpoints are **not** pruned by the consumer. NFS cost is set by
  `GRTECLYN_METRIC_STACK_N_SPACE` (scales as `n_space` cubed): default 33 is
  ~30 MB per t=30 run, HQ's 257 is **~4.2-6.2 GB per run**. The metric stack
  is what post-hoc rescoring reads after plotfiles are purged -- delete it
  only when no further scoring pass (e.g. a later-`t_emit` emission sweep) is
  wanted.
* **Keep every write inside your own space.** Call
  `grteclyn-wrapper/.venv/bin/python -m ...` directly (not `uv run`, which
  writes `~/.cache/uv`) and pin `XDG_CACHE_HOME`, `UV_CACHE_DIR`,
  `MPLCONFIGDIR`, `TMPDIR`, `PYTHONPYCACHEPREFIX` into `$SCRATCH/_cache/`.
  On the shared cluster the complete write set is `/tmp/<scratch>/`
  (transient) and the NFS run directory -- nothing else.

Verify pruning:

```bash
ps aux | grep consume_plotfiles | grep -- '--delete --keep-last 3'
grep -h '\[gc\] deleted' ../runs/grtresna_promote/e118_*/run.log | tail
```

A live run holds at most `keep_last + jobs` plotfiles. If the count keeps
growing, consumers are blocked in NFS I/O (`D` state) -- raise
`plot_interval` or cut extraction cost. Never delete an unprocessed plotfile
to force the count down; that loses the corresponding Psi4/FTL sample.

### Smoke test (all stages)

```bash
bash scripts/campaigns/smoke_test.sh
```

Covers: pytest preflight, QD dry-run, CMA-ES dry-run + warm-start, HQ batch
dry-run + `TOP_K` auto-pick, `replay_eval.py --help`.

---

## Integration with GRTresna and GRTeclyn

Three repos cooperate. The wrapper is the glue; the two C++ trees do the physics.

```
GRTresna (sibling repo, ../GRTresna)        GRTeclyn (this repo, .)
  elliptic constraint solver                  GPU time evolution
  ┌────────────────────┐                     ┌─────────────────────┐
  │ Ham + Mom solve     │  params.txt         │ CCZ4 + matter RHS    │
  │ (MPI, Chombo/AMR)   │ ◄───── wrapper ───► │ (CUDA, AMReX)        │
  │ InitialDataFinal    │   .gridinit         │ ExternalGridInitial  │
  │   .3d.hdf5          │ ──────wrapper────►  │   Data               │
  └────────────────────┘  (+ matter.json)     └─────────────────────┘
                                                        │
                                                  plotfiles (HDF5)
                                                        │
                                                  consume_plotfiles
                                                        ▼
                                             small_data/ + frames/ + score
```

### Per-eval loop (every CMA-ES member and QD candidate)

1. **Sample ansatz parameters → GRTresna MPI solve** (Ham + Mom). Exotic
   (`rho<0`) candidates auto-switch to the K=0 maximal-slicing solver
   (`apply_exotic_safe_solver`).
2. **Reject** if convergence missing, NaN, or above threshold.
3. **Solved-geometry FTL gate** on `.gridinit` (cheap, pre-GPU, ~1 s).
4. **Post-load constraint gate** — short GPU launch (`stop_time=0.01`, **no**
   `consume_plotfiles`) that loads the `.gridinit` and rejects if `L2_Ham`/
   `L2_Mom` exceed thresholds (default `1e-2`). Writes to `postload_gate/`
   only — **no frames here**.
5. **Main** GRTeclyn GPU evolution (`stop_time=16` search / `30–50` HQ) →
   `consume_plotfiles` sidecar → PNG frames in `frames/` → objective score.
6. **Feed fitness** to CMA-ES or MAP-Elites archive cell.

Pre-evolution gates 2–4 are centralized in
`search/grtresna_evaluation_gates.py` (`apply_grtresna_pre_evolution_gates`),
shared by CMA-ES and MAP-Elites.

### Matter–geometry consistency

GRTresna solves the constraints with **independent, signed scalar lumps** (each
lump has its own sign and there are no inter-lump cross terms). The wrapper keeps
GRTeclyn's evolved matter in lock-step end-to-end:

```text
GRTresna signed-lump solve → .gridinit (+ per-lump phi_k/Pi_k channels, matter metadata)
  → GRTeclyn GRTresnaIndependentScalars matter (mass-matched V = ½ m² φ²)
  → post-load gate: L2_Ham / L2_Mom small  ⇒ valid Einstein solution
```

A `gpu_ok` candidate is a genuine, constraint-satisfying solution of GR; whether
it is an *interesting* (transport-capable) one is what the search explores.

**Implementation:** C++ matter `Source/Matter/GRTresnaIndependentScalars.*`,
RadialRecipe dispatch in `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`;
Python bridge `grtresna/lump_fields.py`, `grtresna/matter_wiring.py`,
`grtresna/solver.py`; gate `projection/postload_gate.py`.

### Matter-profile contract — the rail for adding new profiles safely

Wiring a matter profile into the pipeline means the **same physical facts are
repeated across three codes** (Python solver, GRTresna C++ painter, GRTeclyn
evolution): the potential `V = ½m²|Φ|² − ¼λ|Φ|⁴ + ⅙μ₆|Φ|⁶`, the couplings, the
frequency `ω`, the U(1) momentum rule `Π₂ = −(ω/α)φ₁`, the winding phase, the
exotic sign, the profile table format, and the coordinate center. If **any one**
disagrees you don't get an error — you get a *plausible-looking wrong answer*
(silent dispersal, or a t=0 constraint blow-up). Real bugs of this class:
a stale C++ binary painting a torus as `f≡0` (→ "Ham=0%", looks converged), an
initial-data/evolution sign mismatch, a diagnostics column off-by-one.

The **contract** (`grtresna/matter/profile_contract.py`) makes these invariants
executable so they fail in seconds instead of after a wasted GPU run:

- **`MatterProfileSpec`** — the ONE place a profile's physics lives. It emits the
  GRTresna `lump<k>_*` params **and** the Python reference field, so they cannot
  drift.
- **`reference_paint(spec, x, y, z)`** — the Python mirror of the C++ painter
  (`Φ = f(ρ,z)e^{imφ}`, `Π₁=+(ω/α)φ₂`, `Π₂=−(ω/α)φ₁`).
- **`check_gridinit_matches_spec(gridinit, spec)`** — reads the GRTresna-produced
  `.gridinit`, paints the reference at the same cell centres, asserts phi/phi2/
  Pi/Pi2 match (relative tol, default 5%). This is the **t=0 cross-code
  consistency test**.

**To add a NEW matter profile (checklist):**

1. **C++ painter** — add the profile to `../GRTresna/Source/Matter/BosonStarParams.hpp`
   (a loader + a `lump_winding_modulus` / `lump_phi1` branch), give it a `profile`
   enum int, and extend the `read_params` inherit so it picks up its table path.
   **Rebuild GRTresna** (`scripts/wormhole/build/build_grtresna_bosonstar.sh`) —
   forgetting this is the "stale binary" bug the rail catches.
2. **Python mirror** — add the int to `PROFILE_KIND_TO_INT` and one modulus branch
   to `reference_paint` in `profile_contract.py`.
3. **Contract test** — add a case to `tests/grtresna/test_profile_contract.py`
   (a self-contained synthetic-gridinit round-trip needs no GPU/solve).
4. **Wire the gate** into whatever driver emits the ID (see `solve_torus.py`,
   which runs `check_gridinit_matches_spec` after every solve and refuses to hand
   off a mismatched ID).

**Run the rail:**

```bash
cd grteclyn-wrapper
uv run pytest tests/grtresna/test_profile_contract.py -q   # 7 tests: passes good,
#   catches empty-field (stale binary), momentum sign flip, wrong ω, + real torus
```

The gate also runs automatically inside `solve_torus.py` (prints
`MatterProfile consistency (PASS, …)`), so a bad solve never reaches evolution.

### Rotating Q-torus wormhole — genuine stationary eigenstate support + collapse

The wormhole support is a genuine 2D spinning Q-ball eigenstate
`Φ = f(ρ,z)e^{i(mφ−ωt)}`, solved by `grtresna/profiles/qball_torus.py`
(bordered Newton + amplitude continuation) and painted into GRTresna as
`profile == 4`. (The old `(sinθ)^m`-twisted sphere was not stationary and
drained its Noether charge, half-life t~13-16.) Journal:
[`OrbitalPumpPlan.md`](../research/rotatingwormhole/OrbitalPumpPlan.md) Phase 8.

**Solve an isolated torus ID** (throat-free flat background; validates the
eigenstate on its own). `EXOTIC=1` matches the phantom evolution matter;
`TORUS_OMEGA` is the target internal frequency (use the compact/thick-wall band
≈0.2–0.45 — deep thin-wall ω<~0.15 does not fit the box, and the solver errors
with a box-escape message):

```bash
cd /path/to/GRTeclyn
EXOTIC=1 TORUS_OMEGA=0.25 bash grteclyn-wrapper/scripts/wormhole/id/solve_torus.sh 2
# -> runs/rotating_torus_id/torus_m1_om0p250_.../{initial_data.gridinit, qball_torus.dat}
#    prints Ham%/Mom% and the MatterProfile consistency PASS gate
```

**Solve a torus-supported wormhole ID** (adds a bare-mass throat the torus hugs;
`THROAT_MASS` = `bh1_bare_mass`, the positive mass that implodes when support is
cut — the Phase-7 collapse recipe wants b0≈2, i.e. `THROAT_MASS≈1`):

```bash
EXOTIC=1 THROAT_MASS=1.0 TORUS_OMEGA=0.25 \
  bash grteclyn-wrapper/scripts/wormhole/id/solve_torus.sh 2
```

**Evolve** with `wormhole_case.sh` (`--gridinit` override + `--full-box` for a
centered z-symmetric object). Two arms — control (support held) vs collapse
(exotic support ramped to 0, the rotating analogue of the Ellis–Bronnikov
`S_support` cut). Default is one GPU (`--gpu N`); for Arena-OOM cases use
multi-GPU `--np` / `--gpus` — see
[Launch multi-GPU evolutions](#launch-multi-gpu-evolutions-wormhole_case---np----gpus).

```bash
G="$PWD/runs/rotating_torus_id/torus_m1_om0p250_kappa1p00_dx0p5_L64_lam170_mu614450_exotic_throat1/initial_data.gridinit"

# CONTROL — support held
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --gridinit "$G" --full-box \
  --omega 0.25067 --m 1 --dx 0.5 --box-size 64 --max-level 0 --stop-time 20 \
  --mass 0.5 --lambda 170 --mu6 14450 --no-frames --run-suffix torus_wh_control --gpu 0

# COLLAPSE — cut the exotic support over t=8..10 (floor 0)
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --gridinit "$G" --full-box \
  --omega 0.25067 --m 1 --dx 0.5 --box-size 64 --max-level 0 --stop-time 20 \
  --mass 0.5 --lambda 170 --mu6 14450 \
  --support-ramp-t-start 8 --support-ramp-t-end 10 --support-ramp-floor 0 \
  --no-frames --run-suffix torus_wh_collapse --gpu 1
```

Diagnostics land in `runs/rotating_wormhole/evo_*/output/data/collapse_diagnostics.dat`
(1-based: **3** `min_chi`, **8** `max_ah_r`, **9** `min_theta_plus`, **21**
`Q_sphere`, 23 `pump_work`). Collapse ⇒ `min_chi → puncture` + a persistent
`max_ah_r>0`; a clean trapped-surface result needs `--max-level ≥ 4` (AMR; note
the gridinit trilinear-interpolation kink caveat in the wormhole `Debug.md`).

Couplings/ω/exotic/mass **must** match between the `solve_torus` ID and the
`wormhole_case` evolution — the consistency gate enforces the ID side; keep
`--mass/--lambda/--mu6/--omega` equal to the solve's `MASS/LAMBDA/MU6/TORUS_OMEGA`.

#### GW / Ψ₄ extraction is now C++-side (upgrade)

Since 2026-07, Ψ₄ is extracted in-code by GRTeclyn's `WeylExtraction`
(`SphericalExtraction`): `Weyl4_mode_2{0,1,2}.dat` in `output/data/`, one row
per coarse step (dt~0.01), one Re/Im column pair per extraction radius --
dense, multi-radius, independent of `--plot-interval` (which now controls
only frame/movie cadence). The Python `process_wave` sidecar still runs
(coarse cross-check + confinement + frames). For physics plots use
`small_data/psi4_mode_l2m0.dat` -- see the trust note under
[Visualization](#visualization).

```bash
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --gridinit "$G" --full-box \
  --omega 0.25067 --m 1 --dx 0.5 --box-size 64 --max-level 3 --stop-time 40 \
  --mass 0.5 --lambda 170 --mu6 14450 --sponge \
  --extraction-radii 12 18 24 --plot-interval 100 \
  --support-ramp-t-start 8 --support-ramp-t-end 10 --support-ramp-floor 0 \
  --run-suffix torus_wh_collapse_ml3_gw --gpu 0
```

| Flag | Effect |
|------|--------|
| `--extraction-radii R [R ...]` | Weyl4 shell radii (default `12 24`); shells at/behind the sponge (>= L/2-2) auto-dropped |
| `--write-extraction-surfaces` | also dump raw per-step surfaces (off by default -- thousands of files) |
| `--plot-interval N` | frame/movie cadence only (Ψ₄ independent) |

Implementation note: the wormhole `Main` builds a `BHAMR<1>` purely to reuse
its `m_weyl_interpolator` (puncture tracking stays disabled) -- requires the
MPI+CUDA rebuild below.

### One-off GRTresna solve

```python
from pathlib import Path
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve

cfg = GRTresnaConfig(
    mpi_ranks=1, bh1_bare_mass=0.0,
    lump_amp=0.1, lump_width=8.0,
    lump_velocity=(0.2, 0.0, 0.0),
)
gridinit = solve(cfg, work_dir=Path("/tmp/grtresna_run"))
# GRTeclyn: recipe_initial_data_file = <gridinit>
```

Deep bridge docs: [`src/grteclyn_wrapper/grtresna/README.md`](src/grteclyn_wrapper/grtresna/README.md).

### Self-gravitating boson star

ODE solver for genuine self-gravitating boson stars (gravity provides the
binding, not an artificial pump well): `grtresna/profiles/boson_star_ode.py`.
After a four-bug fix a single stable-branch star holds confinement ~0.90 at
t=16 (was 0.58, dispersed); `max_level=3` still develops NaNs at t~6-9 -- an
open numerical-stability issue, not seed/pump physics. Full handoff:
[`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md).

---

## Metrics & scoring

Each episode is measured into an `EpisodeMetrics` dataclass, then mapped to a
scalar `Score` by `score_episode()`. Deep docs:
[`src/grteclyn_wrapper/metrics/README.md`](src/grteclyn_wrapper/metrics/README.md)
(metric groups, FTL stack precedence) and
[`src/grteclyn_wrapper/metrics/score/README.md`](src/grteclyn_wrapper/metrics/score/README.md)
(scoring pipeline).

### Metric groups (`EpisodeMetrics`)

Every group is parsed from a C++ `.dat` diagnostic file or computed from
plotfiles / recipe overrides / `.gridinit`. The `diagnostics/` package holds the
`.dat` parsers; `probes/` holds the computed metrics.

| Group | Source | Key score components |
|-------|--------|----------------------|
| collapse | `collapse_diagnostics.dat` | `numerical_survival` (completion gate), `lapse_health`, `horizon_penalty` |
| constraints | `constraint_norms.dat` | `constraint_health`, `initial_constraint_quality`, `structural_persistence` → survival |
| stability | collapse + `areal_radius.dat` | `stability`, `instability_penalty` |
| growth | derived time series | `constraint_growth` |
| comoving | shell profiles + shift | `comoving_stability` |
| energy_conditions | `energy_conditions.dat` | `energy_condition` |
| curvature | `curvature_invariants.dat` | `curvature_activity`, `nontrivial_geometry` |
| transport | barycenter in collapse dat | `transport_objective` |
| qei | physical + geodesic trajectory | `qei_penalty` |
| confinement | `diagnostics/confinement.py` | `confined_frac`, `rms_radius` (dispersal detectors) |
| psi4 | `small_data/psi4_mode_l2m0.dat` | GW power, beam_ratio, `gw_health_multiplier` |
| ftl (analytic) | recipe overrides t=0 | `ftl_shortcut`, `expansion_asymmetry` |
| general_ftl (t=0) | recipe 2D slice | `ftl_precursor`, `channel_progress` |
| general_ftl_evolved | latest plotfile | **`operational_ftl`** (primary coordinate FTL) |
| ftl_persistence | last N plotfiles | `ftl_persistence` |
| general_ftl_solved | `initial_data.gridinit` | `operational_ftl_solved` |
| geodesic_ftl | null-ray on frozen slice | `operational_ftl_geodesic` (frozen) |
| evolving_geodesic | null-ray on ≥3 plotfile stack (4D metric) | **`ftl_geo_evolving`** (headline gauge-invariant FTL) |
| ftl_timeseries | `small_data/ftl_timeseries.dat` | `ftl_geo_timeavg`, `ftl_geo_peak`, `ftl_lifetime` |
| physical | recipe t=0 proxies | `anec_condition`, `tidal_comfort` |
| boundary_flux | `boundary_flux.dat` / plotfile | `boundary_penalty` |
| effective_ec | ≥3 plotfile stack (`warpfactory.py`) | `exotic_penalty` |

Use `from grteclyn_wrapper.metrics.catalog import list_metrics, get_metric` for
programmatic lookup.

### FTL stack precedence

Scoring rewards **evolved** operational FTL (`general_ftl_evolved`) over t=0
slices. A t=0 shortcut that relaxes away scores nothing on `operational_ftl`.
Geodesic confirmation (`geodesic_ftl` / `evolving_geodesic`) gates the
highest-weight component.

| Probe | Metric field | When | Primary score use |
|-------|--------------|------|-------------------|
| `geodesic_ftl` | single Cauchy slice (∂g/∂t = 0) | per-frame + final | diagnostic when 4D ran; else `operational_ftl_geodesic` (frozen timeavg) |
| `evolving_geodesic` | time-interpolated stack from ≥3 plotfiles | end-of-run when enabled | **`ftl_geo_evolving`** (headline in QD; frozen credit zeroed) |

The evolving probe traces null rays through `EvolvingMetricField` (temporal
linear interpolation between plotfile slices at the ray's coordinate time). It
answers whether a pulse emitted at `t_emit` beats flat-space transit through the
**full evolved geometry**, not a frozen mid-run snapshot. A continuous-emission
sweep fires ray fans at multiple launch times (`t_emit = 0, 2, 4, …`) and
reports the peak `f_geo(t_emit)`.

**Enable:** QD `GRTECLYN_EVOLVING_GEODESIC=1` +
`GRTECLYN_EVOLVING_GEODESIC_MODE=search`; HQ `--evolving-geodesic` or
`GRTECLYN_EVOLVING_GEODESIC_MODE=hq`.

### Scoring `f_geo_evol` for hand-rolled campaigns — the corrected recipe

A campaign queue never reaches the metrics aggregation layer, so
`ftl_timeseries.dat` cols 13/14 stay at their `0.0  0` placeholders and no
`evolving_geodesic.json` exists until the post-hoc pass runs (Debug.md §13).
This is the correct way to run it after the fixes of 2026-07-28 (`4f31f33a`):

```bash
# ONE process per run, all in parallel — the pass is single-core.
# Export the env FIRST, on its own line. Do NOT chain `export … && nohup … &`:
# the trailing `&` backgrounds the whole && list, the exports land in that
# subshell only, and every later launch silently runs in SEARCH mode
# (3 rays, stride-2 slices, 15k steps) while still writing result files.
S=/tmp/grteclyn_scratch/_cache
export GRTECLYN_EVOLVING_GEODESIC_MODE=hq GEODESIC_EMIT_MIN_TIME=0 \
       XDG_CACHE_HOME=$S UV_CACHE_DIR=$S/uv MPLCONFIGDIR=$S/mpl \
       TMPDIR=$S/tmp PYTHONPYCACHEPREFIX=$S/pyc
for k in 0 4 8 16 24 30; do
  nohup grteclyn-wrapper/.venv/bin/python -u \
    grteclyn-wrapper/scripts/campaigns/rl/score_evolving_geodesic.py \
    runs/pump/pump_ladder_m0/lad_m0_tp$k --ftl-l 8 \
    > $S/score_tp$k.log 2>&1 &
done
```

**Verify before trusting output:** the first log line must say `mode=hq` and
the result line must say `rays=5/5` (search mode reports `rays=3/3`). The two
profiles produce different numbers from the same cache with no other visible
difference.

Rules baked into the script (do not work around them):

* **It refuses unfaithful caches.** `cache_fidelity` compares each slice's
  representable `min_chi` against the run's own `collapse_diagnostics.dat`;
  any slice >1.5× too shallow aborts the run's score (Debug.md §15). For a run
  whose *late* slices fail (deep-collapse endgame after the rays have already
  arrived), pass `--max-time <t>` to truncate the stack before the offending
  slices instead of `--force`. Truncation is conservative by construction:
  the strict no-frozen-tail guard fails any ray still in flight past the last
  kept slice rather than letting it coast through frozen geometry.
* **It does not recompute frozen per-slice `f_geo`.** The consumer already
  measured it per plotfile at full AMR fidelity into `ftl_timeseries.dat`
  col 3; the scorer reuses the peak over `geo_trustworthy` rows and records
  provenance in the notes (the old per-slice rebuild cost >60 min/run at
  90-110 GB RSS; the corrected pass is ~3 min/run at ~25 GB).
* **Results are written before they are printed**, so a dead parent shell
  (broken stdout pipe) can no longer discard a finished trace.
* A single launch at `t_emit=0` arrives at t≈12–13 and therefore **cannot
  distinguish rungs that only differ after t=16** — `tp16/tp24/tp30` report
  identical `f_geo_evol` to the last digit because the ray never samples any
  spacetime where they differ. That is the emission protocol, not a bug; use
  the emission sweep (`GRTECLYN_GEO_EMIT_INTERVAL`, `GEODESIC_EMIT_MIN_TIME`)
  to probe late launches.

### Scoring pipeline

`score_episode()` runs phases in fixed order (each mutates
`ScoringContext.components` / `.notes`):

1. `survival.compute_survival_components`
2. `health.compute_health_components`
3. `ftl.compute_ftl_components` — includes stationary, persistence, and
   geodesic-contradiction gates
4. `penalties.compute_penalty_components`
5. `gating.compute_nontriality_gate` (flat-space attractor guard)
6. `objectives.compute_total` → `Score`

The total is then assembled by the chosen [objective mode](#objective-modes).

### Score components and survival

Exact per-component weights live in `metrics/score/objectives.py` (validated
FTL terms 1000x/400x/300x dominate; shaping gradients ~20-100x; health terms
gated by `nontriviality_gate`; penalties: exotic 40x graded 0..-1.6,
stationary 8x, horizon 500x). `survival = numerical_survival x
structural_persistence`, where structural persistence = density retention x
morphological coherence x confined mass fraction -- reaching `stop_time`
alone earns nothing (empty space marches to the end trivially), and a
dissipating or fragmenting structure cannot bank shaping credit.
`SUPERLUMINAL_MARGIN = 0.05` de-saturates the superluminal-fraction
descriptor (only cells genuinely past c=1.05 count).

### Public API

```python
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode, EpisodeMetrics
from grteclyn_wrapper.metrics.score import Score, DEFAULT_WEIGHTS
```

---

## Main results

Run-by-run results live outside this README:

| Where | What |
|-------|------|
| [`MapElites.md`](../research/neuralspacetime/MapElites.md) | FTL search lab journal (SH, shell, wormhole) |
| [`MapElitesDynamics.md`](../research/neuralspacetime/MapElitesDynamics.md) | Trajectory FTL campaign lab journal |
| [`grlab/LabJournal.md`](../research/grlab/LabJournal.md) | GW beam + splash lab journal |
| [`results/`](../results/) | Git-friendly campaign extracts, e.g. [`qball-trajectory-evolving-geodesic-shortcut-search/CAMPAIGN_RESULTS.md`](../results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/CAMPAIGN_RESULTS.md) |
| [`NextSteps.md`](NextSteps.md) | Critical review of the claims' validity + hardening plan |

Three takeaways that shape current work: (1) genuine gauge-invariant FTL
shortcuts exist with exotic matter but are transient -- no configuration yet
holds a shortcut open while keeping matter confined; (2) the validation
pipeline (4D evolving probe, trust flags, HQ resolution ladder) is what
separates physical shortcuts from gauge artifacts -- several would-be headline
results were artifacts; (3) ansatz and matter sector dominate outcomes,
optimizer tuning is secondary.

---

## Operations

### Build GRTresna

Production searches use MPI `mpicxx.gfortran` (`RANKS=8` default). Needs
`CHOMBO_HOME` and `GRTRESNA_ENV` on `PATH` — configure them in
[`.env`](#site-paths-env) first, then:

```bash
cd grteclyn-wrapper && source scripts/lib/env.sh
cd "${GRTRESNA_ROOT}/Examples/ScalarFieldBH"
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

Produces `Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex`.
Serial debug: build with `MPI=FALSE`, run `RANKS=1`.

### Build GRTeclyn (GPU evolution binary)

RadialRecipe supports **both** single-GPU and multi-GPU (MPI+CUDA). Build
**without the conda/grtresna env g++** (nvcc requires gcc ≤ 12).

#### Single-GPU RadialRecipe (default search / HQ)

```bash
cd Examples/RadialRecipe
PATH="/usr/local/cuda/bin:$PATH" NO_MPI_CHECKING=TRUE \
  make USE_MPI=FALSE USE_CUDA=TRUE -j$(nproc)
```

Produces `main3d.gnu.CUDA.ex`. Use for N≲256 / L≲128 HQ promotions that fit one H100.

| Setting | Value | Why |
|---------|-------|-----|
| `USE_MPI` | `FALSE` | Single-GPU runs |
| `USE_CUDA` | `TRUE` | GPU kernels (AMReX + CCZ4 RHS) |
| **Do NOT** put grtresna env on PATH | — | Its `g++ 15.x` breaks nvcc's gcc ≤ 12 check |

#### Multi-GPU RadialRecipe (MPI + CUDA)

For Arena-OOM cases (large base grids, e.g. N=320/L=160 or N=384/L=128 with
`max_level=3`), build the MPI binary the same way as RotatingWormholeCollapse:

```bash
cd Examples/RadialRecipe

export PATH="/usr/local/cuda/bin:$OPENMPI_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$OPENMPI_ROOT/lib:${LD_LIBRARY_PATH:-}"

# USE_PARTICLES=TRUE is required so Weyl/ParticleInterpolator headers resolve.
make -j8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

Produces `main3d.gnu.MPI.CUDA.ex` (~137 MB). HQ replay / validation launchers
drive it with explicit GPU binding:

```bash
# 2 ranks on physical GPUs 4 and 5, reusing an existing gridinit
GRTECLYN_FRAMES=0 GRTECLYN_PSI4=1 \
uv run python grteclyn-wrapper/scripts/campaigns/hq/replay_eval.py \
  runs/grtresna_qd/qball_traj_spiral_v2/eval_000118 \
  --name e118_dl_L160_N320_t30_mpi2_hq_eval000118 \
  --gpu 4,5 --evolution-mpi-ranks 2 \
  --n-full 320 --l-full 160 --stop-time 30 \
  --evolving-geodesic --objective-mode general_ftl \
  --gridinit runs/grtresna_promote/e118_dl_L160_N320_t30_hq_eval000118/initial_data.gridinit

```

| Flag / env | Effect |
|------------|--------|
| `--evolution-mpi-ranks N` / `EVOLUTION_MPI_RANKS` | MPI ranks for GPU evolution (must equal `#` of GPU ids) |
| `--gpu a,b,...` / `GPU_ID` | Explicit physical GPU ids; each rank binds one via `GRTECLYN_GPU_IDS` |

**Binding rule (same as wormhole):** do **not** rely on parent
`CUDA_VISIBLE_DEVICES=0,1,…` + `LOCAL_RANK`. OpenMPI/prterun often drops that
env. The runner wraps each rank with
`CUDA_VISIBLE_DEVICES=${GRTECLYN_GPU_IDS[LOCAL_RANK]}`.

#### MPI status and triage runbook

**Do not trust any "MPI is broken" claim without re-running the checks below.**
That status has flipped twice. It is recorded per date because it is a property
of whichever node the pod currently sits on, not of this repo.

| Layer | Status | Last verified |
|---|---|---|
| `mpirun` itself (local OpenMPI) | **works** — 1 and 4 ranks, correct rank ids | 2026-08-19 |
| GRTeclyn RadialRecipe MPI+CUDA | **works** — 2 ranks, AMR max_level 3, clean past the old crash point | 2026-08-19 |
| GRTresna solver multi-rank | **works** — 8 ranks reproduce the serial residuals digit-for-digit | 2026-08-19 |
| GRTeclyn RotatingWormholeCollapse MPI+CUDA | worked multi-GPU, but only on an **older node** | 2026-06 |

**The July RadialRecipe AMR crash does not reproduce (retested 2026-08-19).**
A 2-rank run on `N=240, L=128, max_level 3` advanced cleanly through all four
levels with zero errors. Memory genuinely splits — **26 GB per card against
49.8 GB on one** — which is the point: multi-GPU buys *headroom*, not speed.
Throughput was 30.7 code units/h on two cards versus 29.4 on one, i.e.
unchanged. So use it to fit a grid that will not fit on one card, and do not
expect a run to finish sooner. Launch with
`--gpu 0,1 --evolution-mpi-ranks 2` (or `GPU_ID="0,1" EVOLUTION_MPI_RANKS=2`
through a promote campaign, which expands a bare `GPU_ID` into consecutive
ids automatically).

**GRTresna multi-rank solves are ~6× faster at no cost in accuracy.** At
`N=256`, 8 ranks ran ~74 s per nonlinear iteration against ~7.7 min
single-rank — a 6-iteration production solve drops from ~46 min to ~8 min —
and iterations 1–4 reproduced the single-rank Ham/Mom residuals to all seven
printed digits. This supersedes the older observation that 8 ranks converged
*worse* than 1 (0.93 % vs 0.63 % Ham); that predates the Chombo rebuild.
Pass `--grtresna-ranks 8`. Load stayed near 11 of 128 cores, so wider is
available if it ever pays.

**Two distinct failures have been mistaken for each other. Keep them apart:**

1. *Node-level.* On one node every MPI job died in PRRTE daemon start-up —
   even `mpirun -np 1 hostname`. Nothing in this repo could fix it; it went
   away when the pod moved. If this is happening, stop and check the node, do
   not rebuild anything.
2. *Toolchain-level.* GRTresna died with SIGILL from mismatched
   `-march=native` objects. Fixed by rebuilding Chombo's MPI libs
   consistently: `scripts/build/rebuild_grtresna_mpi.sh` (it ends with its own
   2-rank smoke test).
3. *Application-level.* RadialRecipe MPI+CUDA segfaults at the **first RK4
   advance under AMR**, inside `amrex::FillPatchIterator::Initialize`
   (`amrex/Src/Amr/AMReX_AmrLevel.cpp:1004`, reached via
   `AmrLevel::FillPatch`). VRAM is fine when it happens. This is specific to
   this example's code path — the wormhole example does not hit it.

**Triage order — run these before concluding anything:**

```bash
export PATH="$OPENMPI_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$OPENMPI_ROOT/lib:${LD_LIBRARY_PATH:-}"

# 1. Is MPI alive at all on this node?  (If this fails, it is the node.)
mpirun -np 1 hostname
mpirun -np 4 bash -c 'echo "rank $OMPI_COMM_WORLD_RANK of $OMPI_COMM_WORLD_SIZE"'

# 2. Does the CPU solver run multi-rank?  (Wins ~40 min per HQ constraint solve.)
bash grteclyn-wrapper/scripts/build/rebuild_grtresna_mpi.sh   # only if step 1 passed

# 3. Does RadialRecipe survive one AMR advance on 2 GPUs?
#    Short stop-time, max_level>0 — the crash is at the FIRST advance, so a
#    run that reaches t>0 has cleared it.  Never launch a long multi-GPU run
#    before this passes.
```

The RadialRecipe MPI+CUDA binary is normally already built and current — check
its timestamp against `RadialRecipeLevel.cpp` before spending a rebuild. A
rebuild does **not** address failure mode 3, which is a code path, not a link.

#### Arena OOM part-way through a long AMR run

A run can start comfortably and still die of GPU memory hours later: as matter
disperses, more cells get tagged and the Arena grows. Observed on a single
H100 at `N=256, L=128, max_level=3`: **62 GB at t=0 climbing to 77 GB, aborting
at t≈37 of 64** while asking for a further 17 MiB. The abort looks like this,
and is *not* the segfault above:

```
amrex::Abort::0::Arena out of memory!!!
Error: cudaMalloc returned 2: out of memory
[The Arena] space allocated (MB): 77208
Free  GPU global memory (MB): 2
```

What to do:

- **Set `checkpoint_interval` > 0 for any multi-hour run.** With the default
  `-1` there is no restart point and the whole run is lost. This is the single
  most valuable change.
- Budget headroom, not just the starting footprint. Measured on one H100:
  N=240 → 49.8 GB, N=256 → 62 GB *initially*, N=288 → OOMs immediately
  (77.8 GB allocated, 2 MB free). A run that starts at 62 GB has far less
  margin than it appears.
- If it recurs at a resolution you need: raise `regrid_threshold` so fewer
  cells are tagged, lower `max_level`, or step down one grid size — in that
  order, since the first two change the least about the physics.

#### Build the RotatingWormholeCollapse binary (MPI + CUDA)

The wormhole evolution binary (`main3d.gnu.MPI.CUDA.ex`, used by
`scripts/wormhole/run/wormhole_case.py`) is built **with MPI** (multi-GPU) and
CUDA. Single-GPU remains the default (`--gpu N`); multi-GPU is
`--np` / `--gpus` (see [Launch multi-GPU evolutions](#launch-multi-gpu-evolutions-wormhole_case---np----gpus)
below). The recurring pain here is that `make` fails immediately unless
**both** `nvcc` **and** the OpenMPI wrapper (`mpicxx`) are on `PATH`. Use this
exact, copy-pasteable recipe:

```bash
cd Examples/RotatingWormholeCollapse

# nvcc (CUDA 12.x) + local OpenMPI must both be on PATH; OpenMPI libs on
# LD_LIBRARY_PATH. Use the SYSTEM g++ (11.4, nvcc-compatible) -- do NOT source
# the grtresna conda env (its g++ 15 fails nvcc's gcc<=12 check).
export PATH="/usr/local/cuda/bin:$OPENMPI_ROOT/bin:$PATH"
export LD_LIBRARY_PATH="$OPENMPI_ROOT/lib:${LD_LIBRARY_PATH:-}"

make -j8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90
```

`CUDA_ARCH=90` targets the H100 (`sm_90`); adjust for other GPUs (A100 = `80`).
The same recipe rebuilds any Example — just `cd` into its dir. After editing
`Source/**` or an Example's `*.cpp/*.hpp`, re-run `make` (incremental).

**Symptom → fix** for the usual "can't find CUDA/MPI" tension:

| Build error | Cause | Fix |
|-------------|-------|-----|
| `/bin/sh: 1: nvcc: not found` | CUDA not on `PATH` | prepend `/usr/local/cuda/bin` (see `export` above) |
| `*** Unknown mpi wrapper. ... Stop.` | `mpicxx` not on `PATH` | prepend `$OPENMPI_ROOT/bin` (from `.env`) |
| `nvcc fatal: unsupported gnu version` / g++ ≥ 13 errors | GRTresna conda env's `g++` too new for nvcc | run in a shell where that env is **not** on `PATH` |
| runtime `libmpi.so not found` | OpenMPI libs missing at run | set `LD_LIBRARY_PATH=$OPENMPI_ROOT/lib` (from `.env`) |

Verify: a fresh `main3d.gnu.MPI.CUDA.ex` (~130 MB) appears in the Example dir.

#### Launch multi-GPU evolutions (`wormhole_case --np` / `--gpus`)

The binary is MPI-linked; **`wormhole_case.py` / `wormhole_case.sh` can drive
multi-GPU runs** (needed when a single H100 OOMs under deep AMR, e.g. fine
convergence rungs). Pass `--np N` and an explicit GPU list:

```bash
# 4 ranks on physical GPUs 0,1,2,4 (OOM-prone fine dx / high max_level)
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh --gridinit "$G" --full-box \
  --omega 0.25067 --m 1 --dx 0.333 --box-size 64 --max-level 3 --stop-time 40 \
  --mass 0.5 --lambda 170 --mu6 14450 --sponge --no-frames \
  --support-ramp-t-start 8 --support-ramp-t-end 10 --support-ramp-floor 0 \
  --run-suffix conv_dx033 --np 4 --gpus 0,1,2,4

# 2 ranks; if --gpus is omitted, uses consecutive devices from --gpu
bash grteclyn-wrapper/scripts/wormhole/run/wormhole_case.sh ... --np 2 --gpu 6
# equivalent: --np 2 --gpus 6,7
```

| Flag | Effect |
|------|--------|
| `--gpu N` | Single-GPU pin (default), or start of a consecutive block when `--np>1` and `--gpus` is omitted |
| `--np N` | MPI ranks = GPUs for this case (default `1`). Uses local OpenMPI `mpirun` |
| `--gpus a,b,...` | Explicit physical GPU ids; length must equal `--np` |

**How binding works:** each rank exports `CUDA_VISIBLE_DEVICES` to **one**
physical id from `--gpus` (via `WORMHOLE_GPU_IDS`). Do **not** rely on parent
`CUDA_VISIBLE_DEVICES=0,1,…` + `LOCAL_RANK` remapping — OpenMPI/prterun often
drops that env, and every job then lands on GPUs `0..N-1`.

`run_rotating_wormhole.sh PARAMS NGPU` is the older params-file multi-GPU
launcher (same binary); prefer `wormhole_case` for the torus / article
campaign so plotfile delete sidecars stay attached.

### Solver-only AMR smoke tests

```bash
cd $GRTRESNA_ROOT/Examples/ScalarFieldBH
EXE=Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
for case in canonical exotic mixed_exotic; do
  PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
    mpirun --oversubscribe -np 8 ./"${EXE}" params_${case}_amr_test.txt
done
```

### Visualization

**Two commands cover every run.** Point them at the run's `output/` directory
(the one holding `data/`, `small_data/`, `frames/`). All figures land in
`<RUN_DIR>/plots/`; all movies in `<RUN_DIR>/movies/`.

```bash
RUN=runs/rotating_wormhole/<episode>/output      # dir containing data/ small_data/ frames/

# 1) All diagnostic + GW figures (radii = extraction shells, e.g. 12 18 24)
bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh "$RUN" 12 18 24

# 2) Field-slice movies (one .mp4 per frames/<field>_z/)
bash grteclyn-wrapper/scripts/plot/make_movies.sh   "$RUN" --framerate 10
```

`plot_diagnostic.sh` writes, into `<RUN>/plots/`:

| Figure | Built from | Shows |
|--------|-----------|-------|
| `constraints_plot.*` | `data/constraint_norms.dat` | L2 Ham/Mom vs `t` (must stay bounded) |
| `collapse_diagnostics_plot.*` | `data/collapse_diagnostics.dat` (23-col) | `min_lapse`, `min_chi`, `max|K|`, `max_ah_r`, `min_theta_plus`, K-decay lifetime |
| `psi4_analysis_M30_D10.*`, `psi4_analysis_M1000_D0.002.*`, … | `small_data/psi4_mode_l2m0.dat` | 6-panel GW: waveform, retarded+QNM fit, PSD, propagation speed, spectrogram, LIGO strain |

`make_movies.sh` writes `movie_<field>_z.mp4` for each of
`chi, chi_minus_1, K, lapse, phi, Pi, Weyl4_Re, Weyl4_Im, Weyl4_Mag`.

Mass/distance configs for the GW panels are baked into `plot_diagnostic.sh`
(`30:10`, `1000:0.002`, `1000:1`); override the first with env `MASS_MSUN` /
`DISTANCE_MPC`. LIGO panel quantity via `LIGO_QUANTITY=asd|hchar`.

> **Which Ψ₄ is trusted?** Use the Python post-hoc extraction
> `small_data/psi4_mode_l2m0.dat` (validated: `m=0` imaginary part ~1e-5, no
> gauge contamination). Do **NOT** feed the dense in-code C++
> `data/Weyl4_mode_2*.dat` into physics plots -- upstream Weyl/CCE extraction
> is still In-progress and carries an O(1) spurious `m=0` imaginary part; it
> is retained for debugging only.

Full module reference: [`src/grteclyn_wrapper/visualisation/README.md`](src/grteclyn_wrapper/visualisation/README.md).

### Falsification tiers

A high score is a proxy, not proof. `validate_tiers.py` records how far each
candidate survives:

| Tier | Name | Gate |
|------|------|------|
| T0 | `constructed` | Constraint-satisfying data exists |
| T1 | `nontrivial` | Non-flat FTL signal |
| T2 | `operational` | Evolved shortcut beats flat control |
| T3 | `persistent` | Stable long evolution |
| T4 | `observer_ec` | Observer-robust EC margin |
| T5 | `converged` | Resolution-ladder replay agrees |
| T6 | `analytic` | Closed-form back-derivation |

```bash
uv run python scripts/search/validate_tiers.py runs/grtresna_qd/<campaign>
uv run python scripts/search/validate_tiers.py runs/grtresna_qd/<campaign> --min-tier 3
```

### Diagnostics

Each run emits under `data/` (parsed by `read_episode_metrics`):

| File | Probes |
|------|--------|
| `constraint_norms.dat` | Ham/Mom L2; `min_rho_req < 0` → exotic needed |
| `energy_conditions.dat` | Evolved matter NEC/WEC/SEC/DEC |
| `curvature_invariants.dat` | Ricci/K invariants |
| `Weyl4_mode_2{0,1,2}.dat` | **In-code C++ Ψ₄** l=2 m=0/1/2 modes, dense (per coarse step), one Re/Im pair per `--extraction-radii` shell. ⚠️ Unreliable — GRTeclyn's Weyl/CCE extraction is 🔧 In progress upstream; use `small_data/psi4_mode_l2m0.dat` for physics. See [GW/Ψ₄ extraction is now C++-side](#gw--ψ₄-extraction-is-now-c-side-upgrade) |

---

## Reference

### What to edit, build, and run

| Goal | Edit here | Rebuild? | Validate with |
|------|-----------|----------|---------------|
| CMA-ES dimensions, ansätze, warm starts | `search/optimize.py`, `__main__.py` | No | `uv run pytest tests/grtresna/test_grtresna_ring_ansatz.py tests/metrics/ftl/test_solved_geometry_ftl.py -q` |
| Launcher defaults / env knobs | `scripts/campaigns/*/run*.sh` | No | `DRY_RUN=1 ... bash scripts/campaigns/cmaes/run.sh` |
| GRTresna invoke / `.gridinit` conversion | `grtresna/solver.py`, `grtresna/io.py` | Usually no | `uv run pytest tests/grtresna/test_grtresna_integration.py -q` |
| Solved-geometry FTL filter | `metrics/ftl_solved_geometry.py`, `search/solved_ftl_gate.py` | No | `uv run pytest tests/metrics/ftl/test_solved_geometry_ftl.py -q` |
| Scoring weights / objectives | `metrics/score/`, `episode_metrics.py` | No | Re-score campaign or metric tests |
| GRTeclyn evolution, plotfiles, gridinit load | `Examples/RadialRecipe/*`, `Source/Matter/*` | Yes (GRTeclyn) | `BUILD=1 bash scripts/radial/run_radialrecipe_gpu_smoke.sh` |
| GRTresna elliptic solver | `../GRTresna/Examples/ScalarFieldBH/*` | Yes (MPI binary) | AMR smoke tests above |
| **New matter profile** (painter + reference + gate) | `../GRTresna/Source/Matter/BosonStarParams.hpp`, `grtresna/matter/profile_contract.py` | Yes (GRTresna) | `uv run pytest tests/grtresna/test_profile_contract.py -q` — see [Matter-profile contract](#matter-profile-contract--the-rail-for-adding-new-profiles-safely) |
| Rotating Q-torus wormhole ID / collapse | `grtresna/profiles/qball_torus.py`, `scripts/wormhole/id/solve_torus.py`, `scripts/wormhole/run/wormhole_case.py` | Yes if C++ painter changes | See [Rotating Q-torus wormhole](#rotating-q-torus-wormhole--genuine-stationary-eigenstate-support--collapse) |

### Tests

```bash
cd grteclyn-wrapper
uv run pytest tests/grtresna/ -q          # bridge, matter, profiles (258 tests)
uv run pytest tests/grtresna/test_profile_contract.py -q  # matter-profile rail (cross-code t=0 gate)
uv run pytest tests/ -q                    # full suite
ruff check src/grteclyn_wrapper            # lint
```

### Campaign triage

```bash
cd runs/grtresna_qd/<campaign>
python3 - <<'PY'
import json, os
rows=[json.loads(l) for l in open("trajectory.jsonl") if l.strip()]
gr=[r for r in rows if r.get("grtresna_rejected")]
sf=[r for r in rows if r.get("solved_ftl_rejected")]
ev=[r for r in rows if not r.get("grtresna_rejected") and not r.get("solved_ftl_rejected") and not r.get("grtresna_failed")]
print("records", len(rows), "grtresna_rej", len(gr), "solved_ftl_rej", len(sf), "gpu/evolved", len(ev))
for r in ev:
    e=f"eval_{r['eval']:06d}"
    print(e, "score", r.get("score"), "exit", r.get("exit_code"), "score.json", os.path.exists(f"{e}/score.json"))
PY
```

### Stop and verify

```bash
pkill -TERM -f 'runs/grtresna_search|runs/grtresna_qd|campaigns/qd/run.sh|campaigns/cmaes/run.sh'
sleep 5
pkill -KILL -f 'runs/grtresna_search|runs/grtresna_qd|campaigns/qd/run.sh|campaigns/cmaes/run.sh'
nvidia-smi   # expect 0 MiB, no running processes
```

### Research manuscript (TikZ / tectonic)

```bash
cd ../research/neuralspacetime/article && tectonic --keep-logs research.tex   # -> research.pdf
```

First run downloads TeX support files (needs network access once). Source:
[`article/research.tex`](../research/neuralspacetime/article/research.tex).

### Related docs

| Doc | Content |
|-----|---------|
| [`.env.example`](.env.example) | Site-path template (`SIM_ROOT`, `GRTRESNA_ENV`, …) — copy to gitignored `.env` |
| [`scripts/campaigns/README.md`](scripts/campaigns/README.md) | Three-stage pipeline detail, env-var reference |
| [`src/grteclyn_wrapper/grtresna/README.md`](src/grteclyn_wrapper/grtresna/README.md) | GRTresna bridge deep docs |
| [`src/grteclyn_wrapper/gw_search/README.md`](src/grteclyn_wrapper/gw_search/README.md) | LIGO matched-filter methodology |
| [`src/grteclyn_wrapper/visualisation/README.md`](src/grteclyn_wrapper/visualisation/README.md) | Plotting module reference |
| [`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md) | Self-grav boson star fix + caveats |
| [`../research/neuralspacetime/article/research.tex`](../research/neuralspacetime/article/research.tex) | Manuscript source (compile with tectonic above) |
| [`../research/neuralspacetime/MapElitesDynamics.md`](../research/neuralspacetime/MapElitesDynamics.md) | FTL trajectory campaign lab journal |
| [`../research/rotatingwormhole/OrbitalPumpPlan.md`](../research/rotatingwormhole/OrbitalPumpPlan.md) | Rotating wormhole: Q-torus eigenstate support, collapse trigger |
| [`../research/grlab/LabJournal.md`](../research/grlab/LabJournal.md) | GW beam + splash lab journal || [`../research/RL/LabJournal.md`](../research/RL/LabJournal.md) | RL chassis handoff |

