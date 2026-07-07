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
see [Operations](#operations).

> **Roadmap:** critical review of the results' validity and the prioritized
> implementation plan (probe calibration, gauge-honest baseline, ANEC/QI,
> convergence, pump accounting, self-grav, wave-zone GW) live in
> [`NextSteps.md`](NextSteps.md).

---

## What is implemented

| Capability | Where | Notes |
|------------|-------|-------|
| **GRTresna bridge** — `params.txt` → MPI solve → Chombo HDF5 → `.gridinit` | `src/grteclyn_wrapper/grtresna/` | Ham + Mom constraint solve; exotic (`rho<0`) auto-switches to K=0 maximal slicing. Deep docs: [`src/grteclyn_wrapper/grtresna/README.md`](src/grteclyn_wrapper/grtresna/README.md) |
| **Matter sectors** — real scalar, complex scalar / U(1) boson star, Q-ball, shell, trajectory | `grtresna/matter/`, `grtresna/fields/`, `grtresna/profiles/` | Per-lump signed scalars (`GRTresnaIndependentScalars` C++); exotic wedge; self-gravitating boson-star ODE solver (isotropic coords) |
| **Plotfile consumer** — streaming `small_data/` + PNG `frames/` + HDF5 deletion | `scripts/lib/`, `src/.../visualisation/` | `consume_plotfiles` sidecar; **required** for every production run |
| **Ψ₄ / GW extraction** — l=2,m=0 mode amplitudes from `Weyl4_Re/Im` | `src/.../visualisation/process_wave/` | Spherical shells at radii from physics `center`; HQ default radii `8 12 24` |
| **Search algorithms** — MAP-Elites (QD) archive, CMA-ES hill-climb | `src/.../search/qd_search/`, `src/.../search/optimize/` | Shared pre-evolution gates; warm-start from any trajectory |
| **Objectives** — `ftl_first`, `robust_ftl`, `general_ftl`, `critical_collapse`, `gw_beam`, `spacetime_shear` | `src/.../metrics/score/objectives.py` | See [Campaigns](#campaigns) for which objective each campaign uses |
| **Descriptors** — `ftl_lifetime`, `speed_horizon`, `wave_focusing`, `spacetime_shear`, `gw_beam` | `src/.../search/qd_search/descriptors.py` | Behavior axes for the MAP-Elites archive |
| **4D null-geodesic probe** — gauge-invariant FTL shortcut measurement | `src/.../metrics/probes/ftl/` | `search` (cheap) and `hq` (full verify) profiles; continuous emission sweep |
| **Falsification tiers** — T0 constructed → T6 analytic | `scripts/search/validate_tiers.py` | Offline ladder; no rerun needed |
| **Geometry-first projection** — motif scout → GRTresna solve | `src/.../initial_data/motif.py`, `grtresna/fit/motif.py`, `projection/` | Additive second stage; never push fitted matter directly into GRTeclyn |
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

Shared search defaults (grid, gates, 4D probe, objective) live in
`scripts/campaigns/lib/search_common.sh`. HQ defaults in
`scripts/campaigns/lib/promote_common.sh`. QD and CMA-ES **must** stay aligned
(same `search_common.sh`) so warm-started scores are comparable. HQ is
intentionally different: higher resolution and longer time stress-test whether
shortcuts survive refinement.

### Campaign catalog

| Campaign | Launcher | Objective | Descriptor | Matter | What it searches for |
|----------|----------|-----------|------------|--------|----------------------|
| **`general_ftl` (wormhole/ring/spin)** | `scripts/campaigns/general_ftl/run_all.sh` | `general_ftl` | `ftl_lifetime` | real scalar (pinned 15-D subspace) | FTL shortcut on a wormhole/ring/spin geometry; current production path |
| **`ftl_4d` (generic QD)** | `scripts/campaigns/qd/run.sh` | `ftl_first` | `ftl_lifetime` | real scalar shell/ring/free | Generic FTL discovery |
| **`qball_trajectory` (spiral)** | `scripts/campaigns/qball_trajectory/run.sh` | `general_ftl` | `ftl_lifetime` | complex scalar Q-ball, 5 per-lump orbits (39-D) | FTL from compact solitons on retrograde spiral orbits |
| **`qball_trajectory` (Lentz)** | `scripts/campaigns/qball_trajectory/run_lentz.sh` | `general_ftl` | `ftl_lifetime` | canonical Q-ball only, v_max=0.5c | Pure positive-energy FTL (no phantom matter) |
| **`qball_trajectory` (shear)** | `scripts/campaigns/qball_trajectory/run_shear.sh` | `spacetime_shear` | `spacetime_shear` | canonical Q-ball | Extreme non-collapsing frame-dragging shear |
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
| `critical_collapse` | GW focusing / collapse | geometric splash (χ-well + Ψ₄ wave + K-crunch) primary; ρ/focus/lapse secondary; FTL ignored |
| `gw_beam` | Directional GW power + Z-beaming | `(1000×quality + 100×peak + health) × gw_health_multiplier + penalties`; collapse → multiplier 0 |
| `spacetime_shear` | Max curvature, avoid horizon | 1000×curvature_activity + confinement; horizon veto (−500); FTL-agnostic |

---

## How to run

### Prerequisites

From the GRTeclyn repo root:

```bash
cd /path/to/GRTeclyn
uv sync   # Python deps including yt, h5py>=3.10 for the Chombo→gridinit bridge
```

First build (single GPU, no MPI):

```bash
BUILD=1 bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

Binaries (GRTresna MPI solver + GRTeclyn GPU binary) — see [Operations](#operations).

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

**Q-ball trajectory spiral** (compact solitons on retrograde orbits, 39-D):

```bash
QD_NAME=qball_traj_spiral_v3 QD_TARGET_EVALS=400 \
GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/campaigns/qball_trajectory/run.sh
```

**Lentz** (pure canonical matter, v_max=0.5c):

```bash
QD_NAME=qball_traj_lentz_v1 QD_TARGET_EVALS=200 \
  bash scripts/campaigns/qball_trajectory/run_lentz.sh
```

**Spacetime shear** (extreme non-collapsing curvature):

```bash
QD_NAME=qball_traj_shear_v1 QD_TARGET_EVALS=200 \
  bash scripts/campaigns/qball_trajectory/run_shear.sh
```

**GW beam** (directional gravitational-wave emission, canonical Q-balls):

```bash
QD_NAME=gw_beam_v1 QD_TARGET_EVALS=200 \
GPU_IDS="0 1 2 3" \
  bash scripts/campaigns/gw_beam/run.sh
```

**Splash** (critical collapse / GW focusing, canonical bosonic shell):

```bash
QD_NAME=spacetime_splash_v13 QD_TARGET_EVALS=100 \
GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/campaigns/splash/run.sh
```

**Boson star** (complex scalar / U(1), 7-D):

```bash
QD_NAME=boson_star_v1 QD_TARGET_EVALS=80 \
GRTRESNA_MATTER_SECTOR=boson_star GRTRESNA_MATTER_COUPLING=canonical \
GPU_IDS="0 1 2 3" \
  bash scripts/campaigns/boson_star/run.sh
```

**Bosonic shell + FTL (RL chassis)** (~18-D, exotic wedge ON):

```bash
uv sync   # h5py>=3.10 required for GRTresna Chombo→gridinit bridge
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
| Population | = GPU count | `POPULATION` defaults to `#GPU_IDS` |
| σ₀ | 0.05–0.08 | local basin width |
| Warm-start | top-K elites | `WARM_START_TOP_K`, `WARM_START_JITTER` |
| Target | eval budget | `TARGET_EVALS` or `MAX_GENERATIONS × pop` |

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
2. **`general_ftl` needs `GRTECLYN_GEO_DIRECTIONS=x y z`** — wormhole shortcuts live on z; x-only scoring replays elites at the wrong fitness.
3. **HQ `CANDIDATES` is eval/gpu pairs** — e.g. `"46 0 39 1"` not a bare eval list.
4. **Search turns frames off, HQ turns them on** — by design (`search_common.sh` vs `promote_common.sh`).
5. **Promotion must use `N > L`** (or same `L` with larger `N`) to refine the grid. `L=N` only enlarges the domain at `dx=1` — no fidelity gain.

### ALWAYS extract frames on the fly (required)

Every GPU evolution run — QD, CMA-ES, HQ promotion, replay — **MUST** stream
plotfiles through `consume_plotfiles` during the simulation. Do not let heavy
`data/plt*` HDF5 directories accumulate; extract PNG frames + `small_data/`
metrics in flight and delete processed plotfiles immediately.

| Requirement | How |
|-------------|-----|
| Sidecar consumer ON | `CONSUME_PLOTFILES=1` (shell) or `consume_plotfiles=True` (Python) |
| Delete heavy HDF5 | `CONSUMER_DELETE=1` with `--keep-last 3` (HQ default) |
| PNG frames written | Set `GRTECLYN_FRAMES_FIELDS` → outputs in `eval_*/frames/` |
| Verify consumer alive | `ps aux \| grep consume_plotfiles` while GPU is busy |

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

### GRTresna → wrapper → GRTeclyn data flow

```text
Search / CLI overrides
        │
        ▼
  GRTresnaConfig  ──►  params.txt  ──►  GRTresna (MPI)  ──►  InitialDataFinal.3d.hdf5
        │                                                        │
        ▼                                                        ▼
  io/conversion  ──►  initial_data.gridinit  (+ optional .matter.json)
        │                                                        │
        ▼                                                        ▼
  post-load gate (short GPU load, L2_Ham/Mom check)  ──►  GRTeclyn ExternalGridInitialData evolution
                                                                 │
                                                                 ▼
                                                       consume_plotfiles → score
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

The wrapper ships an ODE solver for genuine self-gravitating boson stars
(gravity provides the binding, not an artificial pump well). After a four-bug
fix, a single stable-branch star stays mostly confined through a full coarse
evolution. Full handoff: [`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md).

| t | old broken seed | fixed (stable star + pump) |
|---|-----------------|----------------------------|
| 0 | 0.96 | 0.74 (broad seed settling) |
| 6.4 | 0.74 | **0.97** |
| 9.6 | 0.64 | **0.97** |
| 12.8 | 0.61 | **0.97** |
| 16.0 | **0.58** (dispersed) | **0.90** |

**Caveats (read before trusting):** the "confined fraction" is generous; RMS
radius tightens to ~2.1 by t≈9.6 then spreads back to ~4.2 by t=16 (slow leak /
breathing). High resolution (`max_level=3`) still develops NaNs around t≈6–9 —
a separate numerical-relativity stability issue, not the seed/pump physics.

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

### Score components by tier (`ftl_first`)

| Tier | Component | Weight | Meaning |
|------|-----------|-------:|---------|
| **Validated FTL** | `operational_ftl_geodesic` | 1000 | gauge-invariant geodesic shortcut (reliability-gated) |
| | `ftl_geo_evolving` | 1000 | 4D evolving-metric geodesic shortcut |
| | `operational_ftl` | 400 | evolved coordinate-time shortcut vs flat baseline (Dijkstra) |
| | `ftl_persistence` | 300 | shortcut sustained across last retained plotfiles |
| | `operational_ftl_solved` | 50 | constraint-solved t=0 shortcut |
| **Shaping gradients** | `channel_progress` | 100 | `path_closeness × √(ftl_precursor × shift_drive)` |
| | `ftl_precursor` | 30 | local cone-tilt past `c=1` + superluminal area (graded) |
| | `shift_drive` | 20 | frame-drag motor (`max_shift`) |
| **Health/survival** (gated by `nontriviality_gate`) | `survival` | 70 | `numerical_survival × structural_persistence` |
| | `energy_condition` | 40 | evolved NEC/WEC/SEC/DEC |
| | `instability_penalty` | 15 | geometric drift penalty |
| | `stability` | 10 | bounded stability reward |
| | `comoving_stability` | 8 | co-moving-frame drift |
| | `constraint_health` | 6 | evolved Ham/Mom constraint quality |
| **Penalties** | `exotic_penalty` | 40×weight | NEC-violating matter (graded 0..−1.6) |
| | `stationary_artifact_penalty` | 8 | shift-free geometry demotion |
| | `horizon_penalty` | 500 | trapped-surface veto (non-traversable) |

`structural_persistence = density_retention × morphological_coherence` (3D
connected-component count of matter activity). It also gates the Tier-2 shaping
rewards — a structure that dissipates or fragments cannot bank "promising
precursor" credit.

`SUPERLUMINAL_MARGIN = 0.05` de-saturates the superluminal-fraction descriptor:
only cells genuinely past `c=1.05` count, separating cone-tilted lobes from the
broad shift background.

### Survival = numerical_survival × structural_persistence

`numerical_survival` alone (did the integrator reach `stop_time`?) perversely
rewards junk — empty space is the easiest thing to march to the end. It is
gated by **structural persistence**, itself the product of two failure modes:

- **Density retention** — fraction of peak matter energy density still present
  at `stop_time`. A dissipating config sees its peak ρ collapse toward 0.
- **Morphological coherence** — whether surviving matter is still one connected
  structure or has fragmented. Measured in 3D on a level-0 covering grid;
  returns `~1/k` for `k` comparable pieces.

### Public API

```python
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode, EpisodeMetrics
from grteclyn_wrapper.metrics.score import Score, DEFAULT_WEIGHTS
```

---

## Main results

Full per-eval tables, frames, and movies live in three lab journals:
[`research/neuralspacetime/MapElites.md`](../research/neuralspacetime/MapElites.md)
(FTL search — SH, shell, wormhole),
[`research/neuralspacetime/MapElitesDynamics.md`](../research/neuralspacetime/MapElitesDynamics.md)
(trajectory FTL campaigns), and
[`research/grlab/LabJournal.md`](../research/grlab/LabJournal.md) (GW beam +
splash). Headline numbers below, in roughly chronological order.

### Top 3 findings (critical summary)

Across all campaigns — FTL search, wormhole, trajectory, GW beam, splash,
self-grav — the three findings that matter most. A critical review of these
claims' validity (baseline gauge-dependence, missing probe controls, pump
consistency) and the plan to harden them is in [`NextSteps.md`](NextSteps.md):

**1. Genuine gauge-invariant FTL shortcuts exist in GR with exotic matter — but they are transient, not stable warp bubbles.**

The trajectory campaigns produced the strongest evidence: **eval 122**
(`trajectory_5lump_v1`) survived to t=30 at 256³ HQ with a confirmed **9.4%
end-to-end 4D geodesic shortcut** (frozen peak **21%**). **Eval 118**
(`qball_traj_spiral_v2`) peaked at **~23%** mid-run. These are gauge-invariant
(null rays traced through the full evolving metric, not frozen slices), with
5/5 rays reaching their targets, geodesic drift < 0.002, and shift vectors that
*decay* (the opposite of gauge runaway). This is not a coordinate artifact.

The critical caveat: **in every case the matter disperses and the channel
fades.** Eval 118 confinement falls 53%→23% (rms radius 7.6→18.6); the channel
peaks at t≈19 then decays. Eval 122's FTL window lasts ~16.6 code units (55% of
evolution) before closing. The wormhole HQ (eval 046) opens a real throat mid-run
(peak 7.57%) but a **horizon forms at t≈21** and kills it. No campaign found a
configuration that holds a shortcut open while keeping matter confined. The
honest summary: GR permits transient superluminal shortcuts with exotic matter,
but sustaining them requires a confinement mechanism (self-grav boson star, RL
pump) that this work has not yet solved.

**2. The validation pipeline successfully separates physical shortcuts from gauge artifacts — without it, multiple false positives would have been published as FTL.**

The single most important methodological result. Several high-scoring
candidates turned out to be artifacts:

| Candidate | Apparent signal | What it actually was | Caught by |
|-----------|----------------|----------------------|-----------|
| `trajectory_5lump_v1` eval 008 | 24.62% f_geo (Stage 0 leader) | Low-res artifact → 0% at HQ | HQ resolution ladder |
| SH eval 151 | 2.26c max speed | Gauge collapse (geodesics untrusted from t=3.2) | Geodesic trust flag |
| SH eval 101 | f_geo = 0.753 | Single untrusted timestep (shift runaway to 1.01) | Timestep trust + shift monitoring |
| GW beam v3 eval 51 | Score 336 | Numerical bomb: Ham crash → Ψ₄ noise scored as GW | Health multiplier + Ψ₄ truncation |
| Wormhole v22 | operational_ftl = 0 | Real 19% shortcut invisible to coordinate Dijkstra | 4D evolving geodesic vs frozen slice |

The last row is the key insight: a stationary wormhole (β≈0) reads subluminal on
coordinate Dijkstra (`operational_ftl = 0`) because lapse drops below 1, but the
4D null-geodesic tracer measures a real ~19% proper-distance shortcut through
the throat. The pipeline decoupled gauge-dependent coordinate speed from
gauge-invariant traversability. Without the 4D evolving probe + dispersion gate +
geodesic trust flag + HQ ladder, the project would have reported both false
positives (eval 008, 151, 101) and false negatives (the wormhole throat).

**3. Search design — ansatz and matter sector — is the dominant lever; optimizer tuning is secondary.**

The search infrastructure (MAP-Elites, CMA-ES, GRTresna-in-the-loop) is mature,
but what determined success was *what* was searched:

| Comparison | Result | Implication |
|-----------|--------|-------------|
| Trajectory ansatz vs SH | 54% FTL hit rate vs 1.3% (**40×**) | Per-lump independent orbits give the optimizer geometric freedom SH lacks |
| Real scalar vs boson shell | 32/92 FTL vs 0/94 (**zero**) | Complex U(1) boson shells do not open geodesic shortcuts under matched conditions |
| Self-grav boson star | Still disperses at high res (NaN @ t≈6–9) | Confinement is the unsolved bottleneck, not the search |

The GW beam campaign confirmed this from the opposite direction: even with a
working search and hard gates, the best directional GW emission was a **weak
steady hum** (~30% beam ratio at P ~ 10⁻⁵–10⁻⁴), not a laser. The search found
what the physics permits — coherent multi-lump clusters with breathing
quadrupoles — and no amount of optimizer tuning changes that. Future work should
invest in the matter model (self-grav confinement, RL pump actuation) rather
than further search-space refinement.

---

### FTL search — spherical-harmonic `scalar_sh_ftl_v22` (200 evals, general_ftl)

First genuine geodesic FTL: **eval 189** (score 470.6). Source:
[MapElites.md §SH campaign](../research/neuralspacetime/MapElites.md#sh-campaign-results-scalar_sh_ftl_v22-200-evals-2025-06-24).

| Property | Value |
|----------|-------|
| `f_geo_peak` | 3.1% @ t=9.6 (rises dynamically from 0) |
| `ftl_geo_evolving` | 0.101 |
| `ftl_lifetime` | 1.0 (present at every timestep) |
| `max_speed` | 1.33c → 1.21c (shift *decays* 0.36→0.04 — healthy gauge) |
| geodesic trust | 5/5 rays, drift <0.002, all timesteps trusted |
| exotic_fraction | 0.51 (3/5 lumps exotic) |

Two false positives discarded: eval 151 (2.26c = gauge collapse artifact,
geodesics untrusted from t=3.2) and eval 101 (f_geo=0.75 from one untrusted
timestep, shift runaway to 1.01). Pipeline funnel: 200 sampled → 74 GRTresna
rejected → 23 postload rejected → **73 gpu_ok (36%)**.

### `general_ftl` wormhole — `general_ftl_wormhole_v21` (200 evals, 15-D pinned)

Stationary, non-translating wormholes in a 15-D pinned subspace. Source:
[MapElites.md §v22 final results](../research/neuralspacetime/MapElites.md#v22-final-results-top-3--ftl-champions).

| Eval | Score | `ftl_geo_evolving` | `f_geo_peak` | `op_ftl` | Survival | Role |
|------|------:|-------------------:|-------------:|---------:|---------:|------|
| **063** | **165.6** | **19.3%** | 4.2% | 0 | 0.94 | score + 4D record holder |
| **191** | 161.9 | 18.5% | 3.8% | 0 | **1.00** | champion (survival, f_op_peak) |
| **174** | 157.4 | 18.5% | 3.8% | 0 | 1.00 | stable variant |

**Scoring paradox resolved:** `operational_ftl = 0` while `ftl_geo_evolving ≈ 19%`.
In a stationary wormhole (β≈0), lapse-dominated coordinate speed reads subluminal
(Dijkstra → 0), but proper-distance contraction through the throat gives a real
~19% 4D null-geodesic shortcut. The pipeline decouples coordinate gauge artifacts
from physical traversable shortcuts.

**CMA-ES refinement** (eval 063 → eval 046): score 165.6 → **179.8** (+14.2),
`ftl_geo_evolving` 19.3% → **20.3%**, survival → 1.00. Basin-tightening, not a new
mechanism. Source:
[MapElites.md §v22 CMA-ES](../research/neuralspacetime/MapElites.md#v22-cma-es-wormhole-refinement-general_ftl_wormhole_cmaes_v1-2026-06-18).

**HQ promotion** (eval 046, 256³, t=30): throat opens mid-run, peak **7.57%** 4D
@ t≈15.6, then **horizon kills** at t≈21 (score cliff to −546). A mid-run FTL
demonstrator, not a t=30 survivor. Source:
[MapElites.md §HQ eval 046](../research/neuralspacetime/MapElites.md#hq-eval-046-final-results-t30).

### Trajectory FTL — `qball_traj_spiral_v2` (200/200 evals, general_ftl)

Best search candidate **eval 118** — a breathing, retrograde, mostly-exotic
Q-ball shell. Source:
[MapElitesDynamics.md §spiral v2](../research/neuralspacetime/MapElitesDynamics.md#qball_traj_spiral_v2--dispersion-gated-spiral-qd-complete-2026-07-01-200200-evals).

| | QD (128³, t=16) | HQ (256³, t=16) | HQ (256³, t=30) |
|--|----------------:|----------------:|----------------:|
| **Score** | **603.39** | 511.89 | **224.20** |
| `operational_ftl` | 0.347 | 0.346 | 0.099 (dispersal-gated) |
| `ftl_geo_evolving` | 0.306 | 0.225 | 0.150 |
| 4D `f_geo_evol` | peak 17.7% (t_emit≈12) | 13.0% | **13.0%** |
| frozen `f_geo` peak | — | 22.1% @ t≈15.1 | **22.8% @ t≈19.2** |
| `max_local_speed` | 1.47 c | 1.46 c | — |
| confinement | ~53% @ t=0 | ~35% | **23%** (rms 7.6→18.6) |

**Verdict:** A real, gauge-invariant geodesic shortcut (~13% end-to-end, peaking
~23% mid-run). Numerics survive to t=30. But **matter disperses** — confinement
falls 53%→23% and the channel peaks near t≈19 then fades. A transient shortcut
from a dissolving "motor," not a stable warp bubble.

### HQ validation — `trajectory_5lump_v1` (5 elites at 256³, t=30)

Source:
[MapElitesDynamics.md §HQ validation](../research/neuralspacetime/MapElitesDynamics.md#hq-validation-results-trajectory_5lump_v1-only).

| Eval | Stage 0 score | Stage 0 f_geo | HQ f_geo_evol | HQ f_geo_peak | HQ status | Verdict |
|------|--------------:|--------------:|--------------:|--------------:|-----------|---------|
| 122 | 1237.6 | 8.51% | **9.40%** | **20.97%** | survived t=30 | **CONFIRMED** |
| 115 | 1367.9 | 10.63% | 12.5% | 20.3% | crashed t=21 | confirmed (transient) |
| 050 | 1039.5 | 10.82% | 7.4% | 20.3% | crashed t=19 | confirmed (transient) |
| 111 | 1389.6 | 17.37% | 8.6% | 19.8% | crashed t=8.6 | confirmed (short) |
| 008 | 1166.8 | 24.62% | 0.0% | — | survived t=30 | **FALSE POSITIVE** |

All genuinely FTL configs converge to **~20% peak f_geo** at HQ (resolution
ceiling). Eval 008's 24.62% Stage 0 signal was entirely a low-res artifact. Eval
122 is the only eval that both survived to t=30 AND confirmed FTL.

### SH vs trajectory ansatz (Stage 0 head-to-head)

Source:
[MapElitesDynamics.md §SH vs trajectory](../research/neuralspacetime/MapElitesDynamics.md#campaign-comparison-scalar_sh_ftl_v22-vs-trajectory_5lump_v1-2026-06-25).

| Metric | **SH v22** (202 evals) | **Trajectory v1** (130 evals) | Factor |
|--------|----------------------:|---------------------------:|--------|
| Best stable score | 470.6 | **1367.9** | 2.9× |
| Best stable f_geo_peak | 2.12% | **10.63%** | 5.0× |
| Best HQ-confirmed f_geo_evol | — | **9.40%** | — |
| FTL hit rate (per GPU eval) | 1.3% | **54%** | ~40× |

The trajectory ansatz (per-lump independent orbits) decisively outperforms the
spherical-harmonic ansatz for FTL discovery.

### Paired shell — boson vs scalar (200+200 evals, general_ftl)

Source:
[MapElites.md §paired shell](../research/neuralspacetime/MapElites.md#paired-shell-ftl-comparison-boson-vs-scalar-2026-06-23).

| Metric | **Boson** | **Scalar** |
|--------|----------:|-----------:|
| `gpu_ok` | 94 (47%) | 92 (46%) |
| `f_geo_peak > 0` | **0 / 94** | **33 / 92** |
| `ftl_geo_evolving > 0` | **0 / 94** | **32 / 92** |
| Best score | 21.6 (eval 100, no FTL) | **869.3** (eval 166, persist 0.76) |

**Verdict:** Boson static exotic shell **never opens a geodesic shortcut** (0/94);
real scalar exotic shell does (32/92). Boson arm **rejected** for FTL RL;
scalar eval 166/126 promoted for RL chassis Gate 2.

### Wormhole / shell HQ leaderboard (`qd_20260605T155951Z`, t=50)

Source: [MapElites.md §campaign log](../research/neuralspacetime/MapElites.md#campaign-log--runs-analysis).

| Rank | eval | HQ score | `op_ftl` | `channel` | `shift` | Role |
|------|------|--------:|---------:|----------:|--------:|------|
| 1 | **106** | **1423** | **1.000** | 0.423 | 0.179 | HQ winner |
| 2 | **117** | **1346** | 0.920 | 0.436 | 0.190 | Channel backup |
| 3 | **011** | **1274** | 0.885 | 0.302 | 0.091 | Search leader |
| 4 | **094** | **1089** | 0.658 | 0.454 | 0.206 | Best channel/shift |

Resolution ladder (eval 057, `op_ftl`=1.0 holds across all):

| Run | L | N | t | max c | Notes |
|-----|--:|--:|--:|------:|-------|
| `val16hq2` | 128 | 128 | 16 | 1.192 | best t=16 HQ |
| `val30hq` | 128 | 128 | 30 | 1.276 | peak c |
| `val100hq` | 128 | 128 | 100 | 1.205 | long GPU-only |
| `val256hq` | 256 | 256 | 100 | 1.196 | 2× domain; ~3× tighter Ham/Mom |

### GW laser search — `gw_beam_qd100_v4` (100/100 evals, gw_beam objective)

Canonical Q-balls on trajectory orbits, scored for directional Ψ₄ emission
(Z-axis beaming). Source:
[LabJournal.md §gw_beam v4](../research/grlab/LabJournal.md#2026-07-03-gw_beam_qd100_v4-complete-eval-61--88-analysis).

Hard gates held: 77/100 collapse modes crushed to ~−116; 22 healthy survivors;
5 archive elites.

| | eval 88 (best score) | eval 61 (best beam) |
|--|----------------------|---------------------|
| Score | **3.09** | 2.82 |
| mean Ψ₄ power | 4.5×10⁻⁴ | **6.4×10⁻⁴** |
| beam_ratio | 14% | **~30%** (to 40% late) |
| max ‖Ham‖₂ | 0.08 | 0.14 |

**Verdict:** Neither run is a strong GW emitter — both produce a weak steady hum
(P ~ 10⁻⁵–10⁻⁴), not a merger chirp or beamed burst. The t=0 Ψ₄ spike is an
initial/near-zone transient. The ~30% Z-beaming (eval 61) comes from a coherent,
fast, compact multi-lump cluster + breathing quadrupole, not a clean radiative-
zone binary.

**Reward-hacking closure (v3 → v4):** v3's optimizer built a **numerical bomb**
instead of a GW laser — crash the Hamiltonian → grid fills with high-frequency
noise → second-derivative Ψ₄ reports "infinite wave power" (eval 51 @ **336**
vs trustworthy eval 7 @ 3.4). Permanently closed in v4 by three hard gates:
Ψ₄ time-series truncation at the spike, archive admission requiring
`tier ≥ CONSTRUCTED`, and a multiplicative `gw_health_multiplier` (→0 on
collapse). Source:
[LabJournal.md §v3→v4](../research/grlab/LabJournal.md#2026-07-03-gw_beam_qd100_v3--v4-collapse-mode-reward-hacking).

### Splash campaign — `spacetime_splash` (critical_collapse objective)

Canonical bosonic shell, scored for gravitational-wave focusing / critical
collapse. Uses `critical_collapse` objective (geometric splash: χ-well + Ψ₄ wave
+ K-crunch primary; ρ/focus/lapse secondary). Source:
[MapElites.md §handoff to RL](../research/neuralspacetime/MapElites.md#handoff-to-rl).

Interim pump proof: splash boson **`spacetime_splash_v14_moving/eval_000010`**
held as the pump-actuation reference until RL Gate 2 passes on the scalar FTL
chassis. The splash campaign validated the `critical_collapse` scorer and the
`SPLASH_MODE` early-termination (stop once matter disperses after peak, typically
t≈10–12), and feeds the RL handoff path.

### Self-gravitating boson star (single-star smoke)

See [Self-gravitating boson star](#self-gravitating-boson-star) above. After the
four-bug fix, confinement holds ~0.90 at t=16 (was 0.58). Not yet committed;
high-res instability remains open. Source:
[`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md).

---

## Operations

### Build GRTresna

Production searches use MPI `mpicxx.gfortran` (`RANKS=8` default). Needs
`CHOMBO_HOME` and the `grtresna` env on `PATH`.

```bash
GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
CHOMBO_HOME=/home/jovyan/nachevsky/test/simulation/Chombo/lib
cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

Produces `Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex`.
Serial debug: build with `MPI=FALSE`, run `RANKS=1`.

### Build GRTeclyn (GPU evolution binary)

The RadialRecipe GPU binary (`main3d.gnu.CUDA.ex`) is built **without MPI** and
**without the conda/grtresna env g++** (nvcc requires gcc ≤ 12).

```bash
cd /home/jovyan/nachevsky/test/simulation/GRTeclyn/Examples/RadialRecipe
PATH="/usr/local/cuda/bin:$PATH" NO_MPI_CHECKING=TRUE \
  make USE_MPI=FALSE USE_CUDA=TRUE -j$(nproc)
```

| Setting | Value | Why |
|---------|-------|-----|
| `USE_MPI` | `FALSE` | Single-GPU runs |
| `USE_CUDA` | `TRUE` | GPU kernels (AMReX + CCZ4 RHS) |
| **Do NOT** put grtresna env on PATH | — | Its `g++ 15.x` breaks nvcc's gcc ≤ 12 check |

### Solver-only AMR smoke tests

```bash
cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
EXE=Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
for case in canonical exotic mixed_exotic; do
  PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
    mpirun --oversubscribe -np 8 ./"${EXE}" params_${case}_amr_test.txt
done
```

### Visualization

Post-run figures → `<RUN_DIR>/plots/`:

```bash
bash grteclyn-wrapper/scripts/plot/plot_diagnostic.sh \
  runs/grtresna_promote/<name> 8 12 24
```

Writes: `constraints_plot.*` (L2 Ham/Mom), `collapse_diagnostics_plot.*`,
`psi4_analysis_M*_D*.*` (6-panel GW: waveforms, QNM fit, PSD, propagation speed,
spectrogram, LIGO strain). Optional env: `MASS_MSUN`, `DISTANCE_MPC`.

Frame movies → `<RUN_DIR>/movies/`:

```bash
bash grteclyn-wrapper/scripts/plot/make_movies.sh runs/grtresna_promote/<name> --framerate 10
```

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

### Tests

```bash
cd grteclyn-wrapper
uv run pytest tests/grtresna/ -q          # bridge, matter, profiles (258 tests)
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

### Related docs

| Doc | Content |
|-----|---------|
| [`scripts/campaigns/README.md`](scripts/campaigns/README.md) | Three-stage pipeline detail, env-var reference |
| [`src/grteclyn_wrapper/grtresna/README.md`](src/grteclyn_wrapper/grtresna/README.md) | GRTresna bridge deep docs |
| [`src/grteclyn_wrapper/gw_search/README.md`](src/grteclyn_wrapper/gw_search/README.md) | LIGO matched-filter methodology |
| [`src/grteclyn_wrapper/visualisation/README.md`](src/grteclyn_wrapper/visualisation/README.md) | Plotting module reference |
| [`SELFGRAV_HANDOFF.md`](SELFGRAV_HANDOFF.md) | Self-grav boson star fix + caveats |
| [`../research/neuralspacetime/MapElitesDynamics.md`](../research/neuralspacetime/MapElitesDynamics.md) | FTL trajectory campaign lab journal |
| [`../research/grlab/LabJournal.md`](../research/grlab/LabJournal.md) | GW beam + splash lab journal |
| [`../research/RL/LabJournal.md`](../research/RL/LabJournal.md) | RL chassis handoff |
