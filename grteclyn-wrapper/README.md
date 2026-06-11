# grteclyn-wrapper

Run isolated GRTeclyn episodes from the repo root (`GRTeclyn/`). The wrapper can stream plotfiles into `small_data/` during GPU runs and delete heavy HDF5 dirs afterward.

## ALWAYS extract frames on the fly (required)

**Every GPU evolution run — QD, CMA-ES, HQ promotion, replay — MUST stream plotfiles through `consume_plotfiles` during the simulation.** Do not let heavy `data/plt*` HDF5 directories accumulate; extract PNG frames + `small_data/` metrics in flight and delete processed plotfiles immediately.

This is **not optional**. A run without the sidecar consumer is incomplete for physics review.

| Requirement | How |
|-------------|-----|
| Sidecar consumer ON | `consume_plotfiles=True` in Python (`evaluate_overrides` / `run_episode`) or `CONSUME_PLOTFILES=1` in shell launchers |
| Delete heavy HDF5 | `consumer_delete=True` / `CONSUMER_DELETE=1` with `--keep-last 3` (HQ default in `run_promote_qd_batch.sh`) |
| PNG frames written | Set `GRTECLYN_FRAMES_FIELDS` (see `run_promote_qd_batch.sh` / `scripts/lib/promote_common.sh`) — outputs land in `eval_*/frames/` |
| Verify consumer alive | `ps aux \| grep consume_plotfiles` while GPU is busy; `frames/` should grow during evolution |

**Wrong:** launching replay/promotion without frame env vars, or inspecting only `data/plt*` after a long `t=50` run (disk explosion, no reviewable movies).

**Right:** use `run_promote_qd_batch.sh`, or `replay_grtresna_eval.py` with `--consumer-keep-last 3` and the `GRTECLYN_FRAMES_*` exports below.

### Post-load gate vs main evolution (do not confuse them)

HQ promotion has **two** GPU phases. Only the second one produces frames:

| Phase | Dir | `stop_time` | `consume_plotfiles` | Frames? |
|-------|-----|-------------|---------------------|---------|
| **Post-load gate** | `eval_*/postload_gate/postload_gate/` | **`0.01`** (short) | **OFF** | **No** |
| **Main evolution** | `eval_*/` (root `params.txt`) | **`30` or `50`** (HQ script dependent) | **ON** | **Yes** → `eval_*/frames/` |

The post-load gate (`projection/postload_gate.py`) is a **constraint sanity check** after GRTresna — it loads the `.gridinit`, advances a few steps, and reads `L2_Ham` / `L2_Mom`. It must **never** inherit promotion `stop_time=50`; gate-specific params always win over search/promotion overrides.

**Corrupted run symptoms** (GPU busy, zero frames):

- `postload_gate/postload_gate/params.txt` shows `stop_time = 50.0` instead of `0.01`
- `ps aux | grep consume_plotfiles` returns nothing while a GPU is at 100%
- `eval_*/frames/` empty but `postload_gate/postload_gate/run.log` shows hundreds of steps

**Recovery:** stop the run, save `initial_data.gridinit`, wipe the episode dir, relaunch with `--gridinit` (GPU-only, skips GRTresna + postload) or fresh `replay_grtresna_eval.py` after pulling the `postload_gate.py` fix.

**Healthy run check:**

```bash
# Consumer must be alive during main evolution
ps aux | grep 'consume_plotfiles.*l128n256_qd_eval'

# Main evolution uses episode-root params (HQ stop_time), not postload_gate/
grep stop_time runs/grtresna_promote/l128n256_qd_eval074/params.txt
grep stop_time runs/grtresna_promote/l128n256_qd_eval074/postload_gate/postload_gate/params.txt  # only if gate ran
```

```bash
export GRTECLYN_FRAMES_FIELDS="lump_activity scalar_activity phi_lump_sum Pi_lump_sum chi chi_minus_1 local_speed shift1 rho_req"
export GRTECLYN_FRAMES_ZOOM=none   # full-domain slices for HQ
export GRTECLYN_PROJECTION_FIELDS="lump_activity scalar_activity"
export GRTECLYN_PROJECTION_AXES="x y z"
export GRTECLYN_PROJECTION_METHOD=mip
```

## Prerequisites

From the GRTeclyn repo root:

```bash
cd /path/to/GRTeclyn
uv sync   # Python deps including yt for plotfile extraction
```

First build (single GPU, no MPI):

```bash
BUILD=1 bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

## Quick start

Run from `GRTeclyn/grteclyn-wrapper` unless noted.

```bash
# CMA-ES scalar search (see Search pipeline for ansätze and knobs)
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=shell \
  bash scripts/search/run_grtresna_search.sh

# MAP-Elites quality-diversity (diverse FTL families)
# Post-load constraint gate is ON by default (POSTLOAD_GATE=1) so only
# physically self-consistent initial data reaches GPU evolution.
QD_ITERATIONS=8 BINS=8 GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 \
  LUMPS=5 SHELL_PROFILE=compact \
  POSTLOAD_MAX_HAM_L2=1e-2 POSTLOAD_MAX_MOM_L2=1e-2 \
  bash scripts/search/run_grtresna_qd_search.sh

# Cheap smoke: no GRTresna solve, no GPU evolution
DRY_RUN=1 MAX_GENERATIONS=1 GPU_IDS="0 1" GRTRESNA_ANSATZ=ring \
  bash scripts/search/run_grtresna_search.sh
```

**Stop and verify** — kill the whole campaign tree, then check CPU and GPU:

```bash
pkill -TERM -f 'runs/grtresna_search|runs/grtresna_qd|runs/grtresna_refine|run_grtresna_search|run_grtresna_qd_search'
sleep 5
pkill -KILL -f 'runs/grtresna_search|runs/grtresna_qd|runs/grtresna_refine|run_grtresna_search|run_grtresna_qd_search'
ps -eo pid,ppid,pgid,pcpu,pmem,etime,args \
  | awk '/run_grtresna_search|run_grtresna_qd_search|grteclyn_wrapper|Main_ScalarFieldBH3d|main3d.gnu.CUDA.ex|consume_plotfiles|prterun|mpirun|orterun/ && !/awk/ {print}'
nvidia-smi   # expect 0MiB usage, no running processes
```

**Campaign triage** — start with `trajectory.jsonl`:

```bash
cd runs/grtresna_search/optimize_<timestamp>
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

For GPU survivors, inspect `score.json` directly: `metrics.general_ftl_solved.f_op`, `metrics.general_ftl_evolved.f_op`, `max_local_speed`, `max_shift`, and `components.{operational_ftl,ftl_precursor,shift_drive,curvature_activity}`. Treat `F_op ~ 1` or near-degeneracy `max_local_speed` as suspicious; mild `F_op=0.01..0.05` needs HQ replay before any physics claim.

**Falsification tiers** (offline, no rerun):

```bash
uv run python scripts/search/validate_tiers.py runs/grtresna_search/<campaign>
uv run python scripts/search/validate_tiers.py runs/grtresna_search/<campaign> --min-tier 3
```

**Production scripts** — full inventory in [`scripts/README.md`](scripts/README.md):

| Script | Use for |
|--------|---------|
| `run_grtresna_search.sh` | CMA-ES: GRTresna → `.gridinit` → GPU evolution |
| `run_grtresna_qd_search.sh` | MAP-Elites archive over shell space |
| `run_promote_qd_batch.sh` | Batch HQ promotion of QD elites (env-driven candidate list) |
| `replay_grtresna_eval.py` | Single-eval HQ promotion / GPU-only continuation |
| `run_radialrecipe_gpu_smoke.sh` | Single-GPU smoke/build after C++ edits |
| `project_geometry_motif.py` | Geometry-first scout → GRTresna projection → post-load gate |

---

## Geometry-first projection

Geometry-first RadialRecipe searches (`run_ftl_search_cmaes.sh`,
`run_ftl_search_nonspherical.sh`, `run_ftl_search_directional.sh`) are efficient
**scouts**: they find FTL-relevant geometric motifs without solving the full 3D
momentum constraint.  The 1D Hamiltonian reconstruction sets `Pi = 0`, so scouts
find lenses and dipoles, not transport motors by themselves.

The hybrid projection path promotes scout winners to constraint-satisfying initial
data:

```text
geometry-first episode → motif.json → fitted_matter.json → GRTresna solve
  → projection_report_preservation.json → postload_gate.json → optional replay
```

**Critical rule:** fitted matter is never pushed directly into GRTeclyn.  GRTresna
must re-solve both `H` and `M_i` on the scalar-lump ansatz.

| Path | Use when |
|------|----------|
| Geometry-first scout | Cheap motif discovery and FTL shape ranking |
| `project_geometry_motif.py` | Promote top-K scouts to physical initial data |
| Matter-first `ring` / `shell` / `free` | Blind lump-space discovery (see [Search pipeline](#search-pipeline)) |

Matter-first launchers are unchanged.  Projection is an additive second stage.

### Artifacts

| File | Contents |
|------|----------|
| `motif.json` | Support regions, polarity, shortcut metrics, `static_lens_only` |
| `fitted_matter.json` | GRTresna lumps (`amp`, `width`, `center`, `velocity`, `exotic`, …) |
| `momentum_target.json` | Desired motor direction; no velocity when `static_lens_only` |
| `initial_data.gridinit` | GRTresna-projected data |
| `projection_report_preservation.json` | Scout vs solved motif retention |
| `postload_gate.json` | GRTeclyn `L2_Ham` / `L2_Mom` on loaded `.gridinit` |

### Commands

```bash
# Fit only (login node; no GRTresna, no GPU)
uv run python scripts/search/project_geometry_motif.py \
  ../runs/ftl_search/optimize_<ts>/eval_000123 \
  --out-dir ../runs/geometry_projection/eval_000123 \
  --mode fit-only

# GRTresna solve + motif preservation (cluster)
uv run python scripts/search/project_geometry_motif.py \
  ../runs/ftl_search/optimize_<ts>/eval_000123 \
  --out-dir ../runs/geometry_projection/eval_000123 \
  --mode solve-only --mpi-ranks 8 --max-lumps 3

# Solve, post-load gate, short replay
CUDA_VISIBLE_DEVICES=0 uv run python scripts/search/project_geometry_motif.py \
  ../runs/ftl_search/optimize_<ts>/eval_000123 \
  --out-dir ../runs/geometry_projection/eval_000123 \
  --mode solve-and-evolve --cuda-device 0 --stop-time 2.0
```

Reject before long GPU replays when motif preservation, GRTresna convergence, or
post-load constraints fail.  Only claim physics after evolved constraints,
co-moving operational FTL, and resolution replay pass.

**Modules:** `initial_data/motif.py`, `grtresna/motif_fit.py`,
`metrics/motif_preservation.py`, `projection/postload_gate.py`.
**Tests:** `tests/test_hybrid_projection.py`.

---

## Matter–geometry consistency (physically satisfied solutions)

GRTresna solves the constraints with **independent, signed scalar lumps** (each lump
has its own sign and there are no inter-lump cross terms). Earlier, GRTeclyn evolved
a single *summed* scalar, so on overlapping or exotic lumps the evolved stress-energy
disagreed with what GRTresna solved — the loaded initial data violated the constraints
at `t=0`. The matter-first pipeline now keeps the two in lock-step end-to-end:

```text
GRTresna signed-lump solve → .gridinit (+ per-lump phi_k/Pi_k channels, matter metadata)
  → GRTeclyn GRTresnaIndependentScalars matter (mass-matched V = ½ m² φ²)
  → post-load gate: L2_Ham / L2_Mom small  ⇒ valid Einstein solution
```

What this guarantees and what it does **not**:

| "Physically satisfied" means… | Status | Evidence |
|------|--------|----------|
| Hamiltonian & momentum constraints hold at load and through evolution (a valid spacetime) | **Yes** — enforced by the post-load gate | Smokes: post-load `L2_Ham≈4e-3`, `L2_Mom≈2e-4`; survivors show `constraint_health≈0.95`, `survival=1`, `anec_condition=1` |
| The spacetime is a *useful* warp/FTL geometry (`operational_ftl > 0`) | **Separate search goal** — not implied by constraint validity | Found only when the QD/CMA-ES objective surfaces an operational elite |

So a `gpu_ok` candidate is a genuine, constraint-satisfying solution of GR; whether it
is an *interesting* (transport-capable) one is what the search campaign explores.

**Implementation**
- C++ matter: `Source/Matter/GRTresnaIndependentScalars.{hpp,impl.hpp}` (+ `…Vars/D1/D2/Advec`, `GRTresnaScalarLayout.hpp`, `GRTresnaScalarPotential.hpp`); RadialRecipe dispatch in `Examples/RadialRecipe/RadialRecipeMatterDispatch.hpp`, state in `StateVariables.hpp` (`NUM_CCZ4_VARS + 2 + 2·GRTRESNA_MAX_INDEPENDENT_SCALARS`).
- Python bridge: `grtresna/lump_fields.py` (paint per-lump channels), `grtresna/matter_wiring.py` (matter overrides + metadata), `grtresna/solver.py` (`apply_exotic_safe_solver`, `config_has_exotic_lump`).
- Gate: `projection/postload_gate.py`, wired via `search/grtresna_evaluation_gates.py`.
- Tests: `tests/test_matter_geometry_consistency.py`, `tests/test_grtresna_postload_gate_integration.py`.
- End-to-end smokes: `scripts/search/run_matter_geometry_smokes.sh` (canonical, exotic, mixed-shell).

**Knobs** (default the post-load gate **on** in `run_grtresna_qd_search.sh`):

| Env var | Default | Meaning |
|---------|---------|---------|
| `POSTLOAD_GATE` | `1` | Enable the post-load constraint gate |
| `POSTLOAD_MAX_HAM_L2` | `1e-2` | Reject if loaded `L2_Ham` exceeds this |
| `POSTLOAD_MAX_MOM_L2` | `1e-2` | Reject if loaded `L2_Mom` exceeds this |

---

## Search pipeline

GRTresna solves 3D constraints per candidate; survivors evolve on GPU and are scored with `objective_mode=ftl_first`. **CMA-ES** maximizes one scalar and collapses to a basin; **MAP-Elites** keeps diverse families in a behavior archive.

### Stages

```mermaid
flowchart LR
    subgraph shell["16D shell ansatz"]
        V["16D vector"]
        L["LUMPS Fibonacci sites<br/>SHELL_PROFILE=compact"]
        V --> L
    end
    subgraph s1["Stage 1 — Discovery"]
        QD["MAP-Elites<br/>run_grtresna_qd_search.sh"]
        A["archive.json"]
    end
    subgraph s2["Stage 2 — Refine optional"]
        RF["CMA-ES warm start<br/>run_grtresna_search.sh"]
    end
    subgraph s3["Stage 3 — Promote"]
        PR["replay / batch promote"]
        T["validate_tiers.py"]
    end
    L --> QD --> A
    A -->|"precursors only"| RF
    A -->|"operational_ftl winners"| PR
    RF -->|"op_ftl breakthrough"| PR
    PR --> T
```

| Stage | Tool | Runs dir | Typical settings |
|-------|------|----------|------------------|
| 1 — discovery | MAP-Elites | `runs/grtresna_qd/qd_<ts>/` | `L=64 N=64`, `t=8–16`, `BINS=8–10` |
| 2 — refine (optional) | CMA-ES warm start | `runs/grtresna_refine/optimize_<ts>/` | `sigma0=0.10`, `t=16` |
| 3 — promote | replay / batch | `runs/grtresna_promote/` | `L=128 N=256`, `t=50`, fresh GRTresna |

Skip stage 2 when QD already yields `operational_ftl > 0` elites; go straight to stage 3.

**Per-eval loop** (every CMA-ES member and QD candidate):

1. Sample ansatz parameters → GRTresna MPI solve (Ham + Mom). Exotic (`rho<0`) candidates auto-switch to the K=0 maximal-slicing solver (`apply_exotic_safe_solver`).
2. Reject if convergence missing, NaN, or above threshold
3. Solved-geometry FTL gate on `.gridinit` (cheap, pre-GPU)
4. **Post-load constraint gate** — **short** GPU launch (`stop_time=0.01`, **no** `consume_plotfiles`) that loads the `.gridinit` and rejects if `L2_Ham`/`L2_Mom` exceed `--postload-max-{ham,mom}-l2` (default `1e-2`). Writes to `postload_gate/` only — **no frames here**.
5. **Main** GRTeclyn GPU evolution (`stop_time=50` for HQ) → `consume_plotfiles` sidecar → PNG frames in `frames/` → `ftl_first` score
6. Feed fitness to CMA-ES or MAP-Elites archive cell

Pre-evolution gates 2–4 are centralized in `search/grtresna_evaluation_gates.py` (`apply_grtresna_pre_evolution_gates`), shared by CMA-ES and MAP-Elites.

Campaign artifacts: `trajectory.jsonl`, `eval_XXXXXX/{metadata,score}.json`, `initial_data.gridinit`, `initial_data.matter.json` (matter layout), **`frames/`** (required — streamed on the fly by `consume_plotfiles`, not optional).

### Ansätze

| Ansatz | Dims | Use |
|--------|------|-----|
| `ring` | 14D | Equatorial refinement; lumps on one circle (+ small z wobble) |
| `shell` | 16D | Full-sphere discovery (default); Fibonacci lattice + toroidal/poloidal flow |
| `free` | 11×`LUMPS` | Unbiased atlas / audit (`55D` for `LUMPS=5`) |

All modes decode into the same `cfg.lumps` and GRTeclyn path. `LUMPS` sets placement resolution in `ring`/`shell` (not optimizer dims); in `free` it scales search dimension. Ring parameter details: `search/optimize.py`, `tests/test_grtresna_ring_ansatz.py`.

Shell 16D parameters: `amp`, `width`, `radius`, `thickness`, orientation axis `(theta, phi)`, toroidal/poloidal/radial velocities, `omega`, dipole/quadrupole asymmetry, exotic fraction/phase, `mode`, `radial_jitter`.

### Scoring

QD/CMA-ES campaigns score with `objective_mode="ftl_first"` (`metrics/score.py`,
`score_episode`). The total is a weighted sum split into three roles. Scores are
**not comparable** across campaigns with different weights/components — always
re-score with the current code before ranking historical runs.

**Tier 1 — validated FTL** (the actual goal; kept dominant so a real result always wins):

| Component | Weight | Meaning |
|-----------|-------:|---------|
| `operational_ftl_geodesic` | 1500 | Gauge-invariant geodesic shortcut (reliability-gated) |
| `operational_ftl` | 400 | Evolved coordinate-time shortcut vs flat baseline (Dijkstra) |
| `ftl_persistence` | 300 | Shortcut sustained across the last retained plotfiles |
| `operational_ftl_solved` | 180 | Constraint-solved t=0 shortcut, localization + peak-margin gated |

**Tier 2 — shaping gradients** (coordinate cone-tilt heuristics that only *point* toward FTL; cut to ~40% so they guide without out-voting health):

| Component | Weight | Meaning |
|-----------|-------:|---------|
| `channel_progress` | 150 | `path_closeness × √(ftl_precursor × shift_drive)` |
| `ftl_precursor` | 30 | Local cone-tilt past `c=1` + superluminal area (graded, not binary) |
| `shift_drive` | 20 | Frame-drag motor (`max_shift`) |

**Tier 3 — health/survival** (gated by a non-triviality factor so flat Minkowski cannot bank them; boosted so a coherent, persistent, stable structure ranks on par with a strong shaping signal):

| Component | Weight | Meaning |
|-----------|-------:|---------|
| `survival` | 70 | `numerical_survival × structural_persistence` (see below) |
| `instability_penalty` | 15 | Geometric drift penalty |
| `stability` | 10 | Bounded stability reward |
| `comoving_stability` | 8 | Co-moving-frame drift |
| `constraint_health` | 6 | Evolved Ham/Mom constraint quality |

`exotic_penalty` (NEC-violating matter) and a graded `stationary_artifact_penalty` subtract; a trapped-surface `horizon_penalty` (−500 weight) vetoes high-local-speed artifacts that are not traversable.

#### Survival = numerical_survival × structural_persistence

`numerical_survival` alone (did the integrator reach `stop_time`?) perversely
rewards junk — empty/dissipated space is the easiest thing to march to the end.
It is therefore gated by **structural persistence**, itself the product of two
independent failure modes:

- **Density retention** — fraction of the peak matter energy density still
  present at `stop_time` (`final_peak_rho_required / max_rho_required` from
  `constraint_norms.dat`). A configuration that dissipates sees its peak rho
  collapse toward 0.
- **Morphological coherence** — whether the surviving matter is still a single
  connected structure or has fragmented into disconnected pieces. Density
  retention is blind to this (a bubble that splits into two equally-dense lobes
  keeps its peak). `structure_coherence` (in `metrics/probes/ftl/general.py`)
  samples the 3D scalar activity `√(φ² + Π²)` on a level-0 covering grid, labels
  the connected components, and returns `~1/k` for `k` comparable pieces. It is
  measured in **3D** (a single 2D slice is orientation-dependent), so a
  connected dumbbell/ring reads as coherent (1.0) while genuinely dispersed
  matter is penalised. It rides on `general_ftl_evolved` and defaults to 1.0
  when the plotfile is unavailable.

`structural_persistence` also gates the Tier-2 shaping rewards
(`ftl_precursor`, `channel_progress`, `shift_drive`) — a structure that
dissipates or fragments cannot bank "promising precursor" credit for cones
tilting in matter that no longer holds together.

#### Superluminal-fraction margin (de-saturation)

The coordinate light speed exceeds `c=1` wherever there is *any* appreciable
shift, so a broad frame-drag background sits at `c ~ 1.03` across most of the
slice — a gauge feature, not a channel. A bare `> 1.0` test marked ~75% of the
slice superluminal, saturating both the QD `superluminal_fraction` descriptor
(everything piled into one bin) and the precursor area term. `SUPERLUMINAL_MARGIN
= 0.05` (in `metrics/probes/ftl/general.py`) counts only cells genuinely past
`c=1.05`, separating the cone-tilted lobes (`c ~ 1.08–1.18`) from the shift
background. The precursor speed/area scales (`0.15`) and the `speed_super`
descriptor target (`0.15`, in `search/qd_search.py`) are calibrated to the
resulting noise-free range (~`0.04–0.17`).

**`SHELL_PROFILE`** bounds presets (`search/optimize.py`, `SHELL_PROFILE` env):

| Profile | Width | Radius | Use |
|---------|-------|--------|-----|
| `compact` (default) | `1.8–3.0` | `1.5–6.0` | Sharper lumps; less diffuse-blob basin |
| `middle` | `1.8–4.0` | `1.5–6.0` | Reproduce `optimize_20260604T170329Z` |
| `outer_precursor` | `2.8–4.0` | `4.0–6.0` | Bias toward eval-128 outer regime |
| `inner_shift` | `1.8–3.0` | `1.5–3.0` | Bias toward eval-57 compact shift regime |

Matter current scales roughly as `amp² × velocity / width`. Diffuse shells converge cleanly but source tiny shift; the compact profile and `shift_drive` term push away from blob-collapse.

**Blob-collapse red flags** — if top frames all look alike:

| Signal | Good | Red flag |
|--------|------|----------|
| `local_speed_z` | Localized region above `c=1` | Sub-luminal, identical across evals |
| `shift1_z` | Growing shift magnitude | `~10⁻²`, max shift `~0.03` |
| `score.json` | Variation in `shift_drive`, `channel_progress`, `operational_ftl` | `operational_ftl=0`, score from stability only |

Quick triage:

```bash
uv run python3 - <<'PY'
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
root = Path("../runs/grtresna_search/optimize_<timestamp>")
for ep in sorted(root.glob("eval_*")):
    if not (ep / "score.json").exists(): continue
    s = score_episode(read_episode_metrics(ep, ftl_L=8.0), target_stop_time=2.0, objective_mode="ftl_first")
    c = s.components
    print(ep.name, f"score={s.total:.2f}", f"op={c.get('operational_ftl',0):.3f}",
          f"ch={c.get('channel_progress',0):.3f}", f"shift={c.get('shift_drive',0):.3f}")
PY
```

### CMA-ES

```bash
cd /path/to/GRTeclyn/grteclyn-wrapper

# Refinement: 14D ring
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=ring \
  bash scripts/search/run_grtresna_search.sh

# Discovery: 16D shell (default for new families)
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=shell \
  bash scripts/search/run_grtresna_search.sh

# Audit: 55D free
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=free LUMPS=5 \
  bash scripts/search/run_grtresna_search.sh
```

Useful overrides:

```bash
GPU_IDS="0 1 2 3" bash scripts/search/run_grtresna_search.sh   # fewer GPUs
WARM_START_TRAJECTORY=/path/to/trajectory.jsonl GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/search/run_grtresna_search.sh                   # continue prior run
NO_CONSUME=1 GPU_IDS="0 1" MAX_GENERATIONS=1 \
  bash scripts/search/run_grtresna_search.sh                   # keep plotfiles
RANKS=1 bash scripts/search/run_grtresna_search.sh             # serial GRTresna (MPI=FALSE build)
```

| Knob | Default | Meaning |
|------|---------|---------|
| `GRTRESNA_ANSATZ` | `ring` | `ring` / `shell` / `free` |
| `SHELL_PROFILE` | `compact` | Shell bounds preset |
| `GRTRESNA_FULL_Z` | `1` for `shell` | Full-z box (no reflective `z=0` plane) |
| `LUMPS` | `5` | Scalar clouds (55D only when `free`) |
| `RANKS` | `8` | MPI ranks per GRTresna solve |
| `GPU_IDS` | `0 1 2 3` (often all GPUs) | Also sets CMA-ES population |
| `MAX_GENERATIONS` | `50` | CMA-ES generations |
| `OBJECTIVE_MODE` | `ftl_first` | FTL terms dominate health/stability |
| `WARM_START_TRAJECTORY` | unset | Prior `trajectory.jsonl` for warm start |
| `WARM_START_TOP_K` / `JITTER` | `8` / `0.08` | Seeds and jitter fraction |
| `GRTRESNA_TIMEOUT` | `900` | Wall-clock seconds per GRTresna solve |
| `ITERATIONS` | `30` | Max NL iterations per solve (ceiling; early exit usually stops sooner) |
| `GRTRESNA_NL_EXIT_TOLERANCE` | `1.0` | Stop when both Ham and Mom fall below this `%` (`0` disables) |
| `GRTRESNA_NL_STALL_TOLERANCE` | `0.02` | Stop when per-iter Ham/Mom improvement stalls below this fraction (`0` disables) |
| `GRTRESNA_GRIDINIT_WORKERS` | `0` | Parallel Chombo→`.gridinit` threads (`0` = auto, `min(32, cpu_count)`) |
| `STOP_TIME` | `2.0` | GPU evolution stop time |
| `RANDOM_INJECTION_FRACTION` | `0.25` | Per-generation random candidates |

Warm-start example (compact shell + `channel_progress`):

```bash
WARM_START_TRAJECTORY=../runs/grtresna_search/optimize_20260604T170329Z/trajectory.jsonl \
  WARM_START_TOP_K=8 GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 \
  GRTRESNA_ANSATZ=shell LUMPS=5 SHELL_PROFILE=compact \
  SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR=0.9 \
  bash scripts/search/run_grtresna_search.sh
```

Outputs: `runs/grtresna_search/optimize_<timestamp>/`.

### GRTresna solve speed (on by default)

Both `run_grtresna_search.sh` and `run_grtresna_qd_search.sh` enable two complementary
short-cuts so candidates rarely burn the full `ITERATIONS=30` NL loop:

1. **Early NL exit** — `NL_exit_tolerance` stops once Ham **and** Mom are below
   `GRTRESNA_NL_EXIT_TOLERANCE` percent (default `1.0`, well inside the
   post-solve `GRTRESNA_MAX_HAM_PCT=5` gate). `NL_stall_tolerance` (default
   `0.02`) stops diverging or plateaued solves even earlier.
2. **Parallel gridinit conversion** — after each solve, `grtresna/io.py` paints AMR
   boxes with NumPy vectorization plus a `ThreadPoolExecutor`
   (`GRTRESNA_GRIDINIT_WORKERS`, default `0` = auto). Per-lump `phi_lumpK` /
   `Pi_lumpK` channels use vectorized `meshgrid` painting.

Typical converged candidates finish in **~13–19** NL iterations instead of 30.
Disable absolute early exit for debugging: `GRTRESNA_NL_EXIT_TOLERANCE=0`.
Force serial conversion: `GRTRESNA_GRIDINIT_WORKERS=1`.

CLI equivalents (also on `optimize` / `qd`):

```bash
--grtresna-nl-exit-tolerance 1.0 \
--grtresna-nl-stall-tolerance 0.02 \
--grtresna-gridinit-workers 0
```

### MAP-Elites

```bash
QD_ITERATIONS=8 BINS=8 GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 \
  LUMPS=5 SHELL_PROFILE=compact \
  bash scripts/search/run_grtresna_qd_search.sh
```

Long discovery run:

```bash
QD_ITERATIONS=16 BINS=10 STOP_TIME=8 GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 \
  LUMPS=5 SHELL_PROFILE=compact \
  GRTRESNA_MAX_HAM_PCT=10 GRTRESNA_MAX_MOM_PCT=10 \
  bash scripts/search/run_grtresna_qd_search.sh
```

| Env var | Default | Meaning |
|---------|---------|---------|
| `QD_ITERATIONS` | `8` | Archive update rounds |
| `BINS` | `8` | Descriptor grid resolution per axis |
| `SHELL_PROFILE` | `compact` | Same presets as CMA-ES |
| `STOP_TIME` | `2.0` | GPU stop time (use `16` for long runs) |
| `GRTRESNA_TIMEOUT` | `900` | Per-solve wall clock |
| `ITERATIONS` | `30` | Max NL iterations per solve (ceiling; early exit usually stops sooner) |
| `GRTRESNA_NL_EXIT_TOLERANCE` | `1.0` | Stop when both Ham and Mom fall below this `%` (`0` disables) |
| `GRTRESNA_NL_STALL_TOLERANCE` | `0.02` | Stop when per-iter Ham/Mom improvement stalls (`0` disables) |
| `GRTRESNA_GRIDINIT_WORKERS` | `0` | Parallel Chombo→`.gridinit` threads (`0` = auto) |

Outputs: `runs/grtresna_qd/qd_<timestamp>/` with `trajectory.jsonl`, `archive.json`, `eval_XXXXXX/`.

Resume a campaign toward a larger eval budget (replaces the old `continue_center_fixed_qd.sh` wrapper):

```bash
RUNS_DIR=runs/grtresna_qd/ftl_discovery \
QD_NAME=qd_ftl_discovery_20260609T162553Z \
QD_RESUME=1 QD_TARGET_EVALS=400 \
DESCRIPTOR_MODE=speed_horizon \
bash scripts/search/run_grtresna_qd_search.sh
```

| | CMA-ES | MAP-Elites |
|--|--------|------------|
| Goal | One best score | Diverse elites across descriptor cells |
| Descriptors | Score only | `path_closeness`, `mechanism_balance`, tier tags |
| Best for | Fast scalar climb | Operational FTL + path-quality precursors together |

### Post-QD refinement

When QD finds strong precursors but `operational_ftl = 0`, warm-start CMA-ES with small steps:

```bash
cd grteclyn-wrapper/scripts/search

RUNS_DIR=../../runs/grtresna_refine \
WARM_START_TRAJECTORY=../../runs/grtresna_qd/qd_<timestamp>/trajectory.jsonl \
WARM_START_TOP_K=10 WARM_START_JITTER=0.05 SIGMA0=0.10 \
MAX_GENERATIONS=12 POPULATION=8 STOP_TIME=16 PLOT_INTERVAL=48 \
GRTRESNA_ANSATZ=shell SHELL_PROFILE=compact GRTRESNA_FULL_Z=1 \
GRTRESNA_TIMEOUT=900 RANDOM_INJECTION_FRACTION=0.10 \
GPU_IDS="0 1 2 3 4 5 6 7" RANKS=8 \
bash run_grtresna_search.sh
```

CMA-ES writes `trajectory.jsonl` once per generation (empty file during gen-1 GRTresna is normal). If pinned at `op_ftl=0` for several gens: raise `SIGMA0` to `0.15` or return to QD.

### Promotion

> **Resolution rule:** promotion must use **`N > L`** (or same `L` with larger `N`) to refine the grid. `L=N` only enlarges the domain at `dx=1` — no fidelity gain.

**FTL discovery HQ batch** (current production path):

```bash
# Defaults: QD source ftl_discovery campaign, L=128 N=256, t=30, keep-last=3
# Promote specific evals (024, 055, 025 from MapElites.md):
QD_RUN=../runs/grtresna_qd/ftl_discovery/qd_ftl_discovery_20260609T162553Z \
CANDIDATES="024 0 055 1 025 2" \
NAME_PREFIX=ftl_discovery STOP_TIME=30 \
bash scripts/search/run_promote_qd_batch.sh
```

Outputs: `runs/grtresna_promote/l128n256t30_ftl_discovery_qd_eval*/` with **`frames/`** populated during **main** GPU evolution, not during post-load gate.

**Auto-pick top operational elites** from `trajectory.jsonl`:

```bash
QD_RUN=../runs/grtresna_qd/qd_20260605T155951Z \
TOP_K=8 MIN_OPERATIONAL_FTL=0.03 NAME_PREFIX=qd STOP_TIME=50 \
bash scripts/search/run_promote_qd_batch.sh
```

**Center-fixed top-8** (explicit candidate list):

```bash
QD_RUN=../runs/grtresna_qd/center_fixed_search/qd_20260609T094634Z \
CANDIDATES="324 0 114 1 146 2 169 3 358 4 314 5 136 6 228 7" \
NAME_PREFIX=center_fixed STOP_TIME=30 \
bash scripts/search/run_promote_qd_batch.sh
```

**Safe anchor revalidation** (`000114` from MapElites.md):

```bash
QD_RUN=../runs/grtresna_qd/center_fixed_search/qd_20260609T094634Z \
CANDIDATES="114 0" NAME_PREFIX=revalidate STOP_TIME=30 \
bash scripts/search/run_promote_qd_batch.sh
```

Dry-run (print launches without starting GPUs):

```bash
DRY_RUN=1 CANDIDATES="024 0" NAME_PREFIX=ftl_discovery \
bash scripts/search/run_promote_qd_batch.sh
```

**GPU-only continuation** (reuse solved `.gridinit`, skip GRTresna + postload — goes straight to framed evolution):

```bash
cd grteclyn-wrapper
export GRTECLYN_FRAMES_FIELDS="phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req"
export GRTECLYN_FRAMES_ZOOM=none
export GRTECLYN_PROJECTION_FIELDS=scalar_activity
export GRTECLYN_PROJECTION_AXES="x y z"
export GRTECLYN_PROJECTION_METHOD=mip
nohup .venv/bin/python scripts/search/replay_grtresna_eval.py \
  ../runs/grtresna_qd/qd_<ts>/eval_000074 \
  --name l128n256_qd_eval074 --runs-dir ../runs/grtresna_promote --gpu 0 \
  --gridinit /path/to/initial_data.gridinit \
  --n-full 256 --l-full 128 --max-level 3 --stop-time 50 --plot-interval 24 \
  --consumer-keep-last 2 \
  > ../runs/grtresna_promote/l128n256_qd_eval074.log 2>&1 &
```

**Matter-binding check is mandatory for every GPU-only `.gridinit` replay.** GRTresna
solves the geometry together with a matched scalar-matter layout. A replay that
loads the `.gridinit` but drops the matched matter params can still show low
Ham/Mom constraints, but it is not a valid physics validation because GRTeclyn
will evolve the wrong matter model. `replay_grtresna_eval.py` now copies these
keys from the source eval `params.txt` when `--gridinit` is used:

```text
recipe_matter_model = grtresna_independent_scalars
recipe_num_scalar_fields = 5
recipe_scalar_field_signs = ...
recipe_scalar_mass = 0.1
```

Before trusting frames or `score.json`, verify the promoted run root
`params.txt` contains `recipe_initial_data_file = "...gridinit"` **and** the
matched `recipe_matter_model`/scalar-sign lines. If those lines are missing,
stop the run and relaunch after fixing the replay script; do not interpret the
movies as a bound geometry-matter result.

**Top-3 replay** (full path: GRTresna + postload + framed evolution):

```bash
cd grteclyn-wrapper
export GRTECLYN_FRAMES_FIELDS="phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req"
export GRTECLYN_FRAMES_ZOOM=none
export GRTECLYN_PROJECTION_FIELDS=scalar_activity
export GRTECLYN_PROJECTION_AXES="x y z"
export GRTECLYN_PROJECTION_METHOD=mip
QD=../runs/grtresna_qd/qd_20260608T175934Z
for pair in "074 0" "030 1" "023 2"; do
  read -r EVAL GPU <<< "$pair"
  nohup .venv/bin/python scripts/search/replay_grtresna_eval.py \
    "$QD/eval_$(printf '%06d' $((10#${EVAL})))" \
    --name "l128n256_qd_eval${EVAL}" \
    --runs-dir ../runs/grtresna_promote --gpu "$GPU" \
    --n-full 256 --l-full 128 --grtresna-domain-l 128 \
    --max-level 3 --regrid-threshold 0.02 --stop-time 50 --plot-interval 24 \
    --grtresna-ranks 8 --grtresna-timeout 7200 \
    --grtresna-max-ham-pct 10 --grtresna-max-mom-pct 10 \
    --consumer-keep-last 2 \
    > "../runs/grtresna_promote/l128n256_qd_eval${EVAL}.log" 2>&1 &
done
```

| Setting | QD search | HQ promotion |
|---------|-----------|----------------|
| `L_full` | 64 | **128** |
| `N_full` | 64 | **256** |
| `dx = L/N` | 1.0 | **0.5** |
| `stop_time` | 8 | **50** |
| `max_level` | 2 | **3** |

**HQ leaderboard** (`t=50`, `runs/grtresna_promote/l128n256_qd_eval*/score.json`):

| Rank | eval | HQ score | `op_ftl` | `channel` | `shift` | Role |
|------|------|--------:|---------:|----------:|--------:|------|
| 1 | **106** | **1423** | **1.000** | 0.423 | 0.179 | HQ winner |
| 2 | **117** | **1346** | **0.920** | **0.436** | **0.190** | Channel backup |
| 3 | **011** | **1274** | **0.885** | 0.302 | 0.091 | Search leader |
| 4 | **094** | **1089** | 0.658 | **0.454** | **0.206** | Best channel/shift |

**Single-eval replay** (`eval_000057` T5 ladder reference):

```bash
cd grteclyn-wrapper
export GRTECLYN_FRAMES_FIELDS="lump_activity scalar_activity phi_lump_sum Pi_lump_sum chi chi_minus_1 local_speed shift1 rho_req"
export GRTECLYN_FRAMES_ZOOM=none
export GRTECLYN_PROJECTION_FIELDS="lump_activity scalar_activity"
export GRTECLYN_PROJECTION_AXES="x y z"
export GRTECLYN_PROJECTION_METHOD=mip

# Full replay (GRTresna + GPU)
nohup .venv/bin/python scripts/search/replay_grtresna_eval.py \
  ../runs/grtresna_qd/qd_20260605T062448Z/eval_000057 \
  --name val16hq2_qd_eval057 --runs-dir ../runs/grtresna_promote --gpu 0 \
  --n-full 128 --l-full 128 --grtresna-domain-l 128 \
  --max-level 3 --stop-time 16 --plot-interval 48 \
  --grtresna-ranks 8 --grtresna-timeout 7200 \
  --grtresna-max-ham-pct 10 --grtresna-max-mom-pct 10 \
  --consumer-keep-last 3 \
  > ../runs/grtresna_promote/val16hq2_qd_eval057.log 2>&1 &

# GPU-only continuation (reuse gridinit)
nohup .venv/bin/python scripts/search/replay_grtresna_eval.py \
  ../runs/grtresna_qd/qd_20260605T062448Z/eval_000057 \
  --name val30hq_qd_eval057 --runs-dir ../runs/grtresna_promote --gpu 0 \
  --gridinit ../runs/grtresna_promote/val16hq2_qd_eval057/initial_data.gridinit \
  --n-full 128 --l-full 128 --max-level 3 --stop-time 30 --plot-interval 48 \
  --consumer-keep-last 3 \
  > ../runs/grtresna_promote/val30hq_qd_eval057.log 2>&1 &
```

| Run | L | N | t | max c | `op_ftl` | Notes |
|-----|--:|--:|--:|------:|---------:|-------|
| `eval_000057` (QD) | 64 | 64 | 2 | 1.065 | 1.0 | QD breakthrough |
| `val16hq2_qd_eval057` | 128 | 128 | 16 | 1.192 | 1.0 | Best t=16 HQ |
| `val30hq_qd_eval057` | 128 | 128 | 30 | **1.276** | 1.0 | GPU-only; peak c |
| `val100hq_qd_eval057` | 128 | 128 | 100 | 1.205 | 1.0 | Long GPU-only |
| `val256hq_qd_eval057` | 256 | 256 | 100 | 1.196 | 1.0 | 2× domain; ~3× tighter Ham/Mom |

`val100hq` read: central two-lobed structure is real across `shift1`, `rho_req`, `chi`, `local_speed` (scorer reports `max c=1.205`); scalars decay by t=100; outer-domain rings are likely boundary junk on L=128; PNG colorbars auto-scale per frame — compare metrics in `score.json`, not brightness. Larger boxes require a fresh GRTresna solve (misaligned gridinit reuse misplaces the shell).

Stitch movies:

```bash
bash scripts/plot/make_movies.sh runs/grtresna_promote/l128n256_qd_eval{011,106,117} --framerate 8
```

---

## Operations

### Build GRTresna

Production searches use MPI `mpicxx.gfortran` (`RANKS=8` default). Needs `CHOMBO_HOME` and the `grtresna` env on `PATH`.

```bash
GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
CHOMBO_HOME=/home/jovyan/nachevsky/test/simulation/Chombo/lib
cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

Produces `Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex`. Serial debug: build with `MPI=FALSE`, run `RANKS=1`.

After header-only C++ edits, force relink:

```bash
rm -f Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex \
      o/3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI/{Main_ScalarFieldBH,Grids,ScalarField,MyMatterFunctions}.o
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

### Solver-only AMR smoke tests

```bash
cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
EXE=Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
for case in canonical exotic mixed_exotic; do
  PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
    mpirun --oversubscribe -np 8 ./"${EXE}" params_${case}_amr_test.txt
done
```

### Single-GPU run

Pick one initial-data source:

| Env var | Example |
|---------|---------|
| `SEED_NAME` | `ellis_bronnikov` |
| `CANDIDATE_ID` | `bubble_wall_016` |
| `NONSPHERICAL_ID` | `quadrupole_bubble_001` |

```bash
BUILD=0 SEED_NAME=ellis_bronnikov CUDA_VISIBLE_DEVICES_OVERRIDE=0 \
  bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

Outputs: `runs/radialrecipe_gpu_smoke/<name>_gpu_t<stop_time>_<stamp>/`.

### Plotfile consumer (required — default on)

With `CONSUME_PLOTFILES=1`: sidecar `consume_plotfiles` streams `small_data/`, **PNG frames to `frames/`**, deletes processed HDF5 in flight (`--keep-last N`), post-sim drain. **Never disable for production search or promotion runs.**

Runs **without** the sidecar: `postload_gate` (`consume_plotfiles=False` by design — see header callout). Do not mistake a long postload run for main evolution.

| Variable | Default | Meaning |
|----------|---------|---------|
| `CONSUME_PLOTFILES` | `1` | Enable streaming extraction |
| `CONSUMER_DELETE` | `1` | Delete HDF5 plot dirs after extract |
| `CONSUMER_RADII` | `4 8` | Extraction radii |
| `PLOT_INTERVAL` | `10` | Plotfile cadence |
| `STOP_TIME` | `2.0` | Simulation stop time |
| `N_FULL` | `64` | Grid resolution |

Disable: `CONSUME_PLOTFILES=0 bash .../run_radialrecipe_gpu_smoke.sh`

**Frame fixes** (promotion runs): slice zoom now defaults to `L_full` (`GRTECLYN_FRAMES_ZOOM`); stable per-field color limits in `consume_plotfiles.py` (override `GRTECLYN_FRAMES_ZLIM_*` or `GRTECLYN_FRAMES_AUTO_ZLIM=1`).

### Post-run plots

```bash
bash grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode_dir>
bash grteclyn-wrapper/scripts/plot/plot_run_radial.sh runs/<episode_dir> --no-delete   # manual drain
```

### Falsification tiers

A high score is a proxy, not proof. `validation_tiers.py` records how far each candidate survives:

| Tier | Name | Gate |
|------|------|------|
| T0 | `constructed` | Constraint-satisfying data exists |
| T1 | `nontrivial` | Non-flat FTL signal |
| T2 | `operational` | Evolved shortcut beats flat control |
| T3 | `persistent` | Stable long evolution |
| T4 | `observer_ec` | Observer-robust EC margin |
| T5 | `converged` | Resolution-ladder replay agrees |
| T6 | `analytic` | Closed-form back-derivation |

Absent diagnostics report `unavailable`. T5/T6 need promotion runs (`extra={"resolution_converged":...}`).

---

## Reference

### What to edit, build, and run

| Goal | Edit here | Rebuild? | Validate with |
|------|-----------|----------|---------------|
| CMA-ES dimensions, ansätze, warm starts | `search/optimize.py`, `__main__.py` | No | `uv run pytest tests/test_grtresna_ring_ansatz.py tests/test_solved_geometry_ftl.py -q` |
| Launcher defaults / env knobs | `scripts/search/run_grtresna_search.sh` | No | `DRY_RUN=1 MAX_GENERATIONS=1 GPU_IDS="0 1" bash scripts/search/run_grtresna_search.sh` |
| GRTresna invoke / `.gridinit` conversion | `grtresna/solver.py`, `grtresna/io.py` | Usually no | `uv run pytest tests/test_grtresna_integration.py -q` |
| Solved-geometry FTL filter | `metrics/ftl_solved_geometry.py`, `search/solved_ftl_gate.py` | No | `uv run pytest tests/test_solved_geometry_ftl.py -q` |
| Scoring weights | `metrics/score.py`, `episode_metrics.py` | No | Re-score campaign or metric tests |
| GRTeclyn evolution, plotfiles, gridinit load | `Examples/RadialRecipe/*`, `Source/Matter/*` | Yes (GRTeclyn) | `BUILD=1 bash scripts/radial/run_radialrecipe_gpu_smoke.sh` |
| GRTresna elliptic solver | `../GRTresna/Examples/ScalarFieldBH/*` | Yes (MPI binary) | AMR smoke tests above |

### GRTresna bridge

GRTresna solves Hamiltonian + momentum constraints in 3D and hands GRTeclyn constraint-satisfying `.gridinit` data. Deep integration docs: [`src/grteclyn_wrapper/README.md`](src/grteclyn_wrapper/README.md) (*GRTresna integration*).

**Changelog (search pipeline):**

- GRTresna ↔ GRTeclyn bridge via `ExternalGridInitialData` (~2× lower initial constraint error vs radial recipe)
- Momentum-carrying scalar lumps (velocity, rotation, exotic flags) — new frame-dragging search axis
- Exotic-matter AMR on maximal-slicing (`K=0`) path; three AMR smoke fixtures pass
- GRTresna-in-the-loop CMA-ES/QD with graded Ham/Mom rejection before GPU
- Pre-evolution solved-geometry FTL gate (`search/solved_ftl_gate.py`) — cheap filter on `.gridinit` at t=0
- `channel_progress` + `SHELL_PROFILE` presets (2026-06-05)
- GRTresna NL early exit (`NL_exit_tolerance` + `NL_stall_tolerance`) and parallel
  Chombo→`.gridinit` conversion exposed via launcher env vars / CLI (2026-06-09)
- Survival hardened (2026-06-11): `survival = numerical_survival × structural_persistence`,
  with `structural_persistence = density_retention × morphological_coherence` (3D
  connected-component count of the matter activity). Persistence also gates the
  cone-tilt shaping rewards. `ftl_first` weights rebalanced so survival/stability
  matter on par with the (cut) shaping gradients; `ftl_precursor` de-saturated into
  a gradient; `SUPERLUMINAL_MARGIN=0.05` de-saturates the QD superluminal-fraction
  descriptor

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

### Solved-FTL gate (summary)

Scores `.gridinit` at t=0 before GPU (~1 s/candidate). Rejects flat/degenerate slices; graded fitness for near-misses. Policy: `search/solved_ftl_gate.py`. Launcher knobs: `SOLVED_FTL_F_OP_FLOOR`, `SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR`, `SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR`, `SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR`, `SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED`, `SOLVED_FTL_MAX_PHYSICAL_F_OP`, `SOLVED_FTL_REJECTION_SPEED_TARGET`.

Degeneracy guard: slices with `max_local_speed` or `F_op` above physical ceilings are flagged (numerical artifacts near the York existence boundary). Exploration setting: `SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR=0.9` for permissive shell campaigns.

### Exotic matter (summary)

FTL channels need NEC-violating matter. Per-lump exotic flags + maximal-slicing solver handle `rho < 0` in AMR. Search auto-switches exotic candidates to K=0 path. Details: package README + AMR smoke loop above. Rescore prior campaigns: `scripts/search/rescore_grtresna_solved_ftl.py <campaign_dir>`.

### Diagnostics

Each run emits under `data/` (parsed by `read_episode_metrics`):

| File | Probes |
|------|--------|
| `constraint_norms.dat` | Ham/Mom L2; `min_rho_req < 0` → exotic needed |
| `energy_conditions.dat` | Evolved matter NEC/WEC/SEC/DEC |
| `curvature_invariants.dat` | Ricci/K invariants |

Operational FTL (`ftl_general.py`, Dijkstra vs flat baseline) is mechanism-agnostic: `general_ftl` (solved) vs `general_ftl_evolved` (plotfile).

**Two findings:**

1. Canonical `ScalarField` gives `matter_* EC ~ 0` by construction; use `recipe_exotic_matter = 1` (`--phantom`) for genuine evolved exotic violations. Geometry-sourced `T^eff` is post-hoc via `warpfactory.py`.
2. Evolved FTL uses plotfile metric fields; distinguish initial vs evolved channel to spot gauge artifacts.

---

## Campaign snapshots

**Current production:** discovery `runs/grtresna_qd/qd_20260605T155951Z/` (8 operational elites); HQ promotion `runs/grtresna_promote/l128n256_qd_eval{011,016,058,077,094,106,117,121}/`; article montage `research/neuralspacetime/articlepictures/shell_hq_promotion_tiles.png`.

| Campaign | Tool | Path | Best eval | `op_ftl` | Notes |
|----------|------|------|-----------|----------|-------|
| `optimize_20260604T170329Z` | CMA-ES shell | `runs/grtresna_search/optimize_20260604T170329Z/` | 128 / 57 | 0 | Two leader basins (precursor vs shift); pre-`channel_progress` scores |
| `optimize_20260605T051320Z` | CMA-ES compact warm-start | `runs/grtresna_search/optimize_20260605T051320Z/` | 32 | 0 | Compact + warm start; early run |
| `qd_20260605T060902Z` | MAP-Elites | `runs/grtresna_qd/qd_20260605T060902Z/` | 023 | 0 | Probe; exposed QD logger eval-id bug (fixed in next run) |
| `qd_20260605T062448Z` | MAP-Elites | `runs/grtresna_qd/qd_20260605T062448Z/` | **057** | **1.0** | First operational elite; T5 ladder (`val*hq_qd_eval057`) |
| `qd_20260605T133932Z` | MAP-Elites | `runs/grtresna_qd/qd_20260605T133932Z/` | 058 | 0 | Best precursor (`channel_progress≈0.30`); refine source |
| `qd_20260605T155951Z` | MAP-Elites | `runs/grtresna_qd/qd_20260605T155951Z/` | **011** | **0.75** | **Current production**; 8 operational elites at search scale |
| `l128n256_qd_eval*` | HQ promotion | `runs/grtresna_promote/l128n256_qd_eval*/` | **106** | **1.0** | Batch HQ at t=50, dx=0.5; eval 106 saturates `op_ftl` |
| `optimize_20260605T150015Z` | CMA-ES refine | `runs/grtresna_refine/optimize_20260605T150015Z/` | 012 | 0 | Warm-started from `133932Z`; orbiting precursor basin |
