# grteclyn-wrapper

Run isolated GRTeclyn episodes from the repo root (`GRTeclyn/`). For RadialRecipe GPU smoke tests, the wrapper can stream plotfiles into small data during the run and delete heavy HDF5 dirs afterward.

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

## Read This First: What To Edit, Build, And Run

This repository has three cooperating layers. Most mistakes come from editing the
right idea in the wrong layer, or from launching a long search before rebuilding
the C++ side that changed.

| Goal | Edit here | Rebuild needed? | Validate with |
|------|-----------|-----------------|---------------|
| Change CMA-ES search dimensions, the `ring`/`free` ansatz, warm starts, random/exotic injection, or pre-GPU rejection behavior | `grteclyn-wrapper/src/grteclyn_wrapper/search/optimize.py` and CLI wiring in `grteclyn-wrapper/src/grteclyn_wrapper/__main__.py` | No C++ rebuild; Python code is picked up by the wrapper venv | `cd grteclyn-wrapper && uv run pytest tests/test_grtresna_ring_ansatz.py tests/test_solved_geometry_ftl.py -q` |
| Change launcher defaults or add environment knobs for searches | `grteclyn-wrapper/scripts/search/run_grtresna_search.sh` | No | `DRY_RUN=1 MAX_GENERATIONS=1 GPU_IDS="0 1" bash scripts/search/run_grtresna_search.sh` |
| Change how GRTresna is invoked, configured, parsed, or converted to `.gridinit` | `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/solver.py` and `grteclyn-wrapper/src/grteclyn_wrapper/grtresna/io.py` | Usually no C++ rebuild if only wrapper Python changed | `uv run pytest tests/test_grtresna_integration.py tests/test_grtresna_ring_ansatz.py -q` |
| Change the pre-evolution solved-geometry FTL filter or mechanism descriptors | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/ftl_solved_geometry.py` | No | `uv run pytest tests/test_solved_geometry_ftl.py tests/test_ftl_general.py -q` |
| Change scoring weights/components read from episodes | `grteclyn-wrapper/src/grteclyn_wrapper/metrics/score.py`, `episode_metrics.py` | No | Re-score a campaign or run focused metric tests |
| Change GRTeclyn evolution, fields written to plotfiles, external `.gridinit` loading, matter evolution, or diagnostics | `Examples/RadialRecipe/*`, especially `RadialRecipeLevel.cpp`, `SimulationParameters.hpp`, `ExternalGridInitialData.hpp`; shared matter in `Source/Matter/*` | Yes, rebuild GRTeclyn executable | `BUILD=1 bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh` |
| Change the upstream GRTresna elliptic solver, scalar source, AMR tagging, or maximal-slicing solve | `../GRTresna/Examples/ScalarFieldBH/*` and shared Chombo/GRTresna headers | Yes, rebuild the MPI `Main_ScalarFieldBH3d...mpicxx...MPI.ex` used by production searches (`RANKS=8` default); serial build only if you set `RANKS=1` | MPI solver smoke tests below |

### Launch Cookbook

Run these from `GRTeclyn/grteclyn-wrapper` unless noted otherwise.

```bash
# Refinement search: 14D reduced ring ansatz, half-z by default.
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=ring \
  bash scripts/search/run_grtresna_search.sh

# Discovery: 16D full-sphere shell ansatz. Lumps cover the whole 2-sphere with
# an arbitrary orientation axis and poloidal+toroidal currents (reaches 3D
# configurations the planar ring cannot). Defaults to full-z (no reflective
# z=0 plane) and MPI GRTresna (`RANKS=8`). Use to find new families.
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=shell \
  bash scripts/search/run_grtresna_search.sh

# Broader but much harder 55D search: every lump independent.
GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=free LUMPS=5 \
  bash scripts/search/run_grtresna_search.sh

# Cheap launcher/CLI smoke test: no GRTresna solve and no GPU evolution.
DRY_RUN=1 MAX_GENERATIONS=1 GPU_IDS="0 1" GRTRESNA_ANSATZ=ring \
  bash scripts/search/run_grtresna_search.sh

# Keep raw plotfiles for manual debugging (normally they are consumed/deleted).
NO_CONSUME=1 MAX_GENERATIONS=1 GPU_IDS="0 1" \
  bash scripts/search/run_grtresna_search.sh

# Single-GPU RadialRecipe smoke/build run, useful after C++ edits.
cd /path/to/GRTeclyn
BUILD=1 CUDA_VISIBLE_DEVICES_OVERRIDE=0 \
  bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

Search outputs go to `runs/grtresna_search/optimize_<timestamp>/`. Important
files inside a campaign:

| Path | Meaning |
|------|---------|
| `trajectory.jsonl` | One JSON record per CMA-ES evaluation; rejected candidates are included and still teach CMA-ES through their fitness. Lines are in completion order, not sorted by eval id, and begin with `eval` + `status` for quick scanning. |
| `eval_XXXXXX/metadata.json` | GRTresna convergence, solved-geometry FTL, rejection reasons, simulation exit code. |
| `eval_XXXXXX/score.json` | Parsed metrics and score components after GPU evolution. |
| `eval_XXXXXX/initial_data.gridinit` | Constraint-satisfying GRTresna data loaded by GRTeclyn. |
| `eval_XXXXXX/frames/**/frame_*.png` | Consumed plotfile frames (if the candidate reached GPU and frames were enabled). |
| `eval_XXXXXX/data/*.dat` | In-situ GRTeclyn diagnostics. |

### Current Search Modes

`GRTRESNA_ANSATZ=ring` is the default launcher mode. CMA-ES sees a 14D vector:
amplitude, width, radius, phase, tangential/radial/vertical flow, internal
rotation, exotic fraction/phase, dipole/quadrupole asymmetry, and scalar mode.
The wrapper expands that template into `LUMPS` scalar clouds before calling
GRTresna. Use this mode for efficient refinement of momentum-driven FTL
candidates.

`GRTRESNA_ANSATZ=shell` is the **full-sphere discovery** ansatz. CMA-ES sees a
16D vector and the wrapper expands it into `LUMPS` clouds distributed over the
*entire* 2-sphere (a Fibonacci lattice) rather than on one equatorial circle.
It adds an arbitrary orientation axis `(theta, phi)` and decomposes the matter
current into **toroidal** flow (net angular momentum about the axis -- the warp
"motor") and **poloidal** flow (over-the-pole circulation, a current topology
the planar ring cannot represent). Its current bounds are intentionally more
aggressive than the older diffuse shells: compact, higher-amplitude clouds plus
wider toroidal/poloidal velocities are needed because the realized shift scales
through the scalar momentum source, not directly with the requested velocity.
Use it for *discovery* of new FTL families; use `ring` for cheap *refinement*
once a family is found.

`GRTRESNA_ANSATZ=free` exposes the older independent lump basis:
`11 * LUMPS` searched dimensions (`55D` for `LUMPS=5`). Use it for broader atlas
searches, or as a periodic audit: if `free` random injections keep beating the
`ring`/`shell` elites, the reduced ansatz is leaving value on the table.

All modes end at the same GRTresna data structure (`cfg.lumps`) and the same
GRTeclyn evolution path. The difference is only how the optimizer's vector is
decoded into matter.

> **When to use which.** `ring` (14D, planar) and `shell` (16D, full-sphere) are
> reduced *priors*; they make CMA-ES tractable but bias the search toward their
> respective geometries. `ring` cannot leave the equatorial plane, so it is a
> refinement tool for the momentum-loop family. `shell` covers the full sphere
> in 3D and is the right default for *discovery*. `free` (55D) is unbiased but
> expensive -- keep it as an audit/atlas tool, not the workhorse.

### Stop And Verify Runs

Searches launch GRTresna solver processes, GPU processes, and plotfile
consumers. When stopping a run, kill the whole campaign process tree and verify
both CPU and GPU state:

```bash
# Graceful first pass for all wrapper-started GRTresna searches.
pkill -TERM -f 'runs/grtresna_search'
pkill -TERM -f 'run_grtresna_search.sh'
sleep 5

# Force any stubborn solver / consumer leftovers.
pkill -KILL -f 'runs/grtresna_search'
pkill -KILL -f 'run_grtresna_search.sh'

# Verify no search CPU processes remain.
ps -eo pid,ppid,pgid,pcpu,pmem,etime,args \
  | awk '/run_grtresna_search|grteclyn_wrapper|Main_ScalarFieldBH3d|main3d.gnu.CUDA.ex|consume_plotfiles|prterun|mpirun|orterun/ && !/awk/ {print}'

# Verify GPUs are free.
nvidia-smi
```

The GPU table should show `0MiB` usage and `No running processes found`. A high
load average immediately after killing solver processes can be stale decay; rely
on the process scan and `nvidia-smi`.

### How To Validate A Campaign

Start with the trajectory:

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

For any GPU survivor, inspect the raw metrics in `score.json`, not just the
weighted score components:

- `metrics.general_ftl_solved.f_op`: operational FTL on the GRTresna `.gridinit`
  before GPU evolution.
- `metrics.general_ftl_evolved.f_op`: operational FTL from evolved plotfile
  fields.
- `metrics.general_ftl_evolved.max_local_speed`: local coordinate cone opening.
- `metrics.general_ftl_evolved.max_shift`: maximum sampled shift-vector
  magnitude; this is the direct frame-dragging / light-cone-tilt drive.
- `components.operational_ftl`, `components.ftl_precursor`,
  `components.shift_drive`, and `components.curvature_activity`: whether the
  score is rewarding FTL-like structure or just generic survival/cleanliness.
- `metadata.json.grtresna_convergence`: Ham/Mom residuals from the elliptic solve.

Treat `F_op ~ 1` or `max_local_speed` near the degeneracy ceiling as suspicious
until inspected; mild values like `F_op=0.01..0.05` are more plausible first
survivors and should be replayed at higher resolution before any physics claim.

### Diagnosing "all the top frames look the same"

If several high-scoring runs visually collapse to the same weak blob, treat that
as an optimizer failure mode, not as evidence of a discovered geometry. Inspect
these frames first:

| Frame / metric | What a useful candidate should show | Blob-collapse red flag |
|----------------|-------------------------------------|-------------------------|
| `local_speed_z` | A localized region climbing toward or above `c=1`; eventually a connected faster path. | Entire slice stays sub-luminal and nearly identical across evals. |
| `shift1_z` / shift fields | Shift magnitude growing enough to tilt light cones. | Colorbar is `x10^-2` and max shift is only `~0.03`. |
| `Pi_z` and `scalar_activity_*` | Momentum-carrying, asymmetric, coherent current structure. | Smooth two-lobe scalar cloud with little directional structure. |
| `score.json` / `trajectory.jsonl` | Variation in `shift_drive`, `ftl_precursor`, `curvature_activity`, and ideally `operational_ftl`. | `operational_ftl=0`, `ftl_shortcut=0`, identical `ftl_precursor`, default `mechanism_descriptor=0.5`; score differences come from stability/constraint health. |

The physical reason is simple: GRTresna builds the scalar conjugate momentum as
`Pi = -v . grad(phi)`, and the momentum density is `S_i = -Pi d_i phi`, so the
useful source scales roughly like `amp^2 * velocity / width`. Diffuse
`amp ~ 0.15`, `width ~ 3` shells can converge cleanly but only source a tiny
shift. In that regime CMA-ES sees no FTL gradient and learns to optimize a clean,
still scalar blob.

The current shell search is tuned away from that failure mode:

- `grtresna_shell_amp` and velocity bounds are higher, while
  `grtresna_shell_width` is more compact, so the reachable set includes stronger
  matter currents.
- `metrics/ftl_general.py` reports `max_shift`, and `metrics/score.py` rewards
  it through `shift_drive`.
- `curvature_activity` no longer hard-saturates at `1.0` for ordinary shell
  survivors, so `ftl_precursor` has a real gradient instead of a plateau.
- In `objective_mode=ftl_first`, survival/stability/constraint-health bonuses
  are gated by nontriviality, so "healthy but not warping" candidates cannot win
  by cleanliness alone.

Quick triage command for a campaign:

```bash
cd /path/to/GRTeclyn/grteclyn-wrapper
uv run python3 - <<'PY'
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode

root = Path("../runs/grtresna_search/optimize_<timestamp>")
for episode in sorted(root.glob("eval_*")):
    if not (episode / "score.json").exists():
        continue
    score = score_episode(read_episode_metrics(episode, ftl_L=8.0), target_stop_time=2.0, objective_mode="ftl_first")
    c = score.components
    print(
        episode.name,
        f"score={score.total:.2f}",
        f"op={c.get('operational_ftl', 0):.3f}",
        f"prec={c.get('ftl_precursor', 0):.3f}",
        f"shift={c.get('shift_drive', 0):.3f}",
        f"curv={c.get('curvature_activity', 0):.3f}",
    )
PY
```

### Falsification tiers: "did we actually find the solution?"

A high score is only a proxy; it never proves a candidate is *the* FTL geometry.
The validation-tier layer (`src/grteclyn_wrapper/search/validation_tiers.py`)
answers the question honestly by recording, for every candidate, **how far up an
increasingly demanding falsification ladder it survived**:

| Tier | Name | Gate (needs all lower tiers too) |
|------|------|----------------------------------|
| T0 | `constructed` | constraint-satisfying data exists (not rejected) |
| T1 | `nontrivial` | carries some FTL / non-flat signal (not Minkowski) |
| T2 | `operational` | evolved shortcut beats the flat control, constraints bounded, no trapped surface |
| T3 | `persistent` | survives long evolution: stable, non-growing, channel persists |
| T4 | `observer_ec` | observer-robust energy-condition margin available, exotic cost bounded |
| T5 | `converged` | resolution-ladder replay agrees (external evidence) |
| T6 | `analytic` | closed-form back-derivation reproduced the geometry |

Gates whose diagnostics are absent are reported `unavailable` (not silently
passed), so one short episode climbs only as far as its evidence allows; T5/T6
require promotion runs and are injected via `extra={"resolution_converged":...,
"analytic_form":...}`. The MAP-Elites driver annotates each elite with its tier,
writes a per-iteration `validation.json` (tier histogram + the Pareto
**survivor front** = the candidate solutions), and prints an
**archive-convergence** signal: the search has stopped finding new kinds when
coverage stalls, the best score stalls, and the survivor front is stable across
a window. That is the closest thing to "we are done" the search itself can say.

Assess any existing campaign (CMA-ES or QD) offline, no rerun needed:

```bash
# Writes validation.json into the campaign dir and prints the survivor front.
uv run python scripts/search/validate_tiers.py runs/grtresna_search/optimize_20260602T181607Z
uv run python scripts/search/validate_tiers.py runs/grtresna_search/<campaign> --min-tier 3
```

## Single GPU run (one guessed shape)

Pick **one** initial-data source:

| Env var | Example |
|---------|---------|
| `SEED_NAME` | `ellis_bronnikov` |
| `CANDIDATE_ID` | `bubble_wall_016`, `random_000` |
| `NONSPHERICAL_ID` | `quadrupole_bubble_001`, `dipole_lopsided_000` |

```bash
# Known seed
BUILD=0 SEED_NAME=ellis_bronnikov CUDA_VISIBLE_DEVICES_OVERRIDE=0 \
  bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh

# Spherical guesser candidate
BUILD=0 CANDIDATE_ID=bubble_wall_016 CUDA_VISIBLE_DEVICES_OVERRIDE=1 \
  bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh

# Non-spherical guessed shape
BUILD=0 NONSPHERICAL_ID=quadrupole_bubble_001 CUDA_VISIBLE_DEVICES_OVERRIDE=2 \
  bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

Outputs go to `runs/radialrecipe_gpu_smoke/<name>_gpu_t<stop_time>_<stamp>/`.

## Scripts index

All helper scripts live in [`scripts/`](scripts/README.md). The table below is
the practical launch map; prefer these wrappers over calling Python modules by
hand because they set paths, frame fields, plotfile consumers, and GPU lists
consistently.

| Script | Use for | Key env vars |
|--------|---------|--------------|
| `run_grtresna_search.sh` | **Current matter-first production search**. Each candidate: MPI GRTresna solve (`RANKS=8` default) → `.gridinit` → GRTeclyn GPU evolution for survivors; plotfiles streamed/deleted. | `GRTRESNA_ANSATZ`, `GPU_IDS`, `LUMPS`, `RANKS` (MPI ranks per solve), `MAX_GENERATIONS`, `NO_CONSUME`, `WARM_START_TRAJECTORY` |
| `run_radialrecipe_gpu_smoke.sh` | Single GPU smoke/build run for RadialRecipe seeds/candidates. Use after C++ changes. | `BUILD`, `CUDA_VISIBLE_DEVICES_OVERRIDE`, `SEED_NAME`, `CANDIDATE_ID`, `NONSPHERICAL_ID` |
| `run_ftl_search_cmaes.sh` | Older 9D radial CMA-ES geometry-first search. | `GPU_IDS`, `MAX_GENERATIONS` |
| `run_ftl_search_nonspherical.sh` | Gauge-angular non-spherical search. | `GPU_IDS`, `MAX_GENERATIONS` |
| `run_ftl_search_directional.sh` | Full-z directional search variant. | `GPU_IDS`, `MAX_GENERATIONS` |
| `run_tier2_hq_188.sh` / `run_tier2_validation_long16.sh` | Higher-quality replay/validation of selected candidates. | GPU id / candidate-specific args |
| `run_subset.sh`, `run_nonspherical_gpu_batch.sh`, `run_radialrecipe_gpu_promote.sh` | Batch and promotion helpers for fixed candidate lists. | See script headers |
| `validate_campaign.sh` | Post-run campaign validation. | Campaign path |
| `plot/make_movies.sh`, `plot/plot_run_radial.sh` | Post-processing frames/plots (see [`scripts/plot/`](scripts/README.md#plot--visualization-any-example)). | Episode path / field choices |

```bash
# Example: launch the non-spherical FTL campaign
bash grteclyn-wrapper/scripts/search/run_ftl_search_nonspherical.sh
# Validate the winner at high quality (streaming frames, plotfiles deleted)
bash grteclyn-wrapper/scripts/search/run_tier2_hq_188.sh 0 val16hq_nonsph_eval188
```

## In-situ diagnostics & matter sector

Each RadialRecipe run now emits three diagnostic tables under `data/` (read back
automatically by `grteclyn_wrapper.metrics.read_episode_metrics`):

| File | Columns | Probes |
|------|---------|--------|
| `constraint_norms.dat` | `L2_Ham L2_Mom min_rho_req max_rho_req integral_neg_rho` | constraint satisfaction; `min_rho_req < 0` flags geometries needing exotic matter |
| `energy_conditions.dat` | `matter_min_{NEC,WEC,SEC,DEC} matter_integral_NEC_violation` | observer-sampled energy conditions of the **evolved matter** |
| `curvature_invariants.dat` | `max_abs_ricci_scalar max_ricci_tensor_sq max_Kij_sq L2_ricci_scalar` | coordinate-invariant geometry |

A general, mechanism-agnostic operational-FTL measure (`ftl_general.py`,
Dijkstra shortest-coordinate-time vs. flat baseline) is also computed per episode
and exposed as `EpisodeMetrics.general_ftl`. It is not warp-specific: any
geometry whose coordinate light cones open a faster channel scores `f_op > 0`.

### Spacetime → matter: exotic matter is now evolved when needed

A wormhole/warp geometry generally requires exotic (phantom, `rho <= 0`) matter
to satisfy the Hamiltonian constraint. The constrained recipe (`--phantom`)
already solves for the scalar profile under phantom coupling. The C++ level now
evolves the **matching** matter: when `--phantom` is set the wrapper injects
`recipe_exotic_matter = 1`, and `RadialRecipeLevel` evolves an `ExoticScalarField`
(`T_munu = -recipe_support_strength * canonical`) in the RHS, the constraints,
and the energy-condition diagnostic. With a canonical seed it falls back to the
ordinary `ScalarField`. Verified on Ellis–Bronnikov: the evolved matter
NEC/WEC/SEC/DEC go negative (NEC `-0.07`, integrated violation `~2.1`) exactly
where the geometry demands it, instead of the `~0` null result a canonical field
gives.

### Two findings worth keeping in mind

1. **Matter-sector EC is a null result *only* with a canonical field.** A canonical
   `ScalarField` has `rho >= 0`, so its NEC/WEC are `~0` by construction; `--phantom`
   used to shape only the initial data. With `recipe_exotic_matter = 1` the evolved
   matter is genuinely exotic and the `matter_*` columns reveal the violation. The
   curvature invariants and the general FTL measure read the geometry directly and
   are meaningful regardless. The geometry-sourced effective stress energy
   (`T^eff = G / 8pi`) is what the `matter_*` columns *cannot* see — that is
   evaluated post-hoc from plotfiles by `warpfactory.py`.
2. **Evolved-data FTL now uses plotfile metric fields.** The RadialRecipe
   plotfiles used by the wrapper include the metric, shift, scalar, and derived
   fields needed for the evolved operational-FTL probe. The current diagnostics
   distinguish the initial reconstructed channel from the evolved one
   (`general_ftl` vs `general_ftl_evolved`) so short-lived gauge artifacts are
   visible instead of being hidden in the total score.

## GRTresna initial data & momentum-carrying matter

GRTresna is a Chombo-based **elliptic constraint solver**. Instead of guessing a
metric with the 1D radial recipe and tolerating the constraint violation,
GRTresna *solves* the Hamiltonian **and momentum** constraints in full 3D, then
hands GRTeclyn fully constraint-satisfying initial data. The bridge that wires
the two codes together (Python orchestrator + C++ `.gridinit` reader) is
documented in detail in the package README
([`src/grteclyn_wrapper/README.md`](src/grteclyn_wrapper/README.md) →
*GRTresna integration*); this section covers the **search-pipeline** usage and
the **momentum-carrying matter** capability added on top of it.

### What was done

1. **GRTresna ↔ GRTeclyn bridge:** run GRTresna, convert its Chombo HDF5
   output to a `.gridinit` binary, and load it on the GPU via
   `ExternalGridInitialData` (`recipe_initial_data_file` in `params.txt`).
   Round-trip validation shows roughly two orders of magnitude lower initial
   constraint error than the radial recipe.
2. **Momentum-carrying scalar "clouds":** each localized scalar lump carries
   searchable amplitude, width, 3D center, velocity, rotation, azimuthal mode,
   and exotic flag. The conjugate momentum is constructed so the matter has
   net linear and/or angular momentum, then GRTresna solves the momentum
   constraint for the associated geometry.
3. **Exotic-matter AMR robustness:** exotic lumps switch GRTresna onto the
   hardened maximal-slicing (`K=0`) solver with under-relaxed Newton updates,
   a positive `psi` floor, capped matter Jacobian, arithmetic coefficient
   averaging, and `|rho|` refinement tagging. This removes the old AMR `NaN`
   failure for mixed canonical/exotic multi-lump data.
4. **GRTresna-in-the-loop search:** `--grtresna` mode runs the elliptic solve per
   CMA-ES candidate, converts the result, evolves it in GRTeclyn, streams frames
   and metrics, and feeds a score back to CMA-ES.
5. **Search hardening:** non-converged GRTresna candidates are rejected before
   GPU evolution. The rejection is graded by Ham/Mom residuals, so CMA-ES learns
   that mildly bad initial data is less bad than `NaN` or `Ham=100%`.
6. **Solved-geometry FTL pre-filter:** the *constraint-satisfying* `.gridinit`
   is scored for operational FTL at `t = 0` (before any GPU evolution) and
   candidates with no **physical** FTL signal are rejected. A degeneracy guard
   discards numerical artifacts (see *Pre-evolution solved-geometry FTL filter*
   below).

### Why this matters

- **Less junk radiation, more trustworthy scores.** Evolving non-constraint-
  satisfying data emits spurious "junk" gravitational radiation that contaminates
  stability and growth-rate metrics. GRTresna data starts at L2 constraint
  errors ~1e-3 to 1e-4 (vs ~1e-2 for the 1D recipe), so the dynamics you score
  are physical, not numerical transients.
- **A genuinely new FTL axis.** The 1D recipe cannot represent **moving matter**:
  it has `S_i = 0` (no momentum density). GRTresna solves the momentum constraint
  for `S_i = -Π ∂_i φ ≠ 0`, which is the source of **matter-momentum-driven frame
  dragging** — a distinct, physically-grounded operational-FTL ingredient that
  was simply not in the search space before.
- **3D solve, not a 1D ODE.** Constraint-correct *non-spherical* structure in the
  conformal factor becomes possible, instead of the spherically-symmetric radial
  channel the recipe is limited to.

### Building the GRTresna solver (C++)

The closed-loop search shells out to the compiled `ScalarFieldBH` executable, so
build it once before running. Production searches use the **MPI**
`mpicxx.gfortran` binary (`RANKS=8` by default in `run_grtresna_search.sh`);
`grtresna/solver.py` launches it via `mpirun`/`prterun` from the same
`grtresna` conda/micromamba env. GRTresna is a Chombo app: it needs `CHOMBO_HOME`
plus the Fortran/HDF5/MPI toolchain from that env on `PATH`.

```bash
# Toolchain (mpicxx, gfortran, hdf5, OpenMPI) lives in the grtresna env.
GRTRESNA_ENV=/home/jovyan/.mlspace/envs/grtresna
CHOMBO_HOME=/home/jovyan/nachevsky/test/simulation/Chombo/lib

cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
```

This produces `Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex` in
that directory. The launcher passes `--grtresna-ranks "${RANKS}"` (default `8`);
override with `RANKS=16` etc. if you want more MPI parallelism per solve.

> **Serial fallback (`RANKS=1`).** For debugging only, build with `MPI=FALSE`
> (produces `Main_ScalarFieldBH3d...g++.gfortran...ex` without `.MPI.` in the
> name) and run `RANKS=1 bash scripts/search/run_grtresna_search.sh`. The wrapper
> then launches the serial binary directly with no `mpirun`.

> **Rebuilding after editing headers.** Chombo's makefile keys off `.cpp` files,
> so edits to header-only templated code (e.g. `CTTKHybrid.impl.hpp`,
> `MatterParams.hpp`, `RHSTagging.hpp`) can be reported as "up to date" and *not*
> recompiled. Force a clean relink by removing the stale objects + executable
> first, then rebuild:
>
> ```bash
> rm -f Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex \
>       o/3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI/{Main_ScalarFieldBH,Grids,ScalarField,MyMatterFunctions}.o
> PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
>   make all -j4 CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE
> ```

### Solver-only AMR smoke tests

Before launching a full search you can sanity-check the constraint solver on the
three committed AMR fixtures (`max_level=3`), which converge to finite Ham/Mom
residuals with no `NaN` (the Ham/Mom error file is always written; heavy
per-iteration HDF5 is off by default):

```bash
cd /home/jovyan/nachevsky/test/simulation/GRTresna/Examples/ScalarFieldBH
EXE=Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
for case in canonical exotic mixed_exotic; do
  PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
    mpirun --oversubscribe -np 8 ./"${EXE}" params_${case}_amr_test.txt
done
```

### How to use it — closed-loop search

The recommended launcher is:

```bash
cd /path/to/GRTeclyn/grteclyn-wrapper
# GRTresna: MPI solve with RANKS=8 (default); GRTeclyn: one GPU per population slot.
GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/search/run_grtresna_search.sh
```

Useful overrides:

```bash
# Fewer GPUs/candidates per generation.
GPU_IDS="0 1 2 3" bash scripts/search/run_grtresna_search.sh

# Continue from selected prior good candidates.
WARM_START_TRAJECTORY=/path/to/trajectory.jsonl \
  GPU_IDS="0 1 2 3 4 5 6 7" \
  bash scripts/search/run_grtresna_search.sh

# Keep raw plotfiles for manual re-rendering/debugging.
NO_CONSUME=1 GPU_IDS="0 1" MAX_GENERATIONS=1 \
  bash scripts/search/run_grtresna_search.sh

# Serial GRTresna only (requires MPI=FALSE build above).
RANKS=1 bash scripts/search/run_grtresna_search.sh
```

Default production knobs in the launcher:

| Knob | Default | Meaning |
|------|---------|---------|
| `GRTRESNA_ANSATZ` | `ring` | matter parameterization: `ring` (14D, planar refinement) and `shell` (16D, full-sphere discovery) search global template parameters and expand them to `LUMPS` scalar clouds; `free` (55D) searches every lump independently |
| `GRTRESNA_FULL_Z` | `1` for `shell`, `0` otherwise | full-z GRTresna solve and GRTeclyn evolution box (`center=32 32 32`, no reflective `z=0` plane); shell discovery should normally keep this on |
| `LUMPS` | `5` | scalar clouds generated/evolved by GRTresna (`55` searched dimensions only when `GRTRESNA_ANSATZ=free`) |
| `RANKS` | `8` | MPI ranks per GRTresna elliptic solve (`mpirun -np RANKS` on the `.MPI.ex` binary) |
| `GPU_IDS` | `0 1 2 3` in the script, often overridden to all available GPUs | also sets default CMA-ES population |
| `MAX_GENERATIONS` | `50` | CMA-ES generations |
| `OBJECTIVE_MODE` | `ftl_first` | FTL terms dominate health/stability terms |
| `RANDOM_INJECTION_FRACTION` | `0.25` | per-generation random candidates |
| `EXOTIC_INJECTION_FRACTION` | `0.25` | forced exotic-pattern candidates |
| `WARM_START_TRAJECTORY` | unset | comma-separated prior `trajectory.jsonl` files |
| `PROJECTION_FIELDS` | `scalar_activity` | 3D placement diagnostics |
| `PROJECTION_AXES` | `x y z` | max-intensity projections through the 3D cloud |

Per candidate the loop runs:

```
CMA-ES proposes grtresna_lump{k}_*
  → GRTresna 3D elliptic solve under MPI (Hamiltonian + momentum constraints)
  → reject candidate if Ham/Mom convergence is missing, NaN, or above threshold
  → initial_data.gridinit (solved A_ij + phi + Pi)
  → GRTeclyn loads it and evolves on GPU
  → consume plotfiles into metrics, slice frames, and 3D projections
  → ftl_first score → back to CMA-ES
```

> **MPI vs. serial.** Default `RANKS=8`: the wrapper selects
> `Main_ScalarFieldBH3d...mpicxx...MPI.ex` and runs
> `mpirun --oversubscribe -np RANKS ...` (or `prterun` from `$GRTRESNA_ENV`).
> `mpirun` is resolved from `$GRTRESNA_MPIRUN`, then `$GRTRESNA_ENV` /
> `$CONDA_PREFIX`, env roots for `$GRTRESNA_ENV_NAME` (default `grtresna`),
> then `PATH` — see `grtresna/solver.py::_resolve_mpirun`. With `RANKS=1` and a
> serial build present, the non-MPI `g++.gfortran` executable is launched directly.

`--grtresna` **replaces** the radial-recipe search space with a scalar-matter
basis (it is a separate mode from `--nonspherical`). There are two
parameterizations:

- `--grtresna-ansatz ring` (launcher default via `GRTRESNA_ANSATZ=ring`): CMA-ES
  searches **14 global parameters** — ring amplitude/width/radius, phase,
  tangential/radial/vertical flow, rotation, exotic fraction/phase, dipole and
  quadrupole asymmetry, and mode — then deterministically expands them into
  `K` scalar lumps. This keeps the physically useful structure (counterflow,
  angular momentum, exotic support, shear) while reducing the optimizer from
  55D to 14D for `K=5`. **Limitation:** the lumps live on one equatorial circle
  (plus a small `z` wobble), so the ring cannot place matter or currents
  anywhere on the 2-sphere. Treat it as a *refinement* prior for the
  momentum-loop family, not a discovery tool.
- `--grtresna-ansatz shell` (`GRTRESNA_ANSATZ=shell`): the **full-sphere
  discovery** ansatz. CMA-ES searches **16 global parameters** — amplitude,
  width, radius, shell thickness, orientation axis `(theta, phi)`,
  toroidal/poloidal/radial flow, internal rotation, dipole/quadrupole asymmetry
  (Legendre orders along the axis), exotic fraction/phase, mode, and a radial
  jitter — then expands them into `K` lumps spread over the *whole* sphere via a
  Fibonacci lattice. Toroidal flow circulates about the axis (net `L`,
  gravitomagnetic dipole); poloidal flow goes over the poles. This reaches 3D
  configurations the planar ring structurally cannot, at only +2 dimensions.
  The default shell ranges sit between the old diffuse blobs and the too-strong
  compact regime (`amp=0.08..0.28`, `width=1.8..4.0`, toroidal velocity `-2..2`)
  because GRTresna must converge before GPU evolution can score shift/FTL. For
  this mode the launcher defaults to `GRTRESNA_FULL_Z=1`,
- `--grtresna-ansatz free`: the older unconstrained `K`-lump basis. Each lump
  `k` contributes the 11 dimensions below, so `K=5` gives 55 searched
  dimensions. This is maximally expressive but harder for CMA-ES to learn; use
  it as a periodic audit of the reduced ansätze.

Both modes produce the same downstream `cfg.lumps` representation; GRTresna
still solves the full 3D constraints, writes `.gridinit`, and GRTeclyn evolves
the result. Use `--dry-run` to exercise the plumbing without a solve or GPU.

#### Why the `ring` ansatz reduces the search space

The original `free` GRTresna search gives CMA-ES every lump coordinate
independently: `amp`, `width`, `center_{x,y,z}`, `velocity_{x,y,z}`, `omega`,
`mode`, and `exotic`. With `LUMPS=5` that is `5 * 11 = 55` searched dimensions.
It can represent almost any compact scalar-matter cloud, but many dimensions are
symmetry/noise for the physics we are trying to learn: rotating the whole cloud,
permuting equivalent lumps, or nudging one lump independently often describes
nearly the same matter-current pattern while costing CMA-ES extra samples.

The `ring` ansatz searches coherent matter-flow parameters instead, then expands
them into the same five lump dictionaries before calling GRTresna. The searched
14D vector is:

| Ring parameter | Physical role |
|----------------|---------------|
| `grtresna_ring_amp` | common scalar-field amplitude scale |
| `grtresna_ring_width` | common Gaussian cloud width |
| `grtresna_ring_radius` | radius of the generated matter ring |
| `grtresna_ring_z_scale` | sinusoidal out-of-plane thickness/tilt |
| `grtresna_ring_phase` | global rotation angle of the ring |
| `grtresna_ring_tangential_velocity` | circulating matter current / angular momentum source |
| `grtresna_ring_radial_velocity` | inward/outward compression flow |
| `grtresna_ring_vertical_velocity` | alternating vertical motion |
| `grtresna_ring_omega` | shared internal lump rotation |
| `grtresna_ring_exotic_fraction` | fraction of generated lumps that are phantom/exotic |
| `grtresna_ring_exotic_phase` | where the exotic sector sits on the ring |
| `grtresna_ring_dipole_amp` | one-sided asymmetry |
| `grtresna_ring_quadrupole_amp` | two-lobed squeeze/shear asymmetry |
| `grtresna_ring_mode` | shared azimuthal lump mode |

Internally, lump `k` is placed at angle
`theta_k = phase + 2*pi*k/LUMPS`. Its center is on a ring of radius
`grtresna_ring_radius`, with optional sinusoidal `z` displacement; its velocity
is the sum of radial, tangential, and vertical components. Dipole/quadrupole
terms modulate amplitude and width around the ring. The exotic fraction and
phase choose a contiguous sector of the ring to flip to phantom sign. In other
words, CMA-ES searches a **rotating/shearing scalar-current loop** with
localized exotic support, rather than five unrelated blobs.

Use `GRTRESNA_ANSATZ=ring` when you want efficient discovery/refinement of
momentum-driven FTL candidates. Use `GRTRESNA_ANSATZ=free` when you deliberately
want the broader 55D atlas search and can afford many more evaluations.

| Search dimension (per lump `k`) | Meaning |
|---------------------------------|---------|
| `grtresna_lump{k}_amp` | amplitude of the cloud (0 ⇒ disabled) |
| `grtresna_lump{k}_width` | Gaussian width |
| `grtresna_lump{k}_center_{x,y,z}` | cloud position ⇒ where mass/momentum sits |
| `grtresna_lump{k}_velocity_{x,y,z}` | boost velocity ⇒ net **linear** momentum `P_i ~ v_i` |
| `grtresna_lump{k}_omega` | rigid rotation rate about z ⇒ net **angular** momentum `L_z` |
| `grtresna_lump{k}_mode` | azimuthal mode `m` (rounded to int; `m ≥ 1` required for `L_z`) |
| `grtresna_lump{k}_exotic` | phantom/exotic flag (rounded to 0/1) |

`--grtresna-lumps K` sets how many scalar clouds are emitted into GRTresna. In
`ring` mode it changes the generated polygon/ring resolution but not the CMA-ES
dimension count; in `free` mode it changes the optimization dimension
(`11 * K`). Defaults focus on **matter momentum**: black holes are off
(`bh1_bare_mass = 0`) in the search config, so the optimizer explores pure
moving/rotating scalar clouds. Each candidate's heavy HDF5 is deleted right
after conversion to `.gridinit` (`cleanup=True`), so a long conveyor search
stays disk-safe.

The default `free` search bounds are deliberately compact so clouds do not smear
across the whole frame: centers stay near the physical center and width is
bounded. The reduced `shell` ansatz is now compact for a different reason: it
needs enough `amp^2 * velocity / width` to source a visible shift before the
score can climb toward FTL.

For **finding new geometry families** (rather than one optimum), the MAP-Elites
`qd` driver over this same matter basis is the better tool — it keeps a diverse
archive across behavior descriptors instead of collapsing to a single best.

### How to use it — one-off, from Python

```python
from pathlib import Path
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve

cfg = GRTresnaConfig(
    mpi_ranks=1,
    bh1_bare_mass=0.0,            # pure matter-momentum case
    lump_amp=0.1, lump_width=8.0,
    lump_velocity=(0.2, 0.0, 0.0),  # boosted cloud → linear momentum
    # or, for angular momentum:  lump_omega=0.1, lump_mode=1
)
gridinit = solve(cfg, work_dir=Path("/tmp/grtresna_run"))
# then run GRTeclyn with  recipe_initial_data_file = <gridinit>
```

A fast standalone smoke test of the momentum solve lives at
`GRTresna/Examples/ScalarFieldBH/params_momentum_test.txt` (a boosted lump: the
momentum constraint starts at ~98% violation and is solved down to ~0.1%).

### Will GRTeclyn evolve it?

Yes. GRTeclyn's scalar-field matter computes the same momentum density
(`j_i = -Π ∂_i φ`) and evolves the full CCZ4 + Klein–Gordon system; the loaded
solved `A_ij` plus the Γ-driver shift gauge develop the frame-dragging shift
dynamically. No GRTeclyn evolution code changes were needed — only the
`ExternalGridInitialData` loader (already in place) and the upstream solver
profiles.

### Exotic (phantom) matter and the maximal-slicing solver

A genuine FTL channel generally needs **exotic matter** (negative energy density,
`rho < 0`, an NEC violation). We added a per-lump *phantom* capability to the
GRTresna `ScalarFieldBH` example and made the constraint solver able to handle
it robustly in the AMR multi-lump search regime.

**What was added.**

1. **Per-lump exotic flag** (`lump{k}_exotic`). `ScalarField::compute_emtensor`
   uses an independent-field model: each lump contributes its kinetic-energy and
   momentum density with a sign set by its flag (`+1` canonical, `-1` phantom),
   so a configuration can freely mix normal and exotic lumps. The search space
   gained a per-lump `grtresna_lump{k}_exotic` dimension (rounded to a 0/1 flag).
2. **Maximal-slicing (K=0) solve path** (`maximal_slicing` param). The default
   CTTK(Hybrid) method sets `K = sign*sqrt(24 pi G rho + ...)`, which is the
   square root of a **negative** number wherever `rho < 0` — so exotic matter
   produced `NaN` immediately. This is structural, not a tuning bug. The new
   path instead fixes `K = 0` and moves the matter energy into the elliptic
   ψ-solve as a source (a standard York/Lichnerowicz CMC solve:
   `rhs += -2 pi G rho psi^5`, `aCoef += 10 pi G rho psi^4`), which handles
   either sign of `rho`.
3. **Nonlinear-solve robustness** (`psi_relaxation`, `psi_floor`). The K=0 matter
   source makes the linear operator indefinite for `rho < 0`; under-relaxation of
   the Newton step plus a ψ-positivity floor keep it from overshooting into
   `psi <= 0` (where `psi^-7` is `NaN`).
4. **AMR-stable indefinite operator** (`maximal_jacobian_cap`, arithmetic
   coefficient averaging, `|rho|` tagging). On AMR the indefinite K=0 operator
   used to diverge to `NaN` at the coarse-fine boundaries. Three changes fix it:
   (a) the Newton Jacobian contribution from the matter source is **capped**
   (`maximal_jacobian_cap`) and `psi_0` is floored finite/positive in
   `CTTKHybrid` so a single bad cell can no longer poison the multigrid V-cycle;
   (b) coefficient coarsening switches from **harmonic to arithmetic** averaging
   (`coefficient_average_type`), which stays well-defined when the sign-changing
   `aCoef` straddles zero across a coarse-fine face; (c) refinement tagging keys
   on `|rho|` (`RHSTagging`) so negative-energy regions are actually refined
   instead of being cancelled by neighbouring positive matter.
5. **Continuous pre-FTL shaping scores** (`ftl_precursor` and `shift_drive`, in
   `metrics/score.py`). `operational_ftl` is a hard end-to-end-channel gate that
   is flat-zero until a full superluminal channel forms, giving CMA-ES no
   gradient. `ftl_precursor` rewards structured light-cone opening as
   `max_local_speed` climbs toward and beyond `1`; `shift_drive` rewards the
   realized shift magnitude directly, which is the mechanism the scalar-current
   shell is trying to source. Health/stability terms are gated by nontriviality
   in `ftl_first`, so clean static blobs do not dominate.

**What works (validated under AMR).**

- Canonical matter on the K=0 path matches the standard path — the formulation
  is correct.
- A single **isolated** exotic lump converges to clean, finite, constraint-
  satisfying data, where the old solver only produced `NaN`.
- The maximal-slicing path is now **AMR-robust**. Three solver-only smoke
  fixtures at `max_level=3` (in
  `GRTresna/Examples/ScalarFieldBH/params_*_amr_test.txt`) converge cleanly with
  no `NaN`:
  - `params_canonical_amr_test.txt` (positive-energy control, 2 boosted lumps):
    `Converged!` at Ham `0.0056%`, Mom `0.025%`.
  - `params_exotic_amr_test.txt` (one isolated exotic lump, K=0 path):
    `Converged!` at Ham `0.099%`, Mom `0.027%`.
  - `params_mixed_exotic_amr_test.txt` (the production failure mode: 3 lumps,
    mixed canonical/exotic, angular `mode=1`): `Converged!` at Ham `0.058%`,
    Mom `0.029%`.
- There is still a sharp **physical** ceiling at high negative-energy density —
  the Lichnerowicz/York **existence boundary**, beyond which no asymptotically-
  flat constraint-satisfying data exists. `EXOTIC_AMP_SCALE = 0.25` maps the
  search amplitude range into the convergent regime so exotic candidates stay
  solvable.

**How the search uses it.** `build_grtresna_config` auto-detects exotic lumps and
switches that candidate onto the hardened K=0 path (maximal slicing,
`psi_relaxation=0.6`, `psi_floor=0.1`, `maximal_jacobian_cap=25`, arithmetic
coefficient averaging). **Purely-canonical** candidates are left on the standard
CTTK ansatz, which remains the more robust choice for positive multi-lump matter.
The lump-basis search also zeros the legacy radial background scalar
(`dphi=dpi=0`) so the matter is exactly the searched lumps.

**GRTresna quality gate.** A solve that returns missing, non-finite, or large
Ham/Mom residuals is not safe to evolve. The wrapper therefore rejects such
candidates immediately after the GRTresna solve, before launching GRTeclyn. The
fitness penalty is graded by residual size, so CMA-ES can distinguish mildly bad
initial data from `NaN` or `Ham=100%` junk and move away from those regions.

### Pre-evolution solved-geometry FTL filter

**The problem.** The matter-first loop was effectively *blind to FTL* until a
full, expensive GPU evolution finished: the matter ansatz is proposed, GRTresna
solves the constraints, and only the *evolved* spacetime was ever scored for an
operational shortcut. A converged-but-uninteresting solve (a smooth scalar blob
that produces no channel) cost a whole GPU evolution to discover, so the search
wasted most of its budget on candidates that were never going to show FTL. The
quality gate above only checks *constraint residuals*, not whether the resulting
geometry has any FTL potential at all.

**How we solved it.** We now score the constraint-satisfying `.gridinit`
**at `t = 0`**, before the GPU step, using the same mechanism-agnostic Dijkstra
operational-FTL probe used on evolved data
(`metrics/ftl_solved_geometry.py::compute_solved_geometry_ftl`). The result is a
cheap (~1 s/candidate, vectorised x–z slice extraction) signal that the gate
uses to reject candidates with no FTL potential *before* evolving them. The
decision policy lives in `search/solved_ftl_gate.py`; `search/optimize.py` only
applies that policy inside `_objective`. The rejection fitness is graded by how
far the slice is from a signal, so CMA-ES still gets a gradient toward
FTL-promising matter. This gate is **on by default** in
`--grtresna` mode (no flag required); it is what would have pruned the smooth
`eval_000067` blob without spending a GPU evolution on it.

**Gate policy is configurable.** The solved-FTL gate thresholds live in one
policy module (`search/solved_ftl_gate.py`, `SolvedFtlGateConfig`) and are
exposed by the `optimize` CLI plus `scripts/run_grtresna_search.sh`; do **not**
edit `metrics/ftl_solved_geometry.py` just to make a run stricter or more
exploratory. The launcher knobs are:
`SOLVED_FTL_F_OP_FLOOR`, `SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR`,
`SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR`,
`SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR`,
`SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED`, `SOLVED_FTL_MAX_PHYSICAL_F_OP`, and
`SOLVED_FTL_REJECTION_SPEED_TARGET`.

**Exploration threshold.** In production discovery runs the default gate is
deliberately permissive but not flat-space accepting: a clean, non-degenerate
slice passes if it has `F_op > 1e-4`, if it is a subluminal near miss
(`0.99 <= max_local_speed < 1`), or if it has enough superluminal-area precursor.
This lets near-threshold ring candidates reach a short GPU evolution, where the
Gamma-driver shift can amplify a matter-momentum precursor into an evolved
channel. Exactly-flat data (`max_local_speed ~= 1`, `F_op = 0`, no superluminal
area) is still rejected. Candidates far below threshold receive a graded
rejection fitness instead of spending GPU time.

For exploratory `shell` campaigns that otherwise reject nearly everything before
GPU evolution, temporarily lower the near-luminal floor:

```bash
SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR=0.9 \
  GPU_IDS="0 1 2 3 4 5 6 7" GRTRESNA_ANSATZ=shell \
  bash scripts/search/run_grtresna_search.sh
```

This is not a physics claim; it is an exploration setting that lets sub-luminal
but structured candidates reach the GPU long enough for `shift_drive` and the
evolved FTL metrics to sort them.

**Degeneracy guard (the subtle part).** Solves near the Lichnerowicz/York
**existence boundary** (heavy exotic matter) produce isolated near-degenerate
metric cells where `gamma -> 0`. There the coordinate light speed blows up and
the Dijkstra path finds a near-zero-cost edge, so a *numerical artifact* looks
like a spectacular shortcut (`max_local_speed ~ 100`, `F_op ~ 0.99`,
i.e. crossing the whole domain in ~1 % of flat-light time). Letting these pass
would defeat the filter — it would keep garbage and starve the genuinely mild,
physical channels. We therefore flag a slice as `degenerate` and refuse to count
it as a signal when it exceeds configurable physical-plausibility ceilings
(defaults: `SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED=8`,
`SOLVED_FTL_MAX_PHYSICAL_F_OP=0.85`) or is non-finite. A
real constraint-satisfying shortcut on a compact box is mild (`max_c` of order a
few, `F_op` of a few tenths), like `eval_000063` (`F_op = 0.13`,
`max_c = 1.16`), which the filter keeps.

**Mechanism descriptor (for diversity).** Each surviving slice is also tagged
with a continuous `[0, 1]` mechanism axis — shift-warp (`0`), throat-pinch
(`0.5`), portal-compression (`1`) — blended by which proxy is strongest rather
than a hard `argmax` (which collapsed every candidate onto one value). MAP-Elites
(`qd_search.py`) uses this as a behavior descriptor so the archive spreads across
*distinct FTL families* instead of clustering on a single mechanism.

**Validation (offline rescore).** Re-scoring the prior
`runs/grtresna_search/optimize_20260602T181607Z` campaign with the filter cut
the survivors from **20 → 8** of 80 solved candidates: the ~12 dropped were all
degenerate `max_c=100` / `F_op≈0.99` artifacts, while the physically-plausible
candidates (including `eval_000063`) were kept. Replay the rescore with
`scripts/search/rescore_grtresna_solved_ftl.py <campaign_dir>`.

**Net status and current finding.** Exotic matter is implemented end-to-end and
is solvable **inside the full AMR multi-lump regime** the search actually
proposes, validated by the three smoke fixtures above plus
`tests/test_grtresna_integration.py`. A clean gated 8-GPU search
(`runs/grtresna_search/optimize_20260602T181607Z`) produced 54 successful
evolutions and rejected 26 bad GRTresna candidates before GPU evolution. Its best
candidate, `eval_000063`, used three exotic lumps out of five and had clean
GRTresna convergence (`Ham=0.1676%`, `Mom=0.0048%`). It did **not** show an
initial/static shortcut (`ftl_shortcut=0`), but after evolution it produced a
coordinate superluminal channel (`F_op^ev=0.0461`, max local coordinate speed
`1.0888`). This is a promising evolved FTL precursor, not yet a validated
physical warp drive: the next step is deterministic replay, longer evolution,
and resolution convergence.

## Batch: 7 non-spherical shapes on GPUs 0–6

```bash
BUILD=0 bash grteclyn-wrapper/scripts/radial/run_nonspherical_gpu_batch.sh
```

Outputs: `runs/radialrecipe_nonspherical/`. One log per candidate: `<id>_<stamp>.log`.

## Plotfile consumer (default on)

With `CONSUME_PLOTFILES=1` (default), each run:

1. Starts a sidecar `consume_plotfiles` process while the GPU simulation runs.
2. Appends `small_data/shell_profiles.dat`, `small_data/areal_radius.dat`, optional PNG frames.
3. Deletes processed plotfile dirs (`--keep-last 1`).
4. Runs a **post-sim drain** for any backlog.

Useful env vars:

| Variable | Default | Meaning |
|----------|---------|---------|
| `CONSUME_PLOTFILES` | `1` | Enable streaming extraction |
| `CONSUMER_DELETE` | `1` | Delete HDF5 plot dirs after extract |
| `CONSUMER_RADII` | `4 8` | Extraction radii |
| `PLOT_INTERVAL` | `10` if consumer on, else `1` | Plotfile write cadence |
| `STOP_TIME` | `2.0` | Simulation stop time |
| `N_FULL` | `64` | Grid resolution |
| `CUDA_VISIBLE_DEVICES_OVERRIDE` | `0` | GPU index for single run |

Disable consumer (keep all plotfiles):

```bash
CONSUME_PLOTFILES=0 bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh
```

## Post-run plots (no GW / Psi4)

```bash
bash grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh runs/radialrecipe_nonspherical/<episode_dir>
```

Writes EPS/PNG to `grteclyn-wrapper/src/grteclyn_wrapper/visualisation/plots/radial/` (constraints, collapse, shell profiles).

## Manual plotfile drain (if needed)

If a run finished before extraction completed:

```bash
bash grteclyn-wrapper/scripts/plot/plot_run_radial.sh runs/radialrecipe_nonspherical/<episode_dir> --no-delete
# or one-shot batch consume (no watch):
uv run python -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data runs/radialrecipe_nonspherical/<episode_dir> \
  --out runs/radialrecipe_nonspherical/<episode_dir>/small_data \
  --radii 4 8 --no-psi4 --shell-fields chi lapse K --areal-radius \
  --delete --keep-last 1 -j 4
```
