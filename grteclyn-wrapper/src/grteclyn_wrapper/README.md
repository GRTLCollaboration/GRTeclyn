# `grteclyn_wrapper` — Python package reference

Closed-loop search engine for exotic spacetime metrics on top of the
`GRTeclyn` C++/GPU numerical-relativity code. The package proposes initial-data
parameter vectors, runs constrained GPU evolutions, scores the result against
FTL / energy-condition / stability metrics, and drives the search with CMA-ES,
a quality-diversity archive, or a surrogate-assisted loop.

This README documents the **Python modules and the CLI**. For the GPU shell
scripts (`run_radialrecipe_gpu_smoke.sh`, `run_nonspherical_gpu_batch.sh`, …)
see [`../../README.md`](../../README.md).

---

## CLI quick reference

Run as a module from the repo root (`GRTeclyn/`):

```bash
uv run python -m grteclyn_wrapper [GLOBAL OPTIONS] <command> [COMMAND OPTIONS]
```

| Command | Purpose |
|---------|---------|
| `reproduce` | Run one episode from a template + overrides, then score it. |
| `sweep` | Random sweep over wormhole parameters. |
| `atlas` | Low-resolution failure-atlas batch. |
| `optimize` | CMA-ES search over RadialRecipe coefficients (multi-GPU). |
| `qd` | **MAP-Elites quality-diversity** search (Spacetime Failure Atlas). |
| `pareto` | **Multi-objective Pareto-front** extraction from a trajectory. |
| `warpfactory` | **Multi-observer energy conditions** (NEC/WEC/SEC/DEC) of an analytic 4-metric. |
| `validate` | Batch-validate the metric guesser on synthetic candidates. |

Common global options: `--example {RadialRecipe,SupportedWormholeCollapse}`,
`--constrained`, `--phantom`, `--preflight`, `--cuda-devices`, `--runs-dir`,
`--set KEY=VALUE` (params override), `--score-weight KEY=VALUE`, `--ftl-L`.

> **Gotcha:** do **not** `source scripts/lib/env.sh` before `optimize`/`qd`. It sets
> `PYTHONPATH` that shadows `uv`'s resolution of the `cma` package. Single-GPU
> per-process runs use the non-MPI binary and don't need `env.sh`.

### Examples

```bash
# CMA-ES on 8 GPUs with surrogate pre-screening
uv run python -m grteclyn_wrapper --example RadialRecipe \
    --constrained --phantom --preflight \
    optimize --max-generations 50 --population-size 8 \
    --gpu-ids 0 1 2 3 4 5 6 7 --surrogate --surrogate-keep-fraction 0.5

# MAP-Elites quality-diversity campaign (behavior atlas)
uv run python -m grteclyn_wrapper --example RadialRecipe \
    --constrained --phantom --preflight \
    qd --iterations 20 --batch-size 8 --bins 8 --gpu-ids 0 1 2 3 4 5 6 7

# Extract the Pareto front from a finished optimizer run
uv run python -m grteclyn_wrapper pareto \
    --trajectory runs/optimize_.../trajectory.jsonl --output front.json

# Warp Factory-style energy-condition report for an Alcubierre bubble
uv run python -m grteclyn_wrapper warpfactory \
    --metric alcubierre --velocity 0.5 --n-directions 60 --n-speeds 4
```

---

## Module map

Package layout (each area lives in its own subfolder under
`grteclyn-wrapper/src/grteclyn_wrapper/`):

```
grteclyn_wrapper/
  __init__.py
  __main__.py          CLI entry point
  core/
    config.py          repo paths, example/executable resolution
    episode.py         per-run directory layout and metadata
    params.py          render params.txt from template + overrides
    runner.py          launch GRTeclyn binary; plotfile consumer
    evaluation.py      shared candidate → episode → score helper
    plot_consumer.py   sidecar plotfile extraction command builder
  initial_data/
    constrained_recipe.py
    preflight.py
    seeds.py
    candidates.py
    validate_guesser.py
    nonspherical_guesser.py
  metrics/
    episode_metrics.py   parse diagnostics → EpisodeMetrics
    ftl_metrics.py       t=0 FTL shortcut metrics
    ftl_general.py       mechanism-agnostic operational FTL
    physical_metrics.py  ANEC / tidal proxies
    score.py             weighted scalar reward
    warpfactory.py       analytic 4-metric energy conditions
  search/
    optimize.py          CMA-ES driver
    surrogate.py         RBF surrogate screening
    qd_search.py         MAP-Elites quality-diversity
    pareto.py            Pareto-front extraction
    atlas.py             random failure-atlas batch
  grtresna/
    io.py                Chombo HDF5 → .gridinit
    solver.py            GRTresna orchestrator
```

### Import paths

After the refactor, import from the subpackage that owns the module. A few
symbols are re-exported from package `__init__` files for convenience.

| What you need | Import |
|---------------|--------|
| Run/score an episode | `from grteclyn_wrapper.metrics import read_episode_metrics, score_episode` |
| Repo / example config | `from grteclyn_wrapper.core.config import REPO_ROOT, resolve_example` |
| Known seeds | `from grteclyn_wrapper.initial_data.seeds import get_seed` |
| Constrained recipe | `from grteclyn_wrapper.initial_data.constrained_recipe import constrained_overrides` |
| CMA-ES driver | `from grteclyn_wrapper.search.optimize import run_optimize` |
| Warp Factory EC report | `from grteclyn_wrapper.metrics import warpfactory` or `from grteclyn_wrapper import warpfactory` |
| GRTresna bridge | `from grteclyn_wrapper.grtresna import GRTresnaConfig, solve, convert_chombo_to_gridinit` |

Example:

```python
from pathlib import Path

from grteclyn_wrapper.core.config import resolve_example
from grteclyn_wrapper.initial_data.seeds import get_seed
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.search.optimize import run_optimize
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve
```

### Pipeline core (`core/`)
| Module | Role |
|--------|------|
| `__main__.py` | CLI: argument parsing and command dispatch. |
| `core/config.py` | Example/executable resolution, repo paths, default run dir. |
| `core/episode.py` | Per-run directory layout (`Episode`), metadata, JSON writers. |
| `core/params.py` | Render a `params.txt` from a template + overrides. |
| `core/runner.py` | Launch the GRTeclyn binary; optional plotfile consumer. |
| `core/evaluation.py` | **Shared candidate→episode→score helper** used by all drivers. |
| `core/plot_consumer.py` | Build the sidecar `consume_plotfiles` command for streaming extraction. |

### Initial data (`initial_data/`)
| Module | Role |
|--------|------|
| `initial_data/constrained_recipe.py` | Derive `phi` from `chi` via the Hamiltonian constraint; Gaussian basis; Ricci scalar. |
| `initial_data/preflight.py` | Cheap 1D constraint filter; reject bad candidates before the GPU. |
| `initial_data/seeds.py` | Known-solution seeds: flat, Ellis–Bronnikov, Schwarzschild puncture, Alcubierre warp. |
| `initial_data/candidates.py` | Resolve initial-data overrides from seed / candidate / non-spherical IDs. |
| `initial_data/validate_guesser.py` | Batch validation of the metric guesser on synthetic proposals. |
| `initial_data/nonspherical_guesser.py` | Angular-mode non-spherical metric proposals and ray validation. |

### Metrics & scoring (`metrics/`)
| Module | Role |
|--------|------|
| `metrics/episode_metrics.py` | Parse diagnostics into `EpisodeMetrics`; collapse/constraint/stability/comoving + growth-rate (`GrowthMetrics`). |
| `metrics/ftl_metrics.py` | t=0 FTL shortcut metrics: `F_null`, `F_portal`, `F_throat`, `F_asymmetry`, `F_log`, `s_nonflat`. |
| `metrics/ftl_general.py` | Mechanism-agnostic operational FTL on reconstructed or plotfile slices. |
| `metrics/physical_metrics.py` | t=0 gauge-robust proxies: ANEC line integral, curvature/tidal proxy, trapped-surface flag. |
| `metrics/warpfactory.py` | Warp Factory port: 4-metric → Einstein tensor → multi-observer NEC/WEC/SEC/DEC. |
| `metrics/score.py` | Weighted multi-component scalar reward `score_episode()`. |

The `metrics` package re-exports `read_episode_metrics`, `read_growth_metrics`,
`dataclass_to_dict`, and `score_episode` from its `__init__.py`.

### Search drivers (`search/`)
| Module | Role |
|--------|------|
| `search/optimize.py` | CMA-ES loop, multi-GPU parallel generations; optional surrogate screening. |
| `search/surrogate.py` | numpy RBF kernel-ridge surrogate + candidate screening. |
| `search/qd_search.py` | MAP-Elites quality-diversity driver and archive. |
| `search/pareto.py` | non-dominated sorting / Pareto-front extraction. |
| `search/atlas.py` | Random failure-atlas batch runner. |

### GRTresna bridge (`grtresna/`)
| Module | Role |
|--------|------|
| `grtresna/io.py` | Read Chombo AMR checkpoint HDF5 (with ghost cells), flatten to uniform grid, write `.gridinit` binary. |
| `grtresna/solver.py` | Orchestrator: write GRTresna `params.txt`, run via MPI, convert output. `GRTresnaConfig` holds NL early-exit (`nl_exit_tolerance`, `nl_stall_tolerance`) and `gridinit_workers` for parallel conversion. |

---

## New in this work (`feature/interstellar`)

These implement the paper's "Proposed Extensions" (Sec. *Extensions: Implemented
Metrics, Search, and Robustness*). All are numpy-only (no new dependencies) and
GPU-validated; unit tests live in `tests/test_proposed_extensions.py`.

### New physical metrics
- **Constraint growth rate** (`metrics.episode_metrics.GrowthMetrics`, scored as
  `constraint_growth`): log-linear fit of the exponential rate `λ` on the
  `‖H‖`, `|K|_max`, and `1/χ_min` time series. `s_growth = 1/(1 + max(0,λ)/σ_λ)`
  penalizes geometries that merely collapse slowly enough to survive a short
  run — closing the dynamical-stability gap without longer evolutions.
- **ANEC line proxy** (`metrics.physical_metrics`, scored as `anec_condition`): the
  required energy density integrated along the travel axis,
  `A_line = ∫ ρ_req(x) dx`; a directional energy-condition indicator. *t=0,
  χ-sourced — blind to pure-shift warps; the full geodesic ANEC needs the
  evolved probe (future).*
- **Tidal / curvature proxy** (`metrics.physical_metrics`, scored as `tidal_comfort`):
  `max|R| + max|∂²α|` with a companion t=0 trapped-surface flag. Precursor to
  the full Kretschmann / electric-Weyl invariants (future C++ work).
- **Multi-observer energy conditions** (`metrics.warpfactory`, CLI `warpfactory`): a
  numpy port of *Warp Factory* (Helmerich et al., arXiv:2404.03095) with the
  *warpax* refinements (Le, arXiv:2602.18023). Builds the Einstein tensor of a
  full 4-metric by **fourth-order** finite differences, forms `T_{μν}`, then
  verifies the pointwise **NEC/WEC/SEC/DEC** via two pathways:
  (a) **observer sampling** — contract `T` with a sphere of null/timelike
  observers (lower bound); and (b) the **Hawking–Ellis eigenvalue test** — at
  Type I points the conditions reduce to *exact, observer-independent*
  inequalities in the eigenvalues `(-ρ, p_i)` of `T^a_b`. The score uses the
  more-violating "hybrid margin" of the two, so a violation seen by either is
  never falsely certified clean. Reproduces the canonical results (flat
  Minkowski clean; Alcubierre violates NEC/WEC, deepening with velocity).
  `convergence_order(...)` (CLI `warpfactory --convergence`) Richardson-checks
  the finite-difference convergence order (Lousto, arXiv:gr-qc/0503001). Applies
  to the analytic seeds now and to reassembled evolved plotfile metrics later
  (feeding an `energy_conditions` score component). *A PINN constraint-solver
  (Li et al., arXiv:2309.07397) is the planned mesh-free Path B / Level 3.*

### Search optimization
- **Surrogate-assisted CMA-ES** (`search.surrogate` + `optimize --surrogate`):
  fits an RBF regressor on the `(θ→S)` archive each generation and evaluates
  only the top fraction (plus high-uncertainty points) on the GPU; the rest get
  a predicted fitness. ~25–40% fewer GPU evaluations in practice.
- **MAP-Elites quality diversity** (`search.qd_search`, CLI `qd`): keeps the best
  elite per behavior cell over a `(FTL benefit, exoticity)` grid; multi-GPU
  parallel batches; writes `archive.json` + `trajectory.jsonl`.
- **Multi-objective Pareto** (`search.pareto`, CLI `pareto`): non-dominated front
  over `(F_FTL, s_anec, s_growth, s_tidal)` from any optimizer trajectory.

---

## Score components

`metrics.score` maps each violation `v` through a bounded reward `r(v;σ)=1/(1+v/σ)`
and sums weighted components. Defaults (`DEFAULT_WEIGHTS`):

| Component | Weight | Source |
|-----------|--------|--------|
| `ftl_shortcut` (`F_log`) | 5.0 | `metrics.ftl_metrics` |
| `expansion_asymmetry` | 2.0 | `metrics.ftl_metrics` |
| `comoving_stability` | 2.5 | `metrics.episode_metrics` (comoving) |
| `constraint_health` | 2.0 | `metrics.episode_metrics` (constraints) |
| `energy_condition` | 2.0 | `metrics.episode_metrics` (constraints) |
| `survival` | 1.5 | time series |
| `horizon_penalty` | 1.5 | `metrics.episode_metrics` (collapse) |
| `lapse_health` | 1.0 | `metrics.episode_metrics` (collapse) |
| `nonflat_geometry` | 1.0 | `metrics.ftl_metrics` |
| `stability` (Eulerian) | 0.5 | `metrics.episode_metrics` (stability) |
| `initial_constraint_quality` | 0.5 | `metrics.episode_metrics` (constraints) |
| `nontrivial_geometry` | 0.25 | `metrics.episode_metrics` (collapse) |
| **`constraint_growth`** | 2.0 | `metrics.episode_metrics` (growth) |
| **`anec_condition`** | 1.5 | `metrics.physical_metrics` |
| **`tidal_comfort`** | 1.0 | `metrics.physical_metrics` |

New components contribute `0` when their diagnostics are unavailable, so legacy
episodes re-score consistently. Override any weight with `--score-weight key=val`.

---

## GRTresna integration (constraint-satisfying initial data)

GRTresna is a Chombo-based elliptic solver that produces fully
constraint-satisfying 3D initial data (Hamiltonian + momentum constraints)
for numerical relativity. The bridge lets GRTeclyn start evolutions from
GRTresna output instead of the 1D radial recipe, giving ~100x lower
initial constraint violations.

### Architecture

```
theta (BH mass/spin/scalar field params)
  |
  v
grtresna/solver.py  -->  GRTresna (MPI, Chombo)  -->  InitialDataFinal.3d.hdf5
  |                                                          |
  v                                                          v
grtresna/io.py  ----  Chombo HDF5 reader  ---->  initial_data.gridinit
  |                   (strips ghost cells,
  |                    paints coarse-to-fine,
  |                    z-reflection with parity)
  v
GRTeclyn C++:  ExternalGridInitialData  ----  trilinear interp onto AMReX grid
  |            (GPU managed memory)
  v
RadialRecipeLevel::initData()  -->  evolution
```

### Python modules

See **GRTresna bridge (`grtresna/`)** in the module map above. Public API:

```python
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve, convert_chombo_to_gridinit
# or, equivalently:
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, solve
from grteclyn_wrapper.grtresna.io import convert_chombo_to_gridinit
```

`grtresna.io` requires `h5py` for the Chombo→`.gridinit` conversion. It is
imported lazily (so constructing a `GRTresnaConfig`, dry-runs, etc. do not need
it) and is installed in the wrapper env via `uv pip install h5py`.

### C++ side

| File | Role |
|------|------|
| `ExternalGridInitialData.hpp` | Reads `.gridinit` into AMReX managed memory; trilinear-interpolates onto the GPU grid. |
| `SimulationParameters.hpp` | Added `recipe_initial_data_file` parameter. |
| `RadialRecipeLevel.cpp` | `initData()` branches: uses `ExternalGridInitialData` when file is set, else the radial recipe. |

### Usage (runtime)

Add one line to your `params.txt`:

```
recipe_initial_data_file = /path/to/initial_data.gridinit
```

Or from the Python orchestrator:

```python
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve

cfg = GRTresnaConfig(
    mpi_ranks=8,
    bh1_bare_mass=1.0,
    bh1_spin=(0.0, 0.0, 0.5),
    dphi=0.1,
    gridinit_nx=64, gridinit_ny=64, gridinit_nz=64,  # per-axis target resolution
)
gridinit_path = solve(cfg, work_dir=Path("/tmp/grtresna_run"))
# gridinit_path is now ready to pass to GRTeclyn
```

### Momentum-carrying matter (moving / rotating scalar cloud)

The public GRTresna solver assumes conformal flatness, so it cannot represent
arbitrary non-flat metric *shapes* — but it **does** solve the momentum
constraint, which unlocks a real, distinct FTL ingredient the 1D radial recipe
cannot represent: **matter that carries net momentum** (`S_i = -Π ∂_i φ ≠ 0`),
the source of matter-momentum-driven frame dragging.

A localised scalar "lump" is added on top of the legacy spherical profile and
its conjugate momentum is built as the convective derivative of a boosted and/or
rigidly rotating pattern, `Π = -(v · ∇φ) - Ω ∂_φ φ`, so the configuration
carries net linear momentum `P_i ~ v_i` and/or angular momentum `L_z ~ Ω`.
Implemented in the `ScalarFieldBH` example (`MyMatterFunctions.cpp`,
`MatterParams.hpp`); all params default to off, so existing data is unchanged.

| `GRTresnaConfig` field | params.txt key | Meaning |
|------------------------|----------------|---------|
| `lump_amp` | `lump_amp` | amplitude (0 ⇒ disabled) |
| `lump_width` | `lump_width` | Gaussian width |
| `lump_center` | `lump_center` | centre relative to grid centre |
| `lump_velocity` | `lump_velocity` | boost `v` ⇒ linear momentum |
| `lump_omega` | `lump_omega` | rotation rate about z ⇒ `L_z` |
| `lump_mode` | `lump_mode` | azimuthal `m` (`≥ 1` required for `L_z`) |

```python
cfg = GRTresnaConfig(
    mpi_ranks=8,
    bh1_bare_mass=0.0,                # pure matter-momentum case
    lump_amp=0.1, lump_width=8.0,
    lump_velocity=(0.2, 0.0, 0.0),    # boosted cloud → linear momentum
    # angular momentum instead:  lump_omega=0.1, lump_mode=1
)
```

Validated with a boosted lump: the momentum constraint starts at ~98% violation
(real injected momentum) and the solver drives it to ~0.1% by solving for the
longitudinal vector potential `V_i`. A standalone smoke test lives at
`GRTresna/Examples/ScalarFieldBH/params_momentum_test.txt`.

### GRTresna-in-the-loop search (`--grtresna`)

`search/optimize.py` can run the GRTresna solve **per CMA-ES candidate** and
evolve the result on GPU, searching the momentum-cloud parameters:

```bash
uv run python -m grteclyn_wrapper --runs-dir runs \
  optimize --grtresna --gpu-ids 0 1 2 3 \
  --max-generations 30 --grtresna-ranks 8 --grtresna-iterations 50
```

- `grtresna_search_space(K)` defines a **K-lump scalar basis**: per lump `k`,
  the dims `grtresna_lump{k}_{amp,width,center_{x,y,z},velocity_{x,y,z},omega,mode}`
  (10 each; `--grtresna-lumps K`, default 3). `build_grtresna_config()` groups
  the indexed keys into `GRTresnaConfig.lumps` (rounding `mode` to int) while
  everything else flows to GRTeclyn's `params.txt` unchanged. A superposition of
  lumps paints arbitrary `rho`/`S_i`; e.g. two counter-moving lumps give a
  momentum/shear field with no bulk translation.
- The `grtresna_*` keys are intentionally kept **out** of GRTeclyn's
  `params.txt`; only the resulting `recipe_initial_data_file` is set.
- A solver failure for a candidate is penalised (score 0) rather than crashing
  the search; `cleanup=True` deletes the heavy HDF5 after each conversion.
- This mode **replaces** the radial-recipe search space (the initial data comes
  entirely from the solve); it is separate from `--nonspherical`.

### Installing GRTresna from scratch (no root)

The full toolchain is installed via micromamba in a dedicated conda
environment. This section documents every step so it can be reproduced
on a new machine.

**1. Install micromamba** (if not already present):

```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
  | tar -xvj -C /path/to/your/home bin/micromamba
export MAMBA_ROOT_PREFIX=/path/to/your/home/micromamba_root
```

**2. Create the environment with all dependencies:**

```bash
eval "$(/path/to/your/home/micromamba shell hook -s bash)"

micromamba create -n grtresna --override-channels -c conda-forge -y \
    gcc_linux-64 gxx_linux-64 gfortran_linux-64

micromamba install -n grtresna --override-channels -c conda-forge -y \
    openmpi openmpi-mpicxx "hdf5=*=mpi_openmpi*" \
    lapack blas h5py tcsh zlib make perl
```

**3. Create compiler symlinks** (conda-forge uses cross-compiler names):

```bash
micromamba activate grtresna
cd $CONDA_PREFIX/bin
ln -sf x86_64-conda-linux-gnu-gfortran gfortran
ln -sf x86_64-conda-linux-gnu-g++      g++
ln -sf x86_64-conda-linux-gnu-gcc      gcc
ln -sf tcsh                             csh
```

**4. Clone Chombo and GRTresna:**

```bash
cd /path/to/simulation
git clone https://github.com/GRTLCollaboration/Chombo.git
git clone https://github.com/GRTLCollaboration/GRTresna.git
```

**5. Configure Chombo** — create `Chombo/lib/mk/Make.defs.local`:

```makefile
#begin  -- dont change this line

DIM                = 3
DEBUG              = FALSE
OPT                = HIGH
PRECISION          = DOUBLE
OPENMPCC           = FALSE
MPI                = TRUE

CXX                = g++
FC                 = gfortran
MPICXX             = mpicxx

USE_HDF            = TRUE
HDFINCFLAGS        = -I$CONDA_PREFIX/include
HDFLIBFLAGS        = -L$CONDA_PREFIX/lib -lhdf5 -lz -Wl,-rpath,$CONDA_PREFIX/lib
HDFMPIINCFLAGS     = -I$CONDA_PREFIX/include
HDFMPILIBFLAGS     = -L$CONDA_PREFIX/lib -lhdf5 -lz -Wl,-rpath,$CONDA_PREFIX/lib

syslibflags        = -lblas -llapack

cxxoptflags        = -march=native -O3 -fpermissive

#end  -- dont change this line
```

> Replace `$CONDA_PREFIX` with the actual path
> (e.g. `/home/you/.mlspace/envs/grtresna`).

**6. Patch Chombo for csh/tcsh paths:**

Two files hardcode `/bin/csh` which doesn't exist on many systems:

```bash
# Chombo/lib/mk/reverse — change shebang
sed -i "1s|.*|#!$CONDA_PREFIX/bin/tcsh -f|" Chombo/lib/mk/reverse

# Chombo/lib/mk/Make.rules — change CSHELLCMD
sed -i "s|CSHELLCMD=/bin/csh|CSHELLCMD=$CONDA_PREFIX/bin/tcsh|" \
    Chombo/lib/mk/Make.rules
```

**7. Build Chombo:**

```bash
export CHOMBO_HOME=/path/to/simulation/Chombo/lib
cd $CHOMBO_HOME && make lib -j$(nproc)
```

**8. Build GRTresna example:**

```bash
export CHOMBO_HOME=/path/to/simulation/Chombo/lib
cd /path/to/simulation/GRTresna/Examples/ScalarFieldBH
make all -j$(nproc)
```

**9. Test run:**

```bash
cd /path/to/simulation/GRTresna/Examples/ScalarFieldBH
mkdir -p Outputs pout
mpirun --oversubscribe -np 8 \
    ./Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex \
    params.txt max_NL_iterations=25 write_diagnostics=0
# Check convergence:
cat Ham_and_Mom_errors.txt
# Expected: Ham ~ 0.003%, Mom ~ 0.026%
# Output: Outputs/InitialDataFinal.3d.hdf5
```

**10. Convert and use in GRTeclyn:**

```python
from grteclyn_wrapper.grtresna import convert_chombo_to_gridinit

convert_chombo_to_gridinit(
    "GRTresna/Examples/ScalarFieldBH/Outputs/InitialDataFinal.3d.hdf5",
    "initial_data.gridinit",
    nx=64, ny=64, nz=64,  # per-axis target uniform grid cells (non-cubic OK)
    L=128.0               # domain side length (must match GRTresna params)
)
```

Then run GRTeclyn with `recipe_initial_data_file = initial_data.gridinit`.

### Troubleshooting

| Symptom | Fix |
|---------|-----|
| `/bin/csh: not found` during Chombo build | Patch `reverse` shebang and `Make.rules` CSHELLCMD (step 6). |
| `invalid conversion from ch_offset_t*` | Add `-fpermissive` to `cxxoptflags` in `Make.defs.local`. |
| `cannot find -lz` | Install `zlib` via micromamba; add `-Wl,-rpath,$CONDA_PREFIX/lib` to `HDFLIBFLAGS`. |
| GRTresna run hangs with 64 ranks | The grid has few boxes (8 for 64x64x32 / 16^3). Use ~8 MPI ranks to match the box count. |
| `gcc versions later than 12 not supported` (GRTeclyn build) | The micromamba gcc is on PATH; build GRTeclyn with the system gcc (deactivate micromamba first or set `PATH` explicitly). |
| All-zero `.gridinit` data | Ghost cells not stripped. The converter handles this automatically via offset-based ghost inference. |

### Variable mapping

All 27 state variables have **identical names** between GRTresna and
GRTeclyn (both from the GRTL collaboration, CCZ4 formulation):

```
chi  h11 h12 h13 h22 h23 h33  K
A11 A12 A13 A22 A23 A33  Theta
Gamma1 Gamma2 Gamma3  lapse  shift1 shift2 shift3
B1 B2 B3  phi  Pi
```

Component index `i` in the `.gridinit` maps directly to enum value
`c_chi + i` in `CCZ4StateVariables.hpp` + `StateVariables.hpp`.

### Validated results

Round-trip GRTresna ScalarFieldBH → GRTeclyn (64^3, single BH + scalar
field, spinning BH `a=0.5`):

| Stage | Ham constraint | Mom constraint |
|-------|---------------|----------------|
| GRTresna solver output | 0.003% | 0.026% |
| GRTeclyn t=0.04 (first step) | L2 = 1.1e-3 | L2 = 1.6e-4 |

---

## Tests

Run from `grteclyn-wrapper/`:

```bash
cd grteclyn-wrapper

uv run python tests/test_proposed_extensions.py   # growth, ANEC/tidal, surrogate, MAP-Elites, Pareto
uv run python tests/test_ftl_metrics.py
uv run python tests/test_upgraded_scoring.py
uv run python tests/test_stability_score.py
uv run python tests/test_constrained_guesser.py
uv run python tests/test_ftl_general.py           # metrics.ftl_general
uv run python tests/test_warpfactory.py           # metrics.warpfactory
uv run python tests/test_candidates.py            # initial_data.candidates (needs pytest)
```
