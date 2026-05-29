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

> **Gotcha:** do **not** `source scripts/env.sh` before `optimize`/`qd`. It sets
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

### Pipeline core
| Module | Role |
|--------|------|
| `__main__.py` | CLI: argument parsing and command dispatch. |
| `config.py` | Example/executable resolution, repo paths, default run dir. |
| `episode.py` | Per-run directory layout (`Episode`), metadata, JSON writers. |
| `params.py` | Render a `params.txt` from a template + overrides. |
| `runner.py` | Launch the GRTeclyn binary; optional plotfile consumer. |
| `evaluation.py` | **Shared candidate→episode→score helper** used by all drivers. |

### Initial data
| Module | Role |
|--------|------|
| `constrained_recipe.py` | Derive `phi` from `chi` via the Hamiltonian constraint; Gaussian basis; Ricci scalar. |
| `preflight.py` | Cheap 1D constraint filter; reject bad candidates before the GPU. |
| `seeds.py` | Known-solution seeds: flat, Ellis–Bronnikov, Schwarzschild puncture, Alcubierre warp. |
| `candidates.py` | Resolve initial-data overrides from seed / candidate / non-spherical IDs. |
| `validate_guesser.py` | Batch validation of the guesser on synthetic proposals. |

### Metrics & scoring
| Module | Role |
|--------|------|
| `ftl_metrics.py` | t=0 FTL shortcut metrics: `F_null`, `F_portal`, `F_throat`, `F_asymmetry`, `F_log`, `s_nonflat`. |
| `physical_metrics.py` | **NEW** — t=0 gauge-robust proxies: ANEC line integral, curvature/tidal proxy, trapped-surface flag. |
| `warpfactory.py` | **NEW** — Warp Factory port: 4-metric → Einstein tensor → multi-observer NEC/WEC/SEC/DEC. |
| `metrics.py` | Parse diagnostics into `EpisodeMetrics`; collapse/constraint/stability/comoving + **NEW** growth-rate (`GrowthMetrics`). |
| `score.py` | Weighted multi-component scalar reward `score_episode()`. |

### Search drivers
| Module | Role |
|--------|------|
| `optimize.py` | CMA-ES loop, multi-GPU parallel generations; **NEW** optional surrogate screening. |
| `surrogate.py` | **NEW** — numpy RBF kernel-ridge surrogate + candidate screening. |
| `qd_search.py` | **NEW** — MAP-Elites quality-diversity driver and archive. |
| `pareto.py` | **NEW** — non-dominated sorting / Pareto-front extraction. |
| `atlas.py` | Random failure-atlas batch runner. |

---

## New in this work (`feature/interstellar`)

These implement the paper's "Proposed Extensions" (Sec. *Extensions: Implemented
Metrics, Search, and Robustness*). All are numpy-only (no new dependencies) and
GPU-validated; unit tests live in `tests/test_proposed_extensions.py`.

### New physical metrics
- **Constraint growth rate** (`metrics.GrowthMetrics`, scored as
  `constraint_growth`): log-linear fit of the exponential rate `λ` on the
  `‖H‖`, `|K|_max`, and `1/χ_min` time series. `s_growth = 1/(1 + max(0,λ)/σ_λ)`
  penalizes geometries that merely collapse slowly enough to survive a short
  run — closing the dynamical-stability gap without longer evolutions.
- **ANEC line proxy** (`physical_metrics`, scored as `anec_condition`): the
  required energy density integrated along the travel axis,
  `A_line = ∫ ρ_req(x) dx`; a directional energy-condition indicator. *t=0,
  χ-sourced — blind to pure-shift warps; the full geodesic ANEC needs the
  evolved probe (future).*
- **Tidal / curvature proxy** (`physical_metrics`, scored as `tidal_comfort`):
  `max|R| + max|∂²α|` with a companion t=0 trapped-surface flag. Precursor to
  the full Kretschmann / electric-Weyl invariants (future C++ work).
- **Multi-observer energy conditions** (`warpfactory.py`, CLI `warpfactory`): a
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
- **Surrogate-assisted CMA-ES** (`surrogate.py` + `optimize --surrogate`):
  fits an RBF regressor on the `(θ→S)` archive each generation and evaluates
  only the top fraction (plus high-uncertainty points) on the GPU; the rest get
  a predicted fitness. ~25–40% fewer GPU evaluations in practice.
- **MAP-Elites quality diversity** (`qd_search.py`, CLI `qd`): keeps the best
  elite per behavior cell over a `(FTL benefit, exoticity)` grid; multi-GPU
  parallel batches; writes `archive.json` + `trajectory.jsonl`.
- **Multi-objective Pareto** (`pareto.py`, CLI `pareto`): non-dominated front
  over `(F_FTL, s_anec, s_growth, s_tidal)` from any optimizer trajectory.

---

## Score components

`score.py` maps each violation `v` through a bounded reward `r(v;σ)=1/(1+v/σ)`
and sums weighted components. Defaults (`DEFAULT_WEIGHTS`):

| Component | Weight | Source |
|-----------|--------|--------|
| `ftl_shortcut` (`F_log`) | 5.0 | `ftl_metrics` |
| `expansion_asymmetry` | 2.0 | `ftl_metrics` |
| `comoving_stability` | 2.5 | `metrics.comoving` |
| `constraint_health` | 2.0 | `metrics.constraints` |
| `energy_condition` | 2.0 | `metrics.constraints` |
| `survival` | 1.5 | time series |
| `horizon_penalty` | 1.5 | `metrics.collapse` |
| `lapse_health` | 1.0 | `metrics.collapse` |
| `nonflat_geometry` | 1.0 | `ftl_metrics` |
| `stability` (Eulerian) | 0.5 | `metrics.stability` |
| `initial_constraint_quality` | 0.5 | `metrics.constraints` |
| `nontrivial_geometry` | 0.25 | `metrics.collapse` |
| **`constraint_growth`** | 2.0 | `metrics.growth` (NEW) |
| **`anec_condition`** | 1.5 | `physical_metrics` (NEW) |
| **`tidal_comfort`** | 1.0 | `physical_metrics` (NEW) |

New components contribute `0` when their diagnostics are unavailable, so legacy
episodes re-score consistently. Override any weight with `--score-weight key=val`.

---

## Tests

```bash
uv run python tests/test_proposed_extensions.py   # growth, ANEC/tidal, surrogate, MAP-Elites, Pareto
uv run python tests/test_ftl_metrics.py
uv run python tests/test_upgraded_scoring.py
uv run python tests/test_stability_score.py
uv run python tests/test_constrained_guesser.py
```
