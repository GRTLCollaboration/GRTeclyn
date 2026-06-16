# grteclyn-wrapper Tests

Pytest suite organized by domain, mirroring `grteclyn_wrapper/` package layout.

## Layout

| Path | Module under test | Example |
|------|-------------------|---------|
| `conftest.py` | Shared fixtures (`rng`) | — |
| `metrics/scoring/` | `metrics.score` — objectives, tiers, stability, horizon guards | `uv run pytest tests/metrics/scoring/ -q` |
| `metrics/ftl/` | `metrics.probes.ftl`, energy conditions, campaign diagnostics | `uv run pytest tests/metrics/ftl/ -q` |
| `search/` | `search.qd_search`, `search.optimize`, retention, trajectory logging | `uv run pytest tests/search/ -q` |
| `grtresna/` | `grtresna.*` — solver, ansätze, I/O, matter, convergence | `uv run pytest tests/grtresna/ -q` |
| `projection/` | `projection.*` — geometry-first → GRTresna pipeline | `uv run pytest tests/projection/ -q` |
| `initial_data/` | `initial_data.*` — seeds, candidates, constrained guesser | `uv run pytest tests/initial_data/ -q` |
| `analysis/` | `analysis.symbolic_extract` | `uv run pytest tests/analysis/ -q` |
| `scripts/` | Repo script smoke tests (`scripts/wormhole/rollback`) | `uv run pytest tests/scripts/ -q` |
| `tools/` | Manual helpers (plots, validation scripts; not pytest) | `uv run python tests/tools/validate_metric_guesser.py` |

## metrics/scoring/

| File | Focus |
|------|-------|
| `test_robust_ftl_objective.py` | `robust_ftl` vs `ftl_first` scalarization |
| `test_horizon_finder_guard.py` | Apparent-horizon parsing + off-center score guard |
| `test_stability_score.py` | Static vs collapsing episode stability scoring |
| `test_upgraded_scoring.py` | Analytic FTL + full episode scoring |
| `test_validation_tiers.py` | Candidate tier ladder (T0–T3) |

## metrics/ftl/

| File | Focus |
|------|-------|
| `test_ftl_general.py` | Mechanism-agnostic operational FTL |
| `test_ftl_metrics.py` | Analytic FTL on recipe seeds |
| `test_null_geodesic.py` | Gauge-invariant null geodesic ray-tracing |
| `test_alcubierre_validation.py` | Alcubierre positive control |
| `test_solved_geometry_ftl.py` | Operational FTL on `.gridinit` |
| `test_warpfactory.py` | NEC/WEC energy-condition evaluator |
| `test_ftl_peak_metrics.py` | Peak FTL from timeseries |
| `test_ftl_campaign_report.py` | Campaign FTL summaries and ranking |

## search/

| File | Focus |
|------|-------|
| `test_ftl_retention.py` | QD eval-dir retention + FTL champion board |
| `test_optimize_retention.py` | CMA-ES retention + FTL champions |
| `test_trajectory_log.py` | `trajectory.jsonl` field order |
| `test_proposed_extensions.py` | Growth, physical probes, Pareto, QD archive, RBF surrogate |

## grtresna/

| File | Focus |
|------|-------|
| `test_grtresna_integration.py` | End-to-end GRTresna exotic-matter integration |
| `test_grtresna_shell_ansatz.py` | Shell ansatz search space |
| `test_grtresna_ring_ansatz.py` | Ring ansatz search space |
| `test_grtresna_convergence_parse.py` | `parse_convergence` (Ham/Mom NaN) |
| `test_grtresna_postload_gate_integration.py` | Pre-evolution gate path |
| `test_matter_geometry_consistency.py` | Matter–geometry consistency |
| `test_gridinit_conversion.py` | Chombo → gridinit conversion |
| `test_scalar_lambda_potential.py` | Shared λφ⁴ potential convention |

## Run all

```bash
cd grteclyn-wrapper
uv run pytest tests/ -q
```

QD preflight gate (`scripts/campaigns/qd/run.sh` via `search_common.sh`) runs a subset under `grtresna/`, `search/`, and `metrics/ftl/` before launching cluster jobs.
