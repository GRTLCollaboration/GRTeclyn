# Episode Metrics

Measurement and scoring for GRTeclyn simulation episodes. Only code that feeds `EpisodeMetrics` or `score_episode()` lives here.

## Layout

| Path | Role |
|------|------|
| `types/` | Frozen dataclasses (`CollapseMetrics`, `EpisodeMetrics`, …) |
| `diagnostics/` | Parsers for C++ `.dat` diagnostic files |
| `probes/` | Computed metrics from recipe overrides, plotfiles, or `.gridinit` |
| `probes/ftl/` | FTL family: analytic, general (Dijkstra), solved, geodesic |
| `aggregation/` | `read_episode_metrics()` orchestration |
| `io/` | `.dat` parsing, serialization, plotfile thread lock |
| `catalog.py` | `METRIC_REGISTRY` — single source of truth for metric groups |
| `score.py` | `score_episode()` reward engineering |

## Metric groups (in `EpisodeMetrics`)

| Group | Module | Source | Score components |
|-------|--------|--------|------------------|
| collapse | `diagnostics/collapse.py` | `collapse_diagnostics.dat` | survival, lapse_health, horizon_penalty |
| constraints | `diagnostics/constraints.py` | `constraint_norms.dat` | constraint_health, initial_constraint_quality |
| stability | `diagnostics/stability.py` | collapse + `areal_radius.dat` | stability, instability_penalty |
| growth | `diagnostics/growth.py` | derived time series | constraint_growth |
| comoving | `diagnostics/comoving.py` | shell profiles + shift | comoving_stability |
| energy_conditions | `diagnostics/energy_conditions.py` | `energy_conditions.dat` | energy_condition |
| curvature | `diagnostics/curvature.py` | `curvature_invariants.dat` | curvature_activity, nontrivial_geometry |
| transport | `diagnostics/transport.py` | barycenter in collapse dat | transport_objective |
| qei | `diagnostics/qei.py` | physical + geodesic trajectory | qei_penalty |
| ftl (analytic) | `probes/ftl/analytic.py` | recipe overrides t=0 | ftl_shortcut, expansion_asymmetry |
| general_ftl (t=0) | `probes/ftl/general.py` | recipe 2D slice | ftl_precursor, channel_progress |
| general_ftl_evolved | `probes/ftl/general.py` | latest plotfile | **operational_ftl** (primary) |
| ftl_persistence | `aggregation/collector.py` | last N plotfiles | ftl_persistence |
| general_ftl_solved | `probes/ftl/solved.py` | `initial_data.gridinit` | operational_ftl_solved |
| geodesic_ftl | `probes/ftl/geodesic.py` | null-ray on plotfile | operational_ftl_geodesic |
| physical | `probes/physical.py` | recipe t=0 proxies | anec_condition, tidal_comfort |
| boundary_flux | `probes/boundary.py` | `boundary_flux.dat` / plotfile | boundary_penalty |
| effective_ec | `probes/warpfactory.py` | >= 3 plotfile stack | exotic_penalty |

Use `from grteclyn_wrapper.metrics.catalog import list_metrics, get_metric` for programmatic lookup.

## FTL stack precedence

Scoring rewards **evolved** operational FTL (`general_ftl_evolved`) over t=0 slices. A t=0 shortcut that relaxes away scores nothing on `operational_ftl`. Geodesic confirmation (`geodesic_ftl`) gates the highest-weight component.

## Outside this package

| Module | Location | Why |
|--------|----------|-----|
| Motif preservation QA | `projection/motif_preservation.py` | Projection gate, not an episode metric |
| Symbolic regression | `analysis/symbolic_extract.py` | Optional PySR tooling |

## Public API

```python
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode, EpisodeMetrics
```
