# Episode Metrics

Measurement and scoring for GRTeclyn simulation episodes. Only code that feeds `EpisodeMetrics` or `score_episode()` lives here.

## Layout

| Path | Role |
|------|------|
| `types/` | Frozen dataclasses (`CollapseMetrics`, `EpisodeMetrics`, …) |
| `diagnostics/` | Parsers for C++ `.dat` diagnostic files |
| `probes/` | Computed metrics from recipe overrides, plotfiles, or `.gridinit` |
| `probes/ftl/` | FTL family: analytic, general (Dijkstra), solved, frozen geodesic, **evolving geodesic** |
| `aggregation/` | `read_episode_metrics()` orchestration |
| `io/` | `.dat` parsing, serialization, plotfile thread lock |
| `catalog.py` | `METRIC_REGISTRY` — single source of truth for metric groups |
| `score/` | `score_episode()` reward engineering (phased component calculators) |
| `score/scorer.py` | `score_episode()` orchestrator |
| `score/survival.py` | numerical_survival, structural_persistence, survival |
| `score/health.py` | constraint, lapse, stability, comoving, growth, physical, t=0 FTL profile |
| `score/ftl.py` | operational/geodesic FTL, evolving geodesic diagnostic, precursor, channel, shaping gates |
| `score/penalties.py` | exotic, qei, boundary, transport, stationary penalties |
| `score/objectives.py` | `weighted`, `ftl_first`, `robust_ftl` scalarization |

## Metric groups (in `EpisodeMetrics`)

| Group | Module | Source | Score components |
|-------|--------|--------|------------------|
| collapse | `diagnostics/collapse.py` | `collapse_diagnostics.dat` | numerical_survival (completion gate), lapse_health, horizon_penalty |
| constraints | `diagnostics/constraints.py` | `constraint_norms.dat` | constraint_health, initial_constraint_quality, structural_persistence → survival |
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
| geodesic_ftl | `probes/ftl/geodesic.py` | null-ray on latest plotfile (frozen slice) | operational_ftl_geodesic |
| evolving_geodesic | `probes/ftl/evolving_geodesic.py` | null-ray on ≥3 plotfile stack (4D metric) | **ftl_geo_evolving** (weight 0, diagnostic) |
| ftl_timeseries | `diagnostics/ftl_timeseries.py` | `small_data/ftl_timeseries.dat` | ftl_geo_timeavg, ftl_geo_peak, ftl_lifetime |
| physical | `probes/physical.py` | recipe t=0 proxies | anec_condition, tidal_comfort |
| boundary_flux | `probes/boundary.py` | `boundary_flux.dat` / plotfile | boundary_penalty |
| effective_ec | `probes/warpfactory.py` | >= 3 plotfile stack | exotic_penalty |

Use `from grteclyn_wrapper.metrics.catalog import list_metrics, get_metric` for programmatic lookup.

## FTL stack precedence

Scoring rewards **evolved** operational FTL (`general_ftl_evolved`) over t=0 slices. A t=0 shortcut that relaxes away scores nothing on `operational_ftl`. Geodesic confirmation (`geodesic_ftl`) gates the highest-weight component.

### Frozen vs evolving null-geodesic

| Probe | Metric field | When | Primary score use |
|-------|--------------|------|-------------------|
| `geodesic_ftl` | Single Cauchy slice (∂g/∂t = 0) | Per-frame in consumer + final collector | Diagnostic only when 4D ran; else `operational_ftl_geodesic` (frozen timeavg) |
| `evolving_geodesic` | Time-interpolated stack from ≥3 plotfiles | End-of-run when enabled | **`ftl_geo_evolving`** (headline in QD; frozen credit zeroed) |

The evolving probe traces null rays through `EvolvingMetricField` (temporal linear interpolation between plotfile slices at the ray's coordinate time). It answers whether a pulse emitted at `t_emit = times[0]` beats flat-space transit through the **full evolved geometry**, not a frozen mid-run snapshot. Outputs land in `small_data/evolving_geodesic.json` and the last row of `ftl_timeseries.dat` (`f_geo_evol`, `f_geo_evol_ok`).

**Enable:**

- QD search: `GRTECLYN_EVOLVING_GEODESIC=1` + `GRTECLYN_EVOLVING_GEODESIC_MODE=search` (set in `campaigns/lib/search_common.sh`; fast profile)
- HQ promotion: `GRTECLYN_EVOLVING_GEODESIC_MODE=hq` (set in `promote_common.sh`) or `--evolving-geodesic` on `campaigns/hq/replay_eval.py`
- Collector: `read_episode_metrics(..., evolving_geodesic=True)` or env `GRTECLYN_EVOLVING_GEODESIC=1`

Same gate as frozen geodesic (`f_op > 1e-3` or `max_local_speed > 1`) plus `len(plotfiles) >= 3`. Validation: `scripts/validation/run_evolving_geodesic_smoke.sh`.

## Outside this package

| Module | Location | Why |
|--------|----------|-----|
| Motif preservation QA | `projection/motif_preservation.py` | Projection gate, not an episode metric |
| Symbolic regression | `analysis/symbolic_extract.py` | Optional PySR tooling |

## Public API

```python
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode, EpisodeMetrics
```
