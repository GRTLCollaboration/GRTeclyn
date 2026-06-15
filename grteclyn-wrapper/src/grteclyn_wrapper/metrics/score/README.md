# Episode Scoring

Maps `EpisodeMetrics` to a scalar `Score` via weighted component calculators and objective modes. Entry point: `score_episode()`.

## Layout

| Path | Role |
|------|------|
| `__init__.py` | Public re-exports (`score_episode`, `Score`, weights, domain helpers) |
| `scorer.py` | `score_episode()` — builds `ScoringContext`, runs phased pipeline |
| `types.py` | `Score`, `ScoringContext` |
| `weights.py` | `DEFAULT_WEIGHTS`, `HEALTH_COMPONENTS` |
| `helpers.py` | `bounded_reward`, `graded`, `domain_half_width_*` |
| `horizon.py` | Trapped-surface proxy guards (`horizon_penalty_from_collapse`) |
| `survival.py` | `numerical_survival`, `structural_persistence`, `survival` |
| `health.py` | Constraints, lapse, stability, comoving, growth, physical, t=0 FTL profile |
| `ftl.py` | Operational/geodesic FTL, precursor, channel, persistence, shaping gates |
| `penalties.py` | Exotic matter, QEI, boundary, transport, stationary-artifact penalties |
| `gating.py` | `nontriviality_gate` (flat-space attractor guard) |
| `objectives.py` | `weighted`, `ftl_first`, `robust_ftl` total scalarization |

## Pipeline

`score_episode()` runs phases in fixed order (each mutates `ScoringContext.components` / `.notes`):

1. `survival.compute_survival_components`
2. `health.compute_health_components`
3. `ftl.compute_ftl_components` — includes stationary, persistence, and geodesic-contradiction gates
4. `penalties.compute_penalty_components`
5. `gating.compute_nontriviality_gate`
6. `objectives.compute_total` → `Score`

## Objective modes

| Mode | Use |
|------|-----|
| `weighted` | Default; sums `DEFAULT_WEIGHTS` × components, health terms gated by `nontriviality_gate` |
| `ftl_first` | QD/CMA-ES campaigns; geodesic FTL dominates, shaping rewards subordinate |
| `robust_ftl` | Like `ftl_first` but favors persistent, healthy, low-exotic survivors |

## Public API

```python
from grteclyn_wrapper.metrics import score_episode
from grteclyn_wrapper.metrics.score import (
    Score,
    DEFAULT_WEIGHTS,
    domain_half_width_for_episode,
)
```
