# Falsification tier ladder

Records how far each candidate survives an increasingly strict evidence ladder (constructed → nontrivial → operational → persistent → observer EC → converged → analytic). Used by MAP-Elites and offline campaign assessment.

## Layout

| Path | Role |
|------|------|
| `types.py` | `Tier`, `TierConfig`, `TierAssessment`, gate status constants |
| `evaluate.py` | `evaluate_tiers()` — T0–T6 gate chain |
| `survivors.py` | `build_survivors`, `survivor_front` — Pareto front of survivors |
| `convergence.py` | `convergence_signals` — archive stall detection |
| `campaign.py` | `assess_campaign` — offline replay from `eval_*/score.json` |

## Public API

```python
from grteclyn_wrapper.search.validation_tiers import evaluate_tiers, Tier, assess_campaign
```
