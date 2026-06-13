# CMA-ES optimizer

Gradient-free search over RadialRecipe / GRTresna parameter vectors. Each candidate runs a full GRTeclyn episode (optionally with GRTresna initial-data solve) and returns a scalar fitness for CMA-ES.

## Layout

| Path | Role |
|------|------|
| `driver.py` | `run_optimize`, `OptimizeResult` — main loop |
| `eval.py` | Single-candidate and multi-GPU parallel evaluation |
| `spaces.py` | Search-space definitions and `build_search_space()` |
| `candidates.py` | Vector ↔ overrides, warm-start, jitter, exotic injection |
| `config.py` | `build_grtresna_config`, convergence parsing |
| `geometry.py` | GRTresna shell/ring lump layout helpers |
| `dimension.py` | `SearchDimension` dataclass |

## Public API

```python
from grteclyn_wrapper.search.optimize import run_optimize, SearchDimension, build_search_space
```
