# MAP-Elites quality-diversity search

Maintains an archive of the best candidate per behavior-descriptor cell (Spacetime Failure Atlas). Explores the constructibility map instead of converging to a single optimum.

## Layout

| Path | Role |
|------|------|
| `driver.py` | `run_qd_search` — campaign loop, resume, validation snapshots |
| `descriptors.py` | Behavior axes (`legacy`, `channel`, `speed_horizon`, `speed_super`) |
| `archive.py` | `Elite`, `QDArchive` — grid insert, coverage, serialization |
| `sampling.py` | Elite mutation, feasible-box exploration, boundary reflection |
| `io.py` | Trajectory load, eval-dir pruning, target-eval batch math |

## Public API

```python
from grteclyn_wrapper.search.qd_search import run_qd_search, QDArchive, Elite
```
