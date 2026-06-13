# Failure atlas (random batch)

Low-resolution random sampling over parameter ranges. Classifies each episode (constraint blowup, lapse collapse, horizon, trivial geometry) and writes JSONL/CSV for exploratory failure mapping.

## Layout

| Path | Role |
|------|------|
| `driver.py` | `run_atlas` — batch loop |
| `config.py` | Parameter ranges, paths, `sample_overrides` |
| `records.py` | `classify_episode`, `build_record`, `summarize_records` |
| `io.py` | JSONL/CSV append and flatten for export |

## Public API

```python
from grteclyn_wrapper.search.atlas import run_atlas, AtlasPaths
```
