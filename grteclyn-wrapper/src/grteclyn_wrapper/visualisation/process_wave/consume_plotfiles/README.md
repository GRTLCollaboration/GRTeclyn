# consume_plotfiles

Stream AMReX plotfiles during a run: extract small-data time series, render diagnostic PNG frames, optionally delete processed HDF5 to save disk.

Invoked as a sidecar by `plot_consumer.build_consume_command` or directly:

```bash
python -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles --help
```

## Layout

| Path | Role |
|------|------|
| `driver.py` | `main()` — argparse, watch loop, parallel dispatch |
| `worker.py` | `_process_single_plotfile()` — one plotfile → outputs |
| `config.py` | Frame DPI, color limits, default paths |
| `plotfiles.py` | Discover plotfile dirs, readiness, restart detection |
| `fields.py` | Field aliases and yt derived-field registration |
| `sphere.py` | Sphere sampling grid and spin-weighted Y₂₀ |
| `state.py` | `consume_state.json` and `.dat` append helpers |
| `extraction/` | Psi4 mode, shell stats, areal radius, FTL timeseries |
| `frames/` | Slice/projection/embedding PNG rendering and cleanup |

## Public API

```python
from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles import main
from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles import (
    get_extraction_points,
    spin_weighted_sph_harm_s2_l2_m0,
)
```
