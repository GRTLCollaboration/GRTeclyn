# GRTresna Bridge

Python bridge between the external [GRTresna](https://github.com/GRChombo/GRTresna) elliptic initial-data solver and GRTeclyn time evolution. Given a matter/geometry configuration, the package writes `params.txt`, runs GRTresna via MPI, converts the Chombo HDF5 checkpoint to a `.gridinit` file, and wires GRTeclyn evolution parameters so constraint-solved initial data replays correctly on the GPU.

## Layout

| Path | Role | Key entry points |
|------|------|------------------|
| `solver/` | Orchestration: config, params I/O, MPI run, convergence | `GRTresnaConfig`, `solve()`, `write_grtresna_params()`, `parse_convergence()` |
| `domain/` | Solve domain ↔ evolution box mapping, export window | `GRTresnaDomainConfig`, `GridinitExportSpec` |
| `matter/` | Matter registry, GRTeclyn wiring, splash overrides | `matter.models`, `matter.wiring`, `matter.splash` |
| `io/` | Chombo HDF5 → uniform grid → `.gridinit` | `convert_chombo_to_gridinit()`, `read_gridinit()`, `write_gridinit()` |
| `fields/` | Analytic field painting (lumps, boson star, SH ansatz) | `fields.lump`, `fields.boson_star`, `fields.sh` |
| `profiles/` | Reference radial physics (boson star, Q-ball, envelopes) | `profiles.boson_star`, `profiles.qball_ode`, `profiles.envelope` |
| `fit/` | Geometry-first motif → GRTresna config | `fit.motif.build_grtresna_config_from_fitted()` |

## Data flow

```
Search / CLI overrides
        │
        ▼
  GRTresnaConfig  ──►  params.txt
        │
        ▼
  GRTresna (MPI)  ──►  InitialDataFinal.3d.hdf5
        │
        ▼
  io/conversion   ──►  initial_data.gridinit  (+ optional .matter.json)
        │
        ▼
  GRTeclyn ExternalGridInitialData evolution
```

## Public API

Lazy exports from the package root (lightweight import):

```python
from grteclyn_wrapper.grtresna import GRTresnaConfig, solve
from grteclyn_wrapper.grtresna import convert_chombo_to_gridinit, read_chombo_domain
```

Subpackages are the canonical import paths for everything else:

```python
from grteclyn_wrapper.grtresna.matter.models import resolve_matter_selection
from grteclyn_wrapper.grtresna.fields.lump import lump_phi_at
from grteclyn_wrapper.grtresna.fit.motif import fit_matter_from_motif
```

## Related code (outside this package)

GRTresna-specific search and CLI logic lives in sibling modules, not under `grtresna/`:

| Path | Role |
|------|------|
| `search/grtresna_convergence_gate.py` | Ham/Mom convergence gating for search |
| `search/grtresna_evaluation_gates.py` | Pre-GPU evaluation gates |
| `cli/grtresna_args.py`, `cli/grtresna_context.py` | CLI flags and search context |
| `search/optimize/config.py` | Search dimensions → `GRTresnaConfig` |

## Tests

```bash
cd grteclyn-wrapper
uv run pytest tests/grtresna/ -q
```

Cross-package tests that touch this bridge also live under `tests/projection/`, `tests/search/`, and `tests/metrics/ftl/`.
