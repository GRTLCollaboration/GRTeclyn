# CLI command handlers

Each module implements one subcommand from `cli/parser.py`. Called from `cli/main.py` after shared arg parsing.

| Module | Command | Calls into |
|--------|---------|------------|
| `optimize.py` | `optimize` | `search.optimize.run_optimize` |
| `qd.py` | `qd` | `search.qd_search.run_qd_search` |
| `atlas.py` | `atlas` | `search.atlas.run_atlas` |
| `pareto.py` | `pareto` | `search.pareto` |
| `warpfactory.py` | `warpfactory` | `metrics.probes.warpfactory` |

`reproduce`, `sweep`, and `validate` are handled directly in `cli/episode.py` and `initial_data/validate_guesser.py`.
