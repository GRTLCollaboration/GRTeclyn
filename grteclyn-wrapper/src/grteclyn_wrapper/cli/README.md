# CLI

Command-line interface for `python -m grteclyn_wrapper` and the `grteclyn-wrapper` console script.

`__main__.py` is a thin entry point; parsing and dispatch live here.

## Layout

| Path | Role |
|------|------|
| `parser.py` | `build_parser()` — global options and subcommands |
| `main.py` | `main()` — parse args, resolve initial data, dispatch |
| `args.py` | `--set KEY=VALUE` parsing and score-weight helpers |
| `grtresna_args.py` | GRTresna/postload CLI flags shared by optimize and qd |
| `grtresna_context.py` | Build search space + `GRTresnaConfig` from CLI args |
| `episode.py` | `reproduce` and `sweep` single-episode runners |
| `commands/` | One module per subcommand handler |

## Public API

```python
from grteclyn_wrapper.cli import main, build_parser
from grteclyn_wrapper.__main__ import main, build_parser  # same symbols
```
