from __future__ import annotations

import json
from pathlib import Path
from typing import Dict


def _load_state(state_path: Path) -> Dict[str, bool]:
    if not state_path.exists():
        return {}
    try:
        return json.loads(state_path.read_text())
    except Exception:
        return {}


def _save_state(state_path: Path, state: Dict[str, bool]) -> None:
    state_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = state_path.with_suffix(state_path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    tmp.replace(state_path)


def _append_line(path: Path, header: str, line: str) -> None:
    is_new = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        if is_new:
            f.write(header.rstrip() + "\n")
        f.write(line.rstrip() + "\n")
