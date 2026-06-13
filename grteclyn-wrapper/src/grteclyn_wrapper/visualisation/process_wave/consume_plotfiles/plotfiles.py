from __future__ import annotations

import os
import re
import time
from pathlib import Path
from typing import Dict, List


def _iter_plotfile_dirs(data_dir: str) -> List[str]:
    """Return sorted plotfile directories under data_dir."""
    out: List[str] = []
    if not os.path.isdir(data_dir):
        return out
    prefixes = (
        "WormholePlt",
        "SupportedWormholePlt",
        "RotatingWormholePlt",
        "RadialRecipePlt",
        "plt",
    )
    for name in os.listdir(data_dir):
        if not any(name.startswith(prefix) for prefix in prefixes):
            continue
        p = os.path.join(data_dir, name)
        if os.path.isdir(p):
            out.append(p)
    out.sort()
    return out


def _parse_plot_index(plot_dir_basename: str) -> int | None:
    """
    Parse trailing integer index from plotfile directory basename.
    Examples: WormholePlt00010 -> 10, plt000123 -> 123
    """
    m = re.search(r"(\d+)$", plot_dir_basename)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def _should_auto_reset(plot_dirs: List[str], state: Dict[str, bool]) -> bool:
    """
    Heuristic: if the output folder contains a "fresh" plot index (0) but the
    saved state references plotfiles that do not exist anymore, assume the user
    restarted a run in the same directory and reset outputs.
    """
    if not plot_dirs or not state:
        return False
    basenames = [os.path.basename(p) for p in plot_dirs]
    cur_set = set(basenames)
    state_set = set(state.keys())

    cur_idxs = [i for i in (_parse_plot_index(b) for b in basenames) if i is not None]
    if not cur_idxs:
        return False
    min_cur = min(cur_idxs)

    # Restart-like scenario: we see index 0 again, but state refers to old plotfiles.
    if min_cur == 0 and not state_set.issubset(cur_set):
        return True

    # Another restart-like scenario: indices restarted and current max is below
    # what we previously processed.
    state_idxs = [i for i in (_parse_plot_index(b) for b in state_set) if i is not None]
    if min_cur == 0 and state_idxs and max(cur_idxs) < max(state_idxs):
        return True

    return False


def _truncate_if_exists(path: Path) -> None:
    if path.exists() and path.is_file() and path.stat().st_size > 0:
        # Truncate content but keep the file path stable
        path.write_text("", encoding="utf-8")


def _is_plotfile_ready(plot_dir: str, stable_seconds: float) -> bool:
    """Best-effort check that the plotfile is not being written right now."""
    header = os.path.join(plot_dir, "Header")
    if not os.path.isfile(header):
        return False
    try:
        mtime = os.path.getmtime(header)
    except OSError:
        return False
    return (time.time() - mtime) >= stable_seconds
