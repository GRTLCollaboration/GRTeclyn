"""Parse GRTresna Ham/Mom convergence diagnostics."""

from __future__ import annotations

import math
from pathlib import Path


def parse_convergence(work_dir: Path) -> dict[str, float] | None:
    """Read the final Ham/Mom errors from the GRTresna error file."""
    err_file = work_dir / "Ham_and_Mom_errors.txt"
    if not err_file.exists():
        return None
    lines = err_file.read_text().strip().splitlines()
    if len(lines) < 2:
        return None
    last = lines[-1].split()
    if len(last) < 3:
        return None
    try:
        ham_pct = float(last[1])
        mom_pct = float(last[2])
    except (ValueError, IndexError):
        return None
    if math.isfinite(ham_pct) and not math.isfinite(mom_pct):
        mom_pct = 0.0
    return {
        "iteration": int(last[0]),
        "ham_pct": ham_pct,
        "mom_pct": mom_pct,
    }
