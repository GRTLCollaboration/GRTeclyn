"""Curvature invariant diagnostics from ``curvature_invariants.dat``."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import CurvatureInvariantMetrics


def read_curvature_invariant_metrics(path: Path) -> CurvatureInvariantMetrics | None:
    rows = numeric_rows(path, 5)
    if not rows:
        return None
    return CurvatureInvariantMetrics(
        final_time=rows[-1][0],
        max_abs_ricci_scalar=max(row[1] for row in rows),
        max_ricci_tensor_sq=max(row[2] for row in rows),
        max_kij_sq=max(row[3] for row in rows),
        max_l2_ricci_scalar=max(row[4] for row in rows),
    )
