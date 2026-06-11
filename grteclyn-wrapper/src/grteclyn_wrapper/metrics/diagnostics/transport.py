"""Transport objective from barycenter columns in collapse diagnostics."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io.dat import numeric_rows
from ..types.diagnostics import TransportMetrics


def read_transport_metrics(collapse_path: Path) -> TransportMetrics | None:
    rows = numeric_rows(collapse_path, 15)
    if len(rows) < 2:
        return None
    x0 = rows[0][14] if len(rows[0]) >= 15 else None
    xf = rows[-1][14] if len(rows[-1]) >= 15 else None
    if x0 is None or xf is None:
        return None
    translation = float(xf - x0)
    xs = [row[14] for row in rows if len(row) >= 15]
    deformation = float(np.std(xs)) if len(xs) >= 2 else 0.0
    score = max(0.0, translation) / (1.0 + deformation)
    return TransportMetrics(
        initial_barycenter_x=x0,
        final_barycenter_x=xf,
        translation=translation,
        deformation=deformation,
        score=min(score, 1.0),
    )
