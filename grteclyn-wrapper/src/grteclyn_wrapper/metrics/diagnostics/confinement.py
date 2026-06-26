"""Matter-confinement diagnostics from ``confinement.dat``.

The trustworthy "matter dispersed / flew away" detector.  See the consumer-side
extractor (visualisation/.../extraction/confinement.py) for how the per-plotfile
moments are computed; this module just reads the time series back and reduces it
to the aggregates the scorer needs.
"""

from __future__ import annotations

import math
from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import ConfinementMetrics

# confinement.dat columns (0-indexed):
#   0 time  1 total  2 peak  3 rms_radius  4 confined_frac
#   5 bary_x  6 bary_y  7 bary_z  8 r_conf
_RMS_COL = 3
_FRAC_COL = 4
_TOTAL_COL = 1


def _finite(value: float) -> float | None:
    return value if math.isfinite(value) else None


def read_confinement_metrics(path: Path) -> ConfinementMetrics | None:
    rows = numeric_rows(path, 5)
    if not rows:
        return None

    rms = [r[_RMS_COL] for r in rows]
    frac = [r[_FRAC_COL] for r in rows]
    total = [r[_TOTAL_COL] for r in rows]

    initial_rms = _finite(rms[0])
    final_rms = _finite(rms[-1])
    max_rms = _finite(max(rms)) if rms else None
    initial_frac = _finite(frac[0])
    final_frac = _finite(frac[-1])
    min_frac = _finite(min(frac)) if frac else None

    spread_ratio: float | None = None
    if initial_rms is not None and final_rms is not None and initial_rms > 1.0e-12:
        spread_ratio = final_rms / initial_rms

    return ConfinementMetrics(
        n_frames=len(rows),
        final_time=_finite(rows[-1][0]),
        initial_rms_radius=initial_rms,
        final_rms_radius=final_rms,
        max_rms_radius=max_rms,
        initial_confined_frac=initial_frac,
        final_confined_frac=final_frac,
        min_confined_frac=min_frac,
        spread_ratio=spread_ratio,
        initial_total=_finite(total[0]),
        final_total=_finite(total[-1]),
    )
