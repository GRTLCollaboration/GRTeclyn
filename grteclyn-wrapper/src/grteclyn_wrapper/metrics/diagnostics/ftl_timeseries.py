"""Time-resolved FTL feature stream from ``small_data/ftl_timeseries.dat``.

One row per plotfile, written in-flight by the plotfile consumer (general FTL
on every frame, gauge-invariant geodesic gated on a coordinate channel).  This
parser turns that stream into per-frame arrays plus convenience aggregates so
the scorer can average a composite FTL x stability score over the whole run
instead of reading only the diffused final frame.
"""

from __future__ import annotations

import math
from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import FtlTimeSeriesMetrics

# A trustworthy gauge-invariant shortcut below this is integration noise, not a
# certified channel; matches the GEO_FTL_FLOOR used by the scorer.
FTL_LIFETIME_FLOOR: float = 1.0e-3


def read_ftl_timeseries_metrics(path: Path) -> FtlTimeSeriesMetrics | None:
    rows = numeric_rows(path, 12)
    if not rows:
        return None

    t = tuple(float(r[0]) for r in rows)
    f_op = tuple(float(r[1]) for r in rows)
    f_geo = tuple(float(r[2]) for r in rows)
    geo_trustworthy = tuple(bool(round(r[3])) for r in rows)
    max_local_speed = tuple(float(r[4]) for r in rows)
    superluminal_fraction = tuple(float(r[5]) for r in rows)
    structure_coherence = tuple(float(r[7]) for r in rows)
    max_h_rel_drift = tuple(float(r[11]) for r in rows)
    n = len(rows)

    # Peak trustworthy gauge-invariant shortcut + when it occurred.
    f_geo_peak = 0.0
    t_at_f_geo_peak: float | None = None
    for ti, fg, trust in zip(t, f_geo, geo_trustworthy):
        if trust and math.isfinite(fg) and fg > f_geo_peak:
            f_geo_peak = fg
            t_at_f_geo_peak = ti

    f_op_peak = 0.0
    t_at_f_op_peak: float | None = None
    for ti, fo in zip(t, f_op):
        if math.isfinite(fo) and fo > f_op_peak:
            f_op_peak = fo
            t_at_f_op_peak = ti

    geo_alive = sum(
        1
        for fg, trust in zip(f_geo, geo_trustworthy)
        if trust and math.isfinite(fg) and fg > FTL_LIFETIME_FLOOR
    )
    op_alive = sum(1 for fo in f_op if math.isfinite(fo) and fo > FTL_LIFETIME_FLOOR)

    return FtlTimeSeriesMetrics(
        n_frames=n,
        t=t,
        f_op=f_op,
        f_geo=f_geo,
        geo_trustworthy=geo_trustworthy,
        max_local_speed=max_local_speed,
        superluminal_fraction=superluminal_fraction,
        structure_coherence=structure_coherence,
        max_h_rel_drift=max_h_rel_drift,
        f_geo_peak=f_geo_peak,
        t_at_f_geo_peak=t_at_f_geo_peak,
        f_op_peak=f_op_peak,
        t_at_f_op_peak=t_at_f_op_peak,
        ftl_lifetime_fraction=(geo_alive / n) if n else 0.0,
        op_lifetime_fraction=(op_alive / n) if n else 0.0,
    )
