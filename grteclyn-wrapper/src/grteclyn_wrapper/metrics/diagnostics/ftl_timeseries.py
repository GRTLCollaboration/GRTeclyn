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


def _peak_with_time(
    values: tuple[float, ...],
    times: tuple[float, ...],
) -> tuple[float, float | None]:
    peak = float("-inf")
    t_peak: float | None = None
    for value, time in zip(values, times):
        if not math.isfinite(value):
            continue
        if value > peak:
            peak = value
            t_peak = time
    if peak == float("-inf"):
        return 0.0, None
    return peak, t_peak


def _aggregate_ftl_frames(
    *,
    t: tuple[float, ...],
    f_op: tuple[float, ...],
    f_geo: tuple[float, ...],
    geo_trustworthy: tuple[bool, ...],
    max_local_speed: tuple[float, ...],
    superluminal_fraction: tuple[float, ...],
    structure_coherence: tuple[float, ...],
    max_h_rel_drift: tuple[float, ...],
) -> FtlTimeSeriesMetrics:
    """Build metrics by scanning every per-frame sample (peaks are mid-run aware)."""
    n = len(t)

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

    speed_peak, t_speed = _peak_with_time(max_local_speed, t)
    sup_peak, t_sup = _peak_with_time(superluminal_fraction, t)

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
        max_local_speed_peak=speed_peak,
        t_at_max_speed=t_speed,
        superluminal_fraction_peak=sup_peak,
        t_at_superluminal_peak=t_sup,
        ftl_lifetime_fraction=(geo_alive / n) if n else 0.0,
        op_lifetime_fraction=(op_alive / n) if n else 0.0,
    )


def reaggregate_ftl_timeseries(ts: FtlTimeSeriesMetrics) -> FtlTimeSeriesMetrics:
    """Recompute peak/lifetime fields from per-frame arrays (never final-frame only)."""
    return _aggregate_ftl_frames(
        t=ts.t,
        f_op=ts.f_op,
        f_geo=ts.f_geo,
        geo_trustworthy=ts.geo_trustworthy,
        max_local_speed=ts.max_local_speed,
        superluminal_fraction=ts.superluminal_fraction,
        structure_coherence=ts.structure_coherence,
        max_h_rel_drift=ts.max_h_rel_drift,
    )


def read_ftl_timeseries_metrics(path: Path) -> FtlTimeSeriesMetrics | None:
    rows = numeric_rows(path, 12)
    if not rows:
        return None

    f_geo_evol: float | None = None
    f_geo_evol_ok: bool | None = None
    for row in reversed(rows):
        if len(row) >= 14:
            f_geo_evol = float(row[12])
            f_geo_evol_ok = bool(round(row[13]))
            break

    metrics = _aggregate_ftl_frames(
        t=tuple(float(r[0]) for r in rows),
        f_op=tuple(float(r[1]) for r in rows),
        f_geo=tuple(float(r[2]) for r in rows),
        geo_trustworthy=tuple(bool(round(r[3])) for r in rows),
        max_local_speed=tuple(float(r[4]) for r in rows),
        superluminal_fraction=tuple(float(r[5]) for r in rows),
        structure_coherence=tuple(float(r[7]) for r in rows),
        max_h_rel_drift=tuple(float(r[11]) for r in rows),
    )
    if f_geo_evol is None:
        return metrics
    return FtlTimeSeriesMetrics(
        **{
            field.name: getattr(metrics, field.name)
            for field in FtlTimeSeriesMetrics.__dataclass_fields__.values()
            if field.name not in {"f_geo_evol", "f_geo_evol_ok"}
        },
        f_geo_evol=f_geo_evol,
        f_geo_evol_ok=f_geo_evol_ok,
    )
