"""Extract time-resolved FTL peak fields for retention and reporting."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

GEO_FTL_FLOOR = 1.0e-3
GEO_FTL_TARGET = 2.0e-1


@dataclass(frozen=True)
class FtlPeakMetrics:
    """Peak FTL observables used for champion retention and trajectory fallback."""

    f_geo_peak: float
    t_at_f_geo_peak: float | None
    f_op_peak: float
    t_at_f_op_peak: float | None
    max_local_speed_peak: float
    t_at_max_speed: float | None
    superluminal_fraction_peak: float
    t_at_superluminal_peak: float | None
    ftl_lifetime_fraction: float
    ftl_geo_timeavg: float
    ftl_geo_evolving: float = 0.0
    f_geo_evol: float = 0.0
    n_frames: int = 0


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _evolving_geodesic(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    evo = (metrics or {}).get("evolving_geodesic")
    return evo if isinstance(evo, Mapping) else None


def four_d_geodesic_ran(metrics: Mapping[str, Any] | None) -> bool:
    """True when the end-of-run 4D evolving probe produced a report."""
    return _evolving_geodesic(metrics) is not None


def four_d_geodesic_trustworthy(evo: Mapping[str, Any] | None) -> bool:
    if evo is None:
        return False
    n_rays = int(evo.get("n_rays") or 0)
    n_reached = int(evo.get("n_reached") or 0)
    return bool(evo.get("h_quality_ok")) and n_rays > 0 and n_reached == n_rays


def authoritative_geo_strength(
    metrics: Mapping[str, Any] | None,
    components: Mapping[str, float],
) -> float:
    """Normalized MAP-Elites / retention strength from the 4D evolving probe only."""
    evo = _evolving_geodesic(metrics)
    if evo is not None:
        f_geo = evo.get("f_geo")
        if (
            four_d_geodesic_trustworthy(evo)
            and f_geo is not None
            and math.isfinite(float(f_geo))
            and float(f_geo) > GEO_FTL_FLOOR
        ):
            span = GEO_FTL_TARGET - GEO_FTL_FLOOR
            return float(np.clip((float(f_geo) - GEO_FTL_FLOOR) / span, 0.0, 1.0))
        return 0.0

    evo_comp = float(components.get("ftl_geo_evolving", 0.0))
    if evo_comp > 0.0:
        return float(min(max(evo_comp, 0.0), 1.0))
    return 0.0


def authoritative_geo_lifetime(
    metrics: Mapping[str, Any] | None,
    components: Mapping[str, float],
) -> float:
    """MAP-Elites y-axis: end-to-end 4D shortcut present, else zero."""
    return 1.0 if authoritative_geo_strength(metrics, components) > 0.0 else 0.0


def _raw_f_geo_evol(metrics: Mapping[str, Any] | None) -> float:
    evo = _evolving_geodesic(metrics)
    if evo is None:
        return 0.0
    f_geo = _optional_float(evo.get("f_geo"))
    return f_geo if f_geo is not None else 0.0


def from_timeseries_mapping(
    ts: Mapping[str, Any] | None,
    *,
    components: Mapping[str, Any] | None = None,
    metrics: Mapping[str, Any] | None = None,
) -> FtlPeakMetrics | None:
    """Build peak metrics from timeseries + authoritative 4D geodesic fields."""
    if not isinstance(ts, Mapping) or not ts.get("n_frames"):
        if not four_d_geodesic_ran(metrics):
            return None
        n_frames = 0
    else:
        n_frames = int(ts.get("n_frames") or 0)

    components = components or {}
    f_geo_evol = _raw_f_geo_evol(metrics)
    ftl_geo_evolving = float(components.get("ftl_geo_evolving", 0.0))
    if four_d_geodesic_ran(metrics):
        f_geo_peak = f_geo_evol
        ftl_geo_timeavg = 0.0
        ftl_lifetime_fraction = authoritative_geo_lifetime(metrics, components)
    else:
        f_geo_peak = _optional_float(ts.get("f_geo_peak")) if ts else 0.0
        f_geo_peak = f_geo_peak or 0.0
        ftl_geo_timeavg = _optional_float(components.get("ftl_geo_timeavg")) or 0.0
        ftl_lifetime_fraction = (
            _optional_float(ts.get("ftl_lifetime_fraction")) or 0.0 if ts else 0.0
        )

    return FtlPeakMetrics(
        f_geo_peak=f_geo_peak,
        t_at_f_geo_peak=_optional_float(ts.get("t_at_f_geo_peak")) if ts else None,
        f_op_peak=_optional_float(ts.get("f_op_peak")) or 0.0 if ts else 0.0,
        t_at_f_op_peak=_optional_float(ts.get("t_at_f_op_peak")) if ts else None,
        max_local_speed_peak=_optional_float(ts.get("max_local_speed_peak")) or 0.0 if ts else 0.0,
        t_at_max_speed=_optional_float(ts.get("t_at_max_speed")) if ts else None,
        superluminal_fraction_peak=_optional_float(ts.get("superluminal_fraction_peak")) or 0.0 if ts else 0.0,
        t_at_superluminal_peak=_optional_float(ts.get("t_at_superluminal_peak")) if ts else None,
        ftl_lifetime_fraction=ftl_lifetime_fraction,
        ftl_geo_timeavg=ftl_geo_timeavg,
        ftl_geo_evolving=ftl_geo_evolving,
        f_geo_evol=f_geo_evol,
        n_frames=n_frames,
    )


def from_trajectory_record(record: Mapping[str, Any]) -> FtlPeakMetrics | None:
    """Read peak metrics stored on a trajectory record."""
    if record.get("status") != "gpu_ok":
        return None
    details = record.get("descriptor_details") or {}
    components = record.get("components") or {}
    if not details and not components:
        return None
    return FtlPeakMetrics(
        f_geo_peak=_optional_float(details.get("f_geo_peak")) or 0.0,
        t_at_f_geo_peak=_optional_float(details.get("t_at_f_geo_peak")),
        f_op_peak=_optional_float(details.get("f_op_peak")) or 0.0,
        t_at_f_op_peak=_optional_float(details.get("t_at_f_op_peak")),
        max_local_speed_peak=_optional_float(details.get("max_local_speed_peak")) or 0.0,
        t_at_max_speed=_optional_float(details.get("t_at_max_speed")),
        superluminal_fraction_peak=_optional_float(details.get("superluminal_fraction_peak")) or 0.0,
        t_at_superluminal_peak=_optional_float(details.get("t_at_superluminal_peak")),
        ftl_lifetime_fraction=_optional_float(details.get("ftl_lifetime"))
        or _optional_float(details.get("ftl_lifetime_fraction"))
        or 0.0,
        ftl_geo_timeavg=_optional_float(details.get("ftl_geo_timeavg")) or 0.0,
        ftl_geo_evolving=_optional_float(details.get("ftl_geo_evolving"))
        or _optional_float(components.get("ftl_geo_evolving"))
        or 0.0,
        f_geo_evol=_optional_float(details.get("f_geo_evol")) or 0.0,
        n_frames=int(_optional_float(details.get("n_frames")) or 0),
    )


def peak_fields_for_descriptor_details(
    ts: Mapping[str, Any] | None,
    *,
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
) -> dict[str, float]:
    """Serialize peak fields into ``descriptor_details`` for trajectory retention."""
    peaks = from_timeseries_mapping(ts, components=components, metrics=metrics)
    ftl_geo_evolving = float(components.get("ftl_geo_evolving", 0.0))
    f_geo_evol = _raw_f_geo_evol(metrics)

    if peaks is None:
        return {
            "f_geo_peak": 0.0,
            "t_at_f_geo_peak": float("nan"),
            "f_op_peak": float("nan"),
            "t_at_f_op_peak": float("nan"),
            "max_local_speed_peak": float("nan"),
            "t_at_max_speed": float("nan"),
            "superluminal_fraction_peak": float("nan"),
            "t_at_superluminal_peak": float("nan"),
            "n_frames": 0.0,
            "ftl_geo_timeavg": 0.0,
            "ftl_geo_evolving": ftl_geo_evolving,
            "f_geo_evol": f_geo_evol,
        }

    def _time(value: float | None) -> float:
        return float("nan") if value is None else float(value)

    return {
        "f_geo_peak": peaks.f_geo_peak,
        "t_at_f_geo_peak": _time(peaks.t_at_f_geo_peak),
        "f_op_peak": peaks.f_op_peak,
        "t_at_f_op_peak": _time(peaks.t_at_f_op_peak),
        "max_local_speed_peak": peaks.max_local_speed_peak,
        "t_at_max_speed": _time(peaks.t_at_max_speed),
        "superluminal_fraction_peak": peaks.superluminal_fraction_peak,
        "t_at_superluminal_peak": _time(peaks.t_at_superluminal_peak),
        "n_frames": float(peaks.n_frames),
        "ftl_geo_timeavg": peaks.ftl_geo_timeavg,
        "ftl_geo_evolving": peaks.ftl_geo_evolving,
        "f_geo_evol": peaks.f_geo_evol,
    }
