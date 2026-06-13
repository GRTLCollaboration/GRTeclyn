"""Extract time-resolved FTL peak fields for retention and reporting."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any, Mapping

GEO_FTL_FLOOR = 1.0e-3


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
    n_frames: int = 0


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def from_timeseries_mapping(
    ts: Mapping[str, Any] | None,
    *,
    components: Mapping[str, Any] | None = None,
) -> FtlPeakMetrics | None:
    """Build peak metrics from ``metrics.ftl_timeseries`` or equivalent dict."""
    if not isinstance(ts, Mapping) or not ts.get("n_frames"):
        return None

    n_frames = int(ts.get("n_frames") or 0)
    components = components or {}
    return FtlPeakMetrics(
        f_geo_peak=_optional_float(ts.get("f_geo_peak")) or 0.0,
        t_at_f_geo_peak=_optional_float(ts.get("t_at_f_geo_peak")),
        f_op_peak=_optional_float(ts.get("f_op_peak")) or 0.0,
        t_at_f_op_peak=_optional_float(ts.get("t_at_f_op_peak")),
        max_local_speed_peak=_optional_float(ts.get("max_local_speed_peak")) or 0.0,
        t_at_max_speed=_optional_float(ts.get("t_at_max_speed")),
        superluminal_fraction_peak=_optional_float(ts.get("superluminal_fraction_peak")) or 0.0,
        t_at_superluminal_peak=_optional_float(ts.get("t_at_superluminal_peak")),
        ftl_lifetime_fraction=_optional_float(ts.get("ftl_lifetime_fraction")) or 0.0,
        ftl_geo_timeavg=_optional_float(components.get("ftl_geo_timeavg")) or 0.0,
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
        ftl_geo_timeavg=_optional_float(details.get("ftl_geo_timeavg"))
        or _optional_float(components.get("ftl_geo_timeavg"))
        or 0.0,
        n_frames=int(_optional_float(details.get("n_frames")) or 0),
    )


def peak_fields_for_descriptor_details(
    ts: Mapping[str, Any] | None,
    *,
    components: Mapping[str, float],
) -> dict[str, float]:
    """Serialize peak fields into ``descriptor_details`` for trajectory retention."""
    peaks = from_timeseries_mapping(ts, components=components)
    if peaks is None:
        return {
            "f_geo_peak": float("nan"),
            "t_at_f_geo_peak": float("nan"),
            "f_op_peak": float("nan"),
            "t_at_f_op_peak": float("nan"),
            "max_local_speed_peak": float("nan"),
            "t_at_max_speed": float("nan"),
            "superluminal_fraction_peak": float("nan"),
            "t_at_superluminal_peak": float("nan"),
            "n_frames": 0.0,
            "ftl_geo_timeavg": float(components.get("ftl_geo_timeavg", 0.0)),
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
    }
