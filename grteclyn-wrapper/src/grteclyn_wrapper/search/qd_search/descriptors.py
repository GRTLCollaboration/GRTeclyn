"""Behavior descriptors for MAP-Elites archive cells."""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

_SPEED_HORIZON_C_FLOOR = 0.9
_SPEED_HORIZON_C_TARGET = 1.3
_SPEED_HORIZON_THETA_SCALE = 0.5

_SPEED_SUPER_C_FLOOR = 0.95
_SPEED_SUPER_C_TARGET = 1.20
_SPEED_SUPER_FRACTION_TARGET = 0.15

# FTL-lifetime descriptor: separate transient shortcuts from sustained ones.
# x-axis = peak gauge-invariant strength (same floor/target as the scorer),
# y-axis = fraction of the run the shortcut is alive.  An Alcubierre-like
# one-frame spike lands in a low-lifetime cell; a stable warp in a high one.
_GEO_PEAK_FLOOR = 1.0e-3
_GEO_PEAK_TARGET = 2.0e-1


def _ftl_timeseries(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    ts = (metrics or {}).get("ftl_timeseries")
    return ts if isinstance(ts, Mapping) else None


def _geo_peak_axis(ts: Mapping[str, Any] | None) -> float:
    if not ts:
        return 0.0
    peak = ts.get("f_geo_peak")
    if peak is None or not math.isfinite(float(peak)) or float(peak) <= _GEO_PEAK_FLOOR:
        return 0.0
    span = _GEO_PEAK_TARGET - _GEO_PEAK_FLOOR
    return float(np.clip((float(peak) - _GEO_PEAK_FLOOR) / span, 0.0, 1.0))


def _path_closeness_from_report(report: Mapping[str, Any] | None) -> float:
    if not report or not bool(report.get("reachable", False)):
        return 0.0
    t_min = report.get("t_min")
    t_flat = float(report.get("t_flat") or 0.0)
    if t_min is None or not math.isfinite(float(t_min)) or t_flat <= 0.0:
        return 0.0
    t_min_f = float(t_min)
    if t_min_f <= t_flat:
        return 1.0
    excess = (t_min_f - t_flat) / t_flat
    return max(0.0, 1.0 - excess / 0.12)


def _best_ftl_report(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    metrics = metrics or {}
    return (
        metrics.get("general_ftl_evolved")
        or metrics.get("general_ftl_solved")
        or metrics.get("general_ftl")
    )


def _solved_ftl_report(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    metrics = metrics or {}
    return (
        metrics.get("general_ftl_solved")
        or metrics.get("general_ftl_evolved")
        or metrics.get("general_ftl")
    )


def _speed_tilt_axis(
    report: Mapping[str, Any] | None,
    *,
    floor: float,
    target: float,
) -> float:
    if not report:
        return 0.0
    c = float(report.get("max_local_speed") or 0.0)
    if not math.isfinite(c) or c <= floor:
        return 0.0
    span = target - floor
    return float(np.clip((c - floor) / span, 0.0, 1.0))


def _speed_axis_from_report(report: Mapping[str, Any] | None) -> float:
    return _speed_tilt_axis(
        report, floor=_SPEED_HORIZON_C_FLOOR, target=_SPEED_HORIZON_C_TARGET
    )


def _superluminal_axis(
    report: Mapping[str, Any] | None,
    *,
    target: float = 1.0,
) -> float:
    if not report:
        return 0.0
    frac = report.get("superluminal_fraction")
    if frac is None or not math.isfinite(float(frac)):
        return 0.0
    span = target if target > 0.0 else 1.0
    return float(np.clip(float(frac) / span, 0.0, 1.0))


def _horizon_free_axis(metrics: Mapping[str, Any] | None) -> float:
    metrics = metrics or {}
    collapse = metrics.get("collapse") or {}
    theta = collapse.get("min_theta_plus")
    if theta is None or not math.isfinite(float(theta)):
        return 0.5
    scaled = 0.5 + float(theta) / (2.0 * _SPEED_HORIZON_THETA_SCALE)
    return float(np.clip(scaled, 0.0, 1.0))


def _descriptor_details(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
) -> dict[str, float]:
    if mode == "speed_horizon":
        report = _best_ftl_report(metrics)
        speed = _speed_axis_from_report(report)
        horizon_free = _horizon_free_axis(metrics)
        c = float(report.get("max_local_speed")) if report and report.get("max_local_speed") is not None else float("nan")
        collapse = (metrics or {}).get("collapse") or {}
        theta = float(collapse.get("min_theta_plus")) if collapse.get("min_theta_plus") is not None else float("nan")
        return {
            "x": speed,
            "y": horizon_free,
            "speed_tilt": speed,
            "horizon_free": horizon_free,
            "max_local_speed": c,
            "min_theta_plus": theta,
            "ftl_persistence": float(components.get("ftl_persistence", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
        }

    if mode == "speed_super":
        report = _solved_ftl_report(metrics)
        speed = _speed_tilt_axis(
            report, floor=_SPEED_SUPER_C_FLOOR, target=_SPEED_SUPER_C_TARGET
        )
        super_frac = _superluminal_axis(
            report, target=_SPEED_SUPER_FRACTION_TARGET
        )
        raw_frac = (
            float(report.get("superluminal_fraction"))
            if report and report.get("superluminal_fraction") is not None
            else float("nan")
        )
        c = float(report.get("max_local_speed")) if report and report.get("max_local_speed") is not None else float("nan")
        return {
            "x": speed,
            "y": super_frac,
            "speed_tilt": speed,
            "superluminal_fraction": super_frac,
            "superluminal_fraction_raw": raw_frac,
            "max_local_speed": c,
            "ftl_persistence": float(components.get("ftl_persistence", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
        }

    if mode == "channel":
        report = _best_ftl_report(metrics)
        path_closeness = _path_closeness_from_report(report)
        precursor = float(np.clip(components.get("ftl_precursor", 0.0), 0.0, 1.0))
        shift = float(np.clip(components.get("shift_drive", 0.0), 0.0, 1.0))
        mechanism_balance = float(math.sqrt(max(precursor, 0.0) * max(shift, 0.0)))
        t_min = float(report.get("t_min")) if report and report.get("t_min") is not None else float("nan")
        t_flat = float(report.get("t_flat")) if report and report.get("t_flat") is not None else float("nan")
        return {
            "x": path_closeness,
            "y": mechanism_balance,
            "path_closeness": path_closeness,
            "mechanism_balance": mechanism_balance,
            "shift_drive": shift,
            "ftl_precursor": precursor,
            "channel_progress": float(components.get("channel_progress", 0.0)),
            "operational_ftl": float(components.get("operational_ftl", 0.0)),
            "t_min": t_min,
            "t_flat": t_flat,
        }

    if mode == "ftl_lifetime":
        ts = _ftl_timeseries(metrics)
        strength = _geo_peak_axis(ts)
        lifetime = (
            float(np.clip(float(ts.get("ftl_lifetime_fraction") or 0.0), 0.0, 1.0))
            if ts
            else 0.0
        )
        peak = float(ts.get("f_geo_peak")) if ts and ts.get("f_geo_peak") is not None else float("nan")
        t_peak = float(ts.get("t_at_f_geo_peak")) if ts and ts.get("t_at_f_geo_peak") is not None else float("nan")
        n_frames = int(ts.get("n_frames")) if ts and ts.get("n_frames") is not None else 0
        return {
            "x": strength,
            "y": lifetime,
            "ftl_peak_strength": strength,
            "ftl_lifetime": lifetime,
            "f_geo_peak": peak,
            "t_at_f_geo_peak": t_peak,
            "n_frames": float(n_frames),
            "operational_ftl_geodesic": float(components.get("operational_ftl_geodesic", 0.0)),
            "ftl_geo_timeavg": float(components.get("ftl_geo_timeavg", 0.0)),
        }

    ftl_benefit = float(
        np.clip(
            max(
                components.get("operational_ftl_solved", 0.0),
                components.get("operational_ftl", 0.0),
                components.get("ftl_precursor", 0.0),
                components.get("ftl_shortcut", 0.0),
            ),
            0.0,
            1.0,
        )
    )
    mechanism = float(np.clip(components.get("mechanism_descriptor", 0.0), 0.0, 1.0))
    if mechanism <= 0.0:
        mechanism = float(np.clip(1.0 - components.get("anec_condition", 1.0), 0.0, 1.0))
    return {
        "x": ftl_benefit,
        "y": mechanism,
        "ftl_benefit": ftl_benefit,
        "mechanism": mechanism,
    }


def _descriptors(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
) -> tuple[float, float]:
    details = _descriptor_details(components, metrics, mode=mode)
    return details["x"], details["y"]


def _bin_index(value: float, bins: int) -> int:
    return int(min(bins - 1, max(0, math.floor(value * bins))))
