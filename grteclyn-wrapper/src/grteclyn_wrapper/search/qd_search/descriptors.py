"""Behavior descriptors for MAP-Elites archive cells."""

from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

from ..ftl_peak_metrics import (
    authoritative_geo_lifetime,
    authoritative_geo_strength,
    four_d_geodesic_ran,
    peak_fields_for_descriptor_details,
)

_SPEED_HORIZON_C_FLOOR = 0.9
_SPEED_HORIZON_C_TARGET = 1.3
_SPEED_HORIZON_THETA_SCALE = 0.5

_SPEED_SUPER_C_FLOOR = 0.95
_SPEED_SUPER_C_TARGET = 1.20
_SPEED_SUPER_FRACTION_TARGET = 0.15

# FTL-lifetime descriptor: separate transient shortcuts from sustained ones.
# x-axis = 4D end-to-end geodesic strength; y-axis = 4D shortcut present (1/0).


def _ftl_timeseries(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    ts = (metrics or {}).get("ftl_timeseries")
    return ts if isinstance(ts, Mapping) else None


def _geo_strength_axis(
    metrics: Mapping[str, Any] | None,
    components: Mapping[str, float],
) -> float:
    """MAP-Elites x-axis: 4D evolving ``f_geo`` only (never frozen Cauchy peaks)."""
    return authoritative_geo_strength(metrics, components)


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


# Boson profile-width bounds (splash boson-star search space).
_BS_PROFILE_WIDTH_MIN = 2.0
_BS_PROFILE_WIDTH_MAX = 8.0
_BS_OMEGA_MAX = 0.4
# Scalar shell compact profile bounds (splash scalar-shell campaigns).
_SHELL_WIDTH_MIN = 1.8
_SHELL_WIDTH_MAX = 3.0
_SHELL_OMEGA_MAX = 0.5


def _central_metrics(metrics: Mapping[str, Any] | None) -> Mapping[str, Any] | None:
    central = (metrics or {}).get("central")
    return central if isinstance(central, Mapping) else None


def _normalize_profile_width(width: float) -> float:
    span = _BS_PROFILE_WIDTH_MAX - _BS_PROFILE_WIDTH_MIN
    if span <= 0.0:
        return 0.0
    return float(np.clip((width - _BS_PROFILE_WIDTH_MIN) / span, 0.0, 1.0))


def _normalize_shell_width(width: float) -> float:
    span = _SHELL_WIDTH_MAX - _SHELL_WIDTH_MIN
    if span <= 0.0:
        return 0.0
    return float(np.clip((width - _SHELL_WIDTH_MIN) / span, 0.0, 1.0))


def _shell_omega_chromaticity_fallback(overrides: Mapping[str, Any] | None) -> float:
    if not overrides:
        return 0.0
    omega = overrides.get("grtresna_shell_omega")
    if omega is None or not math.isfinite(float(omega)):
        return 0.0
    return float(np.clip(abs(float(omega)) / _SHELL_OMEGA_MAX, 0.0, 1.0))


def _profile_width_from_overrides(overrides: Mapping[str, Any] | None) -> float | None:
    if not overrides:
        return None
    if overrides.get("grtresna_bs_profile_width") is not None:
        return float(overrides["grtresna_bs_profile_width"])
    if overrides.get("grtresna_shell_width") is not None:
        return float(overrides["grtresna_shell_width"])
    return None


def _omega_chromaticity_fallback(overrides: Mapping[str, Any] | None) -> float:
    if not overrides:
        return 0.0
    omega = overrides.get("grtresna_bs_omega")
    if omega is not None and math.isfinite(float(omega)):
        return float(np.clip(float(omega) / _BS_OMEGA_MAX, 0.0, 1.0))
    return _shell_omega_chromaticity_fallback(overrides)


def _descriptor_details(
    components: Mapping[str, float],
    metrics: Mapping[str, Any] | None = None,
    *,
    mode: str = "legacy",
    overrides: Mapping[str, Any] | None = None,
) -> dict[str, float]:
    if mode == "wave_focusing":
        central = _central_metrics(metrics)
        chromaticity = 0.0
        if central and central.get("wave_chromaticity") is not None:
            chromaticity = float(
                np.clip(float(central.get("wave_chromaticity") or 0.0), 0.0, 1.0)
            )
        elif chromaticity <= 0.0:
            chromaticity = _omega_chromaticity_fallback(overrides)
        width = _profile_width_from_overrides(overrides)
        if width is not None:
            if overrides and overrides.get("grtresna_shell_width") is not None:
                profile_axis = _normalize_shell_width(width)
            else:
                profile_axis = _normalize_profile_width(width)
        else:
            profile_axis = 0.5
        return {
            "x": chromaticity,
            "y": profile_axis,
            "wave_chromaticity": chromaticity,
            "profile_width_norm": profile_axis,
            "grtresna_bs_omega": float(
                overrides.get("grtresna_bs_omega", float("nan"))
            )
            if overrides
            else float("nan"),
            "grtresna_shell_omega": float(
                overrides.get("grtresna_shell_omega", float("nan"))
            )
            if overrides
            else float("nan"),
        }

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
        strength = _geo_strength_axis(metrics, components)
        if four_d_geodesic_ran(metrics):
            lifetime = authoritative_geo_lifetime(metrics, components)
        elif ts:
            lifetime = float(np.clip(float(ts.get("ftl_lifetime_fraction") or 0.0), 0.0, 1.0))
        else:
            lifetime = 0.0
        peak_fields = peak_fields_for_descriptor_details(
            ts, components=components, metrics=metrics
        )
        return {
            "x": strength,
            "y": lifetime,
            "ftl_peak_strength": strength,
            "ftl_lifetime": lifetime,
            "ftl_lifetime_fraction": lifetime,
            **peak_fields,
            "ftl_geo_evolving": float(components.get("ftl_geo_evolving", 0.0)),
        }

    ftl_benefit = float(
        np.clip(
            max(
                components.get("operational_ftl_solved", 0.0),
                components.get("operational_ftl", 0.0),
                components.get("ftl_precursor", 0.0),
                components.get("ftl_shortcut", 0.0),
                components.get("ftl_geo_evolving", 0.0),
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
    overrides: Mapping[str, Any] | None = None,
) -> tuple[float, float]:
    details = _descriptor_details(
        components, metrics, mode=mode, overrides=overrides
    )
    return details["x"], details["y"]


def _bin_index(value: float, bins: int) -> int:
    return int(min(bins - 1, max(0, math.floor(value * bins))))
