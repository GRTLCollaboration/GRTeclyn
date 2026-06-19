"""Score components for gravitational splash / critical-collapse campaigns."""

from __future__ import annotations

import math

from .types import ScoringContext

# Typical central rho_req scale on 128³ boson-star runs (order-of-magnitude target).
RHO_PEAK_TARGET = 1.0e-2
INITIAL_RHO_COLLAPSE_FLOOR = 5.0e-3
INITIAL_LAPSE_COLLAPSE_CEILING = 0.2
THRESHOLD_LAPSE_LOW = 0.01
THRESHOLD_LAPSE_HIGH = 0.05
FOCUS_EFFICIENCY_CAP = 5.0
DISPERSION_PEAK_GATE = 0.25
DISPERSION_RETENTION_GATE = 0.7


def _normalize_rho_peak(peak: float) -> float:
    if not math.isfinite(peak) or peak <= 0.0:
        return 0.0
    return float(min(1.0, peak / RHO_PEAK_TARGET))


def _threshold_lapse_reward(min_lapse: float) -> float:
    if not math.isfinite(min_lapse):
        return 0.0
    if min_lapse < THRESHOLD_LAPSE_LOW:
        return -1.0
    if min_lapse <= THRESHOLD_LAPSE_HIGH:
        center = 0.5 * (THRESHOLD_LAPSE_LOW + THRESHOLD_LAPSE_HIGH)
        span = THRESHOLD_LAPSE_HIGH - THRESHOLD_LAPSE_LOW
        return float(max(0.0, 1.0 - abs(min_lapse - center) / (0.5 * span)))
    return 0.0


def _activity_retention(activity: tuple[float, ...]) -> float:
    if len(activity) < 2:
        return 1.0
    initial = float(activity[0])
    final = float(activity[-1])
    if not math.isfinite(initial) or initial <= 0.0:
        return 1.0 if final > 0.0 else 0.0
    return float(min(1.0, final / initial))


def _discovery_lapse_progress(min_lapse: float) -> float:
    """Graded reward as central lapse falls toward a collapse band."""
    if not math.isfinite(min_lapse):
        return 0.0
    if min_lapse >= 0.2:
        return 0.0
    if min_lapse <= 0.05:
        return 1.0
    return float((0.2 - min_lapse) / 0.15)


def _wave_focusing_quality(chromaticity: float, activity: tuple[float, ...]) -> float:
    chroma = float(max(0.0, min(1.0, chromaticity)))
    retention = _activity_retention(activity)
    if retention >= DISPERSION_RETENTION_GATE:
        return chroma
    scale = max(0.0, retention / DISPERSION_RETENTION_GATE)
    return float(chroma * scale * 0.5)


def _dispersion_penalty(peak_norm: float, activity: tuple[float, ...]) -> float:
    retention = _activity_retention(activity)
    if peak_norm < DISPERSION_PEAK_GATE and retention < DISPERSION_RETENTION_GATE:
        return -0.5
    if retention < 0.5:
        return -0.25
    return 0.0


def compute_splash_components(ctx: ScoringContext, *, splash_mode: str = "discovery") -> None:
    central = ctx.metrics.central
    if central is None:
        ctx.components.setdefault("central_energy_peak", 0.0)
        ctx.components.setdefault("focusing_efficiency", 0.0)
        ctx.components.setdefault("wave_focusing_quality", 0.0)
        ctx.components.setdefault("collapse_lapse_progress", 0.0)
        ctx.components.setdefault("dispersion_penalty", 0.0)
        ctx.components.setdefault("central_lapse_collapse", 0.0)
        ctx.components.setdefault("pre_collapsed_penalty", 0.0)
        return

    peak_norm = _normalize_rho_peak(central.peak_rho_req_at_origin)
    ctx.components["central_energy_peak"] = peak_norm
    ctx.components["focusing_efficiency"] = float(
        min(FOCUS_EFFICIENCY_CAP, central.focusing_efficiency)
    )
    ctx.components["wave_focusing_quality"] = _wave_focusing_quality(
        central.wave_chromaticity,
        central.scalar_activity,
    )
    ctx.components["collapse_lapse_progress"] = _discovery_lapse_progress(
        central.min_lapse_at_origin
    )
    ctx.components["dispersion_penalty"] = _dispersion_penalty(
        peak_norm,
        central.scalar_activity,
    )

    pre_penalty = 0.0
    if central.initial_rho_req_at_origin >= INITIAL_RHO_COLLAPSE_FLOOR:
        pre_penalty -= 0.5
    if central.min_lapse_at_origin <= INITIAL_LAPSE_COLLAPSE_CEILING:
        pre_penalty -= 0.25
    ctx.components["pre_collapsed_penalty"] = pre_penalty

    if splash_mode == "threshold":
        ctx.components["central_lapse_collapse"] = _threshold_lapse_reward(
            central.min_lapse_at_origin
        )
    else:
        ctx.components["central_lapse_collapse"] = 0.0

    collapse = ctx.metrics.collapse
    if collapse and collapse.first_corroborated_time is not None:
        t_horizon = float(collapse.first_corroborated_time)
        ctx.components["horizon_formation_time"] = float(
            max(0.0, 1.0 - t_horizon / max(ctx.target_stop_time or 16.0, 1.0))
        )
