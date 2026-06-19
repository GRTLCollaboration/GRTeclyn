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


def compute_splash_components(ctx: ScoringContext, *, splash_mode: str = "discovery") -> None:
    central = ctx.metrics.central
    if central is None:
        ctx.components.setdefault("central_energy_peak", 0.0)
        ctx.components.setdefault("focusing_efficiency", 0.0)
        ctx.components.setdefault("central_lapse_collapse", 0.0)
        ctx.components.setdefault("pre_collapsed_penalty", 0.0)
        return

    peak_norm = _normalize_rho_peak(central.peak_rho_req_at_origin)
    ctx.components["central_energy_peak"] = peak_norm
    ctx.components["focusing_efficiency"] = float(
        min(10.0, central.focusing_efficiency)
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
