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
CUSP_UNRESOLVED_PEAK_SCALE = 0.5

# --- Spacetime-splash (geometric) target scales -----------------------------
# A spacetime splash is a *geometric* event: a converging gravitational wave
# crushes the conformal factor (chi -> 0), spikes the extrinsic-curvature trace
# |K|, and produces a Weyl (Psi4) pulse at the center.  These targets set the
# order of magnitude at which each geometric signature saturates its reward.
GW_WAVE_PEAK_TARGET = 1.0e-2  # |h| proxy or |Re(Psi4)| at center
WEYL4_PEAK_TARGET = GW_WAVE_PEAK_TARGET  # backward-compatible alias
ABS_K_TARGET = 0.35            # |trace K| crunch rate at center (v11-calibrated)
CHI_DROP_TARGET = 0.20         # initial-relative chi drop for a strong well

MIN_CENTER_CHI_WELL_FOR_HORIZON = 0.15  # normalized geometric_curvature_well gate


def _constraint_quality(central) -> float | None:
    ham = central.ham_abs_at_peak
    if ham is None:
        return None
    rho_peak = central.peak_rho_req_at_origin
    if not math.isfinite(rho_peak) or rho_peak <= 0.0:
        return None
    return float(1.0 - min(1.0, float(ham) / rho_peak))


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


def _normalize_to_target(value: float, target: float) -> float:
    if not math.isfinite(value) or value <= 0.0 or target <= 0.0:
        return 0.0
    return float(min(1.0, value / target))


def _compute_geometric_components(ctx: ScoringContext, central) -> None:
    """Spacetime-splash geometry: curvature concentration + GW focusing.

    Independent of matter density (rho).  A static matter blob slowly
    collapsing scores ~0 here; only a genuine converging gravitational wave
    that crushes the geometry at the center earns these rewards.
    """
    if not central.has_geometric_data:
        ctx.components.setdefault("geometric_curvature_well", 0.0)
        ctx.components.setdefault("geometric_wave_arrival", 0.0)
        ctx.components.setdefault("geometric_crunch", 0.0)
        ctx.notes.append("geometric_splash_data_unavailable")
        return

    ctx.components["geometric_curvature_well"] = _normalize_to_target(
        central.chi_drop, CHI_DROP_TARGET
    )
    ctx.components["geometric_wave_arrival"] = _normalize_to_target(
        central.peak_abs_weyl4, GW_WAVE_PEAK_TARGET
    )
    ctx.components["geometric_crunch"] = _normalize_to_target(
        central.peak_abs_K, ABS_K_TARGET
    )


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
        ctx.components.setdefault("pre_collapsed_penalty", 0.0)
        ctx.components.setdefault("horizon_formation_time", 0.0)
        ctx.components.setdefault("constraint_quality", 0.0)
        ctx.components.setdefault("peak_radius", 0.0)
        ctx.components.setdefault("splash_width", 0.0)
        ctx.components.setdefault("compression_ratio", 0.0)
        ctx.components.setdefault("cusp_unresolved", 0.0)
        ctx.components.setdefault("geometric_curvature_well", 0.0)
        ctx.components.setdefault("geometric_wave_arrival", 0.0)
        ctx.components.setdefault("geometric_crunch", 0.0)
        return

    peak_norm = _normalize_rho_peak(central.peak_rho_req_at_origin)
    radial = ctx.metrics.central_radial
    if radial and radial.cusp_unresolved:
        peak_norm *= CUSP_UNRESOLVED_PEAK_SCALE
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
    initial_lapse = central.lapse[0] if central.lapse else 1.0
    if initial_lapse <= INITIAL_LAPSE_COLLAPSE_CEILING:
        pre_penalty -= 0.25
    ctx.components["pre_collapsed_penalty"] = pre_penalty

    if splash_mode == "threshold":
        ctx.components["central_lapse_collapse"] = _threshold_lapse_reward(
            central.min_lapse_at_origin
        )
    else:
        ctx.components["central_lapse_collapse"] = 0.0

    ctx.components["horizon_formation_time"] = 0.0
    _compute_geometric_components(ctx, central)

    collapse = ctx.metrics.collapse
    if collapse and collapse.first_corroborated_time is not None:
        t_horizon = float(collapse.first_corroborated_time)
        horizon_bonus = float(
            max(0.0, 1.0 - t_horizon / max(ctx.target_stop_time or 16.0, 1.0))
        )
        well = float(ctx.components.get("geometric_curvature_well", 0.0))
        if well >= MIN_CENTER_CHI_WELL_FOR_HORIZON:
            ctx.components["horizon_formation_time"] = horizon_bonus
        else:
            ctx.notes.append(
                "horizon_bonus_gated: insufficient_center_chi_drop "
                f"(well={well:.3f} < {MIN_CENTER_CHI_WELL_FOR_HORIZON})"
            )

    constraint = _constraint_quality(central)
    if constraint is None:
        ctx.components["constraint_quality"] = 0.0
        ctx.notes.append("constraint_at_origin_unavailable")
    else:
        ctx.components["constraint_quality"] = constraint

    if radial:
        ctx.components["peak_radius"] = float(radial.peak_radius)
        ctx.components["splash_width"] = float(radial.splash_width)
        ctx.components["compression_ratio"] = float(radial.compression_ratio)
        ctx.components["cusp_unresolved"] = 1.0 if radial.cusp_unresolved else 0.0
    else:
        ctx.components["peak_radius"] = 0.0
        ctx.components["splash_width"] = 0.0
        ctx.components["compression_ratio"] = 0.0
        ctx.components["cusp_unresolved"] = 0.0
