"""Gravitational-wave beam score components for the gw_beam objective mode."""

from __future__ import annotations

import math

from .types import ScoringContext


def compute_gw_beam_components(ctx: ScoringContext) -> None:
    """Populate components used by the ``gw_beam`` objective mode.

    The optimizer is rewarded for strong total GW emission and for a
    preferential Z-axis beam (m=±2 power fraction).  A small regularization
    rewards longer-lived emission so the archive does not overfit to a single
    transient spike.
    """
    metrics = ctx.metrics
    if metrics is None or metrics.psi4 is None or not metrics.psi4.has_data:
        ctx.notes.append("gw_beam: no Psi4 data available")
        return

    psi4 = metrics.psi4

    # Base signal: mean total quadrupole power over the run.
    ctx.components["gw_total_power"] = float(psi4.mean_total_power)

    # Directionality: mean fraction of power in the Z-axis (m=±2) modes.
    ctx.components["gw_beam_ratio"] = float(psi4.mean_beam_ratio)

    # Combined quality: strong signal that is also beamed.
    # Log-scaling the power avoids a single runaway amplitude swamping the
    # beaming term entirely.
    power = float(psi4.mean_total_power)
    if power > 0.0:
        log_power = math.log10(power + 1.0)
    else:
        log_power = 0.0
    ctx.components["gw_beam_quality"] = float(log_power * (1.0 + psi4.mean_beam_ratio))

    # Peak capture so short-lived but intense bursts still get credit.
    ctx.components["gw_peak_power"] = float(psi4.peak_total_power)

    ctx.notes.append(
        f"gw_beam: mean_total_power={psi4.mean_total_power:.3e} "
        f"mean_beam_ratio={psi4.mean_beam_ratio:.3f}"
    )


def gw_beam_total(components: dict[str, float]) -> float:
    """Scalarize gw_beam components into a single MAP-Elites quality score.

    Primary reward is strong directional GW emission.  Health and physics
    penalties are applied on top so the optimizer prefers *sustained, stable,
    non-exotic* emitters over violent single-tick collapses that happen to
    produce a strong transient burst.

    Weight calibration (from observed compact-Q-ball data):
      * GW reward terms total ~3–10 points for a typical emitter.
      * ``exotic_penalty`` is -1.6 across the board (geometric NEC violation
        from the constrained back-solve).  Since it is near-constant it acts
        mainly as a constant offset; a 1x weight (-1.6 pts) is enough to
        matter without swamping the GW reward.
      * ``horizon_penalty`` is -1.0 for ~97% of candidates at this resolution.
        A 1x weight (-1 pt) provides a mild preference for horizon-free
        geometries when they appear, without flooring the population.
      * ``survival`` (0–0.98) and ``stability`` (0–0.11) vary meaningfully;
        5x and 3x weights make a stable, long-lived emitter worth ~3–5 pts
        extra -- a real but not dominant incentive.
      * ``constraint_health`` (0–0.25) gets 2x so numerically trustworthy
        solutions get a small edge.
    """
    return (
        # --- GW emission reward (primary) ---
        1000.0 * components.get("gw_beam_quality", 0.0)
        + 100.0 * components.get("gw_peak_power", 0.0)
        # --- Health rewards (sustained, stable emission) ---
        + 5.0 * components.get("survival", 0.0)
        + 3.0 * components.get("stability", 0.0)
        + 2.0 * components.get("constraint_health", 0.0)
        # --- Physics penalties (avoid pathological geometries) ---
        + 1.0 * components.get("horizon_penalty", 0.0)
        + 1.0 * components.get("exotic_penalty", 0.0)
        + 1.0 * components.get("instability_penalty", 0.0)
    )
