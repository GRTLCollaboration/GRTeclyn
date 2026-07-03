"""Gravitational-wave beam score components for the gw_beam objective mode."""

from __future__ import annotations

import math
import os

from .types import ScoringContext

# Multiplicative veto thresholds (collapse / reward-hacking guard).
GW_FATAL_HAM_L2 = 100.0
GW_ORBIT_RADIUS = 5.0  # typical lump orbit radius (geometric units, c=1)
GW_OBSERVATION_MARGIN = 3.0  # buffer after light-travel time to extraction sphere
GW_DEFAULT_EXTRACTION_RADIUS = 12.0

# Secondary additive demotion (cannot fight unbounded exploits alone).
CONSTRAINT_SPIKE_PENALTY_WEIGHT = 150.0


def gw_min_valid_observation_time(
    extraction_radius: float | None = None,
) -> float:
    """Earliest constraint-spike time that still leaves a clean GW observing window."""
    if extraction_radius is None:
        raw = os.environ.get("CONSUMER_RADII", "").strip().split()
        try:
            extraction_radius = float(raw[0]) if raw else GW_DEFAULT_EXTRACTION_RADIUS
        except ValueError:
            extraction_radius = GW_DEFAULT_EXTRACTION_RADIUS
    return max(0.0, extraction_radius - GW_ORBIT_RADIUS) + GW_OBSERVATION_MARGIN


def gw_health_multiplier(metrics) -> float:
    """Hard gate: 0 when constraints blow up or collapse precedes the wave zone."""
    constraints = metrics.constraints if metrics is not None else None
    if constraints is None:
        return 1.0
    max_ham = constraints.max_hamiltonian_l2 or 0.0
    if max_ham > GW_FATAL_HAM_L2:
        return 0.0
    if constraints.has_constraint_spike:
        spike_t = constraints.constraint_spike_time
        min_t = gw_min_valid_observation_time()
        if spike_t is not None and spike_t < min_t:
            return 0.0
    return 1.0


def compute_gw_beam_components(ctx: ScoringContext) -> None:
    """Populate components used by the ``gw_beam`` objective mode.

    Psi4 aggregates are truncated at the constraint-spike time when present
    (see ``read_psi4_metrics(max_peak_time=...)``).  A multiplicative
    ``gw_health_multiplier`` zeros all GW reward when Ham explodes or the
    spike occurs before the extraction sphere has a clean observing window.
    """
    metrics = ctx.metrics
    if metrics is None or metrics.psi4 is None or not metrics.psi4.has_data:
        ctx.notes.append("gw_beam: no Psi4 data available")
        return

    psi4 = metrics.psi4
    mult = gw_health_multiplier(metrics)
    ctx.components["gw_health_multiplier"] = float(mult)

    ctx.components["gw_total_power"] = float(psi4.mean_total_power)
    ctx.components["gw_beam_ratio"] = float(psi4.mean_beam_ratio)

    power = float(psi4.mean_total_power)
    if power > 0.0:
        log_power = math.log10(power + 1.0)
    else:
        log_power = 0.0
    ctx.components["gw_beam_quality"] = float(log_power * (1.0 + psi4.mean_beam_ratio))
    ctx.components["gw_peak_power"] = float(psi4.peak_total_power)

    if mult == 0.0:
        constraints = metrics.constraints
        if constraints and constraints.has_constraint_spike:
            t = constraints.constraint_spike_time
            ctx.notes.append(
                f"gw_beam: GW reward vetoed (spike t={t:.2f}, "
                f"max_Ham={constraints.max_hamiltonian_l2:.2e})"
            )
        else:
            ctx.notes.append("gw_beam: GW reward vetoed (fatal constraint violation)")
    else:
        ctx.notes.append(
            f"gw_beam: mean_total_power={psi4.mean_total_power:.3e} "
            f"mean_beam_ratio={psi4.mean_beam_ratio:.3f}"
        )


def gw_beam_total(components: dict[str, float]) -> float:
    """Scalarize gw_beam components into a single MAP-Elites quality score.

    GW emission terms are multiplied by ``gw_health_multiplier`` (hard veto on
    collapse / numerical explosions).  Health bonuses share the same gate.
    Physics penalties remain additive.
    """
    mult = components.get("gw_health_multiplier", 1.0)
    gw_reward = (
        1000.0 * components.get("gw_beam_quality", 0.0)
        + 100.0 * components.get("gw_peak_power", 0.0)
    ) * mult
    health_reward = (
        5.0 * components.get("survival", 0.0)
        + 3.0 * components.get("stability", 0.0)
        + 2.0 * components.get("constraint_health", 0.0)
    ) * mult
    penalties = (
        1.0 * components.get("horizon_penalty", 0.0)
        + 1.0 * components.get("exotic_penalty", 0.0)
        + 1.0 * components.get("instability_penalty", 0.0)
        + CONSTRAINT_SPIKE_PENALTY_WEIGHT * components.get("constraint_spike_penalty", 0.0)
    )
    return gw_reward + health_reward + penalties
