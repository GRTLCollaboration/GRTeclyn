from __future__ import annotations

import math

from ..types import STATIONARY_BETA_EPS
from .helpers import graded
from .types import ScoringContext

MATTER_INT_SCALE = 0.5     # integral of negative rho (volume measure)
MATTER_PEAK_SCALE = 5.0e-2  # peak negative energy density
GEO_NEC_SCALE = 5.0e-4      # effective-NEC violation of evolved geometry
# The persistent-FTL winners all live in the int_neg ~ 0.5-2 band, where a
# ceiling of 1.0 saturates and can no longer tell "a bit exotic" from "very
# exotic".  Raise the matter ceiling so the penalty keeps a usable gradient
# right through that band -- this is what lets the search minimise exotic
# content *among* the geometries that actually sustain a channel.
MATTER_CEILING = 1.6
# Pump-work penalty: honest open-system accounting for the PD controller.
# Does NOT restore nabla_mu T^munu=0 by itself; it taxes artificial forcing so
# MAP-Elites prefers configs that need less pump.  Scale is on pump_energy_norm.
PUMP_ENERGY_SCALE = 0.25
PUMP_ENERGY_CEILING = 1.6


def compute_penalty_components(ctx: ScoringContext) -> None:
    metrics = ctx.metrics
    components = ctx.components
    notes = ctx.notes

    # Graded stationary penalty.  A *flat* -1.0 collapses every zero-net-shift
    # geometry to the same score, so the map has no slope and cannot tell
    # "almost propulsive" from "perfectly static" -- the whole stationary basin
    # is a flat floor the search cannot climb out of.  Instead we scale the
    # penalty by how much *coherent net axial shift* the bubble carries
    # (|beta_mean|, the directed frame-drag that an actual warp needs),
    # normalized to the stationarity threshold.  A truly static lens
    # (beta_mean -> 0) still earns the full -1.0, so the eval_083 artifact fix
    # is preserved, but a geometry developing directed shift sees the penalty
    # smoothly relax toward 0 as it approaches the threshold -- a continuous
    # downhill gradient out of the basin toward genuine shift-driven FTL.
    beta_mean_mag = (
        abs(metrics.comoving.beta_mean)
        if metrics.comoving is not None and metrics.comoving.beta_mean is not None
        else 0.0
    )
    shift_coherence = (
        min(beta_mean_mag / STATIONARY_BETA_EPS, 1.0)
        if STATIONARY_BETA_EPS > 0.0
        else 0.0
    )
    components["stationary_artifact_penalty"] = (
        -(1.0 - shift_coherence) if ctx.stationary_artifact else 0.0
    )

    # Exotic-matter penalty.  The goal is FTL *without* exotic matter, so any
    # negative-energy requirement is penalized.  Two independent probes feed it:
    #   1. matter sector at t=0: the constrained solve back-solves the matter
    #      from the proposed geometry; a non-zero integral of negative rho (or a
    #      negative peak min_rho_required) means exotic matter is needed up front.
    #   2. evolved geometry: the effective stress-energy T^eff = G/8pi of the
    #      evolved spacetime violating the NEC (geometric exotic energy that the
    #      matter-sector EC cannot see -- e.g. a shift-driven warp).
    #
    # Each probe is log-scaled against its own reference (the magnitudes differ
    # by orders of magnitude: matter integral ~O(1), geometric NEC ~O(1e-4)) so
    # the penalty is a *smooth gradient* across the observed range rather than a
    # saturated constant.  Steepest near zero, so "almost no exotic" is clearly
    # distinguished from "a little exotic"; heavy exotic saturates at -1.
    matter_exotic = 0.0
    geo_exotic = 0.0
    if metrics.constraints is not None:
        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and int_neg > 0.0:
            matter_exotic = max(
                matter_exotic, graded(int_neg, MATTER_INT_SCALE, MATTER_CEILING)
            )
        if min_rho is not None and min_rho < 0.0:
            matter_exotic = max(
                matter_exotic, graded(-min_rho, MATTER_PEAK_SCALE, MATTER_CEILING)
            )
    if metrics.effective_ec is not None:
        eff = metrics.effective_ec
        eff_candidates = [
            v for v in (eff.nec_slack_min, eff.nec_min, eff.rho_eulerian_min)
            if v is not None
        ]
        eff_worst = min(eff_candidates) if eff_candidates else 0.0
        geo_exotic = graded(-eff_worst, GEO_NEC_SCALE)

    exotic_term = max(matter_exotic, geo_exotic)
    components["exotic_penalty"] = -exotic_term
    if exotic_term > 0.05:
        notes.append(
            "exotic matter required "
            f"(matter={matter_exotic:.2f}, geometric={geo_exotic:.2f} of full penalty)"
        )

    # Pump-energy penalty (RadialRecipe / wormhole pump_work column).  Passive when
    # the diagnostic is absent or the pump was never on.
    pump_term = 0.0
    if metrics.collapse is not None and metrics.collapse.pump_energy_norm is not None:
        pump_term = graded(
            abs(metrics.collapse.pump_energy_norm),
            PUMP_ENERGY_SCALE,
            PUMP_ENERGY_CEILING,
        )
    components["pump_energy_penalty"] = -pump_term
    if pump_term > 0.05:
        pe = (
            metrics.collapse.pump_energy
            if metrics.collapse is not None
            else None
        )
        notes.append(
            f"pump work taxed (norm={pump_term:.2f}"
            + (f", energy={pe:.3g}" if pe is not None else "")
            + "); open-system accounting, not a closed GR solution while pump is on"
        )

    qei_violation = 0.0
    if metrics.qei is not None:
        if metrics.qei.trajectory_violation is not None and metrics.qei.trajectory_violation > 0:
            qei_violation = max(qei_violation, metrics.qei.trajectory_violation)
        if metrics.qei.spatial_proxy is not None and metrics.qei.spatial_proxy > 0:
            qei_violation = max(qei_violation, metrics.qei.spatial_proxy)
    elif metrics.physical is not None and metrics.physical.qei_spatial_proxy is not None:
        qei_violation = max(0.0, metrics.physical.qei_spatial_proxy)
    components["qei_penalty"] = -min(math.log1p(qei_violation / 1.0e-2), 1.0) if qei_violation > 0 else 0.0

    if metrics.boundary_flux is not None and metrics.boundary_flux.reflection_contaminated:
        components["boundary_penalty"] = -1.0
        notes.append("late-time inward scalar flux at boundary (reflection contamination)")
    else:
        components["boundary_penalty"] = 0.0

    if metrics.transport is not None and metrics.transport.score is not None:
        components["transport_objective"] = metrics.transport.score
    else:
        components["transport_objective"] = 0.0
