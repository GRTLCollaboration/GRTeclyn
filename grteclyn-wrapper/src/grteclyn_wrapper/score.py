from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from .metrics import EpisodeMetrics


DEFAULT_WEIGHTS: dict[str, float] = {
    "ftl_shortcut": 5.0,
    "expansion_asymmetry": 2.0,
    "nonflat_geometry": 1.0,
    "comoving_stability": 2.5,
    "survival": 1.5,
    "constraint_health": 2.0,
    "lapse_health": 1.0,
    "horizon_penalty": 1.5,
    "energy_condition": 2.0,
    "stability": 0.5,
    "nontrivial_geometry": 0.25,
    "initial_constraint_quality": 0.5,
    # Proposed extensions (Sec. "Proposed Extensions").
    "constraint_growth": 2.0,
    "anec_condition": 1.5,
    "tidal_comfort": 1.0,
    # Mechanism-agnostic operational FTL (Dijkstra arrival-time advantage on the
    # slice; warp/wormhole/portal/unknown alike) -- this is the primary goal, so
    # it dominates the reward budget.  Curvature activity rewards genuinely
    # non-trivial geometry (and keeps flat space out of the running).
    "operational_ftl": 6.0,
    "curvature_activity": 0.5,
    # Exotic-matter PENALTY.  The objective is FTL *without* exotic matter, so
    # any negative-energy requirement is punished -- both the matter sourced at
    # t=0 by the constrained solve and the geometric exotic energy of the
    # evolved spacetime (effective NEC violation, T^eff = G/8pi).  The component
    # is negative; this weight sets how hard exoticity is penalized.
    "exotic_penalty": 8.0,
}


# "Health/niceness" rewards that a trivial flat spacetime maximises by default.
# They are gated by a non-triviality factor (see ``score_episode``) so flat
# Minkowski cannot bank them and out-rank a genuinely exotic geometry.
HEALTH_COMPONENTS: frozenset[str] = frozenset({
    "survival",
    "constraint_health",
    "energy_condition",
    "initial_constraint_quality",
    "lapse_health",
    "stability",
    "comoving_stability",
    "constraint_growth",
    "anec_condition",
    "tidal_comfort",
})


@dataclass(frozen=True)
class Score:
    total: float
    components: dict[str, float]
    notes: list[str]


def _bounded_reward(value: float, scale: float) -> float:
    if not math.isfinite(value) or value < 0:
        return 0.0
    return 1.0 / (1.0 + value / scale)


def score_episode(
    metrics: EpisodeMetrics,
    *,
    target_stop_time: float | None = None,
    weights: Mapping[str, float] | None = None,
) -> Score:
    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update(weights)

    components: dict[str, float] = {}
    notes: list[str] = []

    final_time = None
    if metrics.collapse and metrics.collapse.final_time is not None:
        final_time = metrics.collapse.final_time
    elif metrics.constraints and metrics.constraints.final_time is not None:
        final_time = metrics.constraints.final_time

    if target_stop_time and target_stop_time > 0 and final_time is not None:
        components["survival"] = min(final_time / target_stop_time, 1.0)
    elif final_time is not None:
        components["survival"] = 1.0
    else:
        components["survival"] = 0.0
        notes.append("no time-series diagnostics were found")

    if metrics.constraints:
        ham = metrics.constraints.max_hamiltonian_l2 or 0.0
        mom = metrics.constraints.max_momentum_l2 or 0.0
        components["constraint_health"] = (
            0.5 * _bounded_reward(ham, 1.0e-2)
            + 0.5 * _bounded_reward(mom, 1.0e-2)
        )

        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and min_rho is not None:
            if min_rho >= 0.0 and int_neg <= 1.0e-12:
                components["energy_condition"] = 1.0
            else:
                components["energy_condition"] = _bounded_reward(
                    int_neg if int_neg > 0 else abs(min_rho), 1.0e-2
                )
        else:
            components["energy_condition"] = 0.0
            notes.append("energy density diagnostics not available")

        if metrics.constraints.final_time is not None and final_time is not None:
            first_ham = metrics.constraints.max_hamiltonian_l2 or 0.0
            components["initial_constraint_quality"] = _bounded_reward(
                first_ham, 1.0e-3
            )
        else:
            components["initial_constraint_quality"] = 0.0
    else:
        components["constraint_health"] = 0.0
        components["energy_condition"] = 0.0
        components["initial_constraint_quality"] = 0.0
        notes.append("constraint_norms.dat missing")

    if metrics.collapse:
        lapse = metrics.collapse.min_lapse or 0.0
        horizon = metrics.collapse.max_horizon_radius or 0.0
        k_activity = metrics.collapse.max_abs_k or 0.0
        scalar_activity = max(
            metrics.collapse.scalar_phi_range or 0.0,
            metrics.collapse.scalar_pi_range or 0.0,
        )

        components["lapse_health"] = min(max(lapse / 1.0e-3, 0.0), 1.0)
        components["horizon_penalty"] = -min(horizon, 1.0)
        components["nontrivial_geometry"] = min(
            math.log1p(k_activity + scalar_activity), 1.0
        )
    else:
        components["lapse_health"] = 0.0
        components["horizon_penalty"] = 0.0
        components["nontrivial_geometry"] = 0.0
        notes.append("collapse_diagnostics.dat missing")

    if metrics.stability and metrics.stability.violation is not None:
        components["stability"] = _bounded_reward(metrics.stability.violation, 1.0)
        if components["stability"] < 0.25:
            notes.append("geometry changes rapidly over the evolution window (Eulerian)")
    else:
        components["stability"] = 0.0
        notes.append("stability diagnostics not available")

    if metrics.comoving and metrics.comoving.stationary:
        components["comoving_stability"] = components["stability"]
        notes.append("co-moving stability uses Eulerian fallback (stationary geometry)")
    elif metrics.comoving and metrics.comoving.score is not None:
        components["comoving_stability"] = metrics.comoving.score
    else:
        components["comoving_stability"] = 0.0
        notes.append("co-moving stability diagnostics not available")

    if metrics.growth is not None and metrics.growth.s_growth is not None:
        components["constraint_growth"] = metrics.growth.s_growth
        if metrics.growth.lambda_effective is not None and metrics.growth.lambda_effective > 0.5:
            notes.append(
                "constraint/collapse series grow exponentially "
                f"(lambda={metrics.growth.lambda_effective:.3f}); slow-collapse penalized"
            )
    else:
        components["constraint_growth"] = 0.0

    if metrics.physical is not None:
        if metrics.physical.s_anec is not None:
            components["anec_condition"] = metrics.physical.s_anec
        else:
            components["anec_condition"] = 0.0
        if metrics.physical.s_tidal is not None:
            components["tidal_comfort"] = metrics.physical.s_tidal
        else:
            components["tidal_comfort"] = 0.0
        for note in metrics.physical.notes:
            notes.append(note)
    else:
        components["anec_condition"] = 0.0
        components["tidal_comfort"] = 0.0

    if metrics.ftl is not None:
        components["ftl_shortcut"] = metrics.ftl.f_log
        components["expansion_asymmetry"] = metrics.ftl.f_asymmetry
        components["nonflat_geometry"] = metrics.ftl.s_nonflat
        if metrics.ftl.f_shortcut <= 0.0 and metrics.ftl.s_nonflat < 0.05:
            notes.append("no FTL shortcut detected in t=0 profile")
        for note in metrics.ftl.notes:
            notes.append(note)
    else:
        components["ftl_shortcut"] = 0.0
        components["expansion_asymmetry"] = 0.0
        components["nonflat_geometry"] = 0.0
        notes.append("FTL profile metrics not available")

    # Mechanism-agnostic operational FTL: reward the larger arrival-time
    # advantage from either the evolved-spacetime search (preferred, when a
    # plotfile survived) or the t=0 reconstructed slice.  f_op is already in
    # [0, 1] (fraction of the flat baseline saved).
    f_op_t0 = (
        metrics.general_ftl.f_op
        if metrics.general_ftl is not None
        else 0.0
    )
    f_op_ev = (
        metrics.general_ftl_evolved.f_op
        if metrics.general_ftl_evolved is not None
        else 0.0
    )
    f_op = max(f_op_t0, f_op_ev)
    # Real arrival-time advantages are small (F_op ~ 1e-3); a raw [0,1] reward is
    # dwarfed by the health basket, so FTL never actually drives the search.
    # Log-amplify against OP_FTL_SCALE so a genuine (small) shortcut is a
    # multi-point reward, while F_op == 0 (no shortcut) stays exactly 0.
    OP_FTL_SCALE = 2.0e-3
    components["operational_ftl"] = (
        min(math.log1p(f_op / OP_FTL_SCALE), 1.0)
        if math.isfinite(f_op) and f_op > 0
        else 0.0
    )
    if metrics.general_ftl_evolved is not None:
        if metrics.general_ftl_evolved.max_local_speed > 1.0:
            notes.append(
                "evolved geometry has a superluminal coordinate channel "
                f"(max c = {metrics.general_ftl_evolved.max_local_speed:.3f})"
            )
    elif metrics.general_ftl is not None and metrics.general_ftl.f_op <= 0.0:
        notes.append("no operational FTL shortcut on the t=0 slice")

    # Coordinate-invariant curvature activity: a saturating reward so flat space
    # scores ~0 while a structured warp/wormhole throat scores up to 1, without
    # letting a near-singular blow-up dominate the objective.
    if metrics.curvature is not None and metrics.curvature.max_l2_ricci_scalar is not None:
        components["curvature_activity"] = min(
            math.log1p(max(0.0, metrics.curvature.max_l2_ricci_scalar)), 1.0
        )
    else:
        components["curvature_activity"] = 0.0

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
    def _graded(x: float, scale: float) -> float:
        if not math.isfinite(x) or x <= 0.0:
            return 0.0
        return min(math.log1p(x / scale), 1.0)

    MATTER_INT_SCALE = 0.5     # integral of negative rho (volume measure)
    MATTER_PEAK_SCALE = 5.0e-2  # peak negative energy density
    GEO_NEC_SCALE = 5.0e-4      # effective-NEC violation of evolved geometry

    matter_exotic = 0.0
    geo_exotic = 0.0
    if metrics.constraints is not None:
        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and int_neg > 0.0:
            matter_exotic = max(matter_exotic, _graded(int_neg, MATTER_INT_SCALE))
        if min_rho is not None and min_rho < 0.0:
            matter_exotic = max(matter_exotic, _graded(-min_rho, MATTER_PEAK_SCALE))
    if metrics.effective_ec is not None:
        eff = metrics.effective_ec
        eff_candidates = [
            v for v in (eff.nec_slack_min, eff.nec_min, eff.rho_eulerian_min)
            if v is not None
        ]
        eff_worst = min(eff_candidates) if eff_candidates else 0.0
        geo_exotic = _graded(-eff_worst, GEO_NEC_SCALE)

    exotic_term = max(matter_exotic, geo_exotic)
    components["exotic_penalty"] = -exotic_term
    if exotic_term > 0.05:
        notes.append(
            "exotic matter required "
            f"(matter={matter_exotic:.2f}, geometric={geo_exotic:.2f} of full penalty)"
        )

    # Non-triviality gate: a trivial flat spacetime aces every "health" reward
    # (survival, stability, clean constraints, no exotic energy) while scoring
    # zero on all exoticity terms, which previously let it out-rank genuine warp
    # / wormhole geometries (the flat-space attractor).  We gate the health
    # rewards by how non-trivial the geometry actually is, so they only count
    # once there is real structure to be healthy about.  The exoticity/FTL terms
    # are never gated -- they supply the gradient out of flatness.
    # Note: exotic energy is deliberately NOT a non-triviality signal anymore --
    # we do not want negative-energy content to "earn" the health rewards.  FTL
    # and curvature structure supply the gradient out of flatness instead.
    nontriviality = max(
        components.get("nonflat_geometry", 0.0),
        components.get("expansion_asymmetry", 0.0),
        components.get("nontrivial_geometry", 0.0),
        components.get("curvature_activity", 0.0),
        components.get("operational_ftl", 0.0),
        components.get("ftl_shortcut", 0.0),
    )
    nontriviality = float(min(max(nontriviality, 0.0), 1.0))
    components["nontriviality_gate"] = nontriviality
    if nontriviality < 0.05:
        notes.append(
            "near-trivial geometry: health rewards gated out (flat-space attractor guard)"
        )

    total = 0.0
    for key, value in components.items():
        if key == "nontriviality_gate":
            continue
        gate = nontriviality if key in HEALTH_COMPONENTS else 1.0
        total += w.get(key, 0.0) * value * gate
    return Score(total=total, components=components, notes=notes)
