from __future__ import annotations

import math

from .helpers import bounded_reward
from .horizon import horizon_penalty_from_collapse
from .types import ScoringContext


def _compute_constraint_spike_penalty(
    metrics,
    components: dict[str, float],
    notes: list[str],
) -> None:
    """Penalize runs where Ham/Mom explode in a single timestep (collapse gaming)."""
    constraints = metrics.constraints
    if constraints is None or not constraints.has_constraint_spike:
        components["constraint_spike_penalty"] = 0.0
        return

    ham_ratio = constraints.ham_spike_ratio or 1.0
    step_ratio = constraints.max_step_ham_ratio or 1.0
    mom_ratio = constraints.mom_spike_ratio or 1.0
    severity = max(
        math.log10(max(ham_ratio, 1.0)) / 6.0,
        math.log10(max(step_ratio, 1.0)) / 4.0,
        math.log10(max(mom_ratio, 1.0)) / 3.0,
    )
    components["constraint_spike_penalty"] = -min(max(severity, 0.25), 0.75)
    spike_t = constraints.constraint_spike_time
    t_note = f" at t={spike_t:.2f}" if spike_t is not None else ""
    notes.append(
        "constraint spike detected"
        f"{t_note}: ham_ratio={ham_ratio:.2e}, step_ratio={step_ratio:.2e}"
    )


def compute_health_components(ctx: ScoringContext) -> None:
    metrics = ctx.metrics
    components = ctx.components
    notes = ctx.notes
    final_time = ctx.final_time

    if metrics.constraints:
        ham = metrics.constraints.max_hamiltonian_l2 or 0.0
        mom = metrics.constraints.max_momentum_l2 or 0.0
        components["constraint_health"] = (
            0.5 * bounded_reward(ham, 1.0e-2)
            + 0.5 * bounded_reward(mom, 1.0e-2)
        )

        int_neg = metrics.constraints.integral_negative_rho
        min_rho = metrics.constraints.min_rho_required
        if int_neg is not None and min_rho is not None:
            if min_rho >= 0.0 and int_neg <= 1.0e-12:
                components["energy_condition"] = 1.0
            else:
                components["energy_condition"] = bounded_reward(
                    int_neg if int_neg > 0 else abs(min_rho), 1.0e-2
                )
        else:
            components["energy_condition"] = 0.0
            notes.append("energy density diagnostics not available")

        if metrics.constraints.final_time is not None and final_time is not None:
            first_ham = metrics.constraints.max_hamiltonian_l2 or 0.0
            components["initial_constraint_quality"] = bounded_reward(
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
        k_activity = metrics.collapse.max_abs_k or 0.0
        scalar_activity = max(
            metrics.collapse.scalar_phi_range or 0.0,
            metrics.collapse.scalar_pi_range or 0.0,
        )

        components["lapse_health"] = min(max(lapse / 1.0e-3, 0.0), 1.0)

        horizon_penalty, horizon_notes = horizon_penalty_from_collapse(
            metrics.collapse,
            domain_half_width=ctx.domain_half_width,
        )
        components["horizon_penalty"] = horizon_penalty
        notes.extend(horizon_notes)
        components["nontrivial_geometry"] = min(
            math.log1p(k_activity + scalar_activity), 1.0
        )
    else:
        components["lapse_health"] = 0.0
        components["horizon_penalty"] = 0.0
        components["nontrivial_geometry"] = 0.0
        notes.append("collapse_diagnostics.dat missing")

    if metrics.stability and metrics.stability.violation is not None:
        components["stability"] = bounded_reward(metrics.stability.violation, 1.0)
        components["instability_penalty"] = -(1.0 - components["stability"])
        if components["stability"] < 0.25:
            notes.append("geometry changes rapidly over the evolution window (Eulerian)")
    else:
        components["stability"] = 0.0
        components["instability_penalty"] = 0.0
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

    _compute_constraint_spike_penalty(metrics, components, notes)

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
