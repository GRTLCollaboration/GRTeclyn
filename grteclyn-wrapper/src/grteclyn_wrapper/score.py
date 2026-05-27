from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from .metrics import EpisodeMetrics


DEFAULT_WEIGHTS: dict[str, float] = {
    "survival": 2.0,
    "constraint_health": 2.0,
    "lapse_health": 1.0,
    "horizon_penalty": 1.0,
    "nontrivial_geometry": 0.5,
    "energy_condition": 1.5,
    "initial_constraint_quality": 1.0,
}


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

    total = sum(w.get(key, 0.0) * value for key, value in components.items())
    return Score(total=total, components=components, notes=notes)
