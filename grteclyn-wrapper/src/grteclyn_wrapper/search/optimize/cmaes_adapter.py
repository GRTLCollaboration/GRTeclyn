"""CMA-ES fitness adapter over shared evaluate_overrides."""

from __future__ import annotations

from ..core.evaluation import Evaluation

DEFAULT_GATE_FITNESS = 100.0


def cmaes_fitness_from_evaluation(
    ev: Evaluation,
    *,
    gate_fitness: float = DEFAULT_GATE_FITNESS,
) -> float:
    """Map an Evaluation to CMA-ES fitness (lower is better)."""
    if ev.preflight_rejected:
        return gate_fitness
    if ev.reason:
        if ev.score < 0:
            return -ev.score
        return gate_fitness
    return -ev.score
