"""Archive convergence detection for MAP-Elites campaigns."""

from __future__ import annotations

from typing import Any, Mapping, Sequence


def _jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 1.0
    union = a | b
    if not union:
        return 1.0
    return len(a & b) / len(union)


def convergence_signals(
    history: Sequence[Mapping[str, Any]],
    *,
    window: int = 3,
    score_eps: float = 1.0e-3,
) -> dict[str, Any]:
    """Decide whether the archive has stopped improving."""
    n = len(history)
    if n == 0:
        return {"converged": False, "reason": "no history"}
    last = history[-1]
    if n < window + 1:
        return {
            "converged": False,
            "reason": f"need {window + 1} iters, have {n}",
            "coverage": last.get("coverage"),
            "best_score": last.get("best_score"),
        }

    recent = history[-(window + 1):]
    coverages = [float(h.get("coverage", 0.0)) for h in recent]
    scores = [float(h.get("best_score", 0.0)) for h in recent]
    coverage_delta = coverages[-1] - coverages[0]
    best_delta = scores[-1] - scores[0]

    front_sets = [set(h.get("front_labels", []) or []) for h in recent]
    pair_stability = [
        _jaccard(front_sets[i], front_sets[i + 1]) for i in range(len(front_sets) - 1)
    ]
    front_stable = bool(pair_stability) and min(pair_stability) >= 0.999

    coverage_stalled = coverage_delta <= 0.0
    score_stalled = best_delta < score_eps
    converged = coverage_stalled and score_stalled and front_stable

    return {
        "converged": bool(converged),
        "coverage": coverages[-1],
        "coverage_delta": coverage_delta,
        "best_score": scores[-1],
        "best_score_delta": best_delta,
        "front_stability": min(pair_stability) if pair_stability else 1.0,
        "coverage_stalled": coverage_stalled,
        "score_stalled": score_stalled,
        "front_stable": front_stable,
        "window": window,
    }
