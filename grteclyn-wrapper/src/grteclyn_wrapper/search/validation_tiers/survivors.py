"""Survivor front and campaign-level survivor collection."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

from ..pareto import ParetoPoint, pareto_front
from .types import SURVIVOR_OBJECTIVES, SurvivorPoint, Tier


def survivor_front(
    survivors: Sequence[SurvivorPoint],
    *,
    objectives: Sequence[str] = SURVIVOR_OBJECTIVES,
) -> list[SurvivorPoint]:
    """Non-dominated set among gauntlet survivors (the candidate solutions)."""
    points = [
        ParetoPoint(label=s.label, objectives=s.objectives, total=s.score, episode=s.episode)
        for s in survivors
    ]
    front_labels = {p.label for p in pareto_front(points, objectives=objectives)}
    front = [s for s in survivors if s.label in front_labels]
    front.sort(key=lambda s: (-s.tier, -s.score))
    return front


def build_survivors(
    items: Sequence[Mapping[str, Any]],
    *,
    min_tier: int = int(Tier.OPERATIONAL),
) -> list[SurvivorPoint]:
    """Collect survivors at or above ``min_tier`` from assessed records."""
    out: list[SurvivorPoint] = []
    for it in items:
        if int(it.get("tier", Tier.REJECTED)) < min_tier:
            continue
        out.append(SurvivorPoint(
            label=str(it.get("label")),
            tier=int(it.get("tier", Tier.REJECTED)),
            score=float(it.get("score", 0.0)),
            objectives=dict(it.get("objectives", {})),
            episode=it.get("episode"),
        ))
    return out
