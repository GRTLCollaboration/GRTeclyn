from __future__ import annotations

from typing import Mapping

from ..types import EpisodeMetrics
from .ftl import compute_ftl_components
from .gating import compute_nontriviality_gate
from .health import compute_health_components
from .objectives import compute_total
from .penalties import compute_penalty_components
from .survival import compute_survival_components
from .types import Score, ScoringContext
from .weights import DEFAULT_WEIGHTS


def score_episode(
    metrics: EpisodeMetrics,
    *,
    target_stop_time: float | None = None,
    weights: Mapping[str, float] | None = None,
    objective_mode: str = "weighted",
    domain_half_width: float | None = None,
) -> Score:
    w = dict(DEFAULT_WEIGHTS)
    if weights:
        w.update(weights)

    ctx = ScoringContext(
        metrics=metrics,
        target_stop_time=target_stop_time,
        domain_half_width=domain_half_width,
        weights=w,
    )

    compute_survival_components(ctx)
    compute_health_components(ctx)
    compute_ftl_components(ctx)
    compute_penalty_components(ctx)
    nontriviality = compute_nontriviality_gate(ctx)
    total = compute_total(ctx, objective_mode, nontriviality)

    return Score(total=total, components=ctx.components, notes=ctx.notes)
