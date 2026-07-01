from __future__ import annotations

import os
from typing import Mapping

from ..types import EpisodeMetrics
from .ftl import compute_ftl_components
from .gating import compute_nontriviality_gate
from .health import compute_health_components
from .objectives import compute_total
from .penalties import compute_penalty_components
from .splash import compute_splash_components
from .survival import compute_survival_components
from .types import Score, ScoringContext
from .weights import DEFAULT_WEIGHTS


def _resolve_splash_mode(splash_mode: str | None) -> str:
    if splash_mode:
        return splash_mode
    env_mode = os.environ.get("SPLASH_MODE", "").strip().lower()
    if env_mode in {"discovery", "threshold"}:
        return env_mode
    return "discovery"


def _resolve_exotic_penalty_weight(
    explicit: float | None = None,
) -> float:
    if explicit is not None:
        return float(explicit)
    env = os.environ.get("SCORE_EXOTIC_PENALTY_WEIGHT", "").strip()
    if env:
        try:
            return float(env)
        except ValueError:
            pass
    return 1.0


def score_episode(
    metrics: EpisodeMetrics,
    *,
    target_stop_time: float | None = None,
    weights: Mapping[str, float] | None = None,
    objective_mode: str = "weighted",
    domain_half_width: float | None = None,
    splash_mode: str | None = None,
    exotic_penalty_weight: float | None = None,
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
    resolved_splash_mode = _resolve_splash_mode(splash_mode)
    if objective_mode == "critical_collapse":
        compute_splash_components(ctx, splash_mode=resolved_splash_mode)
    nontriviality = compute_nontriviality_gate(ctx)
    total = compute_total(
        ctx,
        objective_mode,
        nontriviality,
        splash_mode=resolved_splash_mode,
        exotic_penalty_weight=_resolve_exotic_penalty_weight(exotic_penalty_weight),
    )

    return Score(total=total, components=ctx.components, notes=ctx.notes)
