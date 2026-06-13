"""Falsification-tier layer for the MAP-Elites Spacetime Failure Atlas."""

from __future__ import annotations

from .campaign import assess_campaign
from .convergence import convergence_signals
from .evaluate import evaluate_tiers, survivor_objectives, tier_histogram
from .survivors import build_survivors, survivor_front
from .types import (
    DEFAULT_TIER_CONFIG,
    FAIL,
    PASS,
    SURVIVOR_OBJECTIVES,
    UNAVAILABLE,
    GateResult,
    SurvivorPoint,
    Tier,
    TierAssessment,
    TierConfig,
    TIER_NAMES,
)

__all__ = [
    "DEFAULT_TIER_CONFIG",
    "FAIL",
    "PASS",
    "SURVIVOR_OBJECTIVES",
    "UNAVAILABLE",
    "GateResult",
    "SurvivorPoint",
    "Tier",
    "TierAssessment",
    "TierConfig",
    "TIER_NAMES",
    "assess_campaign",
    "build_survivors",
    "convergence_signals",
    "evaluate_tiers",
    "survivor_front",
    "survivor_objectives",
    "tier_histogram",
]
