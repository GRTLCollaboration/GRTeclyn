"""Episode scoring: component calculators and objective scalarization."""

from .helpers import (
    domain_half_width_for_episode,
    domain_half_width_from_overrides,
    domain_half_width_from_params,
)
from .horizon import (
    HORIZON_OFFCENTER_FRACTION,
    LATE_COLLAPSE_TIME_FRACTION,
    horizon_penalty_from_collapse,
)
from .scorer import score_episode
from .types import Score
from .weights import DEFAULT_WEIGHTS, HEALTH_COMPONENTS

__all__ = [
    "DEFAULT_WEIGHTS",
    "HEALTH_COMPONENTS",
    "HORIZON_OFFCENTER_FRACTION",
    "LATE_COLLAPSE_TIME_FRACTION",
    "Score",
    "domain_half_width_for_episode",
    "domain_half_width_from_overrides",
    "domain_half_width_from_params",
    "horizon_penalty_from_collapse",
    "score_episode",
]
