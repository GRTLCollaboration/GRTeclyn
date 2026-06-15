from __future__ import annotations

from dataclasses import dataclass, field

from ..types import EpisodeMetrics


@dataclass(frozen=True)
class Score:
    total: float
    components: dict[str, float]
    notes: list[str]


@dataclass
class ScoringContext:
    metrics: EpisodeMetrics
    target_stop_time: float | None
    domain_half_width: float | None
    weights: dict[str, float]
    components: dict[str, float] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    stationary_artifact: bool = False
    geo_trustworthy: bool = False
    f_geo: float = 0.0

    @property
    def final_time(self) -> float | None:
        if self.metrics.collapse and self.metrics.collapse.final_time is not None:
            return self.metrics.collapse.final_time
        if self.metrics.constraints and self.metrics.constraints.final_time is not None:
            return self.metrics.constraints.final_time
        return None

    @property
    def structural_persistence(self) -> float:
        return self.components.get("structural_persistence", 1.0)
