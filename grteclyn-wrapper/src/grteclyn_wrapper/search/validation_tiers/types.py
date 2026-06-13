"""Tier ladder types and configuration."""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Any


class Tier(IntEnum):
    """Falsification ladder.  ``REJECTED`` means even T0 failed."""

    REJECTED = -1
    CONSTRUCTED = 0
    NONTRIVIAL = 1
    OPERATIONAL = 2
    PERSISTENT = 3
    OBSERVER_EC = 4
    CONVERGED = 5
    ANALYTIC = 6


TIER_NAMES: dict[int, str] = {
    Tier.REJECTED: "rejected",
    Tier.CONSTRUCTED: "constructed",
    Tier.NONTRIVIAL: "nontrivial",
    Tier.OPERATIONAL: "operational",
    Tier.PERSISTENT: "persistent",
    Tier.OBSERVER_EC: "observer_ec",
    Tier.CONVERGED: "converged",
    Tier.ANALYTIC: "analytic",
}

PASS = "pass"
FAIL = "fail"
UNAVAILABLE = "unavailable"

SURVIVOR_OBJECTIVES: tuple[str, ...] = (
    "operational_ftl",
    "anec_condition",
    "constraint_growth",
    "tidal_comfort",
)


@dataclass(frozen=True)
class TierConfig:
    """Thresholds for the falsification gates."""

    constructed_floor: float = 0.05
    nontrivial_floor: float = 1.0e-3
    operational_floor: float = 1.0e-3
    geodesic_floor: float = 1.0e-3
    constraint_floor: float = 0.20
    max_trapped: float = -0.5
    survival_floor: float = 0.99
    stability_floor: float = 0.25
    growth_floor: float = 0.50
    persistence_floor: float = 1.0e-4
    ec_floor: float = 1.0e-3
    exotic_cap: float = 6.0
    qei_cap: float = 1.0


DEFAULT_TIER_CONFIG = TierConfig()


@dataclass(frozen=True)
class GateResult:
    tier: int
    name: str
    status: str
    detail: str = ""


@dataclass(frozen=True)
class TierAssessment:
    tier: int
    tier_name: str
    gates: list[GateResult]
    objectives: dict[str, float]

    def to_dict(self) -> dict[str, Any]:
        return {
            "tier": int(self.tier),
            "tier_name": self.tier_name,
            "gates": [
                {"tier": g.tier, "name": g.name, "status": g.status, "detail": g.detail}
                for g in self.gates
            ],
            "objectives": dict(self.objectives),
        }


@dataclass(frozen=True)
class SurvivorPoint:
    label: str
    tier: int
    score: float
    objectives: dict[str, float]
    episode: str | None
