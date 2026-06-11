"""Aggregate episode metrics container."""

from __future__ import annotations

from dataclasses import dataclass

from ..probes.boundary import BoundaryFluxMetrics
from ..probes.ftl.analytic import FtlMetrics
from ..probes.ftl.general import GeneralFtlReport
from ..probes.ftl.geodesic import GeodesicFtlReport
from ..probes.physical import PhysicalMetrics
from .diagnostics import (
    CollapseMetrics,
    ComovingMetrics,
    ConstraintMetrics,
    CurvatureInvariantMetrics,
    EffectiveEnergyConditionMetrics,
    EnergyConditionMetrics,
    FtlPersistenceMetrics,
    GrowthMetrics,
    QeiMetrics,
    StabilityMetrics,
    TransportMetrics,
)


@dataclass(frozen=True)
class EpisodeMetrics:
    collapse: CollapseMetrics | None
    constraints: ConstraintMetrics | None
    stability: StabilityMetrics | None
    comoving: ComovingMetrics | None
    ftl: FtlMetrics | None
    termination_reason: str
    growth: GrowthMetrics | None = None
    physical: PhysicalMetrics | None = None
    energy_conditions: EnergyConditionMetrics | None = None
    curvature: CurvatureInvariantMetrics | None = None
    general_ftl: GeneralFtlReport | None = None
    general_ftl_evolved: GeneralFtlReport | None = None
    general_ftl_solved: GeneralFtlReport | None = None
    geodesic_ftl: GeodesicFtlReport | None = None
    mechanism_descriptor: float | None = None
    effective_ec: EffectiveEnergyConditionMetrics | None = None
    boundary_flux: BoundaryFluxMetrics | None = None
    qei: QeiMetrics | None = None
    transport: TransportMetrics | None = None
    ftl_persistence: FtlPersistenceMetrics | None = None
