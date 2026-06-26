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
    EvolvingGeodesicMetrics,
    FtlPersistenceMetrics,
    FtlTimeSeriesMetrics,
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
    evolving_geodesic: EvolvingGeodesicMetrics | None = None
    # When True, only ``ftl_geo_evolving`` may earn geodesic FTL; frozen snapshots
    # and coordinate operational FTL are withheld until the 4D trace completes.
    evolving_geodesic_mode: bool = False
    mechanism_descriptor: float | None = None
    effective_ec: EffectiveEnergyConditionMetrics | None = None
    boundary_flux: BoundaryFluxMetrics | None = None
    qei: QeiMetrics | None = None
    transport: TransportMetrics | None = None
    ftl_persistence: FtlPersistenceMetrics | None = None
    ftl_timeseries: FtlTimeSeriesMetrics | None = None
    confinement: "ConfinementMetrics | None" = None
    central: "CentralFieldMetrics | None" = None
    central_radial: "CentralRadialProfileMetrics | None" = None
