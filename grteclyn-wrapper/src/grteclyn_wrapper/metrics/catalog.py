"""Central registry of episode metric groups."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True)
class MetricSpec:
    group: str
    model: str
    module: str
    source_type: Literal["dat", "recipe", "plotfile", "gridinit", "derived"]
    source_detail: str
    compute_fn: str
    in_episode: bool
    score_components: tuple[str, ...]
    summary: str


METRIC_REGISTRY: dict[str, MetricSpec] = {
    "collapse": MetricSpec(
        group="collapse",
        model="CollapseMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.collapse",
        source_type="dat",
        source_detail="data/collapse_diagnostics.dat",
        compute_fn="read_collapse_metrics",
        in_episode=True,
        score_components=("survival", "lapse_health", "horizon_penalty", "transport_objective"),
        summary="Collapse diagnostics: lapse, chi, K, horizon, barycenter.",
    ),
    "constraints": MetricSpec(
        group="constraints",
        model="ConstraintMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.constraints",
        source_type="dat",
        source_detail="data/constraint_norms.dat",
        compute_fn="read_constraint_metrics",
        in_episode=True,
        score_components=("constraint_health", "initial_constraint_quality", "constraint_growth"),
        summary="Hamiltonian/momentum L2 norms and required energy density.",
    ),
    "stability": MetricSpec(
        group="stability",
        model="StabilityMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.stability",
        source_type="dat",
        source_detail="collapse_diagnostics.dat + areal_radius.dat",
        compute_fn="read_stability_metrics",
        in_episode=True,
        score_components=("stability", "instability_penalty"),
        summary="Fractional growth/drop fractions across the run.",
    ),
    "growth": MetricSpec(
        group="growth",
        model="GrowthMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.growth",
        source_type="derived",
        source_detail="collapse + constraint time series",
        compute_fn="read_growth_metrics",
        in_episode=True,
        score_components=("constraint_growth",),
        summary="Exponential growth rates lambda and s_growth reward.",
    ),
    "comoving": MetricSpec(
        group="comoving",
        model="ComovingMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.comoving",
        source_type="derived",
        source_detail="shell_profiles.dat + recipe/gridinit shift",
        compute_fn="read_comoving_metrics",
        in_episode=True,
        score_components=("comoving_stability", "stationary_artifact_penalty"),
        summary="Co-moving stability from chi drift vs shift magnitude.",
    ),
    "energy_conditions": MetricSpec(
        group="energy_conditions",
        model="EnergyConditionMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.energy_conditions",
        source_type="dat",
        source_detail="data/energy_conditions.dat",
        compute_fn="read_energy_condition_metrics",
        in_episode=True,
        score_components=("energy_condition",),
        summary="In-situ matter-sector NEC/WEC/SEC/DEC from C++.",
    ),
    "curvature": MetricSpec(
        group="curvature",
        model="CurvatureInvariantMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.curvature",
        source_type="dat",
        source_detail="data/curvature_invariants.dat",
        compute_fn="read_curvature_invariant_metrics",
        in_episode=True,
        score_components=("curvature_activity", "nontrivial_geometry"),
        summary="Ricci and extrinsic curvature invariants from C++.",
    ),
    "transport": MetricSpec(
        group="transport",
        model="TransportMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.transport",
        source_type="dat",
        source_detail="collapse_diagnostics.dat barycenter columns",
        compute_fn="read_transport_metrics",
        in_episode=True,
        score_components=("transport_objective",),
        summary="Energy-density barycenter translation and deformation.",
    ),
    "qei": MetricSpec(
        group="qei",
        model="QeiMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.qei",
        source_type="derived",
        source_detail="physical.qei_spatial_proxy + geodesic trajectory",
        compute_fn="read_qei_metrics",
        in_episode=True,
        score_components=("qei_penalty",),
        summary="Quantum energy inequality spatial and trajectory proxies.",
    ),
    "ftl_analytic": MetricSpec(
        group="ftl_analytic",
        model="FtlMetrics",
        module="grteclyn_wrapper.metrics.probes.ftl.analytic",
        source_type="recipe",
        source_detail="episode metadata overrides",
        compute_fn="compute_ftl_metrics",
        in_episode=True,
        score_components=("ftl_shortcut", "expansion_asymmetry", "nonflat_geometry"),
        summary="t=0 1D null-ray, portal, throat, and asymmetry metrics.",
    ),
    "ftl_general_t0": MetricSpec(
        group="ftl_general_t0",
        model="GeneralFtlReport",
        module="grteclyn_wrapper.metrics.probes.ftl.general",
        source_type="recipe",
        source_detail="t=0 2D Dijkstra on recipe slice",
        compute_fn="compute_general_ftl",
        in_episode=True,
        score_components=("ftl_precursor", "channel_progress"),
        summary="Mechanism-agnostic operational FTL at t=0.",
    ),
    "ftl_general_evolved": MetricSpec(
        group="ftl_general_evolved",
        model="GeneralFtlReport",
        module="grteclyn_wrapper.metrics.probes.ftl.general",
        source_type="plotfile",
        source_detail="latest AMReX plotfile",
        compute_fn="compute_general_ftl_from_plotfile",
        in_episode=True,
        score_components=("operational_ftl", "ftl_precursor", "shift_drive"),
        summary="Evolved operational FTL — primary FTL reward signal.",
    ),
    "ftl_persistence": MetricSpec(
        group="ftl_persistence",
        model="FtlPersistenceMetrics",
        module="grteclyn_wrapper.metrics.aggregation.collector",
        source_type="plotfile",
        source_detail="last N retained plotfiles",
        compute_fn="read_episode_metrics",
        in_episode=True,
        score_components=("ftl_persistence",),
        summary="Sustained F_op over recent plotfile window.",
    ),
    "ftl_solved": MetricSpec(
        group="ftl_solved",
        model="GeneralFtlReport",
        module="grteclyn_wrapper.metrics.probes.ftl.solved",
        source_type="gridinit",
        source_detail="initial_data.gridinit",
        compute_fn="compute_solved_geometry_ftl",
        in_episode=True,
        score_components=("operational_ftl_solved", "mechanism_descriptor"),
        summary="Pre-evolution FTL on GRTresna solved geometry.",
    ),
    "geodesic_ftl": MetricSpec(
        group="geodesic_ftl",
        model="GeodesicFtlReport",
        module="grteclyn_wrapper.metrics.probes.ftl.geodesic",
        source_type="plotfile",
        source_detail="latest plotfile null-ray integration",
        compute_fn="compute_geodesic_ftl_from_plotfile",
        in_episode=True,
        score_components=("operational_ftl_geodesic",),
        summary="Gauge-invariant null-geodesic arrival-time advantage.",
    ),
    "evolving_geodesic": MetricSpec(
        group="evolving_geodesic",
        model="EvolvingGeodesicMetrics",
        module="grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic",
        source_type="plotfile",
        source_detail=">=3 plotfile stack, 4D time-interpolated null-ray trace",
        compute_fn="compute_evolving_geodesic_ftl_from_plotfiles",
        in_episode=True,
        score_components=("ftl_geo_evolving",),
        summary="End-to-end gauge-invariant shortcut through evolving geometry (opt-in).",
    ),
    "physical": MetricSpec(
        group="physical",
        model="PhysicalMetrics",
        module="grteclyn_wrapper.metrics.probes.physical",
        source_type="recipe",
        source_detail="t=0 recipe ANEC/tidal proxies",
        compute_fn="compute_physical_metrics",
        in_episode=True,
        score_components=("anec_condition", "tidal_comfort"),
        summary="ANEC line integral, tidal proxy, trapped-surface proxy.",
    ),
    "boundary_flux": MetricSpec(
        group="boundary_flux",
        model="BoundaryFluxMetrics",
        module="grteclyn_wrapper.metrics.probes.boundary",
        source_type="dat",
        source_detail="boundary_flux.dat or live plotfile",
        compute_fn="read_boundary_flux_metrics",
        in_episode=True,
        score_components=("boundary_penalty",),
        summary="Scalar boundary flux and reflection contamination.",
    ),
    "effective_ec": MetricSpec(
        group="effective_ec",
        model="EffectiveEnergyConditionMetrics",
        module="grteclyn_wrapper.metrics.probes.warpfactory",
        source_type="plotfile",
        source_detail="stack of >= 3 plotfiles",
        compute_fn="effective_energy_conditions_from_plotfiles",
        in_episode=True,
        score_components=("exotic_penalty", "energy_condition"),
        summary="Warp Factory T^eff = G/8pi effective energy conditions.",
    ),
    "central": MetricSpec(
        group="central",
        model="CentralFieldMetrics",
        module="grteclyn_wrapper.metrics.diagnostics.central",
        source_type="dat",
        source_detail="small_data/central_timeseries.dat",
        compute_fn="read_central_field_metrics",
        in_episode=True,
        score_components=(
            "central_energy_peak",
            "focusing_efficiency",
            "pre_collapsed_penalty",
            "central_lapse_collapse",
        ),
        summary="Origin-resolved rho/lapse/scalar activity for splash campaigns.",
    ),
}


def list_metrics() -> list[MetricSpec]:
    return list(METRIC_REGISTRY.values())


def get_metric(group: str) -> MetricSpec | None:
    return METRIC_REGISTRY.get(group)
