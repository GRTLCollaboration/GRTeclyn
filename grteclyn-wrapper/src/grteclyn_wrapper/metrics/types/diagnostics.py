"""Diagnostic and derived metric dataclasses parsed from episode outputs."""

from __future__ import annotations

from dataclasses import dataclass

STATIONARY_BETA_EPS: float = 0.05


@dataclass(frozen=True)
class CollapseMetrics:
    final_time: float | None
    min_lapse: float | None
    min_chi: float | None
    max_abs_k: float | None
    max_horizon_radius: float | None
    min_theta_plus: float | None
    scalar_phi_range: float | None
    scalar_pi_range: float | None
    r_at_min_theta_plus: float | None = None
    corroborated_trapped: bool = False
    first_corroborated_time: float | None = None
    barycenter_x: float | None = None
    barycenter_y: float | None = None
    barycenter_z: float | None = None
    rho_sum: float | None = None


@dataclass(frozen=True)
class ConstraintMetrics:
    final_time: float | None
    max_hamiltonian_l2: float | None
    max_momentum_l2: float | None
    final_hamiltonian_l2: float | None
    final_momentum_l2: float | None
    min_rho_required: float | None
    max_rho_required: float | None
    integral_negative_rho: float | None
    final_peak_rho_required: float | None = None
    initial_peak_rho_required: float | None = None
    # Ham/Mom explosion: single-step cliff or runaway max vs early baseline.
    constraint_spike_time: float | None = None
    ham_spike_ratio: float | None = None
    max_step_ham_ratio: float | None = None
    mom_spike_ratio: float | None = None
    has_constraint_spike: bool = False


@dataclass(frozen=True)
class StabilityMetrics:
    final_time: float | None
    k_growth_fraction: float | None
    lapse_drop_fraction: float | None
    chi_drop_fraction: float | None
    horizon_growth_fraction: float | None
    areal_radius_drift_fraction: float | None
    violation: float | None


@dataclass(frozen=True)
class GrowthMetrics:
    """Exponential growth rates of the constraint/collapse time series."""

    lambda_hamiltonian: float | None
    lambda_max_k: float | None
    lambda_inv_chi: float | None
    lambda_effective: float | None
    s_growth: float | None


@dataclass(frozen=True)
class ComovingMetrics:
    beta_mean: float | None
    delta_comoving: float | None
    score: float | None
    stationary: bool = False


@dataclass(frozen=True)
class EnergyConditionMetrics:
    """In-situ observer-sampled energy conditions of the evolved matter sector."""

    final_time: float | None
    min_nec: float | None
    min_wec: float | None
    min_sec: float | None
    min_dec: float | None
    max_integral_nec_violation: float | None


@dataclass(frozen=True)
class CurvatureInvariantMetrics:
    """In-situ coordinate-invariant curvature diagnostics of the evolved geometry."""

    final_time: float | None
    max_abs_ricci_scalar: float | None
    max_ricci_tensor_sq: float | None
    max_kij_sq: float | None
    max_l2_ricci_scalar: float | None


@dataclass(frozen=True)
class EffectiveEnergyConditionMetrics:
    """Effective energy conditions of the evolved geometry from plotfile stacks."""

    nec_min: float | None
    wec_min: float | None
    nec_slack_min: float | None
    rho_eulerian_min: float | None
    wec_violation_fraction: float | None
    s_energy_conditions: float | None
    n_points: int | None


@dataclass(frozen=True)
class TransportMetrics:
    """Co-moving transport objective from energy-density barycenter tracking."""

    initial_barycenter_x: float | None
    final_barycenter_x: float | None
    translation: float | None
    deformation: float | None
    score: float | None


@dataclass(frozen=True)
class QeiMetrics:
    spatial_proxy: float | None
    trajectory_violation: float | None
    s_qei: float | None


@dataclass(frozen=True)
class FtlPersistenceMetrics:
    """Sustained evolved operational FTL across the last retained plotfiles."""

    n_samples: int
    f_op_min: float | None
    f_op_median: float | None
    f_op_last: float | None
    max_local_speed_min: float | None
    max_shift_max: float | None


@dataclass(frozen=True)
class EvolvingGeodesicMetrics:
    """End-to-end null-geodesic shortcut through an evolving metric stack."""

    f_geo: float
    f_geo_frozen_peak: float | None
    t_emit: float
    t_arrival: float | None
    t_flat: float
    n_rays: int
    n_reached: int
    h_quality_ok: bool
    max_h_rel_drift: float


@dataclass(frozen=True)
class FtlTimeSeriesMetrics:
    """Time-resolved FTL features over the WHOLE run (one sample per plotfile).

    The gauge-invariant shortcut peaks mid-run and diffuses, so the final frame
    is half-blind.  This holds the per-frame arrays (so the scorer can average a
    composite FTL x stability score over time, per the avg/sum/divide design)
    plus convenience aggregates (peak + time-of-peak + FTL lifetime fraction)
    used for diagnostics and the MAP-Elites lifetime descriptor.
    """

    n_frames: int
    # Per-frame arrays, ordered by time.
    t: tuple[float, ...]
    f_op: tuple[float, ...]
    f_geo: tuple[float, ...]
    geo_trustworthy: tuple[bool, ...]
    max_local_speed: tuple[float, ...]
    superluminal_fraction: tuple[float, ...]
    structure_coherence: tuple[float, ...]  # nan where the probe omitted it
    max_h_rel_drift: tuple[float, ...]
    # Convenience aggregates.
    f_geo_peak: float  # max trustworthy f_geo over the run
    t_at_f_geo_peak: float | None
    f_op_peak: float
    t_at_f_op_peak: float | None
    max_local_speed_peak: float
    t_at_max_speed: float | None
    superluminal_fraction_peak: float
    t_at_superluminal_peak: float | None
    ftl_lifetime_fraction: float  # frames with trustworthy f_geo>floor / n_frames
    op_lifetime_fraction: float  # frames with f_op>floor / n_frames
    # End-to-end evolving probe (populated on the final timeseries row when enabled).
    f_geo_evol: float | None = None
    f_geo_evol_ok: bool | None = None


@dataclass(frozen=True)
class ConfinementMetrics:
    """Matter spatial-confinement over the run (one sample per plotfile).

    The trustworthy "matter dispersed / flew away" detector.  Built from the
    mass(scalar_activity)-weighted spatial moments in ``confinement.dat``:
    ``rms_radius`` (spread about the true matter barycentre) and
    ``confined_frac`` (fraction of matter still within the lump scale).  Unlike
    peak/total density -- which can RISE under pump injection while the lump
    disperses -- a growing ``rms_radius`` and collapsing ``confined_frac`` are
    unambiguous dispersal signals.
    """

    n_frames: int
    final_time: float | None
    initial_rms_radius: float | None
    final_rms_radius: float | None
    max_rms_radius: float | None
    initial_confined_frac: float | None
    final_confined_frac: float | None
    min_confined_frac: float | None
    #: final_rms / initial_rms -- > 1 means the matter spread out.
    spread_ratio: float | None
    initial_total: float | None
    final_total: float | None
