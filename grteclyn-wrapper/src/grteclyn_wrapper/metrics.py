from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


from .ftl_metrics import FtlMetrics, compute_ftl_metrics, load_overrides_from_episode
from .ftl_general import (
    GeneralFtlReport,
    compute_general_ftl,
    compute_general_ftl_from_plotfile,
    find_latest_plotfile,
)
from .physical_metrics import PhysicalMetrics, compute_physical_metrics


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


@dataclass(frozen=True)
class StabilityMetrics:
    final_time: float | None
    k_growth_fraction: float | None
    lapse_drop_fraction: float | None
    chi_drop_fraction: float | None
    horizon_growth_fraction: float | None
    areal_radius_drift_fraction: float | None
    violation: float | None


STATIONARY_BETA_EPS: float = 0.05


@dataclass(frozen=True)
class GrowthMetrics:
    """Exponential growth rates of the constraint/collapse time series.

    A geometry that merely collapses slowly enough to reach the stop time has
    a positive effective growth rate lambda; a true equilibrium has
    lambda ~ 0.  ``s_growth`` is a bounded reward (higher = closer to
    equilibrium) that closes the dynamical stability gap without requiring
    longer, more expensive runs.
    """

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
    """In-situ observer-sampled energy conditions of the evolved matter sector.

    Written by the C++ ``RadialRecipeLevel`` (``calculate_energy_conditions``)
    to ``data/energy_conditions.dat``.  These are the conditions of the matter
    stress-energy ``T_munu`` itself; the geometry-sourced effective stress
    energy ``T^eff = G/8pi`` (Warp Factory style) is evaluated post-hoc from
    plotfiles by the ``warpfactory`` module.
    """

    final_time: float | None
    min_nec: float | None
    min_wec: float | None
    min_sec: float | None
    min_dec: float | None
    max_integral_nec_violation: float | None


@dataclass(frozen=True)
class CurvatureInvariantMetrics:
    """In-situ coordinate-invariant curvature diagnostics of the evolved
    geometry, written by the C++ ``RadialRecipeLevel``
    (``calculate_curvature_invariants``) to ``data/curvature_invariants.dat``.
    Unlike the matter energy conditions these are sourced purely by the
    geometry and so probe the exotic tidal structure directly.
    """

    final_time: float | None
    max_abs_ricci_scalar: float | None
    max_ricci_tensor_sq: float | None
    max_kij_sq: float | None
    max_l2_ricci_scalar: float | None


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


def _numeric_rows(path: Path, min_columns: int) -> list[list[float]]:
    rows: list[list[float]] = []
    if not path.exists():
        return rows

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                values = [float(item) for item in stripped.split()]
            except ValueError:
                continue
            if len(values) >= min_columns:
                rows.append(values)
    rows.sort(key=lambda row: row[0])
    return rows


def read_collapse_metrics(path: Path) -> CollapseMetrics | None:
    rows = _numeric_rows(path, 4)
    if not rows:
        return None

    final_time = rows[-1][0]
    min_lapse = min(row[1] for row in rows)
    min_chi = min(row[2] for row in rows)
    max_abs_k = max(row[3] for row in rows)
    max_horizon_radius = max((row[7] for row in rows if len(row) >= 8), default=None)
    min_theta_plus = min((row[8] for row in rows if len(row) >= 10), default=None)

    phi_min = min((row[10] for row in rows if len(row) >= 14), default=None)
    phi_max = max((row[11] for row in rows if len(row) >= 14), default=None)
    pi_min = min((row[12] for row in rows if len(row) >= 14), default=None)
    pi_max = max((row[13] for row in rows if len(row) >= 14), default=None)

    return CollapseMetrics(
        final_time=final_time,
        min_lapse=min_lapse,
        min_chi=min_chi,
        max_abs_k=max_abs_k,
        max_horizon_radius=max_horizon_radius,
        min_theta_plus=min_theta_plus,
        scalar_phi_range=(phi_max - phi_min) if phi_min is not None and phi_max is not None else None,
        scalar_pi_range=(pi_max - pi_min) if pi_min is not None and pi_max is not None else None,
    )


def read_constraint_metrics(path: Path) -> ConstraintMetrics | None:
    rows = _numeric_rows(path, 3)
    if not rows:
        return None

    min_rho_req = min((row[3] for row in rows if len(row) >= 6), default=None)
    max_rho_req = max((row[4] for row in rows if len(row) >= 6), default=None)
    max_int_neg = max((row[5] for row in rows if len(row) >= 6), default=None)

    return ConstraintMetrics(
        final_time=rows[-1][0],
        max_hamiltonian_l2=max(row[1] for row in rows),
        max_momentum_l2=max(row[2] for row in rows),
        final_hamiltonian_l2=rows[-1][1],
        final_momentum_l2=rows[-1][2],
        min_rho_required=min_rho_req,
        max_rho_required=max_rho_req,
        integral_negative_rho=max_int_neg,
    )


def _positive_fractional_change(final: float, initial: float, *, floor: float) -> float:
    return max(0.0, final - initial) / max(abs(initial), floor)


def _positive_fractional_drop(initial: float, minimum: float, *, floor: float) -> float:
    return max(0.0, initial - minimum) / max(abs(initial), floor)


def read_stability_metrics(collapse_path: Path, areal_path: Path) -> StabilityMetrics | None:
    collapse_rows = _numeric_rows(collapse_path, 4)
    areal_rows = _numeric_rows(areal_path, 2)
    if not collapse_rows and not areal_rows:
        return None

    final_time = collapse_rows[-1][0] if collapse_rows else (areal_rows[-1][0] if areal_rows else None)
    k_growth_fraction = None
    lapse_drop_fraction = None
    chi_drop_fraction = None
    horizon_growth_fraction = None
    areal_radius_drift_fraction = None

    if collapse_rows:
        initial_lapse = collapse_rows[0][1]
        minimum_lapse = min(row[1] for row in collapse_rows)
        lapse_drop_fraction = _positive_fractional_drop(initial_lapse, minimum_lapse, floor=1.0e-6)

        initial_chi = collapse_rows[0][2]
        minimum_chi = min(row[2] for row in collapse_rows)
        chi_drop_fraction = _positive_fractional_drop(initial_chi, minimum_chi, floor=1.0e-8)

        initial_k = collapse_rows[0][3]
        maximum_k = max(row[3] for row in collapse_rows)
        k_growth_fraction = _positive_fractional_change(maximum_k, initial_k, floor=1.0e-1)

        horizon_values = [row[7] for row in collapse_rows if len(row) >= 8]
        if horizon_values:
            horizon_growth_fraction = _positive_fractional_change(
                max(horizon_values),
                horizon_values[0],
                floor=1.0,
            )

    if areal_rows:
        initial_radius = areal_rows[0][1]
        areal_radius_drift_fraction = max(
            abs(row[1] - initial_radius) / max(abs(initial_radius), 1.0e-8)
            for row in areal_rows
        )

    terms = [
        value for value in (
            k_growth_fraction,
            lapse_drop_fraction,
            chi_drop_fraction,
            horizon_growth_fraction,
            areal_radius_drift_fraction,
        )
        if value is not None
    ]
    violation = sum(terms) if terms else None

    return StabilityMetrics(
        final_time=final_time,
        k_growth_fraction=k_growth_fraction,
        lapse_drop_fraction=lapse_drop_fraction,
        chi_drop_fraction=chi_drop_fraction,
        horizon_growth_fraction=horizon_growth_fraction,
        areal_radius_drift_fraction=areal_radius_drift_fraction,
        violation=violation,
    )


def _parse_shell_profiles(path: Path) -> tuple[list[float], dict[str, list[float]]] | None:
    if not path.exists():
        return None
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        return None
    header = lines[0].lstrip("#").split()
    if not header or header[0] != "time":
        return None
    cols: dict[str, list[float]] = {name: [] for name in header[1:]}
    times: list[float] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(header):
            continue
        times.append(float(parts[0]))
        for name, value in zip(header[1:], parts[1:]):
            cols[name].append(float(value))
    if not times:
        return None
    return times, cols


def _mean_beta_in_bubble(
    overrides: dict[str, object],
    *,
    ftl_L: float | None = None,
) -> float | None:
    from .ftl_metrics import _axis_profiles

    L = ftl_L if ftl_L is not None else float(overrides.get("recipe_basis_radius_max", 8.0))
    _x, _r, _chi, _alpha, beta_x = _axis_profiles(overrides, L=L, n_points=512)
    beta_max = float(max(abs(float(v)) for v in beta_x))
    if beta_max < 1.0e-8:
        return 0.0
    bubble_mask = abs(beta_x) >= 0.1 * beta_max
    if not bubble_mask.any():
        return float(beta_x.mean())
    return float(beta_x[bubble_mask].mean())


def read_comoving_metrics(
    episode_dir: Path,
    overrides: dict[str, object] | None,
    *,
    ftl_L: float | None = None,
) -> ComovingMetrics | None:
    """Estimate co-moving stability from shell chi time series and t=0 shift profile."""
    if overrides is None:
        return None

    beta_mean = _mean_beta_in_bubble(overrides, ftl_L=ftl_L)
    if beta_mean is not None and abs(beta_mean) < STATIONARY_BETA_EPS:
        return ComovingMetrics(
            beta_mean=beta_mean,
            delta_comoving=None,
            score=None,
            stationary=True,
        )

    shell_path = episode_dir / "small_data" / "shell_profiles.dat"
    if not shell_path.exists():
        shell_path = episode_dir / "shell_profiles.dat"
    parsed = _parse_shell_profiles(shell_path)

    chi_series: list[float] | None = None
    times: list[float] | None = None
    if parsed is not None:
        times, cols = parsed
        chi_keys = sorted(key for key in cols if key.startswith("chi_mean_"))
        if chi_keys:
            chi_series = cols[chi_keys[0]]

    if chi_series is None or times is None or len(chi_series) < 2:
        collapse_path = episode_dir / "data" / "collapse_diagnostics.dat"
        if not collapse_path.exists():
            collapse_path = episode_dir / "collapse_diagnostics.dat"
        rows = _numeric_rows(collapse_path, 3)
        if len(rows) < 2:
            return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)
        times = [row[0] for row in rows]
        chi_series = [row[2] for row in rows]

    if beta_mean is None or len(chi_series) < 2 or len(times) < 2:
        return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)

    times_arr = np.asarray(times, dtype=float)
    chi_arr = np.asarray(chi_series, dtype=float)
    if times_arr[-1] <= times_arr[0]:
        return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)

    dchi_dt = np.gradient(chi_arr, times_arr)
    eulerian_rate = float(np.max(np.abs(dchi_dt)))
    shift_scale = max(abs(beta_mean), 1.0e-3)
    delta_comoving = eulerian_rate / shift_scale
    score = 1.0 / (1.0 + delta_comoving)
    return ComovingMetrics(beta_mean=beta_mean, delta_comoving=delta_comoving, score=score)


def _log_growth_rate(
    times: list[float],
    values: list[float],
    *,
    floor: float,
) -> float | None:
    """Least-squares exponential growth rate lambda from ln(value) vs time."""
    t = np.asarray(times, dtype=float)
    v = np.asarray(values, dtype=float)
    if t.size < 2 or float(np.ptp(t)) <= 0.0:
        return None
    v = np.maximum(np.abs(v), floor)
    finite = np.isfinite(v) & np.isfinite(t)
    if int(np.count_nonzero(finite)) < 2:
        return None
    slope = float(np.polyfit(t[finite], np.log(v[finite]), 1)[0])
    return slope


def read_growth_metrics(
    collapse_path: Path,
    constraint_path: Path,
    *,
    sigma_lambda: float = 0.5,
) -> GrowthMetrics | None:
    collapse_rows = _numeric_rows(collapse_path, 4)
    constraint_rows = _numeric_rows(constraint_path, 3)

    lambda_k = lambda_inv_chi = lambda_ham = None

    if len(collapse_rows) >= 2:
        times = [row[0] for row in collapse_rows]
        k_max = [abs(row[3]) for row in collapse_rows]
        chi_min = [row[2] for row in collapse_rows]
        lambda_k = _log_growth_rate(times, k_max, floor=1.0e-3)
        inv_chi = [1.0 / max(abs(c), 1.0e-8) for c in chi_min]
        lambda_inv_chi = _log_growth_rate(times, inv_chi, floor=1.0e-8)

    if len(constraint_rows) >= 2:
        times = [row[0] for row in constraint_rows]
        ham = [row[1] for row in constraint_rows]
        lambda_ham = _log_growth_rate(times, ham, floor=1.0e-12)

    candidates = [v for v in (lambda_ham, lambda_k, lambda_inv_chi) if v is not None]
    if not candidates:
        return None

    lambda_eff = max(candidates)
    s_growth = 1.0 / (1.0 + max(0.0, lambda_eff) / max(sigma_lambda, 1.0e-12))

    return GrowthMetrics(
        lambda_hamiltonian=lambda_ham,
        lambda_max_k=lambda_k,
        lambda_inv_chi=lambda_inv_chi,
        lambda_effective=lambda_eff,
        s_growth=s_growth,
    )


def read_energy_condition_metrics(path: Path) -> EnergyConditionMetrics | None:
    rows = _numeric_rows(path, 6)
    if not rows:
        return None
    return EnergyConditionMetrics(
        final_time=rows[-1][0],
        min_nec=min(row[1] for row in rows),
        min_wec=min(row[2] for row in rows),
        min_sec=min(row[3] for row in rows),
        min_dec=min(row[4] for row in rows),
        max_integral_nec_violation=max(row[5] for row in rows),
    )


def read_curvature_invariant_metrics(path: Path) -> CurvatureInvariantMetrics | None:
    rows = _numeric_rows(path, 5)
    if not rows:
        return None
    return CurvatureInvariantMetrics(
        final_time=rows[-1][0],
        max_abs_ricci_scalar=max(row[1] for row in rows),
        max_ricci_tensor_sq=max(row[2] for row in rows),
        max_kij_sq=max(row[3] for row in rows),
        max_l2_ricci_scalar=max(row[4] for row in rows),
    )


def read_episode_metrics(
    episode_dir: Path,
    *,
    ftl_L: float | None = None,
) -> EpisodeMetrics:
    data_dir = episode_dir / "data"
    small_data_dir = episode_dir / "small_data"
    collapse_path = data_dir / "collapse_diagnostics.dat"
    constraint_path = data_dir / "constraint_norms.dat"
    areal_path = small_data_dir / "areal_radius.dat"
    if not collapse_path.exists():
        collapse_path = episode_dir / "collapse_diagnostics.dat"
    if not constraint_path.exists():
        constraint_path = episode_dir / "constraint_norms.dat"
    if not areal_path.exists():
        areal_path = episode_dir / "areal_radius.dat"

    energy_conditions_path = data_dir / "energy_conditions.dat"
    if not energy_conditions_path.exists():
        energy_conditions_path = episode_dir / "energy_conditions.dat"
    curvature_path = data_dir / "curvature_invariants.dat"
    if not curvature_path.exists():
        curvature_path = episode_dir / "curvature_invariants.dat"

    collapse = read_collapse_metrics(collapse_path)
    constraints = read_constraint_metrics(constraint_path)
    stability = read_stability_metrics(collapse_path, areal_path)
    growth = read_growth_metrics(collapse_path, constraint_path)
    energy_conditions = read_energy_condition_metrics(energy_conditions_path)
    curvature = read_curvature_invariant_metrics(curvature_path)

    # Evolved-spacetime operational FTL: run the search on the latest plotfile
    # if one survived (requires the full metric in amr.plot_vars).  Falls back
    # silently to None when no plotfile is present (e.g. streamed + deleted).
    general_ftl_evolved = None
    try:
        plotfile = find_latest_plotfile(episode_dir)
        if plotfile is not None:
            general_ftl_evolved = compute_general_ftl_from_plotfile(
                plotfile, n=129, L=ftl_L
            )
    except Exception:
        general_ftl_evolved = None

    ftl = None
    comoving = None
    physical = None
    general_ftl = None
    overrides = load_overrides_from_episode(episode_dir)
    if overrides:
        try:
            ftl = compute_ftl_metrics(overrides, L=ftl_L)
        except Exception:
            ftl = None
        try:
            general_ftl = compute_general_ftl(overrides, L=ftl_L, n=97)
        except Exception:
            general_ftl = None
        try:
            comoving = read_comoving_metrics(episode_dir, overrides, ftl_L=ftl_L)
        except Exception:
            comoving = None
        try:
            physical = compute_physical_metrics(overrides, L=ftl_L)
        except Exception:
            physical = None

    if collapse is None and constraints is None:
        reason = "missing_diagnostics"
    else:
        reason = "completed_or_partial"

    return EpisodeMetrics(
        collapse=collapse,
        constraints=constraints,
        stability=stability,
        comoving=comoving,
        ftl=ftl,
        termination_reason=reason,
        growth=growth,
        physical=physical,
        energy_conditions=energy_conditions,
        curvature=curvature,
        general_ftl=general_ftl,
        general_ftl_evolved=general_ftl_evolved,
    )


def dataclass_to_dict(value: object) -> object:
    if hasattr(value, "__dataclass_fields__"):
        return {key: dataclass_to_dict(getattr(value, key)) for key in value.__dataclass_fields__}
    if isinstance(value, dict):
        return {key: dataclass_to_dict(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [dataclass_to_dict(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        # numpy scalar (bool_/integer/floating) -> native Python scalar
        return value.item()
    return value
