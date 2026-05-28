from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


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


@dataclass(frozen=True)
class EpisodeMetrics:
    collapse: CollapseMetrics | None
    constraints: ConstraintMetrics | None
    stability: StabilityMetrics | None
    termination_reason: str


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


def read_episode_metrics(episode_dir: Path) -> EpisodeMetrics:
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

    collapse = read_collapse_metrics(collapse_path)
    constraints = read_constraint_metrics(constraint_path)
    stability = read_stability_metrics(collapse_path, areal_path)

    if collapse is None and constraints is None:
        reason = "missing_diagnostics"
    else:
        reason = "completed_or_partial"

    return EpisodeMetrics(
        collapse=collapse,
        constraints=constraints,
        stability=stability,
        termination_reason=reason,
    )


def dataclass_to_dict(value: object) -> object:
    if hasattr(value, "__dataclass_fields__"):
        return {key: dataclass_to_dict(getattr(value, key)) for key in value.__dataclass_fields__}
    return value
