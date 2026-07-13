"""Collapse diagnostics from ``collapse_diagnostics.dat``."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import CollapseMetrics

# A trapped-surface verdict requires theta+ < threshold AND lapse collapse at
# the same timestep.  theta+ <= 0 with a healthy lapse (common in exotic warp
# channels) is an uncorroborated coordinate artifact, not a traversable horizon.
TRAPPED_THETA_THRESH: float = -0.05
TRAPPED_LAPSE_THRESH: float = 0.2


def _endpoints(rows: list[list[float]], column: int) -> tuple[float | None, float | None]:
    values = [row[column] for row in rows if len(row) > column]
    return (values[0], values[-1]) if values else (None, None)


def _relative_drift(initial: float | None, final: float | None) -> float | None:
    if initial is None or final is None or abs(initial) <= 1.0e-15:
        return None
    return abs(final - initial) / abs(initial)


def _retention(initial: float | None, final: float | None) -> float | None:
    if initial is None or final is None or abs(initial) <= 1.0e-15:
        return None
    return abs(final) / abs(initial)


def _corroborated_trapped(rows: list[list[float]]) -> tuple[bool, float | None]:
    """Return whether any row has corroborated trapped-surface signal and when."""
    first_time: float | None = None
    for row in rows:
        if len(row) < 9:
            continue
        lapse = row[1]
        theta = row[8]
        if theta < TRAPPED_THETA_THRESH and lapse < TRAPPED_LAPSE_THRESH:
            t = row[0]
            if first_time is None or t < first_time:
                first_time = t
    return first_time is not None, first_time


def read_collapse_metrics(path: Path) -> CollapseMetrics | None:
    rows = numeric_rows(path, 4)
    if not rows:
        return None

    final_time = rows[-1][0]
    min_lapse = min(row[1] for row in rows)
    min_chi = min(row[2] for row in rows)
    max_abs_k = max(row[3] for row in rows)
    max_horizon_radius = max((row[7] for row in rows if len(row) >= 8), default=None)
    theta_rows = [row for row in rows if len(row) >= 10]
    min_theta_plus = min((row[8] for row in theta_rows), default=None)
    r_at_min_theta_plus = None
    if theta_rows:
        row_at_min = min(theta_rows, key=lambda row: row[8])
        r_at_min_theta_plus = row_at_min[9]

    corroborated_trapped, first_corroborated_time = _corroborated_trapped(rows)

    phi_min = min((row[10] for row in rows if len(row) >= 14), default=None)
    phi_max = max((row[11] for row in rows if len(row) >= 14), default=None)
    pi_min = min((row[12] for row in rows if len(row) >= 14), default=None)
    pi_max = max((row[13] for row in rows if len(row) >= 14), default=None)
    bary_x = rows[-1][14] if len(rows[-1]) >= 15 else None
    bary_y = rows[-1][15] if len(rows[-1]) >= 16 else None
    bary_z = rows[-1][16] if len(rows[-1]) >= 17 else None
    rho_sum = rows[-1][17] if len(rows[-1]) >= 18 else None
    j_z_initial, j_z_final = _endpoints(rows, 18)
    q_total_initial, q_total_final = _endpoints(rows, 19)
    q_sphere_initial, q_sphere_final = _endpoints(rows, 20)
    rho_sphere_initial, rho_sphere_final = _endpoints(rows, 21)

    return CollapseMetrics(
        final_time=final_time,
        min_lapse=min_lapse,
        min_chi=min_chi,
        max_abs_k=max_abs_k,
        max_horizon_radius=max_horizon_radius,
        min_theta_plus=min_theta_plus,
        r_at_min_theta_plus=r_at_min_theta_plus,
        corroborated_trapped=corroborated_trapped,
        first_corroborated_time=first_corroborated_time,
        scalar_phi_range=(phi_max - phi_min) if phi_min is not None and phi_max is not None else None,
        scalar_pi_range=(pi_max - pi_min) if pi_min is not None and pi_max is not None else None,
        barycenter_x=bary_x,
        barycenter_y=bary_y,
        barycenter_z=bary_z,
        rho_sum=rho_sum,
        j_z_initial=j_z_initial,
        j_z_final=j_z_final,
        q_total_initial=q_total_initial,
        q_total_final=q_total_final,
        q_total_relative_drift=_relative_drift(q_total_initial, q_total_final),
        q_sphere_initial=q_sphere_initial,
        q_sphere_final=q_sphere_final,
        q_sphere_retention=_retention(q_sphere_initial, q_sphere_final),
        rho_sphere_initial=rho_sphere_initial,
        rho_sphere_final=rho_sphere_final,
        rho_sphere_retention=_retention(rho_sphere_initial, rho_sphere_final),
    )
