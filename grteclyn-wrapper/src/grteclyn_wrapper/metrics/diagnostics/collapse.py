"""Collapse diagnostics from ``collapse_diagnostics.dat``."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import CollapseMetrics


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

    phi_min = min((row[10] for row in rows if len(row) >= 14), default=None)
    phi_max = max((row[11] for row in rows if len(row) >= 14), default=None)
    pi_min = min((row[12] for row in rows if len(row) >= 14), default=None)
    pi_max = max((row[13] for row in rows if len(row) >= 14), default=None)
    bary_x = rows[-1][14] if len(rows[-1]) >= 15 else None
    bary_y = rows[-1][15] if len(rows[-1]) >= 16 else None
    bary_z = rows[-1][16] if len(rows[-1]) >= 17 else None
    rho_sum = rows[-1][17] if len(rows[-1]) >= 18 else None

    return CollapseMetrics(
        final_time=final_time,
        min_lapse=min_lapse,
        min_chi=min_chi,
        max_abs_k=max_abs_k,
        max_horizon_radius=max_horizon_radius,
        min_theta_plus=min_theta_plus,
        r_at_min_theta_plus=r_at_min_theta_plus,
        scalar_phi_range=(phi_max - phi_min) if phi_min is not None and phi_max is not None else None,
        scalar_pi_range=(pi_max - pi_min) if pi_min is not None and pi_max is not None else None,
        barycenter_x=bary_x,
        barycenter_y=bary_y,
        barycenter_z=bary_z,
        rho_sum=rho_sum,
    )
