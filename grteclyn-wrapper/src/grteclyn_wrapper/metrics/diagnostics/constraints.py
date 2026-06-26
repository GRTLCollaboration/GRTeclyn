"""Constraint norms from ``constraint_norms.dat``."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import ConstraintMetrics


def read_constraint_metrics(path: Path) -> ConstraintMetrics | None:
    rows = numeric_rows(path, 3)
    if not rows:
        return None

    rho_rows = [row for row in rows if len(row) >= 6]
    min_rho_req = min((row[3] for row in rho_rows), default=None)
    max_rho_req = max((row[4] for row in rho_rows), default=None)
    max_int_neg = max((row[5] for row in rho_rows), default=None)
    final_peak_rho_req = rho_rows[-1][4] if rho_rows else None
    initial_peak_rho_req = rho_rows[0][4] if rho_rows else None

    return ConstraintMetrics(
        final_time=rows[-1][0],
        max_hamiltonian_l2=max(row[1] for row in rows),
        max_momentum_l2=max(row[2] for row in rows),
        final_hamiltonian_l2=rows[-1][1],
        final_momentum_l2=rows[-1][2],
        min_rho_required=min_rho_req,
        max_rho_required=max_rho_req,
        integral_negative_rho=max_int_neg,
        final_peak_rho_required=final_peak_rho_req,
        initial_peak_rho_required=initial_peak_rho_req,
    )
