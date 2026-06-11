"""Constraint/collapse growth-rate metrics."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io.dat import numeric_rows
from ..types.diagnostics import GrowthMetrics


def _log_growth_rate(
    times: list[float],
    values: list[float],
    *,
    floor: float,
) -> float | None:
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
    collapse_rows = numeric_rows(collapse_path, 4)
    constraint_rows = numeric_rows(constraint_path, 3)

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
