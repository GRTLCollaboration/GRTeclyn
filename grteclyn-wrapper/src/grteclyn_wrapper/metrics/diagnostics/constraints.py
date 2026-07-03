"""Constraint norms from ``constraint_norms.dat``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from statistics import median

from ..io.dat import numeric_rows
from ..types.diagnostics import ConstraintMetrics

# Detect collapse-driven constraint cliffs (see LabJournal gw_beam_qd100_v3).
SPIKE_HAM_ABS = 1.0
SPIKE_HAM_RATIO = 100.0
SPIKE_STEP_HAM_RATIO = 50.0
SPIKE_MOM_ABS = 0.1
SPIKE_MOM_RATIO = 20.0


@dataclass(frozen=True)
class ConstraintSpikeInfo:
    constraint_spike_time: float | None
    ham_spike_ratio: float
    max_step_ham_ratio: float
    mom_spike_ratio: float
    has_constraint_spike: bool


def spike_info_from_rows(rows: list[list[float]]) -> ConstraintSpikeInfo:
    """Detect a Ham/Mom explosion from the constraint time series."""
    if len(rows) < 2:
        return ConstraintSpikeInfo(None, 1.0, 1.0, 1.0, False)

    times = [row[0] for row in rows]
    hams = [row[1] for row in rows]
    moms = [row[2] for row in rows]

    baseline_n = max(3, min(len(rows) // 4, 100))
    baseline_ham = max(median(hams[:baseline_n]), 1.0e-6)
    baseline_mom = max(median(moms[:baseline_n]), 1.0e-8)

    max_ham = max(hams)
    max_mom = max(moms)
    ham_ratio = max_ham / baseline_ham
    mom_ratio = max_mom / baseline_mom

    max_step_ham_ratio = 1.0
    spike_time: float | None = None
    for i in range(len(rows) - 1):
        step_ratio = hams[i + 1] / max(hams[i], 1.0e-12)
        if step_ratio > max_step_ham_ratio:
            max_step_ham_ratio = step_ratio
        if step_ratio >= SPIKE_STEP_HAM_RATIO and spike_time is None:
            spike_time = times[i + 1]

    if spike_time is None:
        for i, ham in enumerate(hams):
            if ham >= SPIKE_HAM_ABS or ham >= baseline_ham * SPIKE_HAM_RATIO:
                spike_time = times[i]
                break

    has_spike = (
        max_ham >= SPIKE_HAM_ABS
        or ham_ratio >= SPIKE_HAM_RATIO
        or max_step_ham_ratio >= SPIKE_STEP_HAM_RATIO
        or max_mom >= SPIKE_MOM_ABS
        or mom_ratio >= SPIKE_MOM_RATIO
    )

    return ConstraintSpikeInfo(
        constraint_spike_time=spike_time,
        ham_spike_ratio=ham_ratio,
        max_step_ham_ratio=max_step_ham_ratio,
        mom_spike_ratio=mom_ratio,
        has_constraint_spike=has_spike,
    )


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
    spike = spike_info_from_rows(rows)

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
        constraint_spike_time=spike.constraint_spike_time,
        ham_spike_ratio=spike.ham_spike_ratio,
        max_step_ham_ratio=spike.max_step_ham_ratio,
        mom_spike_ratio=spike.mom_spike_ratio,
        has_constraint_spike=spike.has_constraint_spike,
    )
