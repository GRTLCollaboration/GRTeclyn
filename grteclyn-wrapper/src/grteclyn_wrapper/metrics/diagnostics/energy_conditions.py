"""Matter-sector energy conditions from ``energy_conditions.dat``."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import EnergyConditionMetrics


def read_energy_condition_metrics(path: Path) -> EnergyConditionMetrics | None:
    rows = numeric_rows(path, 6)
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
