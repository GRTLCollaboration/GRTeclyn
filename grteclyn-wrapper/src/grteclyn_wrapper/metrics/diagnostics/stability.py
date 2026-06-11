"""Stability metrics from collapse and areal-radius diagnostics."""

from __future__ import annotations

from pathlib import Path

from ..io.dat import numeric_rows
from ..types.diagnostics import StabilityMetrics
from ._helpers import positive_fractional_change, positive_fractional_drop


def read_stability_metrics(collapse_path: Path, areal_path: Path) -> StabilityMetrics | None:
    collapse_rows = numeric_rows(collapse_path, 4)
    areal_rows = numeric_rows(areal_path, 2)
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
        lapse_drop_fraction = positive_fractional_drop(initial_lapse, minimum_lapse, floor=1.0e-6)

        initial_chi = collapse_rows[0][2]
        minimum_chi = min(row[2] for row in collapse_rows)
        chi_drop_fraction = positive_fractional_drop(initial_chi, minimum_chi, floor=1.0e-8)

        initial_k = collapse_rows[0][3]
        maximum_k = max(row[3] for row in collapse_rows)
        k_growth_fraction = positive_fractional_change(maximum_k, initial_k, floor=1.0e-1)

        horizon_values = [row[7] for row in collapse_rows if len(row) >= 8]
        if horizon_values:
            horizon_growth_fraction = positive_fractional_change(
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
        value
        for value in (
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
