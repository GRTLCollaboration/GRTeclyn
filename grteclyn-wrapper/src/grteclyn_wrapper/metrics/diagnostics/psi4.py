"""Read directional Psi4 metrics produced by the plotfile consumer."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..types.psi4 import Psi4Metrics


_PSI4_DIRECTIONAL_FILENAME = "psi4_directional.dat"


def read_psi4_metrics(small_data_dir: Path) -> Psi4Metrics | None:
    """Aggregate directional Psi4 timeseries from ``small_data/psi4_directional.dat``.

    Columns: time, P_total, P_z_beam, beam_ratio.
    """
    path = small_data_dir / _PSI4_DIRECTIONAL_FILENAME
    if not path.is_file():
        return None

    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                rows.append([float(v) for v in stripped.split()])
            except ValueError:
                continue

    if not rows:
        return None

    rows.sort(key=lambda r: r[0])
    data = np.asarray(rows, dtype=float)
    if data.ndim != 2 or data.shape[1] < 4:
        return None

    p_total = data[:, 1]
    p_z_beam = data[:, 2]
    beam_ratio = data[:, 3]

    def _peak(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite.max()) if finite.size else 0.0

    def _mean(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite.mean()) if finite.size else 0.0

    def _final(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite[-1]) if finite.size else 0.0

    return Psi4Metrics(
        peak_total_power=_peak(p_total),
        peak_z_beam_power=_peak(p_z_beam),
        peak_beam_ratio=_peak(beam_ratio),
        mean_total_power=_mean(p_total),
        mean_z_beam_power=_mean(p_z_beam),
        mean_beam_ratio=_mean(beam_ratio),
        final_total_power=_final(p_total),
        final_z_beam_power=_final(p_z_beam),
        final_beam_ratio=_final(beam_ratio),
        n_samples=int(data.shape[0]),
    )
