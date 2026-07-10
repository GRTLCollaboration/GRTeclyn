"""Read directional Psi4 metrics produced by the plotfile consumer."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..types.psi4 import Psi4Metrics


_PSI4_DIRECTIONAL_FILENAME = "psi4_directional.dat"

# Wave-zone validity: relative std of r*|Psi4| across radii must be below this.
_WAVEZONE_ONE_OVER_R_TOL = 0.3


def read_psi4_metrics(
    small_data_dir: Path,
    *,
    max_peak_time: float | None = None,
    min_valid_time: float | None = None,
) -> Psi4Metrics | None:
    """Aggregate directional Psi4 timeseries from ``small_data/psi4_directional.dat``.

    Columns (v4): time, P_total, P_z_beam, beam_ratio
    Columns (v5): time, P_total, P_z_beam, beam_ratio, beaming_gain, wavezone_std

    When ``max_peak_time`` is set (typically the constraint-spike time), **all**
    aggregates (peak, mean, final) use only samples at or before that time so
    post-collapse numerical noise cannot inflate GW metrics.

    When ``min_valid_time`` is set (typically R_ext + margin, the junk-radiation
    window), samples before that time are excluded from all aggregates.
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

    times = data[:, 0]
    p_total = data[:, 1]
    p_z_beam = data[:, 2]
    beam_ratio = data[:, 3]
    # v5 extended columns (backward compatible — missing columns default to 1.0/0.0)
    has_v5 = data.shape[1] >= 6
    beaming_gain = data[:, 4] if has_v5 else np.ones_like(p_total)
    wavezone_std = data[:, 5] if has_v5 else np.zeros_like(p_total)

    # Apply time-window filters.
    mask = np.ones(len(times), dtype=bool)
    if min_valid_time is not None:
        mask &= times >= min_valid_time - 1.0e-9
    if max_peak_time is not None:
        mask &= times <= max_peak_time + 1.0e-9

    p_total = p_total[mask]
    p_z_beam = p_z_beam[mask]
    beam_ratio = beam_ratio[mask]
    beaming_gain = beaming_gain[mask]
    wavezone_std = wavezone_std[mask]
    times = times[mask]

    if p_total.size == 0:
        return None

    def _peak(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite.max()) if finite.size else 0.0

    def _mean(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite.mean()) if finite.size else 0.0

    def _final(arr: np.ndarray) -> float:
        finite = arr[np.isfinite(arr)]
        return float(finite[-1]) if finite.size else 0.0

    mean_wavezone_std = _mean(wavezone_std)
    wavezone_ok = mean_wavezone_std < _WAVEZONE_ONE_OVER_R_TOL

    # Direction stability: if beaming_gain varies little, direction is stable.
    # Approximate as 1 - coefficient_of_variation of beaming_gain.
    bg_finite = beaming_gain[np.isfinite(beaming_gain)]
    if bg_finite.size > 1 and bg_finite.mean() > 0:
        cv = float(bg_finite.std() / bg_finite.mean())
        direction_stability = float(max(0.0, 1.0 - cv))
    else:
        direction_stability = 0.0

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
        n_samples=int(p_total.size),
        mean_beaming_gain=_mean(beaming_gain),
        peak_beaming_gain=_peak(beaming_gain),
        wavezone_ok=wavezone_ok,
        wavezone_one_over_r_std=mean_wavezone_std,
        direction_stability=direction_stability,
    )
