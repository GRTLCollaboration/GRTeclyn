"""Parser for ``small_data/central_timeseries.dat``."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from ..io.dat import numeric_rows
from ..types.central import CentralFieldMetrics

CENTRAL_TIMESERIES_HEADER = (
    "# time  rho_req  lapse  scalar_activity  phi_re  phi_im  "
    "noether_charge  phase_coherence  ham_abs  mom_abs"
)

_CHROMATICITY_MIN_FRAMES = 4


def _peak_with_time(
    values: tuple[float, ...],
    times: tuple[float, ...],
) -> tuple[float, float | None]:
    peak = float("-inf")
    t_peak: float | None = None
    for value, time in zip(values, times):
        if not math.isfinite(value):
            continue
        if value > peak:
            peak = value
            t_peak = time
    if peak == float("-inf"):
        return 0.0, None
    return peak, t_peak


def _wave_chromaticity(activity: tuple[float, ...], times: tuple[float, ...]) -> float:
    if len(activity) < _CHROMATICITY_MIN_FRAMES:
        return 0.0
    arr = np.asarray(activity, dtype=float)
    t_arr = np.asarray(times, dtype=float)
    mask = np.isfinite(arr) & np.isfinite(t_arr)
    arr = arr[mask]
    t_arr = t_arr[mask]
    if len(arr) < _CHROMATICITY_MIN_FRAMES:
        return 0.0
    dt = float(np.mean(np.diff(t_arr))) if len(t_arr) > 1 else 1.0
    if not math.isfinite(dt) or dt <= 0.0:
        dt = 1.0
    centered = arr - float(np.mean(arr))
    if float(np.max(np.abs(centered))) < 1.0e-14:
        return 0.0
    spectrum = np.abs(np.fft.rfft(centered))
    freqs = np.fft.rfftfreq(len(centered), d=dt)
    if len(spectrum) <= 1:
        return 0.0
    spectrum = spectrum[1:]
    freqs = freqs[1:]
    if len(spectrum) == 0:
        return 0.0
    peak_idx = int(np.argmax(spectrum))
    peak_freq = float(freqs[peak_idx])
    half_max = 0.5 * float(spectrum[peak_idx])
    above = spectrum >= half_max
    if not np.any(above):
        return 0.0
    fwhm = float(freqs[above][-1] - freqs[above][0])
    if fwhm <= 0.0:
        fwhm = float(freqs[1] - freqs[0]) if len(freqs) > 1 else 1.0
    q = peak_freq / max(fwhm, 1.0e-6)
    return float(np.clip(q / (q + 1.0), 0.0, 1.0))


def _optional_column(rows: list[list[float]], index: int) -> tuple[float, ...]:
    out: list[float] = []
    for row in rows:
        if len(row) > index:
            out.append(float(row[index]))
        else:
            out.append(0.0)
    return tuple(out)


def _build_central_metrics(
    rows: list[list[float]],
    *,
    initial_rho_baseline: float | None = None,
) -> CentralFieldMetrics | None:
    if not rows:
        return None

    t = [float(row[0]) for row in rows]
    rho = [float(row[1]) for row in rows]
    lapse = [float(row[2]) for row in rows]
    activity = [float(row[3]) for row in rows]

    t_tuple = tuple(t)
    rho_tuple = tuple(rho)
    lapse_tuple = tuple(lapse)
    activity_tuple = tuple(activity)

    peak_rho, peak_time = _peak_with_time(rho_tuple, t_tuple)
    initial_rho = rho[0] if rho else 0.0
    finite_lapse = [value for value in lapse if math.isfinite(value)]
    min_lapse = min(finite_lapse) if finite_lapse else 1.0
    baseline = initial_rho_baseline if initial_rho_baseline is not None else initial_rho

    return CentralFieldMetrics(
        n_frames=len(t),
        t=t_tuple,
        rho_req=rho_tuple,
        lapse=lapse_tuple,
        scalar_activity=activity_tuple,
        peak_rho_req_at_origin=peak_rho,
        peak_rho_req_time=peak_time,
        initial_rho_req_at_origin=initial_rho,
        min_lapse_at_origin=min_lapse,
        wave_chromaticity=_wave_chromaticity(activity_tuple, t_tuple),
        initial_rho_baseline=max(float(baseline), float(initial_rho)),
        noether_charge=_optional_column(rows, 6),
        phase_coherence=_optional_column(rows, 7),
        ham_abs=_optional_column(rows, 8),
        mom_abs=_optional_column(rows, 9),
    )


def read_central_field_metrics(
    path: Path,
    *,
    initial_rho_baseline: float | None = None,
) -> CentralFieldMetrics | None:
    if not path.is_file():
        return None
    rows = numeric_rows(path, min_columns=4)
    return _build_central_metrics(rows, initial_rho_baseline=initial_rho_baseline)


def read_prefix_central_field_metrics(
    path: Path,
    at_time: float,
    *,
    initial_rho_baseline: float | None = None,
) -> CentralFieldMetrics | None:
    if not path.is_file():
        return None
    rows = [row for row in numeric_rows(path, min_columns=4) if float(row[0]) <= at_time + 1.0e-12]
    return _build_central_metrics(rows, initial_rho_baseline=initial_rho_baseline)
