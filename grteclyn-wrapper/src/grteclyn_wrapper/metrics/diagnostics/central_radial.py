"""Parser for ``small_data/central_radial_profile.dat``."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np

from ..types.central_radial import CentralRadialProfileMetrics

_TIME_RE = re.compile(r"^#\s*time=([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")


def _parse_blocks(text: str) -> list[tuple[float, list[tuple[float, float, float, float]]]]:
    blocks: list[tuple[float, list[tuple[float, float, float, float]]]] = []
    current_time: float | None = None
    current_rows: list[tuple[float, float, float, float]] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        time_match = _TIME_RE.match(line)
        if time_match:
            if current_time is not None and current_rows:
                blocks.append((current_time, current_rows))
            current_time = float(time_match.group(1))
            current_rows = []
            continue
        if line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        current_rows.append(
            (float(parts[0]), float(parts[1]), float(parts[2]), float(parts[3]))
        )

    if current_time is not None and current_rows:
        blocks.append((current_time, current_rows))
    return blocks


def _fwhm(radii: np.ndarray, values: np.ndarray) -> float:
    peak_idx = int(np.argmax(values))
    peak = float(values[peak_idx])
    if peak <= 0.0 or not math.isfinite(peak):
        return 0.0
    half = 0.5 * peak
    above = values >= half
    if not np.any(above):
        return 0.0
    r_above = radii[above]
    return float(r_above[-1] - r_above[0]) if r_above.size > 1 else 0.0


def _frame_metrics(
    rows: list[tuple[float, float, float, float]],
    *,
    dx_finest: float,
) -> tuple[float, float, float, bool, float]:
    radii = np.asarray([row[0] for row in rows], dtype=float)
    rho = np.asarray([row[1] for row in rows], dtype=float)
    mask = np.isfinite(radii) & np.isfinite(rho)
    radii = radii[mask]
    rho = rho[mask]
    if radii.size == 0:
        return 0.0, 0.0, 0.0, False, 0.0
    peak_idx = int(np.argmax(rho))
    peak_radius = float(radii[peak_idx])
    splash_width = _fwhm(radii, rho)
    rho_peak = float(rho[peak_idx])
    rho_outer = float(rho[-1]) if rho.size else 0.0
    compression = rho_peak / max(rho_outer, 1.0e-12)
    cusp_unresolved = splash_width > 0.0 and splash_width < 2.0 * dx_finest
    initial_ring = float(np.max(rho))
    return peak_radius, splash_width, compression, cusp_unresolved, initial_ring


def read_central_radial_profile(
    path: Path,
    *,
    dx_finest: float = 0.25,
    ring_radius: float | None = None,
) -> CentralRadialProfileMetrics | None:
    if not path.is_file():
        return None
    blocks = _parse_blocks(path.read_text(encoding="utf-8"))
    if not blocks:
        return None

    times: list[float] = []
    peak_radii: list[float] = []
    splash_widths: list[float] = []
    compressions: list[float] = []
    cusp_flags: list[bool] = []
    initial_ring = 0.0

    for frame_idx, (time, rows) in enumerate(blocks):
        peak_radius, splash_width, compression, cusp_unresolved, ring_rho = _frame_metrics(
            rows,
            dx_finest=dx_finest,
        )
        times.append(time)
        peak_radii.append(peak_radius)
        splash_widths.append(splash_width)
        compressions.append(compression)
        cusp_flags.append(cusp_unresolved)
        if frame_idx == 0:
            if ring_radius is not None:
                radii = np.asarray([row[0] for row in rows], dtype=float)
                rho = np.asarray([row[1] for row in rows], dtype=float)
                if radii.size:
                    idx = int(np.argmin(np.abs(radii - float(ring_radius))))
                    initial_ring = float(rho[idx]) if math.isfinite(float(rho[idx])) else ring_rho
            else:
                initial_ring = ring_rho

    peak_frame = int(np.argmax(compressions))
    infall_speed: float | None = None
    if len(peak_radii) >= 2:
        dt = times[peak_frame] - times[max(0, peak_frame - 1)]
        if dt > 0.0:
            dr = peak_radii[peak_frame - 1] - peak_radii[peak_frame]
            infall_speed = float(dr / dt)

    return CentralRadialProfileMetrics(
        n_frames=len(times),
        t=tuple(times),
        peak_radius=float(peak_radii[peak_frame]),
        splash_width=float(splash_widths[peak_frame]),
        compression_ratio=float(compressions[peak_frame]),
        infall_speed=infall_speed,
        cusp_unresolved=bool(cusp_flags[peak_frame]),
        dx_finest=float(dx_finest),
        initial_rho_at_ring=float(initial_ring),
    )
