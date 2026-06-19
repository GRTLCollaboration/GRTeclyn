"""Central radial profile metrics for splash diagnostics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CentralRadialProfileMetrics:
    """Aggregates from ``small_data/central_radial_profile.dat``."""

    n_frames: int
    t: tuple[float, ...]
    peak_radius: float
    splash_width: float
    compression_ratio: float
    infall_speed: float | None
    cusp_unresolved: bool
    dx_finest: float
    initial_rho_at_ring: float
