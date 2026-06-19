"""Central-origin field metrics for splash / critical-collapse campaigns."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class CentralFieldMetrics:
    """Aggregates from ``small_data/central_timeseries.dat`` at the grid origin."""

    n_frames: int
    t: tuple[float, ...]
    rho_req: tuple[float, ...]
    lapse: tuple[float, ...]
    scalar_activity: tuple[float, ...]
    peak_rho_req_at_origin: float
    peak_rho_req_time: float | None
    initial_rho_req_at_origin: float
    min_lapse_at_origin: float
    wave_chromaticity: float

    @property
    def focusing_efficiency(self) -> float:
        denom = max(self.initial_rho_req_at_origin, 1.0e-12)
        return float(self.peak_rho_req_at_origin / denom)
