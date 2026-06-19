"""Central ball-averaged field metrics for splash / critical-collapse campaigns."""

from __future__ import annotations

import math

from dataclasses import dataclass


@dataclass(frozen=True)
class CentralFieldMetrics:
    """Aggregates from ``small_data/central_timeseries.dat`` (central ball average)."""

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
    initial_rho_baseline: float = 0.0
    noether_charge: tuple[float, ...] = ()
    phase_coherence: tuple[float, ...] = ()
    ham_abs: tuple[float, ...] = ()
    mom_abs: tuple[float, ...] = ()

    @property
    def focusing_efficiency(self) -> float:
        baseline = self.initial_rho_baseline if self.initial_rho_baseline > 0.0 else self.initial_rho_req_at_origin
        denom = max(baseline, 1.0e-12)
        return float(self.peak_rho_req_at_origin / denom)

    @property
    def ham_abs_at_peak(self) -> float | None:
        if not self.rho_req or not self.ham_abs:
            return None
        peak_idx = 0
        peak_val = float("-inf")
        for idx, value in enumerate(self.rho_req):
            if math.isfinite(value) and value > peak_val:
                peak_val = value
                peak_idx = idx
        if peak_idx >= len(self.ham_abs):
            return None
        value = self.ham_abs[peak_idx]
        return value if math.isfinite(value) else None
