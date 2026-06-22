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
    # Geometric (spacetime-splash) timeseries at the central ball:
    #   chi     -> conformal factor (1 = flat; drops as curvature concentrates)
    #   trace_K -> trace of extrinsic curvature (gravitational "crunch" rate)
    #   weyl4   -> Re(Psi4) Weyl scalar = gravitational-wave content at center
    chi: tuple[float, ...] = ()
    trace_K: tuple[float, ...] = ()
    weyl4: tuple[float, ...] = ()

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

    # --- Spacetime-splash (geometric) properties ----------------------------

    @property
    def has_geometric_data(self) -> bool:
        """True if the timeseries carries chi data (legacy files lack it)."""
        return any(math.isfinite(v) and v != 0.0 for v in self.chi)

    @property
    def min_chi_at_origin(self) -> float:
        """Minimum conformal factor at center over the evolution.

        chi starts near 1 (flat) and drops toward 0 as spacetime curvature
        concentrates at the center -- this is the geometric splash signature.
        """
        finite = [v for v in self.chi if math.isfinite(v) and v > 0.0]
        return min(finite) if finite else 1.0

    @property
    def chi_drop(self) -> float:
        """Fractional deepening of the conformal factor at the origin.

        Uses initial minus minimum chi so static throat profiles (chi < 1 at
        t=0 that relax upward) do not earn false curvature-well credit.
        """
        if not self.has_geometric_data:
            return 0.0
        initial = next((v for v in self.chi if math.isfinite(v) and v > 0.0), 1.0)
        return float(max(0.0, min(1.0, initial - self.min_chi_at_origin)))

    @property
    def peak_abs_K(self) -> float:
        """Peak |trace K| at center: the gravitational 'crunch' rate."""
        finite = [abs(v) for v in self.trace_K if math.isfinite(v)]
        return max(finite) if finite else 0.0

    @property
    def peak_abs_weyl4(self) -> float:
        """Peak |GW signal| at center: Re(Psi4) or A_ij strain proxy magnitude."""
        finite = [abs(v) for v in self.weyl4 if math.isfinite(v)]
        return max(finite) if finite else 0.0

    @property
    def weyl4_peak_time(self) -> float | None:
        """Time of peak |Weyl4| -- when the converging wave focuses."""
        peak = -1.0
        t_peak: float | None = None
        for value, time in zip(self.weyl4, self.t):
            if math.isfinite(value) and abs(value) > peak:
                peak = abs(value)
                t_peak = time
        return t_peak

    def peak_abs_weyl4_after(self, t_min: float) -> float:
        """Peak |GW signal| at center considering only samples at t >= t_min."""
        finite = [
            abs(v)
            for v, time in zip(self.weyl4, self.t)
            if time >= t_min and math.isfinite(v)
        ]
        return max(finite) if finite else 0.0

    def weyl4_peak_time_after(self, t_min: float) -> float | None:
        """Time of peak |Weyl4| among samples at t >= t_min."""
        peak = -1.0
        t_peak: float | None = None
        for value, time in zip(self.weyl4, self.t):
            if time >= t_min and math.isfinite(value) and abs(value) > peak:
                peak = abs(value)
                t_peak = time
        return t_peak
