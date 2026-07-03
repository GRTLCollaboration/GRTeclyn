"""Container for directional Psi4 (gravitational-wave) metrics."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Psi4Metrics:
    """Aggregated directional Psi4 power extracted from small_data/psi4_directional.dat.

    Metrics are computed from the full l=2, s=-2 spin-weighted spherical-harmonic
    decomposition: P_total = sum_m |C_m|^2 and P_z_beam = |C_{-2}|^2 + |C_2|^2.
    """

    peak_total_power: float
    peak_z_beam_power: float
    peak_beam_ratio: float
    mean_total_power: float
    mean_z_beam_power: float
    mean_beam_ratio: float
    final_total_power: float
    final_z_beam_power: float
    final_beam_ratio: float
    n_samples: int

    @property
    def has_data(self) -> bool:
        return self.n_samples > 0

    @property
    def score_signal(self) -> float:
        """Scalar optimization target: total power * (1 + beaming ratio).

        Rewards both strong GW emission and preferential Z-axis beaming.
        """
        if not self.has_data:
            return 0.0
        return float(self.mean_total_power * (1.0 + self.mean_beam_ratio))
