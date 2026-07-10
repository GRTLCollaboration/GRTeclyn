"""Container for directional Psi4 (gravitational-wave) metrics."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Psi4Metrics:
    """Aggregated directional Psi4 power extracted from small_data/psi4_directional.dat.

    Metrics are computed from the full l=2, s=-2 spin-weighted spherical-harmonic
    decomposition: P_total = sum_m |C_m|^2 and P_z_beam = |C_{-2}|^2 + |C_2|^2.

    v5 additions:
      - beaming_gain: max(dP/dΩ) / (P_total / 4π) — directional gain above isotropic.
      - wavezone_ok: 1/r validity check across extraction radii (std/mean of r·Ψ₄ < tol).
      - direction_stability: mean resultant length of beam directions over time [0,1].
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

    # v5 wave-zone / beaming metrics (optional; default to backward-compatible values)
    mean_beaming_gain: float = 1.0
    peak_beaming_gain: float = 1.0
    wavezone_ok: bool = True
    wavezone_one_over_r_std: float = 0.0
    direction_stability: float = 0.0

    @property
    def has_data(self) -> bool:
        return self.n_samples > 0
