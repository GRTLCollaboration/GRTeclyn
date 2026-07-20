#!/usr/bin/env python3
"""Upper bounds on null-timing shortcuts from static Casimir cavities.

Computes the fractional geodesic timing advantage f_geo^bound that a
parallel-plate Casimir cavity of gap *a* and fiber length *L* would
produce via linearised Einstein gravity, and compares it to the ~20 %
shortcut obtained from exotic Q-ball matter in the companion numerical
relativity paper.

All arithmetic is SI; geometric-unit conversions are printed for
cross-reference with the GRTeclyn pipeline (G = c = 1).

References
----------
[1] H.B.G. Casimir, Proc. Kon. Ned. Akad. Wetensch. B51, 793 (1948).
[2] H. White et al., "Warp field mechanics 101", NASA/TF-2011-01, 2011.
[3] C. Wilson et al., Nature 479, 376 (2011)  (dynamical Casimir).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# ---------------------------------------------------------------------------
# Physical constants (SI, CODATA 2018)
# ---------------------------------------------------------------------------
HBAR = 1.054_571_817e-34       # reduced Planck constant  [J s]
C = 2.997_924_58e8             # speed of light           [m/s]
G = 6.674_30e-11               # gravitational constant   [m^3 kg^-1 s^-2]
PI = math.pi


@dataclass
class CasimirBound:
    """Container for one gap / fiber-length scenario."""
    gap_m: float
    fiber_length_m: float
    plate_area_m2: float
    rho_cas_J_m3: float
    rho_cas_kg_m3: float
    E_neg_J: float
    h_metric: float
    delta_t_s: float
    f_geo_bound: float
    amplification_clock: float
    amplification_interesting: float


def casimir_energy_density(a: float) -> float:
    """Casimir energy density rho [J/m^3] for gap *a* [m]."""
    return -(PI**2 * HBAR * C) / (720.0 * a**4)


def casimir_bound(
    gap_m: float,
    fiber_length_m: float = 1_000.0,
    plate_width_m: float = 5e-6,
    plate_area_m2: float | None = None,
) -> CasimirBound:
    """Compute the linearised-GR upper bound on null-timing shortcut.

    Parameters
    ----------
    gap_m : float
        Plate separation (metres).
    fiber_length_m : float
        Length of the Casimir fibre / waveguide (metres).
    plate_width_m : float
        Transverse width of a photonic waveguide face (metres).
        Default 5 μm (silicon-photonics scale).
    plate_area_m2 : float | None
        Face area of one plate (= width × length). Overrides width×length
        when set explicitly.
    """
    if plate_area_m2 is None:
        plate_area_m2 = plate_width_m * fiber_length_m

    rho = casimir_energy_density(gap_m)
    rho_mass = rho / C**2
    # Negative-energy region: gap × face area
    cavity_volume = plate_area_m2 * gap_m
    E_neg = rho * cavity_volume

    # Linearised metric perturbation h ~ (G / c^4) * |E_neg| / L
    kappa_over_8pi = G / C**4          # ~8.26e-45  s^2 kg^-1 m^-1
    h_metric = kappa_over_8pi * abs(E_neg) / fiber_length_m

    # Timing shortcut: delta_t ~ h * (L / c),  f_geo^bound ~ h
    delta_t = h_metric * (fiber_length_m / C)
    f_geo_bound = h_metric

    # Amplification factors needed to reach measurable thresholds
    CLOCK_THRESHOLD = 1e-19            # best atomic-clock precision [s]
    INTERESTING_F = 1e-6               # "interesting" fractional shortcut
    amp_clock = CLOCK_THRESHOLD / delta_t if delta_t > 0 else math.inf
    amp_interesting = INTERESTING_F / f_geo_bound if f_geo_bound > 0 else math.inf

    return CasimirBound(
        gap_m=gap_m,
        fiber_length_m=fiber_length_m,
        plate_area_m2=plate_area_m2,
        rho_cas_J_m3=rho,
        rho_cas_kg_m3=rho_mass,
        E_neg_J=E_neg,
        h_metric=h_metric,
        delta_t_s=delta_t,
        f_geo_bound=f_geo_bound,
        amplification_clock=amp_clock,
        amplification_interesting=amp_interesting,
    )


def print_scenario(b: CasimirBound) -> None:
    print(f"  gap           = {b.gap_m:.0e} m")
    print(f"  fiber length  = {b.fiber_length_m:.0e} m")
    print(f"  plate area    = {b.plate_area_m2:.2e} m^2")
    print(f"  rho_cas       = {b.rho_cas_J_m3:+.4e} J/m^3")
    print(f"  rho_cas(mass) = {b.rho_cas_kg_m3:+.4e} kg/m^3")
    print(f"  E_neg         = {b.E_neg_J:+.4e} J")
    print(f"  h (metric)    = {b.h_metric:.4e}")
    print(f"  delta_t       = {b.delta_t_s:.4e} s")
    print(f"  f_geo^bound   = {b.f_geo_bound:.4e}")
    print(f"  amplification to clock   (1e-19 s) : {b.amplification_clock:.2e}")
    print(f"  amplification to f=1e-6            : {b.amplification_interesting:.2e}")
    print()


def latex_table_row(b: CasimirBound) -> str:
    """One row for the paper's LaTeX table."""
    return (
        f"  {b.gap_m * 1e9:5.0f}"
        f" & ${b.rho_cas_J_m3:+.2e}$"
        f" & ${b.delta_t_s:.1e}$"
        f" & ${b.f_geo_bound:.1e}$"
        f" & ${b.amplification_interesting:.0e}$"
        " \\\\"
    )


# ---------------------------------------------------------------------------
# Main: run canonical scenarios and print results + LaTeX snippet
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Silicon-photonics waveguide: 5 μm width × 1 km length, varying gap.
    scenarios = [
        dict(gap_m=100e-9, fiber_length_m=1_000, plate_width_m=5e-6),
        dict(gap_m=10e-9,  fiber_length_m=1_000, plate_width_m=5e-6),
        dict(gap_m=1e-9,   fiber_length_m=1_000, plate_width_m=5e-6),
    ]

    results: list[CasimirBound] = []
    for kw in scenarios:
        b = casimir_bound(**kw)
        results.append(b)

    print("=" * 64)
    print("Casimir null-timing upper bounds (static parallel plates, SI)")
    print("Waveguide: width = 5 μm, length = 1 km")
    print("=" * 64)
    for b in results:
        print_scenario(b)

    print("-" * 64)
    print("Reference: Q-ball paper headline  f_geo = 0.202  (20.2 %)")
    print(f"Einstein coupling  G/c^4 = {G / C**4:.4e}  s^2 kg^-1 m^-1")
    print("-" * 64)

    print("\nLaTeX table rows (gap [nm] | rho | delta_t | f_geo^bound | amp):")
    print(r"  \begin{tabular}{r l l l l}")
    print(r"  $a$ [nm] & $\rho_{\rm cas}$ [J/m$^3$]"
          r" & $\Delta t$ [s] & $f_{\rm geo}^{\rm bound}$"
          r" & amplification \\")
    print(r"  \hline")
    for b in results:
        print(latex_table_row(b))
    print(r"  \end{tabular}")
