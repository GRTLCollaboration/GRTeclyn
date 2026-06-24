"""Real spherical harmonic expansion for matter field parameterization.

Provides angular modulation of lump amplitudes (and optionally radial
positions) via a truncated real-SH expansion, replacing the fixed
dipole+quadrupole (ell<=2) Legendre modulation in the shell ansatz with
full angular freedom up to configurable ell_max.

Real spherical harmonics (orthonormal, Condon-Shortley phase):

    Y_{l,0}   = N_{l,0} P_l(cos theta)
    Y_{l,m>0} = sqrt(2) N_{l,m} P_l^m(cos theta) cos(m phi)
    Y_{l,m<0} = sqrt(2) N_{l,|m|} P_l^{|m|}(cos theta) sin(|m| phi)

Flat index: idx = l^2 + l + m  (0-based, 0 <= idx < (l_max+1)^2).
"""

from __future__ import annotations

import math
from typing import Sequence


def _factorial(n: int) -> float:
    """Factorial for small n (used in normalization)."""
    result = 1.0
    for k in range(2, n + 1):
        result *= k
    return result


def _normalization(ell: int, m: int) -> float:
    """Normalization factor N_{l,m} for real spherical harmonics."""
    am = abs(m)
    return math.sqrt(
        (2 * ell + 1) / (4.0 * math.pi)
        * _factorial(ell - am)
        / _factorial(ell + am)
    )


def _associated_legendre(ell: int, m: int, x: float) -> float:
    """Associated Legendre polynomial P_l^m(x), m >= 0, |x| <= 1.

    Uses the standard recurrence with Condon-Shortley phase ((-1)^m
    absorbed into the definition).
    """
    if m > ell:
        return 0.0
    # Seed: P_m^m
    pmm = 1.0
    if m > 0:
        sign = -1.0 if m % 2 else 1.0
        somx2 = math.sqrt(max(0.0, 1.0 - x * x))
        fact = 1.0
        for _i in range(1, m + 1):
            pmm *= sign * fact * somx2
            fact += 2.0
            sign = 1.0  # only first iteration carries (-1)^m via seed
        # Correct: P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^{m/2}
        pmm = 1.0
        somx2 = math.sqrt(max(0.0, 1.0 - x * x))
        double_fact = 1.0
        for i in range(1, m + 1):
            double_fact *= 2 * i - 1
        pmm = ((-1) ** m) * double_fact * (somx2 ** m)
    if ell == m:
        return pmm
    # P_{m+1}^m
    pmmp1 = x * (2 * m + 1) * pmm
    if ell == m + 1:
        return pmmp1
    # Recurrence: (l-m) P_l^m = x(2l-1) P_{l-1}^m - (l+m-1) P_{l-2}^m
    pll = 0.0
    for ll in range(m + 2, ell + 1):
        pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    return pll


def sh_flat_index(ell: int, m: int) -> int:
    """Flat index for (ell, m): idx = ell^2 + ell + m."""
    return ell * ell + ell + m


def sh_ell_m(idx: int) -> tuple[int, int]:
    """Recover (ell, m) from flat index."""
    ell = int(math.isqrt(idx))
    m = idx - ell * ell - ell
    return ell, m


def real_sph_harm(ell: int, m: int, theta: float, phi: float) -> float:
    """Evaluate a single real spherical harmonic Y_{l,m}(theta, phi).

    Parameters
    ----------
    ell : non-negative int
    m : int in [-ell, ell]
    theta : polar angle (0 = +z axis)
    phi : azimuthal angle
    """
    cos_theta = math.cos(theta)
    norm = _normalization(ell, m)
    am = abs(m)
    plm = _associated_legendre(ell, am, cos_theta)
    if m == 0:
        return norm * plm
    elif m > 0:
        return math.sqrt(2.0) * norm * plm * math.cos(m * phi)
    else:
        return math.sqrt(2.0) * norm * plm * math.sin(am * phi)


def eval_sh_expansion(
    coeffs: Sequence[float],
    theta: float,
    phi: float,
    ell_max: int,
) -> float:
    """Evaluate sum_{l=0}^{l_max} sum_{m=-l}^{l} c_{idx} Y_{l,m}(theta, phi).

    ``coeffs`` is indexed by the flat index (length ``(ell_max+1)^2``).
    Missing trailing coefficients are treated as zero.
    """
    n_coeffs = (ell_max + 1) ** 2
    total = 0.0
    for idx in range(min(len(coeffs), n_coeffs)):
        c = coeffs[idx]
        if c == 0.0:
            continue
        ell, m = sh_ell_m(idx)
        total += c * real_sph_harm(ell, m, theta, phi)
    return total


def eval_sh_modulation(
    coeffs: Sequence[float],
    theta: float,
    phi: float,
    ell_max: int,
) -> float:
    """Evaluate 1 + sum_{idx>0} c_idx Y_{l,m}(theta, phi).

    The monopole (idx=0) is excluded — it is absorbed into the base
    amplitude.  Returns the multiplicative modulation factor (clamped >= 0).
    """
    total = 1.0
    n_coeffs = (ell_max + 1) ** 2
    for idx in range(1, min(len(coeffs), n_coeffs)):
        c = coeffs[idx]
        if c == 0.0:
            continue
        ell, m = sh_ell_m(idx)
        total += c * real_sph_harm(ell, m, theta, phi)
    return max(0.0, total)


def cartesian_to_spherical(
    x: float, y: float, z: float
) -> tuple[float, float, float]:
    """Convert (x, y, z) to (r, theta, phi)."""
    r = math.sqrt(x * x + y * y + z * z)
    if r < 1e-30:
        return 0.0, 0.0, 0.0
    theta = math.acos(max(-1.0, min(1.0, z / r)))
    phi = math.atan2(y, x)
    return r, theta, phi


def sh_modulated_amplitudes(
    base_amp: float,
    coeffs: Sequence[float],
    positions: Sequence[tuple[float, float, float]],
    ell_max: int,
) -> list[float]:
    """Compute per-position amplitudes from SH-modulated base amplitude.

    For each position (on the unit sphere or any shell), evaluates
    ``base_amp * (1 + sum c_i Y_{l,m})`` where the sum skips the
    monopole.  Result is clamped >= 0.
    """
    amps: list[float] = []
    for pos in positions:
        _, theta, phi = cartesian_to_spherical(*pos)
        mod = eval_sh_modulation(coeffs, theta, phi, ell_max)
        amps.append(base_amp * mod)
    return amps


def sh_modulated_radii(
    base_radius: float,
    coeffs: Sequence[float],
    directions: Sequence[tuple[float, float, float]],
    ell_max: int,
    *,
    max_deviation: float = 0.5,
) -> list[float]:
    """Compute per-direction radii from SH-modulated base radius.

    ``max_deviation`` caps the fractional deviation to avoid
    lumps at negative or extremely large radii.
    """
    radii: list[float] = []
    for d in directions:
        _, theta, phi = cartesian_to_spherical(*d)
        mod = eval_sh_modulation(coeffs, theta, phi, ell_max)
        deviation = max(-max_deviation, min(max_deviation, mod - 1.0))
        radii.append(base_radius * (1.0 + deviation))
    return radii
