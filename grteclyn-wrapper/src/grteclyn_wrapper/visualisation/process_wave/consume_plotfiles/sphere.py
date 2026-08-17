from __future__ import annotations

from math import factorial as _fact, sqrt as _sqrt
from typing import Tuple

import numpy as np


def get_extraction_points(radius: float, n_points: int) -> Tuple[np.ndarray, ...]:
    """Generate sphere points and weights (theta, phi grid)."""
    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")

    X = radius * np.sin(THETA) * np.cos(PHI)
    Y = radius * np.sin(THETA) * np.sin(PHI)
    Z = radius * np.cos(THETA)

    dtheta = np.pi / (n_points - 1)
    dphi = 2.0 * np.pi / n_points
    W = np.sin(THETA) * dtheta * dphi
    return X, Y, Z, THETA, PHI, W


# ---------------------------------------------------------------------------
# l = 2, s = -2 spin-weighted spherical harmonics
# ---------------------------------------------------------------------------

def _spin_weighted_s2_l2_base(theta: np.ndarray, phi: np.ndarray | None) -> dict[int, np.ndarray]:
    """Return the five s=-2, l=2 spin-weighted spherical harmonics.

    Modes are indexed by m = -2, -1, 0, 1, 2.  Values are complex arrays with the
    e^{i m phi} phase factor included.  Normalization follows the standard
    Condon-Shortley phase convention used in NR (matches the existing real m=0
    formula in ``spin_weighted_sph_harm_s2_l2_m0``).
    """
    c = np.cos(theta)
    s = np.sin(theta)
    out: dict[int, np.ndarray] = {}

    # m = 0 (real, matches existing implementation)
    out[0] = np.sqrt(15.0 / (32.0 * np.pi)) * (s ** 2)

    if phi is not None:
        # m = +1
        out[1] = np.sqrt(5.0 / (16.0 * np.pi)) * s * (1.0 + c) * np.exp(1j * phi)
        # m = +2
        out[2] = np.sqrt(5.0 / (64.0 * np.pi)) * ((1.0 + c) ** 2) * np.exp(2j * phi)
        # m = -1
        out[-1] = np.sqrt(5.0 / (16.0 * np.pi)) * s * (1.0 - c) * np.exp(-1j * phi)
        # m = -2
        out[-2] = np.sqrt(5.0 / (64.0 * np.pi)) * ((1.0 - c) ** 2) * np.exp(-2j * phi)
    return out


def spin_weighted_sph_harm_s2_l2_m0(theta: np.ndarray) -> np.ndarray:
    """
    Exact (real) spin-weighted spherical harmonic:  _-2Y_{2,0}
    _-2Y20(theta,phi) = sqrt(15/(32*pi)) * sin^2(theta)
    """
    return np.real(_spin_weighted_s2_l2_base(theta, None)[0])


def spin_weighted_sph_harm_s2_l2(
    theta: np.ndarray,
    phi: np.ndarray,
) -> dict[int, np.ndarray]:
    """All five s=-2, l=2 spin-weighted spherical harmonics keyed by m."""
    return _spin_weighted_s2_l2_base(theta, phi)


# ---------------------------------------------------------------------------
# l = 3, s = -2 spin-weighted spherical harmonics
# ---------------------------------------------------------------------------

def wigner_small_d(l: int, mp: int, m: int, theta: np.ndarray) -> np.ndarray:
    """Wigner small-d matrix element d^l_{mp,m}(theta)."""
    ch = np.cos(theta / 2.0)
    sh = np.sin(theta / 2.0)
    total = np.zeros(np.shape(theta), dtype=float)
    for k in range(max(0, m - mp), min(l + m, l - mp) + 1):
        num = _sqrt(
            _fact(l + m) * _fact(l - m) * _fact(l + mp) * _fact(l - mp)
        )
        den = _fact(l + m - k) * _fact(k) * _fact(mp - m + k) * _fact(l - mp - k)
        total += (
            ((-1.0) ** (mp - m + k))
            * (num / den)
            * ch ** (2 * l + m - mp - 2 * k)
            * sh ** (mp - m + 2 * k)
        )
    return total


def spin_weighted_sph_harm(
    s: int, l: int, m: int, theta: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    """General spin-weighted spherical harmonic sYlm via the Wigner-d formula.

        sYlm(theta,phi) = (-1)^s sqrt((2l+1)/4pi) d^l_{-s,m}(theta) e^{i m phi}

    Verified orthonormal to quadrature accuracy for l = 2,3,4 at s = -2.

    Convention note: this reproduces ``_spin_weighted_s2_l2_base`` exactly for
    m = 0, +-2, and differs by an overall sign for m = +-1.  The l = 2 helpers
    below deliberately keep their original signs because ``psi4_mode_l2_all.dat``
    is a published column contract; an overall per-mode sign does not affect the
    |amplitude| and power fractions quoted from that stream.  Anything that
    compares *across* l (odd-even interference) must use this function for every
    l so the relative phases are consistent.
    """
    return (
        ((-1.0) ** s)
        * _sqrt((2 * l + 1) / (4.0 * np.pi))
        * wigner_small_d(l, -s, m, theta)
        * np.exp(1j * m * phi)
    )


def _spin_weighted_s2_l3_base(theta: np.ndarray, phi: np.ndarray | None) -> dict[int, np.ndarray]:
    """Return the seven s=-2, l=3 spin-weighted spherical harmonics (m=-3..+3).

    Computed from :func:`spin_weighted_sph_harm`.  The previous hand-transcribed
    closed forms here were **not normalised** -- their <|Y|^2> ranged over
    0.27..0.52 instead of 1 -- and were never exercised because only l = 2 was
    wired into the extraction path.
    """
    if phi is None:
        phi = np.zeros(np.shape(theta))
    return {m: spin_weighted_sph_harm(-2, 3, m, theta, phi) for m in range(-3, 4)}


def spin_weighted_sph_harm_s2_l3(
    theta: np.ndarray,
    phi: np.ndarray,
) -> dict[int, np.ndarray]:
    """All seven s=-2, l=3 spin-weighted spherical harmonics keyed by m."""
    return _spin_weighted_s2_l3_base(theta, phi)


# ---------------------------------------------------------------------------
# l = 4, s = -2 spin-weighted spherical harmonics
# ---------------------------------------------------------------------------

def _spin_weighted_s2_l4_base(theta: np.ndarray, phi: np.ndarray | None) -> dict[int, np.ndarray]:
    """Return the nine s=-2, l=4 spin-weighted spherical harmonics (m=-4..+4).

    Computed from :func:`spin_weighted_sph_harm`.  As with l = 3, the previous
    hand-transcribed forms were unnormalised (<|Y|^2> as low as 0.055) and unused.
    """
    if phi is None:
        phi = np.zeros(np.shape(theta))
    return {m: spin_weighted_sph_harm(-2, 4, m, theta, phi) for m in range(-4, 5)}


def spin_weighted_sph_harm_s2_l4(
    theta: np.ndarray,
    phi: np.ndarray,
) -> dict[int, np.ndarray]:
    """All nine s=-2, l=4 spin-weighted spherical harmonics keyed by m."""
    return _spin_weighted_s2_l4_base(theta, phi)


# ---------------------------------------------------------------------------
# Generic multipole access
# ---------------------------------------------------------------------------

def spin_weighted_sph_harm_s2(
    l: int,
    theta: np.ndarray,
    phi: np.ndarray,
) -> dict[int, np.ndarray]:
    """Return s=-2 spin-weighted spherical harmonics for given l, keyed by m."""
    if l == 2:
        return spin_weighted_sph_harm_s2_l2(theta, phi)
    elif l == 3:
        return spin_weighted_sph_harm_s2_l3(theta, phi)
    elif l == 4:
        return spin_weighted_sph_harm_s2_l4(theta, phi)
    else:
        raise ValueError(f"Unsupported l={l}; only l=2,3,4 implemented.")
