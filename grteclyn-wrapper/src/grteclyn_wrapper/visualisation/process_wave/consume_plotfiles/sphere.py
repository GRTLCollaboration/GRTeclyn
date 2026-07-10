from __future__ import annotations

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

def _spin_weighted_s2_l3_base(theta: np.ndarray, phi: np.ndarray | None) -> dict[int, np.ndarray]:
    """Return the seven s=-2, l=3 spin-weighted spherical harmonics (m=-3..+3).

    Standard Condon-Shortley convention; formulas from Goldberg et al. (1967).
    """
    c = np.cos(theta)
    s = np.sin(theta)
    half = theta / 2.0
    ch = np.cos(half)
    sh = np.sin(half)
    out: dict[int, np.ndarray] = {}

    # m = 0
    out[0] = np.sqrt(21.0 / (20.0 * np.pi)) * s * s * np.cos(theta) * (1.0 + 0j)

    if phi is not None:
        # m = +1
        out[1] = np.sqrt(21.0 / (80.0 * np.pi)) * s * (1.0 + c) * (2.0 * c - 1.0) * np.exp(1j * phi)
        # m = +2
        out[2] = np.sqrt(7.0 / (40.0 * np.pi)) * s * ((1.0 + c) ** 2) * 0.5 * np.exp(2j * phi)
        # m = +3
        out[3] = -np.sqrt(7.0 / (80.0 * np.pi)) * ((1.0 + c) ** 3) / (2.0 * np.sqrt(2.0)) * np.exp(3j * phi)
        # m = -1
        out[-1] = np.sqrt(21.0 / (80.0 * np.pi)) * s * (1.0 - c) * (-2.0 * c - 1.0) * np.exp(-1j * phi)
        # m = -2
        out[-2] = np.sqrt(7.0 / (40.0 * np.pi)) * s * ((1.0 - c) ** 2) * 0.5 * np.exp(-2j * phi)
        # m = -3
        out[-3] = np.sqrt(7.0 / (80.0 * np.pi)) * ((1.0 - c) ** 3) / (2.0 * np.sqrt(2.0)) * np.exp(-3j * phi)
    return out


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

    Standard Condon-Shortley convention; formulas from Goldberg et al. (1967).
    """
    c = np.cos(theta)
    s = np.sin(theta)
    s2 = s * s
    out: dict[int, np.ndarray] = {}

    # m = 0
    out[0] = (3.0 / 4.0) * np.sqrt(5.0 / (2.0 * np.pi)) * s2 * (7.0 * c * c - 1.0) / 7.0 * (1.0 + 0j)

    if phi is not None:
        # m = +1
        out[1] = (3.0 / (2.0 * np.sqrt(2.0))) * np.sqrt(5.0 / np.pi) * s * s2 * c * (1.0 + c) / 3.0 * np.exp(1j * phi)
        # m = +2
        out[2] = (3.0 / 4.0) * np.sqrt(5.0 / (2.0 * np.pi)) * s2 * ((1.0 + c) ** 2) * (2.0 * c - 1.0) / 6.0 * np.exp(2j * phi)
        # m = +3
        out[3] = -np.sqrt(5.0 / np.pi) * s * ((1.0 + c) ** 3) / (4.0 * np.sqrt(2.0) * 3.0) * np.exp(3j * phi)
        # m = +4
        out[4] = (3.0 / 16.0) * np.sqrt(5.0 / np.pi) * ((1.0 + c) ** 4) / 6.0 * np.exp(4j * phi)
        # m = -1
        out[-1] = (3.0 / (2.0 * np.sqrt(2.0))) * np.sqrt(5.0 / np.pi) * s * s2 * (-c) * (1.0 - c) / 3.0 * np.exp(-1j * phi)
        # m = -2
        out[-2] = (3.0 / 4.0) * np.sqrt(5.0 / (2.0 * np.pi)) * s2 * ((1.0 - c) ** 2) * (-2.0 * c - 1.0) / 6.0 * np.exp(-2j * phi)
        # m = -3
        out[-3] = np.sqrt(5.0 / np.pi) * s * ((1.0 - c) ** 3) / (4.0 * np.sqrt(2.0) * 3.0) * np.exp(-3j * phi)
        # m = -4
        out[-4] = (3.0 / 16.0) * np.sqrt(5.0 / np.pi) * ((1.0 - c) ** 4) / 6.0 * np.exp(-4j * phi)
    return out


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
