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


def spin_weighted_sph_harm_s2_l2_m0(theta: np.ndarray) -> np.ndarray:
    """
    Exact (real) spin-weighted spherical harmonic:  _-2Y_{2,0}
    _-2Y20(theta,phi) = sqrt(15/(32*pi)) * sin^2(theta)
    """
    return np.sqrt(15.0 / (32.0 * np.pi)) * (np.sin(theta) ** 2)
