"""Unit tests for the spin-weighted spherical harmonics used in Psi4 extraction."""

from __future__ import annotations

import numpy as np
import pytest

from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.sphere import (
    spin_weighted_sph_harm_s2_l2,
    spin_weighted_sph_harm_s2_l2_m0,
)


def _make_grid(n: int = 128):
    theta = np.linspace(0.0, np.pi, n)
    phi = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    dtheta = np.pi / (n - 1)
    dphi = 2.0 * np.pi / n
    W = np.sin(THETA) * dtheta * dphi
    return THETA, PHI, W


def test_m0_matches_real_l2m0():
    """The new m=0 mode must equal the existing real implementation."""
    THETA, PHI, _ = _make_grid(64)
    ylm = spin_weighted_sph_harm_s2_l2(THETA, PHI)
    np.testing.assert_allclose(np.real(ylm[0]), spin_weighted_sph_harm_s2_l2_m0(THETA))
    np.testing.assert_allclose(np.imag(ylm[0]), 0.0, atol=1e-15)


def test_m0_peaks_at_equator():
    THETA, PHI, _ = _make_grid(64)
    ylm = spin_weighted_sph_harm_s2_l2(THETA, PHI)
    vals = np.abs(ylm[0])
    assert vals.max() > 0.0
    # m=0 is ~sin^2(theta), so maximum is at equator.
    equator_idx = np.unravel_index(vals.argmax(), vals.shape)
    assert np.isclose(THETA[equator_idx], np.pi / 2, atol=0.1)


def test_m2_peaks_at_poles():
    THETA, PHI, _ = _make_grid(64)
    ylm = spin_weighted_sph_harm_s2_l2(THETA, PHI)
    vals = np.abs(ylm[2])
    assert vals.max() > 0.0
    # m=+2 peaks at north pole, m=-2 at south pole.
    north_idx = np.unravel_index(vals.argmax(), vals.shape)
    assert np.isclose(THETA[north_idx], 0.0, atol=0.1)

    vals_minus = np.abs(ylm[-2])
    south_idx = np.unravel_index(vals_minus.argmax(), vals.shape)
    assert np.isclose(THETA[south_idx], np.pi, atol=0.1)


def test_orthonormality():
    """The five s=-2, l=2 modes must be orthonormal on the sphere."""
    THETA, PHI, W = _make_grid(128)
    ylm = spin_weighted_sph_harm_s2_l2(THETA, PHI)
    for m in (-2, -1, 0, 1, 2):
        norm = np.sum(np.abs(ylm[m]) ** 2 * W)
        np.testing.assert_allclose(norm, 1.0, rtol=1e-2, err_msg=f"mode {m} norm")
    for m1 in (-2, -1, 0, 1, 2):
        for m2 in (-2, -1, 0, 1, 2):
            if m1 == m2:
                continue
            overlap = np.sum(ylm[m1] * np.conj(ylm[m2]) * W)
            np.testing.assert_allclose(
                overlap, 0.0, atol=1e-2, err_msg=f"overlap {m1},{m2}"
            )
