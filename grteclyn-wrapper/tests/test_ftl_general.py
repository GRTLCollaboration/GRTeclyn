"""Tests for the general, mechanism-agnostic operational FTL diagnostics."""

from __future__ import annotations

import numpy as np

from grteclyn_wrapper.metrics.probes.ftl.general import (
    build_tslice_fields_xz,
    compute_general_ftl,
    coordinate_light_speed,
    operational_ftl_on_grid,
)


def test_flat_space_light_speed_is_unity():
    gamma = np.eye(2)
    beta = np.zeros(2)
    for direction in [(1.0, 0.0), (0.0, 1.0), (1.0, 1.0), (-1.0, 0.3)]:
        s = coordinate_light_speed(1.0, beta, gamma, direction)
        assert abs(s - 1.0) < 1e-12


def test_shift_creates_superluminal_coordinate_channel():
    # Standard ADM sign: the forward coordinate light speed along n is
    # s = -beta.n + alpha (flat gamma), so a +x shortcut needs beta_x < 0
    # (Alcubierre's convention for a bubble travelling in +x).
    gamma = np.eye(2)
    beta = np.array([-0.5, 0.0])  # drags light in +x
    forward = coordinate_light_speed(1.0, beta, gamma, (1.0, 0.0))
    backward = coordinate_light_speed(1.0, beta, gamma, (-1.0, 0.0))
    assert forward > 1.0
    assert backward < 1.0
    assert abs(forward - 1.5) < 1e-12
    assert abs(backward - 0.5) < 1e-12


def test_proper_shortcut_from_conformal_compression():
    # gamma = delta / chi with chi > 1 shrinks coordinate proper distance and
    # raises the coordinate light speed above 1 -> a portal-type shortcut with
    # no shift at all (mechanism-agnostic detection).
    gamma = np.eye(2) / 4.0  # chi = 4
    beta = np.zeros(2)
    s = coordinate_light_speed(1.0, beta, gamma, (1.0, 0.0))
    assert abs(s - 2.0) < 1e-12  # sqrt(chi)


def test_flat_grid_has_no_shortcut():
    n = 41
    alpha = np.ones((n, n))
    beta = np.zeros((n, n, 2))
    gamma = np.zeros((n, n, 2, 2))
    gamma[:, :, 0, 0] = 1.0
    gamma[:, :, 1, 1] = 1.0
    report = operational_ftl_on_grid(
        alpha, beta, gamma, spacing=(1.0, 1.0), source=(1, n // 2), target=(n - 2, n // 2)
    )
    assert report.reachable
    assert report.f_op < 1e-9
    assert abs(report.max_local_speed - 1.0) < 1e-9


def test_grid_shortcut_detected_with_central_speedup():
    # A central band with chi > 1 (faster light) along the travel axis should
    # yield t_min < t_flat -> f_op > 0.
    n = 61
    alpha = np.ones((n, n))
    beta = np.zeros((n, n, 2))
    gamma = np.zeros((n, n, 2, 2))
    gamma[:, :, 0, 0] = 1.0
    gamma[:, :, 1, 1] = 1.0
    # chi = 4 in a horizontal channel around z = mid -> gamma = 1/4 there.
    mid = n // 2
    gamma[:, mid - 3 : mid + 4, 0, 0] = 0.25
    gamma[:, mid - 3 : mid + 4, 1, 1] = 0.25
    report = operational_ftl_on_grid(
        alpha, beta, gamma, spacing=(1.0, 1.0), source=(1, mid), target=(n - 2, mid)
    )
    assert report.reachable
    assert report.f_op > 0.0
    assert report.max_local_speed > 1.5


def test_compute_general_ftl_flat_recipe_is_zero():
    overrides = {
        "recipe_num_bases": 4,
        "recipe_basis_width": 1.0,
        "recipe_basis_radius_max": 8.0,
        "recipe_chi_asymptotic": 1.0,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_beta_asymptotic": 0.0,
    }
    report = compute_general_ftl(overrides, L=8.0, n=61)
    assert report.reachable
    assert report.f_op < 1e-6


def test_compute_general_ftl_warp_recipe_positive():
    # Negative shift bump -> superluminal +x channel -> f_op > 0.
    overrides = {
        "recipe_num_bases": 4,
        "recipe_basis_width": 2.0,
        "recipe_basis_radius_max": 8.0,
        "recipe_chi_asymptotic": 1.0,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_beta_asymptotic": 0.0,
        "recipe_beta_coeff_0": -0.6,
        "recipe_beta_coeff_1": -0.6,
        "recipe_beta_coeff_2": -0.4,
        "recipe_beta_coeff_3": -0.2,
    }
    report = compute_general_ftl(overrides, L=8.0, n=81)
    assert report.reachable
    assert report.max_local_speed > 1.0
    assert report.f_op > 0.0
