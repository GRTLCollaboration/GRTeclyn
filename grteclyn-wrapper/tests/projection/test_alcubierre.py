"""Tests for the Alcubierre warp-drive matter reconstruction pipeline.

Covers:
  - motif_from_alcubierre: flat chi, axial shift, toroidal exotic support
  - Toroidal support regions: negative rho, ring at bubble wall
  - Toroidal exotic ansatz: lumps in yz-plane ring, exotic flag, axial velocity
  - Warp-target mismatch: flat chi target, non-zero A_ij target, K=0
  - Warp-factory cross-check: reconstructed rho matches T=G/8pi, NEC violated
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.fit.motif import (
    fit_matter_from_motif,
    fit_momentum_from_motif,
)
from grteclyn_wrapper.initial_data.motif import (
    GeometryMotif,
    MomentumTarget,
    alcubierre_shape_function,
    motif_from_alcubierre,
)
from grteclyn_wrapper.projection.mismatch import (
    _is_warp_motif,
    _target_2d_slice_warp,
    _warp_params_from_motif,
    warp_factory_cross_check,
)


# ---------------------------------------------------------------------------
# Shape function
# ---------------------------------------------------------------------------

class TestAlcubierreShapeFunction:
    def test_inside_bubble_is_one(self):
        r = np.array([0.0, 0.5, 1.0])
        f = alcubierre_shape_function(r, bubble_radius=2.0, sigma=2.0)
        assert all(f > 0.98), f"Expected f~1 inside bubble, got {f}"

    def test_outside_bubble_is_zero(self):
        r = np.array([5.0, 8.0, 10.0])
        f = alcubierre_shape_function(r, bubble_radius=2.0, sigma=2.0)
        assert all(f < 0.01), f"Expected f~0 outside bubble, got {f}"

    def test_wall_is_half(self):
        r = np.array([2.0])
        f = alcubierre_shape_function(r, bubble_radius=2.0, sigma=2.0)
        assert abs(f[0] - 0.5) < 0.01, f"Expected f~0.5 at wall, got {f[0]}"

    def test_scalar_input(self):
        f = alcubierre_shape_function(0.0, bubble_radius=2.0, sigma=2.0)
        assert abs(float(f) - 1.0) < 0.01


# ---------------------------------------------------------------------------
# Motif generation
# ---------------------------------------------------------------------------

class TestMotifFromAlcubierre:
    def test_flat_chi(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        num_bases = int(m.overrides["recipe_num_bases"])
        for n in range(num_bases):
            assert float(m.overrides[f"recipe_chi_coeff_{n}"]) == 0.0
        assert float(m.overrides["recipe_chi_asymptotic"]) == 1.0

    def test_axial_shift_present(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        num_bases = int(m.overrides["recipe_num_bases"])
        has_shift = any(
            abs(float(m.overrides[f"recipe_beta_coeff_{n}"])) > 1.0e-6
            for n in range(num_bases)
        )
        assert has_shift, "Expected non-zero beta coefficients for axial shift"

    def test_exotic_needed(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        assert m.exotic_needed is True

    def test_toroidal_momentum_template(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        assert m.momentum_target.template == "toroidal"
        assert m.momentum_target.credible is True

    def test_support_regions_are_exotic(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        assert len(m.support_regions) >= 1
        for region in m.support_regions:
            assert region.exotic is True
            assert region.peak_rho < 0.0, "Expected negative rho (exotic)"

    def test_support_at_bubble_wall(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        for region in m.support_regions:
            # The ring should be at cylindrical radius ~ bubble_radius
            assert abs(region.radial_center - 2.0) < 1.5, (
                f"Expected ring at ~R=2.0, got radial_center={region.radial_center}"
            )

    def test_episode_path_encodes_params(self):
        m = motif_from_alcubierre(velocity=0.7, bubble_radius=3.0, sigma=1.5)
        assert "v=0.7" in m.episode_path
        assert "R=3.0" in m.episode_path
        assert "sigma=1.5" in m.episode_path

    def test_min_rho_is_negative(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        assert m.min_rho_required is not None
        assert m.min_rho_required < 0.0


# ---------------------------------------------------------------------------
# Toroidal matter fitting
# ---------------------------------------------------------------------------

class TestToroidalFitting:
    def test_toroidal_ring_in_yz_plane(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        assert len(fitted.lumps) == 5
        # All lumps should be at x~0 (yz-plane ring)
        for lump in fitted.lumps:
            assert abs(lump["center"][0]) < 0.1, (
                f"Expected x~0 for toroidal ring, got {lump['center']}"
            )

    def test_all_lumps_exotic(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        for lump in fitted.lumps:
            assert lump["exotic"] == 1

    def test_maximal_slicing_required(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        assert fitted.maximal_slicing is True

    def test_axial_velocity_after_momentum_fit(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        fitted = fit_momentum_from_motif(m, fitted)
        for lump in fitted.lumps:
            assert lump["velocity"][0] > 0.0, "Expected axial (x) velocity"
            assert abs(lump["velocity"][1]) < 1.0e-10
            assert abs(lump["velocity"][2]) < 1.0e-10

    def test_toroidal_omega_and_mode(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        fitted = fit_momentum_from_motif(m, fitted)
        for lump in fitted.lumps:
            assert lump["omega"] > 0.0, "Expected non-zero toroidal omega"
            assert lump["mode"] >= 1, "Expected mode >= 1 for toroidal"

    def test_ring_radius_matches_bubble(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        fitted = fit_matter_from_motif(m, max_lumps=5)
        # Check that lumps are at cylindrical radius ~ bubble_radius from x-axis
        for lump in fitted.lumps:
            cy, cz = lump["center"][1], lump["center"][2]
            r_cyl = math.sqrt(cy**2 + cz**2)
            assert abs(r_cyl - 2.0) < 1.0, (
                f"Expected ring at r~2.0, got r_cyl={r_cyl}"
            )


# ---------------------------------------------------------------------------
# Warp-target mismatch
# ---------------------------------------------------------------------------

class TestWarpTargetMismatch:
    def test_is_warp_motif(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        assert _is_warp_motif(m) is True

    def test_not_warp_motif_for_spherical(self):
        m = GeometryMotif(
            episode_path="/tmp/fake",
            overrides={"recipe_num_bases": 4, "recipe_basis_width": 1.0,
                       "recipe_basis_radius_max": 8.0, "recipe_chi_asymptotic": 1.0},
            transport_axis=(1.0, 0.0, 0.0),
            polarity=0.5, f_shortcut=0.1, f_op=0.08, f_null=0.05,
            f_portal=0.0, f_throat=0.0, f_asymmetry=0.3,
            beta_max=0.3, beta_mean=0.1, exotic_needed=True,
            min_rho_required=-0.01, integral_negative_rho=-0.005,
            static_lens_only=False,
            support_regions=(),
            momentum_target=MomentumTarget(
                direction=(1.0, 0.0, 0.0), support_center=(0.0, 0.0, 0.0),
                strength=0.5, template="axial_boost", credible=True,
            ),
        )
        assert _is_warp_motif(m) is False

    def test_warp_params_extraction(self):
        m = motif_from_alcubierre(velocity=0.7, bubble_radius=3.0, sigma=1.5)
        params = _warp_params_from_motif(m)
        assert params["v"] == 0.7
        assert params["R"] == 3.0
        assert params["sigma"] == 1.5

    def test_warp_target_chi_is_flat(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        _, _, chi, _, _, _, _ = _target_2d_slice_warp(m, n=32, L=8.0)
        assert np.allclose(chi, 1.0), "Expected flat chi=1 for Alcubierre"

    def test_warp_target_beta_is_axial(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        _, _, _, bx, bz, _, _ = _target_2d_slice_warp(m, n=32, L=8.0)
        # beta_x should be ~-v inside bubble, ~0 outside
        center = 16
        assert abs(bx[center, center] - (-0.5)) < 0.01, (
            f"Expected beta_x~-0.5 at center, got {bx[center, center]}"
        )
        # beta_z should be zero everywhere
        assert np.allclose(bz, 0.0), "Expected beta_z=0 for axial shift"

    def test_warp_target_K_is_zero(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        _, _, _, _, _, K, _ = _target_2d_slice_warp(m, n=32, L=8.0)
        assert np.allclose(K, 0.0), "Expected K=0 for Alcubierre"

    def test_warp_target_A_ij_nonzero(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        _, _, _, _, _, _, A = _target_2d_slice_warp(m, n=32, L=8.0)
        assert A.max() > 0.01, (
            f"Expected non-zero A_ij at bubble wall, got max={A.max()}"
        )
        # A_ij should be near zero far from the bubble
        corner = 0
        assert A[corner, corner] < 0.01


# ---------------------------------------------------------------------------
# Warp-factory cross-check
# ---------------------------------------------------------------------------

class TestWarpFactoryCrossCheck:
    def test_rho_matches_warpfactory(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        cc = warp_factory_cross_check(m, n=32, L=8.0)
        # The analytic formula and T=G/8pi should agree to within FD noise
        assert cc.rho_l2 < 0.001, (
            f"Expected rho_l2 < 0.001, got {cc.rho_l2}"
        )

    def test_rho_min_is_negative(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        cc = warp_factory_cross_check(m, n=32, L=8.0)
        assert cc.rho_min_motif < 0.0, "Expected negative rho (exotic)"
        assert cc.rho_min_warpfactory < 0.0

    def test_nec_violated(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        cc = warp_factory_cross_check(m, n=32, L=8.0)
        assert cc.nec_violation_fraction > 0.0, "Expected NEC violation"
        assert any("NEC" in note for note in cc.notes)

    def test_exotic_energy_positive(self):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0)
        cc = warp_factory_cross_check(m, n=32, L=8.0)
        assert cc.exotic_energy > 0.0, "Expected positive exotic energy budget"

    def test_higher_velocity_more_exotic(self):
        m_low = motif_from_alcubierre(velocity=0.1, bubble_radius=2.0, sigma=2.0)
        m_high = motif_from_alcubierre(velocity=0.9, bubble_radius=2.0, sigma=2.0)
        cc_low = warp_factory_cross_check(m_low, n=32, L=8.0)
        cc_high = warp_factory_cross_check(m_high, n=32, L=8.0)
        assert cc_high.exotic_energy > cc_low.exotic_energy, (
            f"Expected higher velocity to need more exotic energy: "
            f"low={cc_low.exotic_energy}, high={cc_high.exotic_energy}"
        )
