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
    _compute_si_mismatch,
    warp_factory_cross_check,
)
from grteclyn_wrapper.projection.warp_gridinit import (
    alcubierre_analytic_fields,
    alcubierre_analytic_Si,
    paint_alcubierre_warp_on_gridinit,
    write_alcubierre_gridinit,
)
from grteclyn_wrapper.grtresna.io.gridinit import read_gridinit, write_gridinit


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


# ---------------------------------------------------------------------------
# Analytic gridinit writer (Phase A0)
# ---------------------------------------------------------------------------

class TestAnalyticGridinit:
    def test_round_trip_through_read_gridinit(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        assert g.data.shape == (16, 16, 16, 37)
        assert "shift1" in g.comp_names

    def test_chi_is_one_everywhere(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        chi = g.data[:, :, :, g.comp_names.index("chi")]
        assert np.allclose(chi, 1.0)

    def test_lapse_is_one_everywhere(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        lapse = g.data[:, :, :, g.comp_names.index("lapse")]
        assert np.allclose(lapse, 1.0)

    def test_K_is_zero_everywhere(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        K = g.data[:, :, :, g.comp_names.index("K")]
        assert np.allclose(K, 0.0)

    def test_shift_at_bubble_centre(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        shift1 = g.data[:, :, :, g.comp_names.index("shift1")]
        # Centre of 16^3 grid, bubble at (8, 8, 8)
        assert abs(shift1[8, 8, 8] - (-0.5)) < 0.02, (
            f"Expected shift1~ -0.5 at centre, got {shift1[8, 8, 8]}"
        )

    def test_shift_far_field_is_zero(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        shift1 = g.data[:, :, :, g.comp_names.index("shift1")]
        assert abs(shift1[0, 0, 0]) < 0.01, (
            f"Expected shift1~0 at corner, got {shift1[0, 0, 0]}"
        )

    def test_A_ij_traceless(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        trace = (
            g.data[:, :, :, idx["A11"]]
            + g.data[:, :, :, idx["A22"]]
            + g.data[:, :, :, idx["A33"]]
        )
        assert np.allclose(trace, 0.0, atol=1e-10), (
            f"Expected traceless A_ij, got max|trace|={np.abs(trace).max()}"
        )

    def test_A_ij_peaks_at_wall(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=32, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        A_mag = np.sqrt(
            g.data[:, :, :, idx["A11"]]**2
            + 2 * g.data[:, :, :, idx["A12"]]**2
            + g.data[:, :, :, idx["A33"]]**2
        )
        # A_ij should peak at the bubble wall (r_s ~ R=2.0), not at the centre
        centre = A_mag[16, 16, 16]
        assert A_mag.max() > centre, (
            f"Expected A_ij to peak at wall, not centre (max={A_mag.max()}, centre={centre})"
        )


# ---------------------------------------------------------------------------
# Warp-shift painting (Phase A1)
# ---------------------------------------------------------------------------

class TestWarpPainting:
    def test_painter_overwrites_shift(self, tmp_path):
        p = tmp_path / "test.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        # Zero out the shift to simulate a solved gridinit
        g = read_gridinit(p)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        g.data[:, :, :, idx["shift1"]] = 0.0
        write_gridinit(g.data, g.comp_names, g.dx_xyz, g.origin, p)
        # Paint
        paint_alcubierre_warp_on_gridinit(
            p, velocity=0.5, bubble_radius=2.0, sigma=2.0, paint_aij=True,
        )
        g2 = read_gridinit(p)
        shift1 = g2.data[:, :, :, idx["shift1"]]
        assert abs(shift1[8, 8, 8] - (-0.5)) < 0.02

    def test_painter_preserves_header(self, tmp_path):
        p = tmp_path / "test.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g_before = read_gridinit(p)
        paint_alcubierre_warp_on_gridinit(
            p, velocity=0.5, bubble_radius=2.0, sigma=2.0, paint_aij=True,
        )
        g_after = read_gridinit(p)
        assert g_after.comp_names == g_before.comp_names
        assert np.allclose(g_after.dx_xyz, g_before.dx_xyz)
        assert np.allclose(g_after.origin, g_before.origin)
        assert g_after.data.shape == g_before.data.shape

    def test_paint_aij_false_leaves_A_ij_untouched(self, tmp_path):
        p = tmp_path / "test.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        # Set A_ij to a sentinel value
        g.data[:, :, :, idx["A11"]] = 42.0
        write_gridinit(g.data, g.comp_names, g.dx_xyz, g.origin, p)
        paint_alcubierre_warp_on_gridinit(
            p, velocity=0.5, bubble_radius=2.0, sigma=2.0, paint_aij=False,
        )
        g2 = read_gridinit(p)
        A11 = g2.data[:, :, :, idx["A11"]]
        assert np.allclose(A11, 42.0), "A_ij should be unchanged with paint_aij=False"

    def test_paint_aij_true_overwrites_A_ij(self, tmp_path):
        p = tmp_path / "test.gridinit"
        write_alcubierre_gridinit(p, n=16, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        g.data[:, :, :, idx["A11"]] = 42.0
        write_gridinit(g.data, g.comp_names, g.dx_xyz, g.origin, p)
        paint_alcubierre_warp_on_gridinit(
            p, velocity=0.5, bubble_radius=2.0, sigma=2.0, paint_aij=True,
        )
        g2 = read_gridinit(p)
        A11 = g2.data[:, :, :, idx["A11"]]
        assert not np.allclose(A11, 42.0), "A_ij should be overwritten with paint_aij=True"


# ---------------------------------------------------------------------------
# Analytic S_i momentum density (Phase A2)
# ---------------------------------------------------------------------------

class TestAlcubierreSi:
    def test_Si_returns_correct_shapes(self):
        x, z, Sx, Sy, Sz = alcubierre_analytic_Si(
            velocity=0.5, bubble_radius=2.0, sigma=2.0, L=8.0, n=32,
        )
        assert len(x) == 32 and len(z) == 32
        assert Sx.shape == (32, 32)
        assert Sy.shape == (32, 32)
        assert Sz.shape == (32, 32)

    def test_Si_dominant_x_component(self):
        _, _, Sx, Sy, Sz = alcubierre_analytic_Si(
            velocity=0.5, bubble_radius=2.0, sigma=2.0, L=8.0, n=32,
        )
        # S_x should be the dominant component (axial transport)
        assert np.abs(Sx).max() > np.abs(Sz).max()
        # S_y should be exactly zero (axial symmetry on xz-slice)
        assert np.allclose(Sy, 0.0)

    def test_Si_peaks_at_wall(self):
        _, _, Sx, _, _ = alcubierre_analytic_Si(
            velocity=0.5, bubble_radius=2.0, sigma=2.0, L=8.0, n=64,
        )
        # S_x should peak near the bubble wall, not at the centre
        centre = Sx[32, 32]
        assert np.abs(Sx).max() > np.abs(centre)

    def test_Si_higher_velocity_more_momentum(self):
        _, _, Sx_low, _, _ = alcubierre_analytic_Si(
            velocity=0.1, bubble_radius=2.0, sigma=2.0, L=8.0, n=32,
        )
        _, _, Sx_high, _, _ = alcubierre_analytic_Si(
            velocity=0.9, bubble_radius=2.0, sigma=2.0, L=8.0, n=32,
        )
        assert np.abs(Sx_high).max() > np.abs(Sx_low).max()


# ---------------------------------------------------------------------------
# S_i mismatch (Phase A2)
# ---------------------------------------------------------------------------

class TestSiMismatch:
    def test_si_mismatch_self_is_small(self, tmp_path):
        p = tmp_path / "analytic.gridinit"
        write_alcubierre_gridinit(p, n=32, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0, ftl_L=8.0)
        si = _compute_si_mismatch(m, p, n=32, L=8.0)
        assert si < 0.05, f"Expected small si_l2 for self-comparison, got {si}"

    def test_si_mismatch_increases_with_perturbation(self, tmp_path):
        m = motif_from_alcubierre(velocity=0.5, bubble_radius=2.0, sigma=2.0, ftl_L=8.0)
        # Correct A_ij
        p_good = tmp_path / "good.gridinit"
        write_alcubierre_gridinit(p_good, n=32, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        si_good = _compute_si_mismatch(m, p_good, n=32, L=8.0)
        # Doubled A_ij
        p_bad = tmp_path / "bad.gridinit"
        write_alcubierre_gridinit(p_bad, n=32, L=8.0, velocity=0.5, bubble_radius=2.0, sigma=2.0)
        g = read_gridinit(p_bad)
        idx = {n: i for i, n in enumerate(g.comp_names)}
        for name in ("A11", "A12", "A13", "A22", "A23", "A33"):
            g.data[:, :, :, idx[name]] *= 3.0
        write_gridinit(g.data, g.comp_names, g.dx_xyz, g.origin, p_bad)
        si_bad = _compute_si_mismatch(m, p_bad, n=32, L=8.0)
        assert si_bad > si_good, (
            f"Expected si_l2 to increase with perturbed A_ij: "
            f"good={si_good}, bad={si_bad}"
        )
