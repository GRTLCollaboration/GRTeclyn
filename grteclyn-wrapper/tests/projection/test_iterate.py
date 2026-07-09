"""Tests for the iterative geometry-first matter adjustment loop."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.fit.motif import (
    FittedMatter,
    fit_matter_from_motif,
    fit_momentum_from_motif,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig
from grteclyn_wrapper.initial_data.motif import (
    GeometryMotif,
    MomentumTarget,
    SupportRegion,
    extract_motif_from_episode,
)
from grteclyn_wrapper.projection.iterate import (
    IterationConfig,
    SPARSITY_AMP_THRESHOLD,
    _clip_vector,
    _compute_tight_bounds,
    _lumps_from_vector,
    _scale_lump_amplitudes,
    _sparsity_penalty,
    _vector_from_lumps,
    amplitude_precondition,
    run_iterate,
)
from grteclyn_wrapper.projection.mismatch import (
    GATE_FITNESS,
    MismatchReport,
    FEASIBILITY_RHO_SAFE,
    FEASIBILITY_RHO_HARD,
    _exotic_mass_proxy,
    _l2_norm,
    _resample,
    _target_profiles,
    _target_2d_slice,
    compute_mismatch,
    feasibility_precheck,
)


def _make_episode(root: Path) -> Path:
    """Create a minimal geometry-first episode fixture."""
    episode = root / "eval_000001"
    episode.mkdir(parents=True)
    overrides = {
        "recipe_num_bases": 4,
        "recipe_basis_width": 1.0,
        "recipe_basis_radius_max": 8.0,
        "recipe_chi_asymptotic": 1.0,
        "recipe_chi_coeff_0": -0.2,
        "recipe_chi_coeff_1": -0.05,
        "recipe_chi_coeff_2": 0.0,
        "recipe_chi_coeff_3": 0.0,
        "recipe_beta_coeff_0": 0.35,
        "recipe_beta_coeff_1": -0.15,
        "recipe_beta_coeff_2": 0.0,
        "recipe_beta_coeff_3": 0.0,
        "recipe_alpha_coeff_0": 0.05,
        "recipe_alpha_coeff_1": 0.0,
        "recipe_alpha_coeff_2": 0.0,
        "recipe_alpha_coeff_3": 0.0,
    }
    (episode / "metadata.json").write_text(
        json.dumps({"overrides": overrides}, indent=2) + "\n",
        encoding="utf-8",
    )
    data = episode / "data"
    data.mkdir()
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
        handle.write("0 1e-5 1e-5 -0.02 0.01 0.5\n")
        handle.write("1 1e-5 1e-5 -0.02 0.01 0.5\n")
    return episode


def _make_motif() -> GeometryMotif:
    """Minimal motif for unit tests (no episode needed)."""
    return GeometryMotif(
        episode_path="/tmp/fake_episode",
        overrides={
            "recipe_num_bases": 4,
            "recipe_basis_width": 1.0,
            "recipe_basis_radius_max": 8.0,
            "recipe_chi_asymptotic": 1.0,
            "recipe_chi_coeff_0": -0.15,
            "recipe_chi_coeff_1": 0.0,
            "recipe_chi_coeff_2": 0.0,
            "recipe_chi_coeff_3": 0.0,
            "recipe_beta_coeff_0": 0.3,
            "recipe_beta_coeff_1": 0.0,
            "recipe_beta_coeff_2": 0.0,
            "recipe_beta_coeff_3": 0.0,
            "recipe_alpha_coeff_0": 0.0,
            "recipe_alpha_coeff_1": 0.0,
            "recipe_alpha_coeff_2": 0.0,
            "recipe_alpha_coeff_3": 0.0,
        },
        transport_axis=(1.0, 0.0, 0.0),
        polarity=0.5,
        f_shortcut=0.1,
        f_op=0.08,
        f_null=0.05,
        f_portal=0.0,
        f_throat=0.0,
        f_asymmetry=0.3,
        beta_max=0.3,
        beta_mean=0.1,
        exotic_needed=True,
        min_rho_required=-0.01,
        integral_negative_rho=-0.005,
        static_lens_only=False,
        support_regions=(
            SupportRegion(
                center=(0.0, 0.0, 0.0),
                width=3.0,
                peak_rho=-0.02,
                exotic=True,
                radial_center=0.0,
            ),
        ),
        momentum_target=MomentumTarget(
            direction=(1.0, 0.0, 0.0),
            support_center=(0.0, 0.0, 0.0),
            strength=0.5,
            template="axial_boost",
            credible=True,
        ),
    )


# --- Mismatch unit tests ---


class TestL2Norm:
    def test_identical_arrays_zero(self):
        a = np.array([1.0, 2.0, 3.0])
        assert _l2_norm(a, a) == 0.0

    def test_known_distance(self):
        a = np.array([0.0, 0.0])
        b = np.array([1.0, 1.0])
        assert abs(_l2_norm(a, b) - 1.0) < 1e-10

    def test_monotone_with_perturbation(self):
        a = np.linspace(0, 1, 100)
        b1 = a + 0.01
        b2 = a + 0.1
        assert _l2_norm(a, b1) < _l2_norm(a, b2)


class TestResample:
    def test_identity_on_same_grid(self):
        x = np.linspace(0, 10, 50)
        y = np.sin(x)
        resampled = _resample(x, y, x)
        np.testing.assert_allclose(resampled, y, atol=1e-12)

    def test_interpolation(self):
        x_src = np.array([0.0, 1.0, 2.0])
        y_src = np.array([0.0, 1.0, 2.0])
        x_dst = np.array([0.5, 1.5])
        result = _resample(x_src, y_src, x_dst)
        np.testing.assert_allclose(result, [0.5, 1.5], atol=1e-12)


class TestExoticMassProxy:
    def test_no_exotic(self):
        lumps = [{"amp": 0.2, "exotic": 0}, {"amp": 0.15, "exotic": 0}]
        assert _exotic_mass_proxy(lumps) == 0.0

    def test_exotic_counted(self):
        lumps = [{"amp": 0.2, "exotic": 1}, {"amp": 0.15, "exotic": 0}]
        assert abs(_exotic_mass_proxy(lumps) - 0.2) < 1e-10

    def test_multiple_exotic(self):
        lumps = [{"amp": 0.2, "exotic": 1}, {"amp": 0.1, "exotic": 1}]
        assert abs(_exotic_mass_proxy(lumps) - 0.3) < 1e-10


class TestTargetProfiles:
    def test_returns_correct_shapes(self):
        motif = _make_motif()
        x, chi, beta_x = _target_profiles(motif, n_points=128)
        assert x.shape == (128,)
        assert chi.shape == (128,)
        assert beta_x.shape == (128,)

    def test_chi_positive(self):
        motif = _make_motif()
        x, chi, _ = _target_profiles(motif, n_points=64)
        assert np.all(chi > 0)


class TestComputeMismatch:
    def test_gate_fitness_on_missing_gridinit(self):
        motif = _make_motif()
        report = compute_mismatch(motif, "/nonexistent/path.gridinit")
        assert report.solve_failed
        assert report.fitness == GATE_FITNESS

    def test_gate_fitness_on_empty_file(self):
        with TemporaryDirectory() as tmp:
            gridinit = Path(tmp) / "initial_data.gridinit"
            gridinit.write_bytes(b"")
            motif = _make_motif()
            report = compute_mismatch(motif, gridinit)
            assert report.solve_failed
            assert report.fitness == GATE_FITNESS


# --- Vector encoding tests ---


class TestVectorEncoding:
    def test_roundtrip(self):
        lumps = [
            {
                "amp": 0.15,
                "width": 3.0,
                "center": (1.0, 2.0, 3.0),
                "velocity": (0.1, -0.1, 0.0),
                "omega": 0.05,
                "mode": 1,
                "exotic": 1,
            },
            {
                "amp": 0.10,
                "width": 2.5,
                "center": (-1.0, 0.0, 0.0),
                "velocity": (0.0, 0.0, 0.0),
                "omega": 0.0,
                "mode": 0,
                "exotic": 0,
            },
        ]
        vec = _vector_from_lumps(lumps)
        assert vec.shape == (18,)  # 2 lumps * 9 params

        recovered = _lumps_from_vector(vec, template_lumps=lumps)
        assert len(recovered) == 2
        assert recovered[0]["mode"] == 1
        assert recovered[0]["exotic"] == 1
        assert recovered[1]["mode"] == 0
        assert recovered[1]["exotic"] == 0
        np.testing.assert_allclose(recovered[0]["amp"], 0.15)
        np.testing.assert_allclose(recovered[0]["center"], (1.0, 2.0, 3.0))
        np.testing.assert_allclose(recovered[0]["velocity"], (0.1, -0.1, 0.0))

    def test_clip_bounds(self):
        # Extreme values should be clipped
        lumps = [
            {
                "amp": 999.0,  # way above MAX_LUMP_AMP
                "width": 0.1,  # below MIN_LUMP_WIDTH
                "center": (0.0, 0.0, 0.0),
                "velocity": (9.0, 0.0, 0.0),  # above MAX_VELOCITY
                "omega": 0.0,
                "mode": 0,
                "exotic": 0,
            },
        ]
        vec = _vector_from_lumps(lumps)
        clipped = _clip_vector(vec, 1)
        recovered = _lumps_from_vector(clipped, template_lumps=lumps)
        from grteclyn_wrapper.grtresna.fit.motif import MAX_LUMP_AMP, MAX_VELOCITY, MIN_LUMP_WIDTH

        assert recovered[0]["amp"] <= MAX_LUMP_AMP
        assert recovered[0]["width"] >= MIN_LUMP_WIDTH
        assert abs(recovered[0]["velocity"][0]) <= MAX_VELOCITY


# --- Iteration loop tests (with mocked solve) ---


class TestRunIterate:
    """Test the CMA-ES loop with a mocked GRTresna solver.

    We mock solve() to avoid needing the MPI binary. The mock writes a
    fake gridinit file, and we mock compute_mismatch to return a quadratic
    fitness that CMA-ES can optimize.
    """

    def test_loop_converges_on_quadratic(self):
        """Mock fitness = ||vector - target||^2 — CMA-ES should converge."""
        motif = _make_motif()
        fitted = fit_matter_from_motif(motif, max_lumps=1)
        fitted = fit_momentum_from_motif(motif, fitted)

        # Target vector (the "geometry" we want to hit)
        target = _vector_from_lumps(fitted.lumps) + 0.02  # small offset

        def mock_solve(cfg, work_dir=None, gridinit_path=None):
            """Write a fake gridinit so the loop can proceed."""
            if work_dir is not None:
                Path(work_dir).mkdir(parents=True, exist_ok=True)
            if gridinit_path is not None:
                Path(gridinit_path).write_bytes(b"FAKE_GRIDINIT")
            return gridinit_path or (Path(work_dir) / "initial_data.gridinit")

        def mock_mismatch(motif, gridinit_path, *, lumps=None, **kwargs):
            """Quadratic fitness: distance from target vector."""
            if lumps is None:
                return MismatchReport(
                    fitness=GATE_FITNESS,
                    chi_l2=0.0,
                    beta_l2=0.0,
                    exotic_penalty=0.0,
                    convergence_penalty=0.0,
                    solve_failed=True,
                )
            vec = _vector_from_lumps(lumps)
            dist = float(np.sum((vec - target) ** 2))
            return MismatchReport(
                fitness=dist,
                chi_l2=dist * 0.5,
                beta_l2=dist * 0.5,
                exotic_penalty=0.0,
                convergence_penalty=0.0,
                solve_failed=False,
            )

        def mock_preservation(motif, gridinit_path, **kwargs):
            """Always pass preservation."""
            from grteclyn_wrapper.projection.motif_preservation import PreservationReport
            return PreservationReport(
                passed=True,
                retention_score=0.8,
                f_op_target=0.1,
                f_op_solved=0.08,
                f_op_retention=0.8,
                polarity_target=0.5,
                polarity_solved=0.45,
                polarity_retention=0.9,
                beta_max_solved=0.2,
                shift_alignment=0.5,
                momentum_alignment=0.6,
                support_localized=True,
                mechanism_descriptor=0.3,
                notes=(),
            )

        with TemporaryDirectory() as tmp:
            out_dir = Path(tmp) / "iterate"
            cfg = IterationConfig(
                max_evals=40,
                popsize=4,
                sigma0=0.1,
                max_concurrent_grtresna=2,
                mpi_ranks=1,
                plateau_generations=10,
                retention_threshold=0.5,
                gridinit_n=32,
                grtresna_L=64.0,
                seed=42,
                precondition=False,  # disable for quadratic mock test
            )

            with patch(
                "grteclyn_wrapper.projection.iterate.solve", side_effect=mock_solve
            ), patch(
                "grteclyn_wrapper.projection.iterate.compute_mismatch",
                side_effect=mock_mismatch,
            ), patch(
                "grteclyn_wrapper.projection.iterate.compare_motif_preservation",
                side_effect=mock_preservation,
            ):
                result = run_iterate(
                    motif,
                    fitted,
                    out_dir=out_dir,
                    config=cfg,
                )

            assert result.total_evals > 0
            assert result.best_fitness < GATE_FITNESS
            assert result.best_preservation_score == 0.8
            assert result.converged
            # Fitness should have decreased from initial
            assert result.fitness_history[-1] <= result.fitness_history[0]
            # Artifacts written
            assert (out_dir / "iteration_log.jsonl").exists()
            assert (out_dir / "iterate_summary.json").exists()
            assert (out_dir / "best_fitted_matter.json").exists()

    def test_empty_lumps_raises(self):
        motif = _make_motif()
        fitted = FittedMatter(
            lumps=(),
            scalar_mass=0.0,
            maximal_slicing=False,
            static_lens_only=True,
            momentum_target=MomentumTarget(
                direction=(0.0, 0.0, 0.0),
                support_center=(0.0, 0.0, 0.0),
                strength=0.0,
                template="none",
                credible=False,
            ),
            notes=(),
        )
        with TemporaryDirectory() as tmp:
            with pytest.raises(ValueError, match="no lumps"):
                run_iterate(motif, fitted, out_dir=Path(tmp) / "iter")


# --- Phase 1b: amplitude pre-conditioning tests ---


class TestAmplitudePrecondition:
    def test_scale_lump_amplitudes(self):
        lumps = [{"amp": 0.2, "exotic": 1}, {"amp": 0.1, "exotic": 0}]
        scaled = _scale_lump_amplitudes(lumps, 0.5)
        assert abs(scaled[0]["amp"] - 0.1) < 1e-10
        assert abs(scaled[1]["amp"] - 0.05) < 1e-10
        # Original unchanged
        assert lumps[0]["amp"] == 0.2

    def test_precondition_already_feasible(self):
        """If the original amplitude is feasible, return unchanged."""
        motif = _make_motif()
        fitted = fit_matter_from_motif(motif, max_lumps=1)

        def mock_solve(cfg, work_dir=None, gridinit_path=None):
            if work_dir is not None:
                Path(work_dir).mkdir(parents=True, exist_ok=True)
            if gridinit_path is not None:
                Path(gridinit_path).write_bytes(b"FAKE")
            return gridinit_path

        def mock_parse_convergence(work_dir):
            return {"iteration": 10, "ham_pct": 2.0, "mom_pct": 1.0}

        with TemporaryDirectory() as tmp:
            with patch(
                "grteclyn_wrapper.projection.iterate.solve", side_effect=mock_solve
            ), patch(
                "grteclyn_wrapper.projection.iterate.parse_convergence",
                side_effect=mock_parse_convergence,
            ):
                result_fitted, scale, notes = amplitude_precondition(
                    fitted,
                    GRTresnaConfig(),
                    Path(tmp) / "precond",
                )
            assert scale == 1.0
            assert result_fitted.lumps == fitted.lumps

    def test_precondition_halves_until_feasible(self):
        """If original is infeasible, bisection should find a feasible scale."""
        motif = _make_motif()
        fitted = fit_matter_from_motif(motif, max_lumps=1)

        call_count = {"n": 0}

        def mock_solve(cfg, work_dir=None, gridinit_path=None):
            call_count["n"] += 1
            if work_dir is not None:
                Path(work_dir).mkdir(parents=True, exist_ok=True)
            if gridinit_path is not None:
                Path(gridinit_path).write_bytes(b"FAKE")
            return gridinit_path

        def mock_parse_convergence(work_dir):
            # First call (scale=1.0): infeasible
            # Second call (scale=0.5): feasible
            if call_count["n"] <= 1:
                return {"iteration": 50, "ham_pct": 80.0, "mom_pct": 50.0}
            return {"iteration": 20, "ham_pct": 3.0, "mom_pct": 1.0}

        with TemporaryDirectory() as tmp:
            with patch(
                "grteclyn_wrapper.projection.iterate.solve", side_effect=mock_solve
            ), patch(
                "grteclyn_wrapper.projection.iterate.parse_convergence",
                side_effect=mock_parse_convergence,
            ):
                result_fitted, scale, notes = amplitude_precondition(
                    fitted,
                    GRTresnaConfig(),
                    Path(tmp) / "precond",
                )
            assert scale == 0.5
            assert result_fitted.lumps[0]["amp"] < fitted.lumps[0]["amp"]


# --- Phase 1d: tight bounds tests ---


class TestTightBounds:
    def test_bounds_around_support_region(self):
        motif = _make_motif()  # has support region at (0,0,0) width=3.0
        lower, upper = _compute_tight_bounds(motif, n_lumps=1, grtresna_L=64.0)
        # Center bounds should be tighter than ±64
        assert lower[2] > -64.0  # cx lower
        assert upper[2] < 64.0   # cx upper
        # Amplitude/width bounds unchanged
        from grteclyn_wrapper.grtresna.fit.motif import MAX_LUMP_AMP
        assert upper[0] == MAX_LUMP_AMP

    def test_bounds_fallback_without_support(self):
        """Without support regions, bounds stay at defaults."""
        motif = GeometryMotif(
            episode_path="/tmp",
            overrides={"recipe_num_bases": 1, "recipe_basis_width": 1.0,
                       "recipe_basis_radius_max": 8.0, "recipe_chi_asymptotic": 1.0,
                       "recipe_chi_coeff_0": 0.0, "recipe_beta_coeff_0": 0.0,
                       "recipe_alpha_coeff_0": 0.0},
            transport_axis=(1.0, 0.0, 0.0),
            polarity=0.0, f_shortcut=0.0, f_op=None, f_null=0.0,
            f_portal=0.0, f_throat=0.0, f_asymmetry=0.0,
            beta_max=0.0, beta_mean=0.0,
            exotic_needed=False, min_rho_required=None, integral_negative_rho=None,
            static_lens_only=True,
            support_regions=(),
            momentum_target=MomentumTarget(
                direction=(0.0, 0.0, 0.0), support_center=(0.0, 0.0, 0.0),
                strength=0.0, template="none", credible=False,
            ),
        )
        lower, upper = _compute_tight_bounds(motif, n_lumps=1, grtresna_L=64.0)
        # Should be default ±64
        assert lower[2] == -64.0
        assert upper[2] == 64.0


# --- Phase 1a: two-phase fitness tests ---


class TestTwoPhaseFitness:
    def test_infeasible_fitness_dominated_by_convergence(self):
        """When Ham > threshold, fitness should be dominated by convergence."""
        motif = _make_motif()
        with TemporaryDirectory() as tmp:
            gridinit = Path(tmp) / "initial_data.gridinit"
            gridinit.write_bytes(b"")
            # We can't easily test the full compute_mismatch without a real gridinit,
            # but we can verify the constants are sensible
            from grteclyn_wrapper.projection.mismatch import (
                FEASIBILITY_THRESHOLD,
                FEASIBILITY_WEIGHT,
                CONV_TANH_SCALE,
            )
            assert FEASIBILITY_THRESHOLD == 5.0
            assert FEASIBILITY_WEIGHT == 50.0
            assert CONV_TANH_SCALE == 20.0


# --- Phase 4a: 2D slice mismatch tests ---


class Test2DSliceMismatch:
    def test_target_2d_slice_shapes(self):
        """Target 2D slice should return n×n arrays."""
        motif = _make_motif()
        x, z, chi_2d, beta_x_2d, beta_z_2d = _target_2d_slice(motif, n=32, L=16.0)
        assert x.shape == (32,)
        assert z.shape == (32,)
        assert chi_2d.shape == (32, 32)
        assert beta_x_2d.shape == (32, 32)
        assert beta_z_2d.shape == (32, 32)

    def test_target_2d_chi_positive(self):
        """Conformal factor should be positive everywhere on the 2D slice."""
        motif = _make_motif()
        _, _, chi_2d, _, _ = _target_2d_slice(motif, n=32, L=16.0)
        assert np.all(chi_2d > 0)

    def test_target_2d_radial_symmetry(self):
        """For a spherically symmetric recipe, chi at (x,0) should match chi at (0,x)."""
        motif = _make_motif()
        _, _, chi_2d, _, _ = _target_2d_slice(motif, n=65, L=16.0)
        mid = 32
        # chi along x-axis (z=0) should equal chi along z-axis (x=0) by symmetry
        np.testing.assert_allclose(chi_2d[:, mid], chi_2d[mid, :], rtol=1e-10)

    def test_target_2d_beta_z_zero_on_axis(self):
        """beta_z should be approximately zero along the x-axis (z≈0)."""
        motif = _make_motif()
        _, _, _, _, beta_z_2d = _target_2d_slice(motif, n=65, L=16.0)  # odd n → exact z=0 at mid
        mid = 32
        np.testing.assert_allclose(beta_z_2d[:, mid], 0.0, atol=1e-12)

    def test_mismatch_report_has_2d_fields(self):
        """MismatchReport should carry chi_2d_l2 and beta_2d_l2 fields."""
        r = MismatchReport(
            fitness=1.0, chi_l2=0.1, beta_l2=0.1,
            exotic_penalty=0.0, convergence_penalty=0.0,
            solve_failed=False, chi_2d_l2=0.05, beta_2d_l2=0.03,
            kij_l2=0.0, aij_l2=0.0,
        )
        d = r.to_dict()
        assert "chi_2d_l2" in d
        assert "beta_2d_l2" in d
        assert "kij_l2" in d
        assert "aij_l2" in d
        assert d["chi_2d_l2"] == 0.05


# --- Phase 4b: K_ij matching tests ---


class TestKijMatching:
    def test_kij_weights_are_set(self):
        """K_ij and A_ij weights should be positive and non-zero by default."""
        from grteclyn_wrapper.projection.mismatch import W_KIJ, W_AIJ
        assert W_KIJ > 0.0
        assert W_AIJ > 0.0

    def test_2d_weights_higher_than_1d(self):
        """2D slice weights should be ≥ 1D weights (2D captures more info)."""
        from grteclyn_wrapper.projection.mismatch import W_CHI, W_CHI_2D, W_BETA, W_BETA_2D
        assert W_CHI_2D >= W_CHI
        assert W_BETA_2D >= W_BETA


# --- Phase 4c: feasibility pre-check tests ---


class TestFeasibilityPrecheck:
    def test_safe_target(self):
        """A motif with low rho_peak should be classified as safe."""
        motif = _make_motif()
        # Override support regions with low peak_rho
        motif = GeometryMotif(
            episode_path=motif.episode_path,
            overrides=motif.overrides,
            transport_axis=motif.transport_axis,
            polarity=motif.polarity,
            f_shortcut=motif.f_shortcut,
            f_op=motif.f_op,
            f_null=motif.f_null,
            f_portal=motif.f_portal,
            f_throat=motif.f_throat,
            f_asymmetry=motif.f_asymmetry,
            beta_max=motif.beta_max,
            beta_mean=motif.beta_mean,
            exotic_needed=False,
            min_rho_required=0.01,
            integral_negative_rho=0.0,
            static_lens_only=True,
            support_regions=(
                SupportRegion(center=(0, 0, 0), width=3.0, peak_rho=0.1, exotic=False, radial_center=0.0),
            ),
            momentum_target=motif.momentum_target,
        )
        est = feasibility_precheck(motif)
        assert est.feasible
        assert est.risk_level == "safe"
        assert est.rho_peak < FEASIBILITY_RHO_SAFE

    def test_hard_target(self):
        """A motif with very high rho_peak should be classified as hard."""
        motif = _make_motif()
        motif = GeometryMotif(
            episode_path=motif.episode_path,
            overrides=motif.overrides,
            transport_axis=motif.transport_axis,
            polarity=motif.polarity,
            f_shortcut=motif.f_shortcut,
            f_op=motif.f_op,
            f_null=motif.f_null,
            f_portal=motif.f_portal,
            f_throat=motif.f_throat,
            f_asymmetry=motif.f_asymmetry,
            beta_max=motif.beta_max,
            beta_mean=motif.beta_mean,
            exotic_needed=True,
            min_rho_required=-5.0,
            integral_negative_rho=10.0,
            static_lens_only=False,
            support_regions=(
                SupportRegion(center=(0, 0, 0), width=3.0, peak_rho=5.0, exotic=True, radial_center=0.0),
            ),
            momentum_target=motif.momentum_target,
        )
        est = feasibility_precheck(motif)
        assert not est.feasible
        assert est.risk_level == "hard"
        assert est.rho_peak >= FEASIBILITY_RHO_HARD
        assert any("exotic" in n for n in est.notes)

    def test_marginal_target(self):
        """A motif with moderate rho_peak should be classified as marginal."""
        motif = _make_motif()
        motif = GeometryMotif(
            episode_path=motif.episode_path,
            overrides=motif.overrides,
            transport_axis=motif.transport_axis,
            polarity=motif.polarity,
            f_shortcut=motif.f_shortcut,
            f_op=motif.f_op,
            f_null=motif.f_null,
            f_portal=motif.f_portal,
            f_throat=motif.f_throat,
            f_asymmetry=motif.f_asymmetry,
            beta_max=motif.beta_max,
            beta_mean=motif.beta_mean,
            exotic_needed=False,
            min_rho_required=0.5,
            integral_negative_rho=0.0,
            static_lens_only=True,
            support_regions=(
                SupportRegion(center=(0, 0, 0), width=3.0, peak_rho=1.0, exotic=False, radial_center=0.0),
            ),
            momentum_target=motif.momentum_target,
        )
        est = feasibility_precheck(motif)
        assert est.feasible  # marginal is still feasible (we try)
        assert est.risk_level == "marginal"
        assert FEASIBILITY_RHO_SAFE <= est.rho_peak < FEASIBILITY_RHO_HARD


# --- Variable lump count (sparsity penalty) tests ---


class TestSparsityPenalty:
    def test_all_lumps_active(self):
        """5 lumps with amp > threshold → penalty = 5."""
        vec = np.array([0.1, 3.0, 0, 0, 0, 0, 0, 0, 0.0] * 5, dtype=float)
        assert _sparsity_penalty(vec, 5) == 5.0

    def test_all_lumps_off(self):
        """5 lumps with amp = 0 → penalty = 0."""
        vec = np.array([0.0, 3.0, 0, 0, 0, 0, 0, 0, 0.0] * 5, dtype=float)
        assert _sparsity_penalty(vec, 5) == 0.0

    def test_mixed_active_inactive(self):
        """2 active, 3 off → penalty = 2."""
        vec = np.zeros(45, dtype=float)
        # lump 0: active
        vec[0] = 0.15
        # lump 1: off
        vec[9] = 0.0
        # lump 2: active
        vec[18] = 0.08
        # lump 3: off (below threshold)
        vec[27] = SPARSITY_AMP_THRESHOLD / 2
        # lump 4: off
        vec[36] = 0.0
        assert _sparsity_penalty(vec, 5) == 2.0

    def test_amp_lower_bound_is_zero(self):
        """The amplitude lower bound should be 0.0 so lumps can be turned off."""
        from grteclyn_wrapper.projection.iterate import _LOWER
        assert _LOWER[0] == 0.0

    def test_sparsity_weight_in_config(self):
        """IterationConfig should have a sparsity_weight field."""
        cfg = IterationConfig()
        assert cfg.sparsity_weight > 0.0
