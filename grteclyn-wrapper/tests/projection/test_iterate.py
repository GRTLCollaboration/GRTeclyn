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
from grteclyn_wrapper.initial_data.motif import (
    GeometryMotif,
    MomentumTarget,
    SupportRegion,
    extract_motif_from_episode,
)
from grteclyn_wrapper.projection.iterate import (
    IterationConfig,
    _clip_vector,
    _lumps_from_vector,
    _vector_from_lumps,
    run_iterate,
)
from grteclyn_wrapper.projection.mismatch import (
    GATE_FITNESS,
    MismatchReport,
    _exotic_mass_proxy,
    _l2_norm,
    _resample,
    _target_profiles,
    compute_mismatch,
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

        call_count = {"n": 0}

        def mock_solve(cfg, work_dir=None, gridinit_path=None):
            """Write a fake gridinit so the loop can proceed."""
            if work_dir is not None:
                Path(work_dir).mkdir(parents=True, exist_ok=True)
            if gridinit_path is not None:
                Path(gridinit_path).write_bytes(b"FAKE_GRIDINIT")
            return gridinit_path or (Path(work_dir) / "initial_data.gridinit")

        def mock_mismatch(motif, gridinit_path, *, lumps=None, **kwargs):
            """Quadratic fitness: distance from target vector."""
            call_count["n"] += 1
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
