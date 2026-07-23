"""Tests for the pure-geometry MAP-Elites atlas."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from grteclyn_wrapper.grtresna.io.gridinit import read_gridinit
from grteclyn_wrapper.initial_data.adm_stationary import pack_ccz4_grid
from grteclyn_wrapper.initial_data.teo import TeoWormholeConfig, build_adm, build_grid
from grteclyn_wrapper.metrics.probes.ftl.geodesic import build_metric_3d_from_gridinit
from grteclyn_wrapper.search.geometry_atlas import (
    GeometryAtlasConfig,
    GeometryGenomeConfig,
    RenderConfig,
    calibrate_atlas_probe,
    evaluate_genome,
    mutate_genome,
    render_and_write,
    render_genome,
    run_geometry_atlas,
    run_geometry_cmaes,
    sample_genome,
    seed_alcubierre_genome,
    zero_genome,
)
from grteclyn_wrapper.search.geometry_atlas.genome import (
    PARAMS_PER_CENTER,
    compact_envelope,
    decode_fields_at_points,
    fibonacci_centers,
    unpack_alcubierre,
)
from grteclyn_wrapper.search.geometry_atlas.score import (
    cell_from_descriptors,
    descriptor_axes,
    probe_half_length,
)
from grteclyn_wrapper.search.qd_search.archive import QDArchive


def test_zero_genome_is_minkowski():
    genome = zero_genome(GeometryGenomeConfig(n_centers=5, support_radius=10.0))
    pts = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [9.0, 0.0, 0.0]])
    alpha, beta, gamma, kij = decode_fields_at_points(genome, pts)
    assert np.allclose(alpha, 1.0)
    assert np.allclose(beta, 0.0)
    assert np.allclose(gamma, np.eye(3))
    assert np.allclose(kij, 0.0)


def test_sample_and_mutate_deterministic_and_bounded():
    cfg = GeometryGenomeConfig(n_centers=4, alpha_amp=0.1, beta_amp=0.2, log_metric_amp=0.05)
    rng1 = np.random.default_rng(0)
    rng2 = np.random.default_rng(0)
    g1 = sample_genome(rng1, cfg)
    g2 = sample_genome(rng2, cfg)
    assert np.allclose(g1.coeffs, g2.coeffs)
    assert np.allclose(g1.centers, g2.centers)
    alpha_slots = g1.coeffs[: cfg.n_centers * PARAMS_PER_CENTER : PARAMS_PER_CENTER]
    assert np.all(np.abs(alpha_slots) <= cfg.alpha_amp + 1e-12)

    child = mutate_genome(g1, np.random.default_rng(1), sigma=0.2)
    assert child.coeffs.shape == g1.coeffs.shape
    child_alpha = child.coeffs[: cfg.n_centers * PARAMS_PER_CENTER : PARAMS_PER_CENTER]
    assert np.all(np.abs(child_alpha) <= cfg.alpha_amp + 1e-12)


def test_compact_envelope_and_centers():
    r = np.array([0.0, 5.0, 10.0, 12.0])
    env = compact_envelope(r, 10.0)
    assert env[0] == pytest.approx(1.0)
    assert env[-2] == pytest.approx(0.0)
    assert env[-1] == pytest.approx(0.0)
    centers = fibonacci_centers(7, 4.0)
    assert centers.shape == (7, 3)
    assert np.allclose(centers[0], 0.0)


def test_rendered_metric_is_spd_and_asymptotically_flat(tmp_path: Path):
    cfg = GeometryGenomeConfig(n_centers=5, support_radius=10.0, rbf_width=3.0)
    genome = sample_genome(np.random.default_rng(3), cfg)
    rendered = render_genome(genome, RenderConfig(n=12, L=48.0))
    eig = np.linalg.eigvalsh(rendered.gamma)
    assert float(np.min(eig)) > 0.0
    assert float(np.min(rendered.alpha)) > 0.0
    assert rendered.diagnostics["boundary_max_abs_alpha_m1"] < 0.05
    assert rendered.diagnostics["boundary_max_abs_beta"] < 0.05

    path = tmp_path / "g.gridinit"
    render_and_write(genome, path, RenderConfig(n=12, L=48.0))
    grid = read_gridinit(path)
    assert grid.data.shape[-1] == 37
    assert "chi" in grid.comp_names


def test_minkowski_render_has_near_zero_exotic_cost():
    genome = zero_genome(GeometryGenomeConfig(n_centers=3, support_radius=8.0))
    rendered = render_genome(genome, RenderConfig(n=10, L=40.0))
    assert rendered.diagnostics["integral_negative_rho"] == pytest.approx(0.0, abs=1e-8)
    assert rendered.diagnostics["max_abs_alpha_m1"] == pytest.approx(0.0, abs=1e-12)
    # Stationary Minkowski constraints should be near truncation noise.
    assert rendered.diagnostics["ham_l2"] < 1e-6
    assert rendered.diagnostics["mom_l2"] < 1e-6


def test_gridinit_metric_bridge_roundtrip(tmp_path: Path):
    genome = zero_genome(GeometryGenomeConfig(n_centers=3, support_radius=8.0))
    path = tmp_path / "flat.gridinit"
    render_and_write(genome, path, RenderConfig(n=10, L=40.0))
    g, origin, spacing = build_metric_3d_from_gridinit(path)
    assert g.shape[-2:] == (4, 4)
    # Minkowski: g_tt ≈ -1, spatial ≈ I.
    assert g[5, 5, 5, 0, 0] == pytest.approx(-1.0, abs=1e-8)
    assert g[5, 5, 5, 1, 1] == pytest.approx(1.0, abs=1e-8)
    assert origin.shape == (3,)
    assert all(s > 0 for s in spacing)


def test_teo_still_builds_after_adm_extract():
    cfg = TeoWormholeConfig(nx=16, ny=16, nz=16, Lx=32, Ly=32, Lz=32, center=(16, 16, 16))
    alpha, beta, gamma = build_adm(cfg)
    assert alpha.shape == (16, 16, 16)
    grid = build_grid(cfg)
    assert grid.data.shape[-1] == 37
    packed = pack_ccz4_grid(alpha, beta, gamma, cfg.dx, cfg.origin)
    assert packed.data.shape == grid.data.shape


def test_descriptors_and_cells():
    x, y = descriptor_axes(0.25, 1.0e-3)
    assert 0.0 <= x <= 1.0
    assert 0.0 <= y <= 1.0
    cell = cell_from_descriptors((0.0, 0.0), bins=8)
    assert cell == (0, 0)
    cell2 = cell_from_descriptors((0.999, 0.999), bins=8)
    assert cell2 == (7, 7)


def test_shift_fraction_separates_families():
    """Morphology descriptor: Alcubierre (shift) ≫ conformal/lapse deformation."""
    cfg = GeometryGenomeConfig(n_centers=3, support_radius=8.0)
    render = RenderConfig(n=16, L=48.0)

    alc = render_genome(seed_alcubierre_genome(cfg), render)
    alc_shift = alc.diagnostics["shift_fraction"]

    lapse = zero_genome(cfg)
    # Push a pure lapse bump (alpha slot of the first center) with zero shift.
    lapse.coeffs[0] = 0.2
    lens = render_genome(lapse, render)
    lens_shift = lens.diagnostics["shift_fraction"]

    assert 0.0 <= lens_shift <= 1.0
    assert alc_shift > 0.5
    assert alc_shift > lens_shift


def test_evaluate_minkowski_genome():
    genome = zero_genome(GeometryGenomeConfig(n_centers=3, support_radius=8.0))
    ev = evaluate_genome(
        genome,
        eval_id=0,
        render_cfg=RenderConfig(n=10, L=40.0),
        bins=8,
        n_rays=3,
        compute_ff=False,
    )
    assert not ev.rejected
    assert ev.f_geo == pytest.approx(0.0, abs=1e-8)
    assert ev.integral_negative_rho == pytest.approx(0.0, abs=1e-8)


def test_alcubierre_seed_produces_shift():
    cfg = GeometryGenomeConfig(n_centers=3, support_radius=10.0, enable_alcubierre=True)
    genome = seed_alcubierre_genome(cfg)
    v, r, s, axis = unpack_alcubierre(genome)
    assert v > 0.0
    assert axis == 0
    pts = np.array([[0.0, 0.0, 0.0]])
    _a, beta, _g, _k = decode_fields_at_points(genome, pts)
    assert abs(float(beta[0, axis])) > 0.1


def test_probe_half_length_caps():
    assert probe_half_length(support_radius=12.0, box_length=64.0) == pytest.approx(18.0)
    assert probe_half_length(support_radius=40.0, box_length=64.0) == pytest.approx(25.6)


def test_calibrate_alcubierre_positive_control(tmp_path: Path):
    report = calibrate_atlas_probe(
        tmp_path / "cal",
        n_rays=3,
        alc_n=32,
        alc_L=12.0,
        alc_velocity=0.9,
        alc_radius=3.5,
        alc_sigma=1.2,
        cand146_path=None,
        localise=True,
    )
    assert report["verdict"]["alcubierre_localised_ok"]
    local = next(c for c in report["cases"] if c["name"] == "alcubierre_analytic")
    full = next(c for c in report["cases"] if c["name"] == "alcubierre_analytic_fullbox")
    assert local["f_geo"] >= 0.10
    assert local["f_geo"] >= full["f_geo"] - 1e-12


def test_cmaes_smoke(tmp_path: Path):
    cfg = GeometryAtlasConfig(
        runs_dir=tmp_path / "runs",
        name="cma_smoke",
        target_evals=4,
        bins=4,
        seed=3,
        n_rays=3,
        compute_ff=False,
        genome=GeometryGenomeConfig(
            n_centers=3,
            support_radius=10.0,
            rbf_width=3.0,
            enable_alcubierre=True,
            alc_velocity_max=1.2,
        ),
        render=RenderConfig(n=12, L=40.0),
    )
    root, best, summary = run_geometry_cmaes(
        cfg,
        max_evals=4,
        population_size=4,
        sigma0=0.2,
        objective="f_geo",
        alc_only=True,
    )
    assert root.exists()
    assert (root / "best.gridinit").exists()
    assert summary["evals"] == 4
    assert summary["best_f_geo"] is not None
    assert best.coeffs.shape[0] == 3 * PARAMS_PER_CENTER + 4


def test_archive_resume_and_smoke(tmp_path: Path):
    cfg = GeometryAtlasConfig(
        runs_dir=tmp_path / "runs",
        name="smoke",
        target_evals=4,
        bins=4,
        seed=11,
        batch_size=1,
        n_rays=3,
        compute_ff=False,
        resume=False,
        genome=GeometryGenomeConfig(
            n_centers=3, support_radius=8.0, rbf_width=3.0, enable_alcubierre=True
        ),
        render=RenderConfig(n=10, L=40.0),
    )
    root, archive, summary = run_geometry_atlas(cfg)
    assert root.exists()
    assert (root / "archive.json").exists()
    assert (root / "state.json").exists()
    assert summary["evals"] == 4
    assert isinstance(archive, QDArchive)
    assert len(archive.cells) >= 1
    # Elite gridinit written for at least the Minkowski seed.
    elite_files = list((root / "elites").glob("*.gridinit"))
    assert elite_files
    # Accepted elites must not be mislabeled as rejected, and must point at
    # the Stage-2 handoff gridinit path.
    for elite in archive.cells.values():
        assert elite.tier_name != "rejected"
        assert elite.episode is not None
        assert Path(elite.episode).exists()

    # Resume should not reset eval counter.
    cfg2 = GeometryAtlasConfig(
        runs_dir=tmp_path / "runs",
        name="smoke",
        target_evals=5,
        bins=4,
        seed=11,
        batch_size=1,
        n_rays=3,
        compute_ff=False,
        resume=True,
        genome=GeometryGenomeConfig(
            n_centers=3, support_radius=8.0, rbf_width=3.0, enable_alcubierre=True
        ),
        render=RenderConfig(n=10, L=40.0),
    )
    root2, archive2, summary2 = run_geometry_atlas(cfg2)
    assert root2 == root
    assert summary2["evals"] == 5
    assert (root / "evals" / "eval_000004.json").exists()
