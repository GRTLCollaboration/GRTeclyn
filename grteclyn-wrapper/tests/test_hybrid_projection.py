"""Tests for the geometry-first -> GRTresna projection pipeline."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from grteclyn_wrapper.grtresna.motif_fit import (
    build_grtresna_config_from_fitted,
    estimate_momentum_source,
    fit_matter_from_motif,
    fit_momentum_from_motif,
    read_fitted_matter_json,
    write_fitted_matter_json,
)
from grteclyn_wrapper.grtresna.solver import GRTresnaConfig, write_grtresna_params
from grteclyn_wrapper.initial_data.motif import (
    extract_motif_from_episode,
    read_motif_json,
    write_motif_json,
)
from grteclyn_wrapper.projection.motif_preservation import compare_motif_preservation
from grteclyn_wrapper.projection.postload_gate import (
    PostLoadGateConfig,
    evaluate_constraint_gate,
    run_postload_gate,
)
from grteclyn_wrapper.search.optimize import build_grtresna_config


def _write_geometry_episode(root: Path, *, with_motor: bool) -> Path:
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
        "recipe_beta_coeff_0": 0.35 if with_motor else 0.0,
        "recipe_beta_coeff_1": -0.15 if with_motor else 0.0,
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


def test_motif_extraction_from_fixture() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=True)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        assert motif.support_regions
        assert motif.exotic_needed
        assert not motif.static_lens_only
        assert motif.momentum_target.credible
        assert motif.momentum_target.template == "axial_boost"


def test_static_lens_motif_has_no_motor() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=False)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        assert motif.static_lens_only
        assert not motif.momentum_target.credible


def test_fit_matter_serializes_to_grtresna_config() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=True)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        fitted = fit_matter_from_motif(motif, max_lumps=2)
        fitted = fit_momentum_from_motif(motif, fitted)
        out = Path(tmp) / "fitted_matter.json"
        write_fitted_matter_json(fitted, out)
        loaded = read_fitted_matter_json(out)
        cfg = build_grtresna_config_from_fitted(loaded)
        assert cfg.scalar_mass == 0.0
        assert len(cfg.lumps) >= 1
        assert cfg.maximal_slicing
        write_grtresna_params(cfg, Path(tmp) / "grtresna_params.txt")
        text = (Path(tmp) / "grtresna_params.txt").read_text(encoding="utf-8")
        assert "scalar_mass = 0.0" in text
        assert "num_lumps" in text


def test_boosted_lump_produces_momentum_source() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=True)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        fitted = fit_momentum_from_motif(
            motif,
            fit_matter_from_motif(motif, max_lumps=1),
        )
        center = fitted.lumps[0]["center"]
        sample = (center[0] + 0.2, center[1], center[2])
        source = np.asarray(
            estimate_momentum_source(fitted.lumps, sample_point=sample),
            dtype=float,
        )
        assert float(np.linalg.norm(source)) > 1.0e-6
        assert source[0] * motif.momentum_target.direction[0] > 0.0


def test_static_lens_skips_momentum_fit() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=False)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        fitted = fit_momentum_from_motif(
            motif,
            fit_matter_from_motif(motif, max_lumps=1),
        )
        assert fitted.static_lens_only
        assert all(lump["velocity"] == (0.0, 0.0, 0.0) for lump in fitted.lumps)


def test_project_geometry_motif_fit_only() -> None:
    script = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "search"
        / "project_geometry_motif.py"
    )
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=True)
        out_dir = Path(tmp) / "projection"
        proc = subprocess.run(
            [
                sys.executable,
                str(script),
                str(episode),
                "--out-dir",
                str(out_dir),
                "--mode",
                "fit-only",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert proc.returncode == 0
        assert (out_dir / "motif.json").exists()
        assert (out_dir / "fitted_matter.json").exists()
        assert (out_dir / "momentum_target.json").exists()
        motif = read_motif_json(out_dir / "motif.json")
        assert motif.episode_path == str(episode.resolve())


def test_postload_gate_on_known_good_constraints() -> None:
    with TemporaryDirectory() as tmp:
        episode = Path(tmp) / "episode"
        data = episode / "data"
        data.mkdir(parents=True)
        with (data / "constraint_norms.dat").open("w", encoding="utf-8") as handle:
            handle.write("0 1e-6 1e-6 -0.01 0.01 0.1\n")
            handle.write("0.01 1e-6 1e-6 -0.01 0.01 0.1\n")
        gate = evaluate_constraint_gate(data / "constraint_norms.dat")
        assert gate.passed
        gate = evaluate_constraint_gate(
            data / "constraint_norms.dat",
            config=PostLoadGateConfig(max_hamiltonian_l2=1.0e-8),
        )
        assert not gate.passed


def test_postload_gate_dry_run_writes_episode() -> None:
    with TemporaryDirectory() as tmp:
        gridinit = Path(tmp) / "initial_data.gridinit"
        gridinit.write_bytes(b"")
        result = run_postload_gate(
            gridinit,
            out_dir=Path(tmp) / "gate",
            dry_run=True,
        )
        assert result.passed
        assert result.episode_path is not None
        assert "dry_run" in result.notes[0]


def test_matter_first_build_grtresna_config_unchanged() -> None:
    overrides = {
        "grtresna_ring_lumps": 5,
        "grtresna_ring_amp": 0.15,
        "grtresna_ring_width": 3.0,
        "grtresna_ring_radius": 4.0,
        "grtresna_ring_phase": 0.0,
        "grtresna_ring_tangential_velocity": 0.4,
        "grtresna_ring_radial_velocity": 0.0,
        "grtresna_ring_vertical_velocity": 0.0,
        "grtresna_ring_exotic_fraction": 0.4,
        "grtresna_ring_exotic_phase": 0.0,
        "grtresna_ring_mode": 1.0,
    }
    cfg = build_grtresna_config(overrides)
    assert len(cfg.lumps) == 5
    assert sum(lump["exotic"] for lump in cfg.lumps) == 2
    assert cfg.lumps[0]["center"][0] == 4.0


def test_motif_preservation_rejects_missing_gridinit() -> None:
    with TemporaryDirectory() as tmp:
        episode = _write_geometry_episode(Path(tmp), with_motor=True)
        motif = extract_motif_from_episode(episode, phantom=True, ftl_L=8.0)
        report = compare_motif_preservation(
            motif,
            Path(tmp) / "missing.gridinit",
            ftl_L=8.0,
        )
        assert not report.passed
        assert report.f_op_solved == 0.0


if __name__ == "__main__":
    for name in sorted(n for n in dir() if n.startswith("test_")):
        globals()[name]()
        print(f"{name} passed")
