"""Unit tests for HQ FTL Richardson / completeness analyzer."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "campaigns"
    / "hq"
    / "ftl_convergence_analyze.py"
)


def _load_mod():
    spec = importlib.util.spec_from_file_location("ftl_convergence_analyze", _SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_richardson_monotonic_ok():
    mod = _load_mod()
    spacings = {"coarse": 0.5, "medium": 0.25, "fine": 0.125}
    values = {k: 1.0 + 2.0 * h * h for k, h in spacings.items()}
    result = mod.unequal_richardson(values, spacings)
    assert result["status"] == "ok"
    assert result["order"] == pytest.approx(2.0, rel=1e-3)
    assert result["continuum"] == pytest.approx(1.0, abs=1e-3)


def test_richardson_non_monotonic_rejected():
    mod = _load_mod()
    spacings = {"coarse": 0.5, "medium": 0.25, "fine": 0.125}
    values = {"coarse": 0.10, "medium": 0.13, "fine": 0.11}
    result = mod.unequal_richardson(values, spacings)
    assert result["status"] == "rejected_extrapolation"
    assert result["order"] is None
    assert result["continuum"] is None
    assert result["reason"] == "non_monotonic_or_zero_dmf"


def test_relative_diff_and_trusted_time():
    mod = _load_mod()
    assert mod.relative_diff(0.13, 0.117) == pytest.approx(0.1, rel=1e-6)
    t = mod.trusted_time(box_half=64.0, extraction_radius=24.0)
    assert t == pytest.approx(40.0 / (2.0**0.5), rel=1e-9)


def test_existing_hq_incomplete_for_validation(tmp_path: Path):
    mod = _load_mod()
    real = (
        Path(__file__).resolve().parents[3]
        / "runs"
        / "grtresna_promote"
        / "qball_traj_spiral_v2_hq_eval000118"
    )
    if real.is_dir() and (real / "small_data" / "evolving_geodesic.json").is_file():
        summary = mod.extract_run_summary(real)
        missing = mod.assess_completeness_for_validation(summary)
        assert "emit_sweep (>=2 launches)" in missing
        assert summary["f_geo"] == pytest.approx(0.13039151448076688, rel=1e-6)
        assert summary["h_quality_ok"] is True
        return

    ep = tmp_path / "hq"
    (ep / "small_data").mkdir(parents=True)
    (ep / "data").mkdir()
    geo = {
        "f_geo": 0.13,
        "h_quality_ok": True,
        "n_rays": 5,
        "n_reached": 5,
        "emit_sweep": [],
    }
    (ep / "small_data" / "evolving_geodesic.json").write_text(json.dumps(geo))
    (ep / "score.json").write_text(json.dumps({"score": {"total": 100.0, "components": {}}}))
    (ep / "metadata.json").write_text(json.dumps({"overrides": {"L_full": 128, "N_full": 256}}))
    summary = mod.extract_run_summary(ep)
    missing = mod.assess_completeness_for_validation(summary)
    assert "emit_sweep (>=2 launches)" in missing
    assert "psi4_mode_l2m0.dat with data rows" in missing


def test_complete_fixture_passes(tmp_path: Path):
    mod = _load_mod()
    ep = tmp_path / "hq"
    (ep / "small_data").mkdir(parents=True)
    (ep / "data").mkdir()
    geo = {
        "f_geo": 0.12,
        "h_quality_ok": True,
        "n_rays": 5,
        "n_reached": 5,
        "emit_sweep": [[0.0, 0.1, 5], [2.0, 0.12, 5]],
    }
    (ep / "small_data" / "evolving_geodesic.json").write_text(json.dumps(geo))
    (ep / "score.json").write_text(
        json.dumps({"score": {"total": 100.0, "components": {"numerical_survival": 1.0}}})
    )
    (ep / "metadata.json").write_text(
        json.dumps({"overrides": {"L_full": 128.0, "N_full": 256, "stop_time": 30.0}})
    )
    (ep / "data" / "constraint_norms.dat").write_text(
        "0.0 1e-3 1e-4 0 0 0\n1.0 2e-3 2e-4 0 0 0\n"
    )
    (ep / "small_data" / "psi4_mode_l2m0.dat").write_text(
        "# time Re Im\n0.0 1e-5 0.0\n1.0 2e-5 0.0\n"
    )
    (ep / "small_data" / "psi4_mode_l2_all.dat").write_text("# placeholder\n")
    summary = mod.extract_run_summary(ep)
    assert mod.assess_completeness_for_validation(summary) == []
