"""Tests for the f_geo_depth scoring scalarization (raw uncapped depth)."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.score.objectives import _f_geo_depth_total


def _write_minimal_episode(root: Path) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1.0 0.98 0.02 0 0 0 5.0 0.05 5.0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_depth_dominates_all_shaping() -> None:
    # A 5% raw shortcut must beat a candidate with every shaping component
    # maxed out but no shortcut at all.
    shaped_only = _f_geo_depth_total(
        {
            "ftl_geo_depth": 0.0,
            "operational_ftl_geodesic": 1.0,
            "curvature_activity": 1.0,
        },
        [],
    )
    geo_only = _f_geo_depth_total({"ftl_geo_depth": 0.05}, [])
    assert geo_only > shaped_only


def test_depth_is_uncapped_past_saturation_target() -> None:
    # Under f_geo_max both of these saturate to the same reward; here the
    # deeper corridor must strictly win, linearly.
    at_target = _f_geo_depth_total({"ftl_geo_depth": 0.20}, [])
    deeper = _f_geo_depth_total({"ftl_geo_depth": 0.38}, [])
    assert deeper - at_target == 10000.0 * (0.38 - 0.20)


def test_survival_stability_confinement_have_no_effect() -> None:
    base = {"ftl_geo_depth": 0.3, "nontriviality_gate": 1.0}
    clean = _f_geo_depth_total(dict(base), [])
    with_health = _f_geo_depth_total(
        {
            **base,
            "survival": 1.0,
            "stability": 1.0,
            "constraint_health": 1.0,
            "confinement_final_frac": 1.0,
            "structural_persistence": 1.0,
            "ftl_persistence": 1.0,
            "exotic_penalty": -5.0,
        },
        [],
    )
    assert with_health == clean


def test_horizon_and_pump_penalties_still_cost() -> None:
    base = {"ftl_geo_depth": 0.3}
    clean = _f_geo_depth_total(dict(base), [])
    collapsed = _f_geo_depth_total({**base, "horizon_penalty": -1.0}, [])
    pumped = _f_geo_depth_total({**base, "pump_energy_penalty": -1.0}, [])
    assert collapsed < clean
    assert pumped < clean


def test_shaping_gives_gradient_at_zero_depth() -> None:
    flat = _f_geo_depth_total({}, [])
    bending = _f_geo_depth_total({"operational_ftl_geodesic": 0.1}, [])
    assert bending > flat


def test_geo_magnitude_is_uncapped() -> None:
    # The old min(..., 1.0) saturation at GEO_FTL_TARGET made every mode blind
    # to depth beyond 20% (qball_traj_fgeo_v1: the deepest corridor on record
    # scored 4x below a shallower survivor).  Depth must keep paying linearly.
    from grteclyn_wrapper.metrics.score.ftl import (
        GEO_FTL_FLOOR,
        GEO_FTL_TARGET,
        _geo_magnitude,
    )

    assert _geo_magnitude(GEO_FTL_TARGET) == 1.0
    assert _geo_magnitude(0.38) > _geo_magnitude(0.20) > _geo_magnitude(0.10)
    span = GEO_FTL_TARGET - GEO_FTL_FLOOR
    assert _geo_magnitude(0.4) - _geo_magnitude(0.3) == pytest.approx(
        (0.4 - 0.3) / span
    )
    assert _geo_magnitude(0.0) == 0.0  # floored, never negative


def test_f_geo_depth_runs_end_to_end() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(
            metrics, target_stop_time=2.0, objective_mode="f_geo_depth"
        )
        assert score.total == score.total  # finite
        assert any("f_geo_depth" in n for n in score.notes)
