"""Tests for the f_geo_max scoring scalarization."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.score.objectives import _f_geo_max_total


def _write_minimal_episode(root: Path) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1.0 0.98 0.02 0 0 0 5.0 0.05 5.0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_f_geo_dominates_all_shaping() -> None:
    # A 5% measured evolving shortcut must beat a candidate with every shaping
    # component maxed out but no shortcut at all.
    shaped_only = _f_geo_max_total(
        {
            "nontriviality_gate": 1.0,
            "ftl_geo_evolving": 0.0,
            "operational_ftl_geodesic": 1.0,
            "ftl_persistence": 1.0,
            "curvature_activity": 1.0,
            "survival": 1.0,
            "stability": 1.0,
            "constraint_health": 1.0,
        },
        [],
    )
    geo_only = _f_geo_max_total({"ftl_geo_evolving": 0.05}, [])
    assert geo_only > shaped_only


def test_exotic_penalty_has_no_effect() -> None:
    base = {
        "nontriviality_gate": 1.0,
        "ftl_geo_evolving": 0.1,
        "survival": 0.5,
    }
    clean = _f_geo_max_total(dict(base), [])
    exotic = _f_geo_max_total({**base, "exotic_penalty": -5.0}, [])
    assert exotic == clean


def test_horizon_and_pump_penalties_still_cost() -> None:
    base = {"nontriviality_gate": 1.0, "ftl_geo_evolving": 0.1}
    clean = _f_geo_max_total(dict(base), [])
    collapsed = _f_geo_max_total({**base, "horizon_penalty": -1.0}, [])
    pumped = _f_geo_max_total({**base, "pump_energy_penalty": -1.0}, [])
    assert collapsed < clean
    assert pumped < clean


def test_shaping_gives_gradient_at_zero_f_geo() -> None:
    flat = _f_geo_max_total({"nontriviality_gate": 1.0}, [])
    bending = _f_geo_max_total(
        {"nontriviality_gate": 1.0, "operational_ftl_geodesic": 0.1},
        [],
    )
    assert bending > flat


def test_f_geo_max_runs_end_to_end() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, objective_mode="f_geo_max")
        assert score.total == score.total  # finite
        assert any("f_geo_max" in n for n in score.notes)
