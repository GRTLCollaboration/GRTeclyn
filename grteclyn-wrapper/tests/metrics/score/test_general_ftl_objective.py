"""Tests for the general_ftl scoring scalarization."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.score.objectives import _general_ftl_total


def _write_minimal_episode(root: Path) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1.0 0.98 0.02 0 0 0 5.0 0.05 5.0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_general_ftl_ignores_warp_shaping_terms() -> None:
    notes: list[str] = []
    warp_only = _general_ftl_total(
        {
            "nontriviality_gate": 1.0,
            "shift_drive": 1.0,
            "channel_progress": 1.0,
            "operational_ftl_solved": 1.0,
            "ftl_precursor": 1.0,
            "stationary_artifact_penalty": -1.0,
            "operational_ftl_geodesic": 0.0,
            "ftl_geo_evolving": 0.0,
        },
        notes,
    )
    static_shortcut = _general_ftl_total(
        {
            "nontriviality_gate": 1.0,
            "shift_drive": 0.0,
            "channel_progress": 0.0,
            "operational_ftl_geodesic": 0.08,
            "ftl_geo_evolving": 0.0,
        },
        notes,
    )
    assert static_shortcut > warp_only
    assert any("general_ftl" in n for n in notes)


def test_general_ftl_ranks_geodesic_above_zero_shortcut() -> None:
    with_geo = _general_ftl_total(
        {"operational_ftl_geodesic": 0.1, "nontriviality_gate": 1.0},
        [],
    )
    without_geo = _general_ftl_total(
        {"operational_ftl_geodesic": 0.0, "nontriviality_gate": 1.0},
        [],
    )
    assert with_geo > without_geo


def test_general_ftl_runs_end_to_end() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, objective_mode="general_ftl")
        assert score.total == score.total  # finite
        assert any("general_ftl" in n for n in score.notes)
