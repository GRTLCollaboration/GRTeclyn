"""Tests for the spacetime_shear scoring scalarization."""

from __future__ import annotations

import argparse
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from grteclyn_wrapper.cli.grtresna_args import grtresna_solved_ftl_gate_enabled
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.score.objectives import _spacetime_shear_total
from grteclyn_wrapper.search.qd_search.descriptors import _descriptor_details


def _write_minimal_episode(root: Path) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1.0 0.98 0.02 0 0 0 5.0 0.05 5.0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_spacetime_shear_prefers_high_curvature() -> None:
    notes: list[str] = []
    low = _spacetime_shear_total(
        {"nontriviality_gate": 1.0, "curvature_activity": 0.1, "horizon_penalty": 0.0},
        notes,
    )
    high = _spacetime_shear_total(
        {"nontriviality_gate": 1.0, "curvature_activity": 0.8, "horizon_penalty": 0.0},
        notes,
    )
    assert high > low
    assert any("spacetime_shear" in n for n in notes)


def test_spacetime_shear_penalizes_horizon() -> None:
    notes: list[str] = []
    no_horizon = _spacetime_shear_total(
        {"nontriviality_gate": 1.0, "curvature_activity": 0.5, "horizon_penalty": 0.0},
        notes,
    )
    horizon = _spacetime_shear_total(
        {"nontriviality_gate": 1.0, "curvature_activity": 0.5, "horizon_penalty": -1.0},
        notes,
    )
    assert no_horizon > horizon


def test_spacetime_shear_gates_health_by_nontriviality() -> None:
    notes: list[str] = []
    trivial = _spacetime_shear_total(
        {
            "nontriviality_gate": 0.0,
            "curvature_activity": 0.0,
            "survival": 1.0,
            "horizon_penalty": 0.0,
        },
        notes,
    )
    non_trivial = _spacetime_shear_total(
        {
            "nontriviality_gate": 1.0,
            "curvature_activity": 0.5,
            "survival": 1.0,
            "horizon_penalty": 0.0,
        },
        notes,
    )
    assert non_trivial > trivial


def test_spacetime_shear_runs_end_to_end() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(
            metrics, target_stop_time=2.0, objective_mode="spacetime_shear"
        )
        assert score.total == score.total  # finite
        assert any("spacetime_shear" in n for n in score.notes)


def test_spacetime_shear_disables_solved_ftl_gate() -> None:
    args = argparse.Namespace(grtresna_solved_ftl_gate=None)
    assert not grtresna_solved_ftl_gate_enabled(
        args, use_grtresna=True, objective_mode="spacetime_shear"
    )


def test_spacetime_shear_exotic_penalty_weight_reduces_cost() -> None:
    components = {
        "nontriviality_gate": 1.0,
        "curvature_activity": 0.5,
        "horizon_penalty": 0.0,
        "exotic_penalty": -1.0,
        "survival": 0.5,
    }
    full = _spacetime_shear_total(dict(components), [])
    reduced = _spacetime_shear_total(
        dict(components), [], exotic_penalty_weight=0.2
    )
    assert reduced > full
    # Exotic penalty contribution: 40 * weight * (-1.0)
    assert reduced == full + 40.0 * (1.0 - 0.2) * 1.0


def test_score_episode_reads_exotic_penalty_weight_from_env() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)

        def _force_exotic_penalty(ctx) -> None:
            ctx.components["exotic_penalty"] = -1.0
            ctx.components.setdefault("nontriviality_gate", 1.0)

        with patch(
            "grteclyn_wrapper.metrics.score.scorer.compute_penalty_components",
            side_effect=_force_exotic_penalty,
        ):
            score_full = score_episode(
                metrics,
                target_stop_time=2.0,
                objective_mode="spacetime_shear",
                exotic_penalty_weight=1.0,
            )
            with patch.dict(
                "os.environ",
                {"SCORE_EXOTIC_PENALTY_WEIGHT": "0.2"},
            ):
                score_env = score_episode(
                    metrics,
                    target_stop_time=2.0,
                    objective_mode="spacetime_shear",
                )
        assert score_env.total > score_full.total
        assert score_env.total == score_full.total + 40.0 * (1.0 - 0.2) * 1.0


def test_spacetime_shear_descriptor() -> None:
    details = _descriptor_details(
        {"curvature_activity": 0.7, "horizon_penalty": 0.0, "nontrivial_geometry": 0.5},
        None,
        mode="spacetime_shear",
    )
    assert details["x"] == 0.7
    assert details["y"] == 1.0
    assert details["curvature_activity"] == 0.7
    assert details["horizon_free"] == 1.0

    collapsed = _descriptor_details(
        {"curvature_activity": 0.7, "horizon_penalty": -1.0},
        None,
        mode="spacetime_shear",
    )
    assert collapsed["y"] == 0.0
