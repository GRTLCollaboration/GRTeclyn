"""Tests for the general_ftl scoring scalarization."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from grteclyn_wrapper.cli.grtresna_args import grtresna_solved_ftl_gate_enabled
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode
from grteclyn_wrapper.metrics.aggregation.collector import _operational_precursor_satisfied
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import EvolvingGeodesicFtlReport
from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
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


def test_general_ftl_disables_solved_ftl_gate_by_default() -> None:
    args = argparse.Namespace(grtresna_solved_ftl_gate=None)
    assert not grtresna_solved_ftl_gate_enabled(
        args, use_grtresna=True, objective_mode="general_ftl"
    )
    assert grtresna_solved_ftl_gate_enabled(
        args, use_grtresna=True, objective_mode="ftl_first"
    )


def test_general_ftl_cli_can_force_solved_ftl_gate() -> None:
    args = argparse.Namespace(grtresna_solved_ftl_gate=True)
    assert grtresna_solved_ftl_gate_enabled(
        args, use_grtresna=True, objective_mode="general_ftl"
    )


def test_general_ftl_skips_operational_precursor_for_evolving_trace() -> None:
    evolved = GeneralFtlReport(
        f_op=0.0,
        t_min=18.0,
        t_flat=16.0,
        reachable=True,
        path_offaxis=False,
        notes=(),
        max_local_speed=0.94,
        superluminal_fraction=0.0,
        max_shift=0.0,
        structure_coherence=None,
    )
    assert not _operational_precursor_satisfied(evolved, objective_mode="ftl_first")
    assert _operational_precursor_satisfied(evolved, objective_mode="general_ftl")


def test_collector_runs_evolving_trace_without_coordinate_precursor() -> None:
    with TemporaryDirectory() as tmp:
        episode = Path(tmp) / "eval"
        data = episode / "data"
        small = episode / "small_data"
        data.mkdir(parents=True)
        small.mkdir(parents=True)
        (data / "collapse_diagnostics.dat").write_text(
            "0.0 1.0 1.0 0.0 0 0 0 0.0 0.1\n", encoding="utf-8"
        )
        (data / "constraint_norms.dat").write_text(
            "0.0 0.001 0.0001 -0.001 0.01 0.0\n", encoding="utf-8"
        )
        (small / "ftl_timeseries.dat").write_text(
            "# time  f_op  f_geo  geo_trustworthy  max_local_speed  "
            "superluminal_fraction  max_shift  structure_coherence  reachable  "
            "n_rays  n_reached  max_h_rel_drift\n"
            "8.0 0.0 0.0 1 0.94 0.0 0.0 nan 1 5 5 0.001\n",
            encoding="utf-8",
        )
        evolved = GeneralFtlReport(
            f_op=0.0,
            t_min=18.0,
            t_flat=16.0,
            reachable=True,
            path_offaxis=False,
            notes=(),
            max_local_speed=0.94,
            superluminal_fraction=0.0,
            max_shift=0.0,
            structure_coherence=None,
        )
        fake_report = EvolvingGeodesicFtlReport(
            f_geo=0.12,
            f_geo_frozen_peak=None,
            t_emit=8.0,
            t_arrival=14.0,
            t_flat=16.0,
            n_rays=5,
            n_reached=5,
            max_h_drift=1.0e-8,
            h_quality_ok=True,
            max_h_rel_drift=1.0e-9,
            notes=("test",),
        )
        with (
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.wait_for_plotfile_complete",
            ),
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.matter_coherence_from_plotfile",
                return_value=None,
            ),
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.find_latest_plotfile",
                return_value=episode / "plotfile_000",
            ),
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.find_recent_plotfiles",
                return_value=[episode / f"plot_{i}" for i in range(3)],
            ),
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.compute_general_ftl_from_plotfile",
                return_value=evolved,
            ),
            patch(
                "grteclyn_wrapper.metrics.aggregation.collector.compute_evolving_geodesic_ftl_from_plotfiles",
                return_value=fake_report,
            ),
        ):
            metrics = read_episode_metrics(
                episode,
                ftl_L=8.0,
                evolving_geodesic=True,
                objective_mode="general_ftl",
            )
        assert metrics.evolving_geodesic is not None
        assert metrics.evolving_geodesic.f_geo == 0.12


def test_general_ftl_exotic_penalty_weight_reduces_cost() -> None:
    components = {
        "nontriviality_gate": 1.0,
        "operational_ftl_geodesic": 0.0,
        "ftl_geo_evolving": 0.0,
        "exotic_penalty": -1.0,
        "survival": 0.5,
    }
    full = _general_ftl_total(dict(components), [])
    reduced = _general_ftl_total(
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
                objective_mode="general_ftl",
                exotic_penalty_weight=1.0,
            )
            with patch.dict(
                "os.environ", {"SCORE_EXOTIC_PENALTY_WEIGHT": "0.1"}
            ):
                score_env = score_episode(
                    metrics,
                    target_stop_time=2.0,
                    objective_mode="general_ftl",
                )
            score_explicit = score_episode(
                metrics,
                target_stop_time=2.0,
                objective_mode="general_ftl",
                exotic_penalty_weight=0.1,
            )

        assert score_env.total == score_explicit.total
        assert score_env.total > score_full.total
        # general_ftl exotic term: 40 * weight * (-1.0)
        assert math.isclose(
            score_env.total - score_full.total, 40.0 * (1.0 - 0.1), rel_tol=0.0, abs_tol=1e-9
        )

