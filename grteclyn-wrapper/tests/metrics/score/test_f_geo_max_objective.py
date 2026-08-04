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


def test_objective_mode_whitelists_are_single_sourced() -> None:
    # A mode missing from any parser's --objective-mode whitelist crashes that
    # entry point at argparse time.  For the plotfile consumer the pipeline
    # swallows the crash, which silently disables the metric-stack cache: the
    # evolving trace falls back to the last ~5 surviving plotfiles, a window
    # too short for a full ray crossing, and every evolving f_geo reads 0.
    # (This is exactly how f_geo_max shipped broken.)  So: no parser may hold
    # its own literal copy of the mode list -- everyone imports the canonical
    # tuples from grteclyn_wrapper.objective_modes.
    import inspect

    import grteclyn_wrapper.cli.parser as cli_parser
    import grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.driver as consumer_driver
    from grteclyn_wrapper.metrics.score import objectives
    from grteclyn_wrapper.objective_modes import OBJECTIVE_MODES, QD_OBJECTIVE_MODES

    assert "f_geo_max" in OBJECTIVE_MODES

    from pathlib import Path

    for module in (cli_parser, consumer_driver):
        src = Path(module.__file__).read_text(encoding="utf-8")
        assert 'choices=["weighted"' not in src, (
            f"{module.__name__} holds a private copy of the objective-mode "
            "whitelist; import OBJECTIVE_MODES / QD_OBJECTIVE_MODES instead"
        )
        assert "OBJECTIVE_MODES" in src

    # The whitelist is only honest if compute_total actually dispatches every
    # advertised mode ("weighted" is the default fall-through).
    dispatch_src = inspect.getsource(objectives.compute_total)
    for mode in OBJECTIVE_MODES:
        if mode == "weighted":
            continue
        assert f'"{mode}"' in dispatch_src, (
            f"objective mode {mode!r} is advertised but not dispatched in "
            "compute_total"
        )

    # The main parser must accept every canonical mode end-to-end.
    parser = cli_parser.build_parser()
    for mode in QD_OBJECTIVE_MODES:
        args = parser.parse_args(
            ["qd", "--objective-mode", mode, "--iterations", "1"]
        )
        assert args.objective_mode == mode


def test_f_geo_max_runs_end_to_end() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, objective_mode="f_geo_max")
        assert score.total == score.total  # finite
        assert any("f_geo_max" in n for n in score.notes)
