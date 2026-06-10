"""Tests for the proposed-extension modules: growth metrics, ANEC/tidal
proxies, the RBF surrogate, MAP-Elites archive, and Pareto fronts."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from grteclyn_wrapper.metrics import read_episode_metrics, read_growth_metrics, score_episode
from grteclyn_wrapper.metrics.physical_metrics import compute_physical_metrics
from grteclyn_wrapper.search.pareto import ParetoPoint, dominates, pareto_front
from grteclyn_wrapper.search.qd_search import (
    QDArchive,
    Elite,
    _bin_index,
    _descriptor_details,
    _descriptors,
    _iterations_for_target_evals,
    run_qd_search,
)
from grteclyn_wrapper.core.evaluation import Evaluation
from grteclyn_wrapper.initial_data.seeds import get_seed
from grteclyn_wrapper.search import qd_search as qd_module
from grteclyn_wrapper.search.optimize import SearchDimension
from grteclyn_wrapper.search.surrogate import RBFSurrogate, screen_candidates


# --------------------------------------------------------------------------
# Growth-rate metric (closes the stability gap).
# --------------------------------------------------------------------------

def _write_series(root: Path, rows: list[tuple[float, float, float, float]], ham: list[float]) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t, lapse, chi, k in rows:
            fh.write(f"{t:g} {lapse:g} {chi:g} {k:g} 0 0 0 0 0 0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for (t, *_), h in zip(rows, ham):
            fh.write(f"{t:g} {h:g} 1e-4 0 0 0\n")


def test_growth_rate_distinguishes_equilibrium_from_slow_collapse() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        stable = root / "stable"
        collapsing = root / "collapsing"
        _write_series(
            stable,
            rows=[(0.0, 1.0, 0.90, 0.05), (1.0, 0.99, 0.899, 0.051), (2.0, 0.985, 0.898, 0.052)],
            ham=[1e-4, 1.0e-4, 1.1e-4],
        )
        _write_series(
            collapsing,
            rows=[(0.0, 1.0, 0.90, 0.05), (1.0, 0.6, 0.55, 0.4), (2.0, 0.2, 0.25, 1.4)],
            ham=[1e-4, 5e-3, 4e-2],
        )
        g_stable = read_growth_metrics(stable / "data" / "collapse_diagnostics.dat", stable / "data" / "constraint_norms.dat")
        g_coll = read_growth_metrics(collapsing / "data" / "collapse_diagnostics.dat", collapsing / "data" / "constraint_norms.dat")
        assert g_stable is not None and g_coll is not None
        assert g_stable.s_growth > 0.7
        assert g_coll.s_growth < 0.4
        assert g_coll.lambda_effective > g_stable.lambda_effective


# --------------------------------------------------------------------------
# ANEC line proxy and tidal proxy.
# --------------------------------------------------------------------------

def test_flat_space_satisfies_anec_and_has_low_tidal() -> None:
    seed = get_seed("flat_minkowski")
    phys = compute_physical_metrics(seed.overrides, L=8.0)
    assert abs(phys.anec_line) < 1e-6
    assert phys.s_anec > 0.99
    assert phys.s_tidal > 0.99
    assert not phys.has_trapped_proxy


def test_ellis_bronnikov_curvature_and_throat() -> None:
    seed = get_seed("ellis_bronnikov", b0=0.5)
    phys = compute_physical_metrics(seed.overrides, L=8.0)
    # A throat has strong curvature -> lower tidal comfort than flat space.
    assert phys.tidal_proxy > 0.0
    assert phys.s_tidal < 1.0


def test_physical_metrics_flow_into_score() -> None:
    seed = get_seed("ellis_bronnikov")
    with TemporaryDirectory() as tmp:
        episode = Path(tmp) / "ep"
        episode.mkdir()
        (episode / "metadata.json").write_text(json.dumps({"overrides": dict(seed.overrides)}), encoding="utf-8")
        data = episode / "data"
        data.mkdir()
        with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
            for t in (0.0, 1.0, 2.0):
                fh.write(f"{t:g} 0.9 {0.88 - 0.005 * t:g} 0.05 0 0 0 0 0 0 0 0 0 0\n")
        with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
            for t in (0.0, 1.0, 2.0):
                fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")
        metrics = read_episode_metrics(episode, ftl_L=8.0)
        assert metrics.growth is not None
        assert metrics.physical is not None
        score = score_episode(metrics, target_stop_time=2.0)
        for key in ("constraint_growth", "anec_condition", "tidal_comfort"):
            assert key in score.components


# --------------------------------------------------------------------------
# RBF surrogate.
# --------------------------------------------------------------------------

def test_surrogate_interpolates_and_screens() -> None:
    rng = np.random.default_rng(0)
    lower = np.array([-1.0, -1.0])
    upper = np.array([1.0, 1.0])

    def truth(x: np.ndarray) -> np.ndarray:
        return -(x[:, 0] ** 2 + x[:, 1] ** 2)  # peak at origin

    x_train = rng.uniform(-1, 1, size=(60, 2))
    y_train = truth(x_train)
    model = RBFSurrogate(lower=lower, upper=upper).fit(x_train, y_train)

    # Prediction near origin should beat prediction near a corner.
    near = model.predict(np.array([[0.05, -0.05]]))[0]
    far = model.predict(np.array([[0.95, 0.95]]))[0]
    assert near > far

    candidates = rng.uniform(-1, 1, size=(8, 2))
    eval_idx, predicted = screen_candidates(model, candidates, keep_fraction=0.5, min_eval=1)
    assert len(eval_idx) >= 4
    assert len(predicted) == 8
    # The globally-best (closest to origin) candidate must be selected.
    best = int(np.argmax(predicted))
    assert best in eval_idx


# --------------------------------------------------------------------------
# MAP-Elites archive.
# --------------------------------------------------------------------------

def test_map_elites_archive_keeps_best_per_cell() -> None:
    arch = QDArchive(bins=4)
    a = Elite(cell=(1, 2), score=3.0, descriptors=(0.3, 0.6), params={}, episode=None)
    b = Elite(cell=(1, 2), score=5.0, descriptors=(0.3, 0.6), params={}, episode=None)
    c = Elite(cell=(0, 0), score=1.0, descriptors=(0.0, 0.0), params={}, episode=None)
    assert arch.insert(a)
    assert arch.insert(b)  # higher score replaces
    assert not arch.insert(a)  # lower score rejected
    assert arch.insert(c)
    assert arch.cells[(1, 2)].score == 5.0
    assert len(arch.cells) == 2
    assert arch.best.score == 5.0


def test_descriptors_and_bins() -> None:
    d1, d2 = _descriptors({"ftl_shortcut": 0.5, "anec_condition": 0.25})
    assert abs(d1 - 0.5) < 1e-9
    assert abs(d2 - 0.75) < 1e-9
    assert _bin_index(0.0, 8) == 0
    assert _bin_index(0.999, 8) == 7
    assert _bin_index(0.5, 8) == 4


def test_channel_descriptors_separate_path_from_shift_blob() -> None:
    components = {
        "ftl_precursor": 0.25,
        "shift_drive": 0.36,
        "channel_progress": 0.12,
        "operational_ftl": 0.0,
    }
    metrics = {
        "general_ftl_evolved": {
            "reachable": True,
            "t_min": 16.38,
            "t_flat": 15.75,
        }
    }

    details = _descriptor_details(components, metrics, mode="channel")
    d1, d2 = _descriptors(components, metrics, mode="channel")

    assert details["path_closeness"] > 0.6
    assert abs(details["mechanism_balance"] - 0.3) < 1e-9
    assert d1 == details["path_closeness"]
    assert d2 == details["mechanism_balance"]

    slow = {
        "general_ftl_evolved": {
            "reachable": True,
            "t_min": 25.52,
            "t_flat": 15.75,
        }
    }
    assert _descriptor_details(components, slow, mode="channel")["path_closeness"] == 0.0


def test_speed_super_y_axis_tracks_superluminal_fraction() -> None:
    components = {"operational_ftl": 0.0, "ftl_persistence": 0.0}

    # Same cone-tilt (max_local_speed) but different superluminal_fraction must
    # land in different y-bins, so the new axis carries diversity signal even
    # when the horizon-free axis would not.  The descriptor reads the *solved*
    # report (general_ftl_solved), where the superluminal region survives.
    localized = {
        "general_ftl_solved": {
            "max_local_speed": 1.05,
            "superluminal_fraction": 0.06,
        }
    }
    widespread = {
        "general_ftl_solved": {
            "max_local_speed": 1.05,
            "superluminal_fraction": 0.30,
        }
    }

    d_loc = _descriptor_details(components, localized, mode="speed_super")
    d_wide = _descriptor_details(components, widespread, mode="speed_super")

    # x-axis (recalibrated cone-tilt) identical for equal speeds...
    assert abs(d_loc["x"] - d_wide["x"]) < 1e-9
    assert d_loc["speed_tilt"] == d_loc["x"]
    # ...but y-axis separates localized from widespread superluminal regions.
    # y is rescaled by the observed solved ceiling (~0.30): 0.06 -> 0.20, 0.30 -> 1.0,
    # so the realistic 0-0.30 range fills the grid instead of clustering in bin 0.
    assert d_wide["y"] > d_loc["y"]
    assert abs(d_loc["y"] - 0.20) < 1e-9
    assert abs(d_wide["y"] - 1.0) < 1e-9
    assert _bin_index(d_loc["y"], 8) < _bin_index(d_wide["y"], 8)
    # Raw (un-rescaled) fraction is preserved for diagnostics.
    assert abs(d_loc["superluminal_fraction_raw"] - 0.06) < 1e-9

    # The solved report is preferred over the (collapsed) evolved report.
    mixed = {
        "general_ftl_evolved": {"max_local_speed": 0.96, "superluminal_fraction": 0.0},
        "general_ftl_solved": {"max_local_speed": 1.15, "superluminal_fraction": 0.30},
    }
    d_mixed = _descriptor_details(components, mixed, mode="speed_super")
    assert abs(d_mixed["max_local_speed"] - 1.15) < 1e-9
    assert abs(d_mixed["y"] - 1.0) < 1e-9

    # Recalibrated x-axis spreads realistic speeds across the grid instead of
    # saturating: 0.95 -> bin 0, 1.20 -> top bins.
    floor = _descriptor_details(
        components,
        {"general_ftl_solved": {"max_local_speed": 0.95, "superluminal_fraction": 0.0}},
        mode="speed_super",
    )
    target = _descriptor_details(
        components,
        {"general_ftl_solved": {"max_local_speed": 1.20, "superluminal_fraction": 0.30}},
        mode="speed_super",
    )
    assert floor["x"] == 0.0
    assert target["x"] == 1.0


def test_qd_search_flushes_initial_trajectory(tmp_path, monkeypatch) -> None:
    def fake_evaluate_overrides(overrides, *, out_dir, name, **_kwargs):
        idx = int(name.rsplit("_", 1)[1])
        episode = out_dir / name
        episode.mkdir()
        return Evaluation(
            score=float(idx),
            components={
                "constraint_growth": 1.0,
                "ftl_precursor": 0.05 * idx,
                "shift_drive": 0.2,
                "channel_progress": 0.01 * idx,
                "operational_ftl": 0.0,
                "anec_condition": 1.0,
                "tidal_comfort": 1.0,
            },
            notes=[],
            episode_path=str(episode),
            exit_code=0,
            preflight_rejected=False,
            reason=None,
            metrics={
                "general_ftl_evolved": {
                    "reachable": True,
                    "t_flat": 10.0,
                    "t_min": 10.0 + 0.1 * idx,
                },
            },
        )

    monkeypatch.setattr(qd_module, "evaluate_overrides", fake_evaluate_overrides)

    run_qd_search(
        runs_dir=tmp_path,
        name="qd_live",
        iterations=0,
        init_random=2,
        batch_size=1,
        bins=4,
        seed=7,
        search_space=[SearchDimension("toy_param", 0.0, 1.0)],
        consume_plotfiles=False,
        check_params=False,
    )

    lines = (tmp_path / "qd_live" / "trajectory.jsonl").read_text(encoding="utf-8").splitlines()
    rows = [json.loads(line) for line in lines]

    assert [row["eval"] for row in rows] == [1, 2]
    assert all(row["episode"].endswith(f"eval_{row['eval']:06d}") for row in rows)
    assert all(row["status"] == "gpu_ok" for row in rows)
    assert (tmp_path / "qd_live" / "archive.json").exists()


def test_iterations_for_target_evals() -> None:
    assert _iterations_for_target_evals(target_evals=88, batch=8) == 10
    assert _iterations_for_target_evals(target_evals=400, batch=8) == 49
    assert _iterations_for_target_evals(
        target_evals=400, batch=8, completed_evals=88,
    ) == 39


def test_qd_search_resume_continues_eval_counter(tmp_path, monkeypatch) -> None:
    def fake_evaluate_overrides(overrides, *, out_dir, name, **_kwargs):
        idx = int(name.rsplit("_", 1)[1])
        episode = out_dir / name
        episode.mkdir()
        return Evaluation(
            score=float(idx),
            components={
                "constraint_growth": 1.0,
                "ftl_precursor": 0.05,
                "shift_drive": 0.2,
                "channel_progress": 0.01,
                "operational_ftl": 0.0,
                "anec_condition": 1.0,
                "tidal_comfort": 1.0,
            },
            notes=[],
            episode_path=str(episode),
            exit_code=0,
            preflight_rejected=False,
            reason=None,
            metrics={
                "general_ftl_evolved": {
                    "reachable": True,
                    "t_flat": 10.0,
                    "t_min": 10.2,
                },
            },
        )

    monkeypatch.setattr(qd_module, "evaluate_overrides", fake_evaluate_overrides)

    run_qd_search(
        runs_dir=tmp_path,
        name="qd_resume",
        iterations=1,
        init_random=2,
        batch_size=1,
        bins=4,
        seed=7,
        search_space=[SearchDimension("toy_param", 0.0, 1.0)],
        consume_plotfiles=False,
        check_params=False,
    )

    run_qd_search(
        runs_dir=tmp_path,
        name="qd_resume",
        resume=True,
        target_evals=4,
        batch_size=1,
        bins=4,
        seed=7,
        search_space=[SearchDimension("toy_param", 0.0, 1.0)],
        consume_plotfiles=False,
        check_params=False,
    )

    lines = (tmp_path / "qd_resume" / "trajectory.jsonl").read_text(encoding="utf-8").splitlines()
    evals = [json.loads(line)["eval"] for line in lines]
    assert evals == [1, 2, 3, 4]
    archive = QDArchive.from_dict(
        json.loads((tmp_path / "qd_resume" / "archive.json").read_text(encoding="utf-8"))
    )
    assert len(archive.cells) >= 1


# --------------------------------------------------------------------------
# Pareto front.
# --------------------------------------------------------------------------

def test_pareto_dominance_and_front() -> None:
    keys = ("ftl_shortcut", "anec_condition")
    assert dominates({"ftl_shortcut": 1.0, "anec_condition": 1.0}, {"ftl_shortcut": 0.5, "anec_condition": 0.5}, keys)
    assert not dominates({"ftl_shortcut": 1.0, "anec_condition": 0.1}, {"ftl_shortcut": 0.5, "anec_condition": 0.9}, keys)

    points = [
        ParetoPoint("a", {"ftl_shortcut": 0.9, "anec_condition": 0.1, "constraint_growth": 0.0, "tidal_comfort": 0.0}, 1.0, None),
        ParetoPoint("b", {"ftl_shortcut": 0.1, "anec_condition": 0.9, "constraint_growth": 0.0, "tidal_comfort": 0.0}, 1.0, None),
        ParetoPoint("c", {"ftl_shortcut": 0.05, "anec_condition": 0.05, "constraint_growth": 0.0, "tidal_comfort": 0.0}, 0.1, None),
    ]
    front = pareto_front(points)
    labels = {p.label for p in front}
    assert "a" in labels and "b" in labels  # both trade-off extremes
    assert "c" not in labels  # dominated by a and b


if __name__ == "__main__":
    test_growth_rate_distinguishes_equilibrium_from_slow_collapse()
    test_flat_space_satisfies_anec_and_has_low_tidal()
    test_ellis_bronnikov_curvature_and_throat()
    test_physical_metrics_flow_into_score()
    test_surrogate_interpolates_and_screens()
    test_map_elites_archive_keeps_best_per_cell()
    test_descriptors_and_bins()
    test_pareto_dominance_and_front()
    print("proposed-extension tests passed")
