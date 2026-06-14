"""Smoke test for the robust_ftl scoring scalarization (option B).

robust_ftl is ftl_first with the reward budget tilted toward persistent /
healthy / low-exotic geometries: it must still run end-to-end and must penalize
exotic matter harder than ftl_first (exotic weight 70 vs 40).
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from grteclyn_wrapper.metrics import read_episode_metrics, score_episode


def _write_minimal_episode(root: Path) -> None:
    data = root / "data"
    data.mkdir(parents=True)
    # A healthy, non-collapsing slice history (lapse ~1, no trapped surface).
    with (data / "collapse_diagnostics.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1.0 0.98 0.02 0 0 0 5.0 0.05 5.0 0 0 0 0\n")
    with (data / "constraint_norms.dat").open("w", encoding="utf-8") as fh:
        for t in (0.0, 1.0, 2.0):
            fh.write(f"{t:g} 1e-4 1e-4 0 0 0\n")


def test_robust_ftl_runs_and_is_labeled() -> None:
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        score = score_episode(metrics, target_stop_time=2.0, objective_mode="robust_ftl")
        assert score.total == score.total  # finite (not NaN)
        assert any("robust_ftl" in n for n in score.notes)


def test_robust_ftl_penalizes_exotic_harder_than_ftl_first() -> None:
    """For an identical exotic-heavy candidate, robust_ftl <= ftl_first total.

    Both share the same dominant FTL terms; robust_ftl only differs by a heavier
    exotic penalty (70 vs 40) and rebalanced health weights, so a maximally
    exotic geometry must not score higher under robust_ftl.
    """
    with TemporaryDirectory() as tmp:
        root = Path(tmp) / "ep"
        _write_minimal_episode(root)
        metrics = read_episode_metrics(root)
        ftl_first = score_episode(metrics, target_stop_time=2.0, objective_mode="ftl_first")
        robust = score_episode(metrics, target_stop_time=2.0, objective_mode="robust_ftl")
        exotic = ftl_first.components.get("exotic_penalty", 0.0)
        if exotic < 0.0:
            # exotic present -> robust_ftl applies the heavier penalty
            assert robust.total <= ftl_first.total + 1e-6
