"""Tests for artifact cleanup tiers."""

from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.core.artifact_cleanup import cleanup_episode_artifacts


def _episode_with_score(tmp_path: Path) -> Path:
    episode = tmp_path / "eval_000001"
    episode.mkdir()
    plot_dir = episode / "RadialRecipePlt00000"
    plot_dir.mkdir()
    (plot_dir / "Header").write_text("plot")
    (episode / "initial_data.gridinit").write_bytes(b"grid")
    (episode / "score.json").write_text(
        json.dumps({"score": {"general_ftl_solved": {"f_op": 1.0}}}),
        encoding="utf-8",
    )
    return episode


def test_cleanup_removes_plotfiles_immediately_after_score(tmp_path: Path) -> None:
    episode = _episode_with_score(tmp_path)
    report = cleanup_episode_artifacts(episode, tier="plotfiles_only")
    assert not any(episode.glob("RadialRecipePlt*"))
    assert (episode / "initial_data.gridinit").is_file()
    assert (episode / "score.json").is_file()
    assert report.removed_paths


def test_gridinit_not_removed_mid_campaign(tmp_path: Path) -> None:
    episode = _episode_with_score(tmp_path)
    cleanup_episode_artifacts(episode, tier="plotfiles_only")
    assert (episode / "initial_data.gridinit").is_file()


def test_gridinit_removed_on_final_prune_non_hq(tmp_path: Path) -> None:
    episode = _episode_with_score(tmp_path)
    report = cleanup_episode_artifacts(episode, tier="full_non_hq", keep_for_hq=False)
    assert not (episode / "initial_data.gridinit").exists()
    assert (episode / "score.json").is_file()
    assert report.removed_paths


def test_cleanup_warns_when_gridinit_had_solved_ftl_in_score(tmp_path: Path) -> None:
    episode = _episode_with_score(tmp_path)
    report = cleanup_episode_artifacts(episode, tier="full_non_hq", keep_for_hq=False)
    assert any("general_ftl_solved" in warning for warning in report.warnings)


def test_cleanup_never_touches_score_json(tmp_path: Path) -> None:
    episode = _episode_with_score(tmp_path)
    cleanup_episode_artifacts(episode, tier="full_non_hq", keep_for_hq=False)
    assert (episode / "score.json").is_file()
