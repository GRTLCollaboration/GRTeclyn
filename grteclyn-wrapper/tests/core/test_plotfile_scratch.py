"""Plotfiles belong on node-local disk; results belong on the shared one.

A plotfile is written once, read once, and deleted.  Routing that traffic over
NFS is what held concurrency at two runs: consumers blocked in I/O wait and
backlogs grew.  The promotion queues moved their own transients to ``/tmp`` by
rewriting ``amr.plot_file`` in their params; these tests cover the same policy
applied where QD and CMA-ES pass through, so no launcher has to remember it.

The failure this guards against is quiet in both directions: send too little to
scratch and NFS stalls again; send too much and the scoring step finds no
plotfiles and reports a healthy run that measured nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from grteclyn_wrapper.core import scratch
from grteclyn_wrapper.core.artifact_cleanup import cleanup_episode_artifacts
from grteclyn_wrapper.core.episode import Episode
from grteclyn_wrapper.core.params import episode_path_overrides
from grteclyn_wrapper.core.plot_consumer import build_consume_command
from grteclyn_wrapper.core.scratch import (
    cache_env,
    plotfile_dir,
    purge_orphan_scratch,
    purge_plotfile_scratch,
    scratch_root,
)
from grteclyn_wrapper.metrics.probes.ftl.general import find_latest_plotfile
from grteclyn_wrapper.search.qd_search.io import _prune_eval_dirs


@pytest.fixture(autouse=True)
def _fresh_scratch_resolution():
    """The root is memoized per env value; tests change it between cases."""
    scratch._root_cache.clear()
    yield
    scratch._root_cache.clear()


def _episode(tmp_path: Path, name: str = "eval_000001") -> Path:
    episode = tmp_path / "runs" / "campaign_v1" / name
    (episode / "small_data").mkdir(parents=True)
    (episode / "data").mkdir()
    return episode


def _make_plotfile(parent: Path, name: str) -> Path:
    plotfile = parent / name
    (plotfile / "Level_0").mkdir(parents=True)
    (plotfile / "Header").write_text("HyperCLaw-V1.1\n", encoding="utf-8")
    (plotfile / "Level_0" / "Cell_D_00000").write_text("x", encoding="utf-8")
    return plotfile


def test_transients_go_to_scratch_and_results_stay_put(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)

    overrides = episode_path_overrides(episode, example="RadialRecipe")

    assert overrides["output_path"] == episode.resolve()
    for key in ("amr.plot_file", "amr.check_file"):
        assert (tmp_path / "scratch") in Path(overrides[key]).parents
        assert episode.resolve() not in Path(overrides[key]).parents


def test_the_scratch_directory_is_created_before_the_run(tmp_path, monkeypatch):
    """AMReX will not make the parent for us, and the run dies at first write."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)

    overrides = episode_path_overrides(episode, example="RadialRecipe")

    assert Path(overrides["amr.plot_file"]).parent.is_dir()


def test_turning_scratch_off_restores_the_old_layout(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", "0")
    episode = _episode(tmp_path)

    overrides = episode_path_overrides(episode, example="RadialRecipe")

    assert Path(overrides["amr.plot_file"]).parent == episode.resolve()
    assert plotfile_dir(episode) == episode.resolve()


def test_scratch_is_on_by_default(tmp_path, monkeypatch):
    """Nothing in a launcher should have to ask for this."""
    monkeypatch.delenv("GRTECLYN_SCRATCH", raising=False)
    monkeypatch.setattr(scratch, "DEFAULT_SCRATCH_ROOT", tmp_path / "default")

    assert scratch_root() == tmp_path / "default"


def test_two_campaigns_with_the_same_eval_number_do_not_collide(tmp_path, monkeypatch):
    """Every campaign has an ``eval_000001``; they run at the same time."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    first = tmp_path / "runs" / "qd_alpha" / "eval_000001"
    second = tmp_path / "runs" / "qd_beta" / "eval_000001"
    first.mkdir(parents=True)
    second.mkdir(parents=True)

    assert plotfile_dir(first) != plotfile_dir(second)


def test_an_unusable_scratch_root_falls_back_rather_than_failing(tmp_path, monkeypatch):
    """A missing directory is not worth killing a campaign over."""
    unwritable = tmp_path / "not-a-dir"
    unwritable.write_text("this is a file", encoding="utf-8")
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(unwritable))
    episode = _episode(tmp_path)

    assert scratch_root() is None
    assert plotfile_dir(episode) == episode.resolve()


def test_the_consumer_reads_scratch_and_writes_the_episode(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode_dir = _episode(tmp_path)
    (episode_dir / "params.txt").write_text("stop_time = 16\n", encoding="utf-8")
    episode = Episode(path=episode_dir)

    command = build_consume_command(episode, profile="radial")

    data = command[command.index("--data") + 1]
    out = command[command.index("--out") + 1]
    stop_sim = command[command.index("--stop-sim-path") + 1]
    assert Path(data) == plotfile_dir(episode_dir)
    assert Path(out) == episode_dir / "small_data"
    # The runner watches the episode for this file; it must not follow the
    # plotfiles onto a disk the runner is not looking at.
    assert Path(stop_sim) == episode_dir / ".stop_sim"


def test_scoring_still_finds_a_plotfile_that_moved(tmp_path, monkeypatch):
    """Without this the scorer sees an empty episode and scores it as fine."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00000")
    expected = _make_plotfile(target, "RadialRecipePlt00024")

    found = find_latest_plotfile(episode, complete_only=False)

    assert found == expected


def test_purge_keeps_what_the_consumer_never_recorded(tmp_path, monkeypatch):
    """Deleting from scratch is final: an unextracted sample is not recoverable."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00000")
    kept = _make_plotfile(target, "RadialRecipePlt00024")
    (episode / "small_data" / "consume_state.json").write_text(
        json.dumps({"RadialRecipePlt00000": True}), encoding="utf-8"
    )

    removed, warning = purge_plotfile_scratch(episode)

    assert kept.is_dir()
    assert not (target / "RadialRecipePlt00000").exists()
    assert removed
    assert warning is not None and "RadialRecipePlt00024" in warning


def test_purge_reclaims_everything_once_the_ledger_is_complete(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00000")
    _make_plotfile(target, "RadialRecipePlt00024")
    (episode / "small_data" / "consume_state.json").write_text(
        json.dumps({"RadialRecipePlt00000": True, "RadialRecipePlt00024": True}),
        encoding="utf-8",
    )

    _removed, warning = purge_plotfile_scratch(episode)

    assert not target.exists()
    assert warning is None


def test_a_discarded_eval_takes_its_scratch_with_it(tmp_path, monkeypatch):
    """Nothing will come back for these, and the local disk is the small one."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00024")  # never extracted

    _removed, _warning = purge_plotfile_scratch(episode, force=True)

    assert not target.exists()


def test_pruning_an_eval_dir_reclaims_its_scratch(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    campaign = tmp_path / "runs" / "campaign_v1"
    episode = _episode(tmp_path)
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00024")
    records = [{"eval": 1, "score": 0.5, "episode": str(episode)}]

    deleted = _prune_eval_dirs(campaign, records, keep_eval_ids=set())

    assert deleted == 1
    assert not episode.exists()
    assert not target.exists()


def test_cleanup_after_scoring_reclaims_the_scratch(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    (episode / "score.json").write_text("{}", encoding="utf-8")
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00000")
    (episode / "small_data" / "consume_state.json").write_text(
        json.dumps({"RadialRecipePlt00000": True}), encoding="utf-8"
    )

    report = cleanup_episode_artifacts(episode, tier="plotfiles_only")

    assert not target.exists()
    assert report.removed_paths
    assert not report.warnings


def test_cleanup_reports_what_it_had_to_keep(tmp_path, monkeypatch):
    """Silence here would read as 'reclaimed' while the disk kept filling."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    episode = _episode(tmp_path)
    (episode / "score.json").write_text("{}", encoding="utf-8")
    target = plotfile_dir(episode, create=True)
    _make_plotfile(target, "RadialRecipePlt00024")

    report = cleanup_episode_artifacts(episode, tier="plotfiles_only")

    assert target.is_dir()
    assert any("unextracted" in warning for warning in report.warnings)


def test_a_killed_campaign_does_not_hold_the_disk_forever(tmp_path, monkeypatch):
    """No cleanup path runs when a campaign is killed; the next one pays."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    abandoned = _episode(tmp_path, "eval_000001")
    live = _episode(tmp_path, "eval_000002")
    dead_scratch = plotfile_dir(abandoned, create=True)
    live_scratch = plotfile_dir(live, create=True)
    _make_plotfile(dead_scratch, "RadialRecipePlt00024")
    _make_plotfile(live_scratch, "RadialRecipePlt00024")
    import shutil as _shutil

    _shutil.rmtree(abandoned)  # the campaign was killed and its dirs cleared

    removed = purge_orphan_scratch()

    assert removed == [str(dead_scratch)]
    assert live_scratch.is_dir()


def test_the_sweep_leaves_directories_it_does_not_own(tmp_path, monkeypatch):
    """``_cache`` and anything hand-made are not ours to delete."""
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    root = tmp_path / "scratch"
    (root / "_cache" / "mpl").mkdir(parents=True)
    (root / "someone_elses_work").mkdir(parents=True)

    assert purge_orphan_scratch() == []
    assert (root / "_cache" / "mpl").is_dir()
    assert (root / "someone_elses_work").is_dir()


def test_library_caches_are_pinned_into_scratch(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    for name in ("XDG_CACHE_HOME", "MPLCONFIGDIR", "TMPDIR", "PYTHONPYCACHEPREFIX"):
        monkeypatch.delenv(name, raising=False)

    env = cache_env()

    assert env
    for value in env.values():
        assert (tmp_path / "scratch") in Path(value).parents or Path(value) == (
            tmp_path / "scratch" / "_cache"
        )


def test_an_operator_pinned_cache_wins(tmp_path, monkeypatch):
    monkeypatch.setenv("GRTECLYN_SCRATCH", str(tmp_path / "scratch"))
    monkeypatch.setenv("MPLCONFIGDIR", "/somewhere/deliberate")

    assert "MPLCONFIGDIR" not in cache_env()
