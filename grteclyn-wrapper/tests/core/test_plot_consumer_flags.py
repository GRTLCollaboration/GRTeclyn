"""Tests for plot consumer splash/FTL flag isolation."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from grteclyn_wrapper.core.config import RADIAL_RECIPE_DIR
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.core.params import write_params
from grteclyn_wrapper.core.plot_consumer import (
    build_consume_command,
    consumer_delete_plotfiles_enabled,
    consumer_radii_from_env,
    post_run_frames_enabled,
)


@pytest.fixture()
def episode(tmp_path: Path):
    ep = create_episode(tmp_path, name="consumer-flags")
    template = RADIAL_RECIPE_DIR / "params.txt"
    if not template.is_file():
        pytest.skip("RadialRecipe template unavailable")
    write_params(template, ep.params_path, episode_dir=ep.path, overrides={"stop_time": 16.0})
    return ep


def test_ftl_consumer_has_no_central_by_default(episode, monkeypatch) -> None:
    monkeypatch.delenv("GRTECLYN_CENTRAL_TIMESERIES", raising=False)
    monkeypatch.delenv("GRTECLYN_CONSUMER_DRAIN_MINIMAL", raising=False)
    command = build_consume_command(
        episode,
        ftl_timeseries=True,
        incremental_score=True,
        objective_mode="ftl_first",
    )
    assert "--ftl-timeseries" in command
    assert "--central-timeseries" not in command


def test_splash_consumer_wires_central_and_incremental(episode, monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_CENTRAL_TIMESERIES", "1")
    monkeypatch.setenv("GRTECLYN_CENTRAL_BALL", "1")
    monkeypatch.setenv("GRTECLYN_CENTRAL_RADIAL", "1")
    monkeypatch.setenv("GRTECLYN_INCREMENTAL_SCORE", "1")
    monkeypatch.setenv("GRTECLYN_SPLASH_EARLY_TERM", "1")
    monkeypatch.delenv("GRTECLYN_CONSUMER_DRAIN_MINIMAL", raising=False)
    command = build_consume_command(
        episode,
        ftl_timeseries=False,
        incremental_score=True,
        objective_mode="critical_collapse",
    )
    assert "--central-timeseries" in command
    assert "--central-ball" in command
    assert "--central-radial-profile" in command
    assert "--incremental-score" in command
    assert "--splash-early-term" in command
    assert "--objective-mode" in command
    assert "critical_collapse" in command


def test_radial_psi4_disabled_by_default(episode, monkeypatch) -> None:
    monkeypatch.delenv("GRTECLYN_PSI4", raising=False)
    command = build_consume_command(episode, profile="radial")
    assert "--no-psi4" in command
    assert "--psi4" not in command


def test_radial_psi4_enabled_when_env_set(episode, monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_PSI4", "1")
    command = build_consume_command(episode, profile="radial")
    assert "--psi4" in command
    assert "--no-psi4" not in command
    assert "--frames-corner" not in command


def test_psi4_n_points_bump_when_enabled(episode, monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_PSI4", "1")
    monkeypatch.setenv("GRTECLYN_PSI4_N_POINTS", "96")
    command = build_consume_command(episode, profile="radial", n_points=64)
    idx = command.index("--n-points")
    assert command[idx + 1] == "96"


def test_consumer_radii_from_env(monkeypatch) -> None:
    monkeypatch.setenv("CONSUMER_RADII", "8 12 24")
    assert consumer_radii_from_env() == (8.0, 12.0, 24.0)


def test_consumer_radii_default_when_env_empty(monkeypatch) -> None:
    monkeypatch.delenv("CONSUMER_RADII", raising=False)
    assert consumer_radii_from_env() == (8.0, 12.0, 24.0)


def test_metrics_only_consumer_does_not_delete_plotfiles(episode, monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_FRAMES", "0")
    command = build_consume_command(episode, profile="radial", frames=False, delete=True)
    assert "--delete" not in command
    assert "--frames-fields" not in command


def test_metrics_only_consumer_can_delete_when_forced(episode, monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_FRAMES", "0")
    monkeypatch.setenv("GRTECLYN_DELETE_WITHOUT_FRAMES", "1")
    command = build_consume_command(episode, profile="radial", frames=False, delete=True)
    assert "--delete" in command


def test_post_run_frames_default_on(monkeypatch) -> None:
    monkeypatch.delenv("GRTECLYN_CONSUMER_DRAIN", raising=False)
    assert post_run_frames_enabled() is True


def test_post_run_frames_off_for_fast_drain(monkeypatch) -> None:
    monkeypatch.setenv("GRTECLYN_CONSUMER_DRAIN", "1")
    assert post_run_frames_enabled() is False


def test_consumer_delete_plotfiles_enabled_matrix() -> None:
    assert consumer_delete_plotfiles_enabled(frames=True, delete_requested=True) is True
    assert consumer_delete_plotfiles_enabled(frames=False, delete_requested=True) is False
    assert consumer_delete_plotfiles_enabled(frames=False, delete_requested=False) is False
