"""Tests for plot consumer splash/FTL flag isolation."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from grteclyn_wrapper.core.config import RADIAL_RECIPE_DIR
from grteclyn_wrapper.core.episode import create_episode
from grteclyn_wrapper.core.params import write_params
from grteclyn_wrapper.core.plot_consumer import build_consume_command


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
