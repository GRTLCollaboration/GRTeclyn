"""Integration tests for critical_collapse scoring."""

from __future__ import annotations

import os
from unittest.mock import patch

from grteclyn_wrapper.metrics.score.scorer import score_episode
from grteclyn_wrapper.metrics.types.central import CentralFieldMetrics
from grteclyn_wrapper.metrics.types.diagnostics import CollapseMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _episode(**kwargs) -> EpisodeMetrics:
    central = CentralFieldMetrics(
        n_frames=3,
        t=(0.0, 1.0, 2.0),
        rho_req=(1.0e-4, 8.0e-3, 3.0e-3),
        lapse=(1.0, 0.4, 0.5),
        scalar_activity=(0.01, 0.05, 0.02),
        peak_rho_req_at_origin=8.0e-3,
        peak_rho_req_time=1.0,
        initial_rho_req_at_origin=1.0e-4,
        min_lapse_at_origin=0.4,
        wave_chromaticity=0.5,
    )
    collapse = CollapseMetrics(
        final_time=2.0,
        min_lapse=0.4,
        min_chi=0.5,
        max_abs_k=0.1,
        max_horizon_radius=5.0,
        min_theta_plus=0.1,
        scalar_phi_range=0.2,
        scalar_pi_range=0.1,
        first_corroborated_time=5.0,
    )
    defaults = dict(
        collapse=collapse,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="completed_or_partial",
        central=central,
    )
    defaults.update(kwargs)
    return EpisodeMetrics(**defaults)


def test_score_episode_critical_collapse_end_to_end() -> None:
    score = score_episode(
        _episode(),
        target_stop_time=16.0,
        objective_mode="critical_collapse",
    )
    assert score.total > 0.0
    assert score.components["central_energy_peak"] > 0.0
    assert any("critical_collapse" in note for note in score.notes)


@patch.dict(os.environ, {"SPLASH_MODE": "threshold"}, clear=False)
def test_splash_mode_env_respected() -> None:
    score = score_episode(
        _episode(),
        target_stop_time=16.0,
        objective_mode="critical_collapse",
    )
    assert "central_lapse_collapse" in score.components
