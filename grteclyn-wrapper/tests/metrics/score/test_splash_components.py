"""Tests for splash score components."""

from __future__ import annotations

from grteclyn_wrapper.metrics.score.splash import compute_splash_components
from grteclyn_wrapper.metrics.score.types import ScoringContext
from grteclyn_wrapper.metrics.types.central import CentralFieldMetrics
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _central(**kwargs) -> CentralFieldMetrics:
    defaults = dict(
        n_frames=3,
        t=(0.0, 1.0, 2.0),
        rho_req=(1.0e-4, 5.0e-3, 2.0e-3),
        lapse=(1.0, 0.5, 0.4),
        scalar_activity=(0.01, 0.02, 0.03),
        peak_rho_req_at_origin=5.0e-3,
        peak_rho_req_time=1.0,
        initial_rho_req_at_origin=1.0e-4,
        min_lapse_at_origin=0.4,
        wave_chromaticity=0.6,
    )
    defaults.update(kwargs)
    return CentralFieldMetrics(**defaults)


def _ctx(central: CentralFieldMetrics | None) -> ScoringContext:
    return ScoringContext(
        metrics=EpisodeMetrics(
            collapse=None,
            constraints=None,
            stability=None,
            comoving=None,
            ftl=None,
            termination_reason="completed_or_partial",
            central=central,
        ),
        target_stop_time=16.0,
        domain_half_width=8.0,
        weights={},
    )


def test_zero_central_metrics_yields_zero_splash_components() -> None:
    ctx = _ctx(None)
    compute_splash_components(ctx)
    assert ctx.components["central_energy_peak"] == 0.0
    assert ctx.components["focusing_efficiency"] == 0.0
    assert ctx.components["pre_collapsed_penalty"] == 0.0


def test_central_energy_peak_normalization() -> None:
    ctx = _ctx(_central(peak_rho_req_at_origin=1.0e-2))
    compute_splash_components(ctx)
    assert ctx.components["central_energy_peak"] == 1.0


def test_focusing_efficiency_capped() -> None:
    ctx = _ctx(_central(initial_rho_req_at_origin=1.0e-6, peak_rho_req_at_origin=1.0))
    compute_splash_components(ctx)
    assert ctx.components["focusing_efficiency"] == 10.0


def test_pre_collapsed_penalty_when_initial_rho_high() -> None:
    ctx = _ctx(_central(initial_rho_req_at_origin=1.0e-2, min_lapse_at_origin=0.15))
    compute_splash_components(ctx, splash_mode="discovery")
    assert ctx.components["pre_collapsed_penalty"] < 0.0


def test_threshold_mode_lapse_band_reward() -> None:
    good = _ctx(_central(min_lapse_at_origin=0.03))
    compute_splash_components(good, splash_mode="threshold")
    assert good.components["central_lapse_collapse"] > 0.0

    bad = _ctx(_central(min_lapse_at_origin=0.005))
    compute_splash_components(bad, splash_mode="threshold")
    assert bad.components["central_lapse_collapse"] == -1.0
