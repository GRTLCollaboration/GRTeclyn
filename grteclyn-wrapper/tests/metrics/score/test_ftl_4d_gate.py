"""4D evolving geodesic authoritative scoring gates."""

from __future__ import annotations

from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
from grteclyn_wrapper.metrics.score import score_episode
from grteclyn_wrapper.metrics.types.diagnostics import (
    EvolvingGeodesicMetrics,
    FtlTimeSeriesMetrics,
)
from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics


def _frozen_timeseries(*, peak: float = 0.05) -> FtlTimeSeriesMetrics:
    return FtlTimeSeriesMetrics(
        n_frames=3,
        t=(0.0, 8.0, 16.0),
        f_op=(0.0, 0.1, 0.0),
        f_geo=(0.0, peak, 0.0),
        geo_trustworthy=(True, True, True),
        max_local_speed=(1.0, 1.2, 1.0),
        superluminal_fraction=(0.0, 0.5, 0.0),
        structure_coherence=(float("nan"), float("nan"), float("nan")),
        max_h_rel_drift=(0.0, 0.001, 0.0),
        f_geo_peak=peak,
        t_at_f_geo_peak=8.0,
        f_op_peak=0.1,
        t_at_f_op_peak=8.0,
        max_local_speed_peak=1.2,
        t_at_max_speed=8.0,
        superluminal_fraction_peak=0.5,
        t_at_superluminal_peak=8.0,
        ftl_lifetime_fraction=0.33,
        op_lifetime_fraction=0.33,
    )


def _base_metrics(**kwargs) -> EpisodeMetrics:
    return EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="completed",
        **kwargs,
    )


def test_4d_success_zeros_frozen_and_credits_evolving(monkeypatch) -> None:
    # Search-tier unit test: no igniter filter. A leaked RL_PUMP_STOP_TIME from
    # a paper/HQ shell would reject t_emit=0 and zero ftl_geo_evolving.
    monkeypatch.delenv("RL_PUMP_STOP_TIME", raising=False)
    metrics = _base_metrics(
        ftl_timeseries=_frozen_timeseries(peak=0.05),
        evolving_geodesic=EvolvingGeodesicMetrics(
            f_geo=0.04,
            f_geo_frozen_peak=0.05,
            t_emit=0.0,
            t_arrival=12.0,
            t_flat=14.0,
            n_rays=3,
            n_reached=3,
            h_quality_ok=True,
            max_h_rel_drift=0.001,
        ),
    )
    score = score_episode(metrics, objective_mode="ftl_first")
    assert score.components["ftl_geo_evolving"] > 0.0
    assert score.components["operational_ftl_geodesic"] == 0.0
    assert score.components["ftl_geo_timeavg"] == 0.0
    assert score.components["ftl_geo_peak"] == 0.0


def test_4d_failure_zeros_frozen_credit() -> None:
    metrics = _base_metrics(
        ftl_timeseries=_frozen_timeseries(peak=0.05),
        evolving_geodesic=EvolvingGeodesicMetrics(
            f_geo=0.0,
            f_geo_frozen_peak=0.05,
            t_emit=0.0,
            t_arrival=None,
            t_flat=14.0,
            n_rays=3,
            n_reached=0,
            h_quality_ok=False,
            max_h_rel_drift=3.0,
        ),
    )
    score = score_episode(metrics, objective_mode="ftl_first")
    assert score.components["ftl_geo_evolving"] == 0.0
    assert score.components["operational_ftl_geodesic"] == 0.0
    assert score.components["ftl_geo_timeavg"] == 0.0
    assert score.components["ftl_geo_peak"] == 0.0


def test_frozen_fallback_when_4d_not_run() -> None:
    metrics = _base_metrics(
        ftl_timeseries=_frozen_timeseries(peak=0.05),
        evolving_geodesic=None,
    )
    score = score_episode(metrics, objective_mode="ftl_first")
    assert score.components["operational_ftl_geodesic"] > 0.0
    assert score.components["ftl_geo_evolving"] == 0.0


def test_4d_mode_suppresses_frozen_when_trace_not_run() -> None:
    metrics = _base_metrics(
        ftl_timeseries=_frozen_timeseries(peak=0.05),
        evolving_geodesic=None,
        evolving_geodesic_mode=True,
        general_ftl_evolved=GeneralFtlReport(
            f_op=0.08,
            t_min=None,
            t_flat=1.0,
            max_local_speed=1.2,
            superluminal_fraction=0.1,
            path_offaxis=False,
            reachable=True,
            notes=(),
        ),
    )
    score = score_episode(metrics, objective_mode="ftl_first")
    assert score.components["operational_ftl_geodesic"] == 0.0
    assert score.components["ftl_geo_evolving"] == 0.0
    assert score.components["operational_ftl"] == 0.0
    assert score.components["ftl_geo_timeavg"] == 0.0
    assert score.components["ftl_geo_peak"] == 0.0
