"""One ray-bundle trust bar, applied identically by every consumer.

Captured rays (fell into a puncture throat / horizon) are physics, not
integration failures, so they leave the denominator:
``n_reached == n_rays - n_captured``.

The probes adopted that bar; the consumers of their reports kept the older
``n_reached == n_rays``.  A run where any ray was captured was therefore
certified by the probe and simultaneously marked untrusted by the scorer, the
``f_geo_evol_ok`` column, the ``geo_trustworthy`` column, and the QD search
gate.  These tests pin every one of those to the shared predicate so the bars
cannot drift apart again.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import (
    EvolvingGeodesicFtlReport,
    evolving_report_trustworthy,
)
from grteclyn_wrapper.metrics.probes.ftl.geodesic import rays_complete
from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport
from grteclyn_wrapper.search.ftl_peak_metrics import four_d_geodesic_trustworthy


def _report(**kw) -> EvolvingGeodesicFtlReport:
    base = dict(
        f_geo=0.22,
        f_geo_frozen_peak=0.15,
        t_emit=0.0,
        t_arrival=7.5,
        t_flat=10.0,
        n_rays=5,
        n_reached=5,
        max_h_drift=1.0e-8,
        h_quality_ok=True,
        max_h_rel_drift=1.0e-9,
        notes=("test",),
        n_captured=0,
    )
    base.update(kw)
    return EvolvingGeodesicFtlReport(**base)


# --------------------------------------------------------------- the predicate


def test_reduces_to_old_bar_when_nothing_captured() -> None:
    assert rays_complete(5, 5, 0)
    assert not rays_complete(5, 4, 0)


def test_captured_rays_leave_the_denominator() -> None:
    # 2 of 5 fell into the throat; the remaining 3 all reached.  Certified.
    assert rays_complete(5, 3, 2)
    # 2 captured but only 2 of the remaining 3 reached — one genuinely lost.
    assert not rays_complete(5, 2, 2)


def test_degenerate_bundles_are_not_certified() -> None:
    assert not rays_complete(0, 0, 0)      # no rays fired
    assert not rays_complete(5, 0, 5)      # every ray captured, nothing left


# ------------------------------------------------------- consumers agree with it


def test_probe_predicate_matches_shared_bar() -> None:
    assert evolving_report_trustworthy(_report(n_reached=3, n_captured=2))
    assert not evolving_report_trustworthy(_report(n_reached=2, n_captured=2))
    # h-quality still vetoes independently of the bundle count.
    assert not evolving_report_trustworthy(
        _report(n_reached=3, n_captured=2, h_quality_ok=False)
    )


def test_qd_search_gate_matches_shared_bar() -> None:
    evo = {"h_quality_ok": True, "n_rays": 5, "n_reached": 3, "n_captured": 2}
    assert four_d_geodesic_trustworthy(evo)
    assert not four_d_geodesic_trustworthy({**evo, "n_reached": 2})
    # Payloads written before n_captured existed must keep the old meaning.
    assert four_d_geodesic_trustworthy(
        {"h_quality_ok": True, "n_rays": 5, "n_reached": 5}
    )


def test_scorer_credits_a_captured_ray_bundle(monkeypatch) -> None:
    from grteclyn_wrapper.metrics.score import score_episode
    from grteclyn_wrapper.metrics.types.diagnostics import (
        EvolvingGeodesicMetrics,
        FtlTimeSeriesMetrics,
    )
    from grteclyn_wrapper.metrics.types.episode import EpisodeMetrics

    # Search-tier: a leaked RL_PUMP_STOP_TIME would reject t_emit=0 for its own
    # reason and mask what this test is checking.
    monkeypatch.delenv("RL_PUMP_STOP_TIME", raising=False)
    metrics = EpisodeMetrics(
        collapse=None,
        constraints=None,
        stability=None,
        comoving=None,
        ftl=None,
        termination_reason="completed",
        ftl_timeseries=FtlTimeSeriesMetrics(
            n_frames=3,
            t=(0.0, 8.0, 16.0),
            f_op=(0.0, 0.1, 0.0),
            f_geo=(0.0, 0.05, 0.0),
            geo_trustworthy=(True, True, True),
            max_local_speed=(1.0, 1.2, 1.0),
            superluminal_fraction=(0.0, 0.5, 0.0),
            structure_coherence=(float("nan"),) * 3,
            max_h_rel_drift=(0.0, 0.001, 0.0),
            f_geo_peak=0.05,
            t_at_f_geo_peak=8.0,
            f_op_peak=0.1,
            t_at_f_op_peak=8.0,
            max_local_speed_peak=1.2,
            t_at_max_speed=8.0,
            superluminal_fraction_peak=0.5,
            t_at_superluminal_peak=8.0,
            ftl_lifetime_fraction=0.33,
            op_lifetime_fraction=0.33,
        ),
        evolving_geodesic=EvolvingGeodesicMetrics(
            f_geo=0.04,
            f_geo_frozen_peak=0.05,
            t_emit=0.0,
            t_arrival=12.0,
            t_flat=14.0,
            n_rays=5,
            n_reached=3,
            h_quality_ok=True,
            max_h_rel_drift=0.001,
            n_captured=2,
        ),
    )
    score = score_episode(metrics, objective_mode="ftl_first")
    # Under the stale bar this was 0.0 — a certified shortcut earning nothing.
    assert score.components["ftl_geo_evolving"] > 0.0


# ----------------------------------------------------- the f_geo_evol_ok column


def _minimal_episode(tmp_path: Path) -> Path:
    episode = tmp_path / "eval"
    (episode / "data").mkdir(parents=True)
    (episode / "small_data").mkdir(parents=True)
    (episode / "data" / "collapse_diagnostics.dat").write_text(
        "0.0 1.0 1.0 0.0 0 0 0 0.0 0.1\n", encoding="utf-8"
    )
    (episode / "data" / "constraint_norms.dat").write_text(
        "0.0 0.001 0.0001 -0.001 0.01 0.0\n", encoding="utf-8"
    )
    (episode / "small_data" / "ftl_timeseries.dat").write_text(
        "# time  f_op  f_geo  geo_trustworthy  max_local_speed  "
        "superluminal_fraction  max_shift  structure_coherence  reachable  "
        "n_rays  n_reached  max_h_rel_drift\n"
        "8.0 0.05 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001\n",
        encoding="utf-8",
    )
    return episode


def test_collector_marks_a_captured_ray_bundle_ok(tmp_path: Path) -> None:
    episode = _minimal_episode(tmp_path)
    fake_report = _report(n_reached=3, n_captured=2)
    evolved = GeneralFtlReport(
        f_op=0.05,
        t_min=9.0,
        t_flat=10.0,
        reachable=True,
        path_offaxis=False,
        notes=(),
        max_local_speed=1.2,
        superluminal_fraction=0.1,
        max_shift=0.05,
        structure_coherence=None,
    )
    c = "grteclyn_wrapper.metrics.aggregation.collector"
    with (
        patch(f"{c}.wait_for_plotfile_complete"),
        patch(f"{c}.matter_coherence_from_plotfile", return_value=None),
        patch(f"{c}.find_latest_plotfile", return_value=episode / "plotfile_000"),
        patch(f"{c}.find_recent_plotfiles", return_value=[]),
        patch(f"{c}.compute_general_ftl_from_plotfile", return_value=evolved),
        patch(f"{c}.compute_geodesic_ftl_from_plotfile", return_value=None),
        patch(
            f"{c}.compute_evolving_geodesic_ftl_from_metric_stack_cache",
            return_value=fake_report,
        ),
        patch(f"{c}.slice_count", return_value=3),
    ):
        from grteclyn_wrapper.metrics.aggregation.collector import read_episode_metrics

        metrics = read_episode_metrics(episode, ftl_L=8.0, evolving_geodesic=True)

    assert metrics.evolving_geodesic is not None
    assert metrics.evolving_geodesic.n_captured == 2

    row = [
        line
        for line in (episode / "small_data" / "ftl_timeseries.dat")
        .read_text(encoding="utf-8")
        .splitlines()
        if line and not line.startswith("#")
    ][-1]
    # Last column is f_geo_evol_ok.  Stale bar wrote "0" here.
    assert row.split()[-1] == "1"

    payload = json.loads(
        (episode / "small_data" / "evolving_geodesic.json").read_text(encoding="utf-8")
    )
    assert payload["n_captured"] == 2
