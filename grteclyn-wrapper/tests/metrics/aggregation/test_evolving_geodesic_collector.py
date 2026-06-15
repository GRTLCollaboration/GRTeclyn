"""Collector hook for opt-in 4D evolving null-geodesic scoring."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

from grteclyn_wrapper.metrics.aggregation.collector import read_episode_metrics
from grteclyn_wrapper.metrics.probes.ftl.evolving_geodesic import EvolvingGeodesicFtlReport
from grteclyn_wrapper.metrics.probes.ftl.general import GeneralFtlReport


def _minimal_episode(tmp_path: Path) -> Path:
    episode = tmp_path / "eval"
    data = episode / "data"
    small = episode / "small_data"
    data.mkdir(parents=True)
    small.mkdir(parents=True)
    (data / "collapse_diagnostics.dat").write_text(
        "0.0 1.0 1.0 0.0 0 0 0 0.0 0.1\n", encoding="utf-8"
    )
    (data / "constraint_norms.dat").write_text(
        "0.0 0.001 0.0001 -0.001 0.01 0.0\n", encoding="utf-8"
    )
    (small / "ftl_timeseries.dat").write_text(
        "# time  f_op  f_geo  geo_trustworthy  max_local_speed  "
        "superluminal_fraction  max_shift  structure_coherence  reachable  "
        "n_rays  n_reached  max_h_rel_drift\n"
        "8.0 0.05 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001\n",
        encoding="utf-8",
    )
    return episode


def test_collector_writes_evolving_geodesic_when_flag_on(tmp_path: Path) -> None:
    episode = _minimal_episode(tmp_path)
    fake_report = EvolvingGeodesicFtlReport(
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
    )
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

    with (
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.wait_for_plotfile_complete",
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.matter_coherence_from_plotfile",
            return_value=None,
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.find_latest_plotfile",
            return_value=episode / "plotfile_000",
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.find_recent_plotfiles",
            return_value=[episode / f"plot_{i}" for i in range(3)],
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.compute_general_ftl_from_plotfile",
            return_value=evolved,
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.compute_geodesic_ftl_from_plotfile",
            return_value=None,
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.compute_evolving_geodesic_ftl_from_metric_stack_cache",
            return_value=fake_report,
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.slice_count",
            return_value=3,
        ),
        patch(
            "grteclyn_wrapper.metrics.aggregation.collector.compute_evolving_geodesic_ftl_from_plotfiles",
        ) as plotfiles_mock,
    ):
        metrics = read_episode_metrics(episode, ftl_L=8.0, evolving_geodesic=True)

    plotfiles_mock.assert_not_called()

    assert metrics.evolving_geodesic is not None
    assert metrics.evolving_geodesic.f_geo == 0.22
    json_path = episode / "small_data" / "evolving_geodesic.json"
    assert json_path.is_file()
    payload = json.loads(json_path.read_text(encoding="utf-8"))
    assert payload["f_geo"] == 0.22

    ts_text = (episode / "small_data" / "ftl_timeseries.dat").read_text(encoding="utf-8")
    assert "f_geo_evol" in ts_text
    assert "2.2000000000000000e-01" in ts_text


def test_collector_skips_evolving_geodesic_when_flag_off(tmp_path: Path) -> None:
    episode = _minimal_episode(tmp_path)
    with patch(
        "grteclyn_wrapper.metrics.aggregation.collector.compute_evolving_geodesic_ftl_from_plotfiles",
    ) as mocked:
        read_episode_metrics(episode, evolving_geodesic=False)
    mocked.assert_not_called()
