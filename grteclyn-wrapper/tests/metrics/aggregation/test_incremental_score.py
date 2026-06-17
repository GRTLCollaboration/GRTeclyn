"""Tests for prefix incremental scoring."""

from __future__ import annotations

import json
from pathlib import Path

from grteclyn_wrapper.metrics.aggregation.incremental import IncrementalScoreWriter
from grteclyn_wrapper.visualisation.process_wave.consume_plotfiles.extraction.ftl import (
    FTL_TIMESERIES_HEADER as HEADER,
)


def _write_collapse(path: Path, rows: list[tuple[float, float, float, float, float, float, float, float, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            # time lapse chi K ... horizon_r ... theta+
            handle.write(
                f"{row[0]:.6e} {row[1]:.6e} {row[2]:.6e} {row[3]:.6e} "
                f"0 0 0 {row[7]:.6e} {row[8]:.6e}\n"
            )


def _write_constraints(path: Path, rows: list[tuple[float, float, float, float, float, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(
                f"{row[0]:.6e} {row[1]:.6e} {row[2]:.6e} {row[3]:.6e} {row[4]:.6e} {row[5]:.6e}\n"
            )


def test_incremental_score_writer_appends_jsonl(tmp_path: Path) -> None:
    episode = tmp_path / "eval"
    small = episode / "small_data"
    data = episode / "data"
    small.mkdir(parents=True)
    data.mkdir(parents=True)

    ftl_path = small / "ftl_timeseries.dat"
    ftl_path.write_text(HEADER + "\n", encoding="utf-8")
    row1 = "8.0 0.05 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001"
    row2 = "16.0 0.04 0.03 1 1.08 0.08 0.01 nan 1 5 5 0.001"
    with ftl_path.open("a", encoding="utf-8") as handle:
        handle.write(row1 + "\n")
        handle.write(row2 + "\n")

    _write_collapse(
        data / "collapse_diagnostics.dat",
        [
            (0.0, 1.0, 1.0, 0.0, 0, 0, 0, 0.0, 0.1),
            (8.0, 0.9, 0.8, 0.1, 0, 0, 0, 1.0, 0.05),
            (16.0, 0.85, 0.75, 0.1, 0, 0, 0, 1.0, 0.02),
        ],
    )
    _write_constraints(
        data / "constraint_norms.dat",
        [
            (0.0, 0.001, 0.0001, -0.001, 0.01, 0.0),
            (8.0, 0.002, 0.0002, -0.001, 0.009, 0.0),
            (16.0, 0.003, 0.0003, -0.001, 0.008, 0.0),
        ],
    )

    writer = IncrementalScoreWriter(
        episode,
        objective_mode="ftl_first",
        target_stop_time=30.0,
        ftl_L=8.0,
    )
    first = writer.append(8.0)
    second = writer.append(16.0)
    assert first is not None
    assert second is not None
    assert second["score"] != first["score"]

    out = small / "score_timeseries.jsonl"
    lines = out.read_text(encoding="utf-8").strip().splitlines()
    assert len(lines) == 2
    payload = json.loads(lines[-1])
    assert payload["t"] == 16.0
    assert "score" in payload
    assert "components" in payload
    assert payload["f_geo_peak"] > 0.0


def test_incremental_4d_mode_withholds_frozen_ftl(tmp_path: Path) -> None:
    episode = tmp_path / "eval"
    small = episode / "small_data"
    data = episode / "data"
    small.mkdir(parents=True)
    data.mkdir(parents=True)

    ftl_path = small / "ftl_timeseries.dat"
    ftl_path.write_text(HEADER + "\n", encoding="utf-8")
    row = "8.0 0.08 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001"
    with ftl_path.open("a", encoding="utf-8") as handle:
        handle.write(row + "\n")

    _write_collapse(
        data / "collapse_diagnostics.dat",
        [(0.0, 1.0, 1.0, 0.0, 0, 0, 0, 0.0, 0.1), (8.0, 0.9, 0.8, 0.1, 0, 0, 0, 1.0, 0.05)],
    )
    _write_constraints(
        data / "constraint_norms.dat",
        [(0.0, 0.001, 0.0001, -0.001, 0.01, 0.0), (8.0, 0.002, 0.0002, -0.001, 0.009, 0.0)],
    )

    writer = IncrementalScoreWriter(
        episode,
        objective_mode="ftl_first",
        target_stop_time=30.0,
        ftl_L=8.0,
        evolving_geodesic_mode=True,
    )
    record = writer.append(8.0)
    assert record is not None
    assert record["operational_ftl_geodesic"] == 0.0
    assert record["ftl_geo_evolving"] == 0.0
    assert record["score"] < 500.0


def test_ftl_timeseries_parser_backward_compatible_12_and_14_columns(tmp_path: Path) -> None:
    from grteclyn_wrapper.metrics.diagnostics.ftl_timeseries import read_ftl_timeseries_metrics

    path12 = tmp_path / "ftl12.dat"
    path12.write_text(
        HEADER + "\n"
        "8.0 0.05 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001\n"
        "16.0 0.04 0.03 1 1.08 0.08 0.01 nan 1 5 5 0.001\n",
        encoding="utf-8",
    )
    ts12 = read_ftl_timeseries_metrics(path12)
    assert ts12 is not None
    assert ts12.n_frames == 2
    assert ts12.f_geo_evol is None
    assert ts12.f_geo_evol_ok is None

    path14 = tmp_path / "ftl14.dat"
    path14.write_text(
        HEADER + "  f_geo_evol  f_geo_evol_ok\n"
        "8.0 0.05 0.04 1 1.10 0.10 0.01 nan 1 5 5 0.001 0.0 0\n"
        "16.0 0.04 0.03 1 1.08 0.08 0.01 nan 1 5 5 0.001 0.18 1\n",
        encoding="utf-8",
    )
    ts14 = read_ftl_timeseries_metrics(path14)
    assert ts14 is not None
    assert ts14.f_geo_evol == 0.18
    assert ts14.f_geo_evol_ok is True
