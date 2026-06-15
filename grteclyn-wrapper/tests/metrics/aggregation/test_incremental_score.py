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
