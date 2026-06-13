"""JSONL/CSV persistence for atlas records."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping

CSV_FIELDS = [
    "index",
    "episode",
    "exit_code",
    "score_total",
    "labels",
    "final_time",
    "target_stop_time",
    "min_lapse",
    "min_chi",
    "max_abs_k",
    "max_horizon_radius",
    "min_theta_plus",
    "max_hamiltonian_l2",
    "max_momentum_l2",
    "min_rho_required",
    "integral_negative_rho",
    "overrides_json",
]


def append_jsonl(path: Path, record: Mapping[str, Any]) -> None:
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def flatten_record(record: Mapping[str, Any]) -> dict[str, Any]:
    metrics = record.get("metrics") or {}
    collapse = metrics.get("collapse") or {}
    constraints = metrics.get("constraints") or {}
    score = record.get("score") or {}
    return {
        "index": record.get("index"),
        "episode": record.get("episode"),
        "exit_code": record.get("exit_code"),
        "score_total": score.get("total"),
        "labels": ",".join(record.get("labels") or []),
        "final_time": record.get("final_time"),
        "target_stop_time": record.get("target_stop_time"),
        "min_lapse": collapse.get("min_lapse"),
        "min_chi": collapse.get("min_chi"),
        "max_abs_k": collapse.get("max_abs_k"),
        "max_horizon_radius": collapse.get("max_horizon_radius"),
        "min_theta_plus": collapse.get("min_theta_plus"),
        "max_hamiltonian_l2": constraints.get("max_hamiltonian_l2"),
        "max_momentum_l2": constraints.get("max_momentum_l2"),
        "min_rho_required": constraints.get("min_rho_required"),
        "integral_negative_rho": constraints.get("integral_negative_rho"),
        "overrides_json": json.dumps(record.get("overrides") or {}, sort_keys=True),
    }


def append_csv(path: Path, record: Mapping[str, Any]) -> None:
    row = flatten_record(record)
    should_write_header = not path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        if should_write_header:
            writer.writeheader()
        writer.writerow(row)
