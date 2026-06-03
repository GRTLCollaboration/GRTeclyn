#!/usr/bin/env python3
"""Plot rotating wormhole support-scan outcomes."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _omega_from_case(case: str, manifest: dict) -> float:
    cfg = manifest.get("config", {})
    if "lump_omega" in cfg:
        return float(cfg["lump_omega"])
    match = re.search(r"omega_([pm])([0-9p]+)", case)
    if not match:
        return float("nan")
    sign = 1.0 if match.group(1) == "p" else -1.0
    return sign * float(match.group(2).replace("p", "."))


def _load_rows(summary: Path) -> list[dict]:
    rows = []
    with summary.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))
    return rows


def _diagnostic_path(episode: Path) -> Path | None:
    for candidate in (
        episode / "data" / "collapse_diagnostics.dat",
        episode / "collapse_diagnostics.dat",
    ):
        if candidate.exists():
            return candidate
    return None


def _load_diagnostic(path: Path) -> np.ndarray | None:
    rows: list[list[float]] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            try:
                rows.append([float(part) for part in stripped.split()])
            except ValueError:
                continue
    if not rows:
        return None
    return np.asarray(rows, dtype=float)


def _classify_episode(episode: Path) -> tuple[str, float, float]:
    diag_path = _diagnostic_path(episode)
    if diag_path is None:
        return "missing", float("nan"), float("nan")
    data = _load_diagnostic(diag_path)
    if data is None:
        return "missing", float("nan"), float("nan")

    min_lapse = float(np.nanmin(data[:, 1]))
    max_ah_r = float(np.nanmax(data[:, 7])) if data.shape[1] >= 8 else 0.0
    if max_ah_r > 0.0 or min_lapse < 0.2:
        return "collapse", min_lapse, max_ah_r
    return "open", min_lapse, max_ah_r


def build_points(rows: list[dict]) -> list[dict]:
    points: list[dict] = []
    for row in rows:
        episode = Path(row["episode"])
        outcome, min_lapse, max_ah_r = _classify_episode(episode)
        points.append({
            "case": row["case"],
            "omega": _omega_from_case(row["case"], row.get("manifest", {})),
            "support": float(row["support_strength"]),
            "outcome": outcome,
            "min_lapse": min_lapse,
            "max_ah_r": max_ah_r,
        })
    return points


def plot(points: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.0, 4.5))

    markers = {"collapse": "o", "open": "s", "missing": "x"}
    colors = {"collapse": "tab:red", "open": "tab:blue", "missing": "tab:gray"}
    labels_seen: set[str] = set()
    for point in points:
        label = point["outcome"]
        ax.scatter(
            point["omega"],
            point["support"],
            marker=markers[label],
            color=colors[label],
            s=70,
            label=label if label not in labels_seen else None,
        )
        labels_seen.add(label)

    ax.set_xlabel("GRTresna rotation rate omega")
    ax.set_ylabel("Runtime phantom support strength")
    ax.set_title("Rotating wormhole collapse threshold scan")
    ax.grid(True, alpha=0.3)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=200)


def write_points(points: list[dict], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        "\n".join(json.dumps(point, sort_keys=True) for point in points) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path, help="scan_rotating_wormhole_support.py JSONL summary.")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("src/visualisation/process_wave/rotating_wormhole_support_scan.png"),
    )
    parser.add_argument(
        "--points-jsonl",
        type=Path,
        default=Path("src/visualisation/process_wave/rotating_wormhole_support_points.jsonl"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    rows = _load_rows(args.summary)
    points = build_points(rows)
    write_points(points, args.points_jsonl)
    plot(points, args.output)
    print(f"Wrote {args.output}")
    print(f"Wrote {args.points_jsonl}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
