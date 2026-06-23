#!/usr/bin/env python3
"""Compare Gate 0A baseline vs RL-neutral constraint / collapse trajectories."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _load_small_data(path: Path) -> tuple[list[float], list[list[float]]]:
    times: list[float] = []
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        times.append(float(parts[0]))
        rows.append([float(x) for x in parts[1:]])
    return times, rows


def _max_rel_err(
    ref_t: list[float],
    ref_v: list[float],
    other_t: list[float],
    other_v: list[float],
    *,
    atol: float = 1.0e-30,
) -> tuple[float, float, int]:
    other_map = {round(t, 12): v for t, v in zip(other_t, other_v)}
    max_err = 0.0
    worst_t = 0.0
    n = 0
    for t, rv in zip(ref_t, ref_v):
        key = round(t, 12)
        if key not in other_map:
            continue
        ov = other_map[key]
        denom = max(abs(rv), atol)
        err = abs(rv - ov) / denom
        if err > max_err:
            max_err = err
            worst_t = t
        n += 1
    return max_err, worst_t, n


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline_dir", type=Path)
    parser.add_argument("rl_dir", type=Path)
    parser.add_argument("--rtol", type=float, default=1.0e-10)
    args = parser.parse_args()

    base_data = args.baseline_dir / "data"
    rl_data = args.rl_dir / "data"
    checks: list[tuple[str, Path, int, int]] = [
        ("L2_Ham", base_data / "constraint_norms.dat", rl_data / "constraint_norms.dat", 0),
        ("min_lapse", base_data / "collapse_diagnostics.dat", rl_data / "collapse_diagnostics.dat", 0),
    ]

    failed = False
    for name, ref_path, other_path, col in checks:
        if not ref_path.exists() or not other_path.exists():
            print(f"FAIL missing file for {name}: {ref_path} / {other_path}")
            failed = True
            continue
        ref_t, ref_rows = _load_small_data(ref_path)
        other_t, other_rows = _load_small_data(other_path)
        ref_v = [row[col] for row in ref_rows]
        other_v = [row[col] for row in other_rows]
        max_err, worst_t, n = _max_rel_err(ref_t, ref_v, other_t, other_v)
        ok = max_err <= args.rtol
        status = "PASS" if ok else "FAIL"
        print(
            f"{status} {name}: max_rel_err={max_err:.3e} at t={worst_t:.6f} "
            f"(n={n}, rtol={args.rtol:.0e})"
        )
        if not ok:
            failed = True

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
