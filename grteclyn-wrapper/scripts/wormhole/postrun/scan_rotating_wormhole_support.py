#!/usr/bin/env python3
"""Sweep runtime phantom support for rotating GRTresna wormhole data."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

WRAPPER_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(WRAPPER_ROOT / "src"))
from grteclyn_wrapper.core.site_paths import grtresna_env, openmpi_root, sim_root

SIM_ROOT = sim_root()
LOCAL_OPENMPI = openmpi_root()
_GRTRESNA_ENV = grtresna_env()
GRTRESNA_ENV = _GRTRESNA_ENV if _GRTRESNA_ENV is not None else Path()


def _parse_float_list(raw: str) -> list[float]:
    values = [float(part) for part in raw.replace(",", " ").split()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one value")
    return values


def _case_dirs(paths: list[Path]) -> list[Path]:
    cases: list[Path] = []
    for path in paths:
        if (path / "initial_data.gridinit").exists():
            cases.append(path)
            continue
        cases.extend(
            child for child in sorted(path.iterdir())
            if child.is_dir() and (child / "initial_data.gridinit").exists()
        )
    return cases


def _load_manifest(case_dir: Path) -> dict:
    manifest_path = case_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    return json.loads(manifest_path.read_text(encoding="utf-8"))


def _episode_name(case_dir: Path, support: float) -> str:
    token = f"{support:.4g}".replace(".", "p")
    return f"{case_dir.name}_support_{token}"


def run_case(args: argparse.Namespace, case_dir: Path, support: float) -> dict:
    episode_name = _episode_name(case_dir, support)
    cmd = [
        sys.executable,
        "-m",
        "grteclyn_wrapper",
        "--example",
        "RotatingWormholeCollapse",
        "--runs-dir",
        str(args.runs_dir),
        "--mpi-ranks",
        str(args.mpi_ranks),
        "--cuda-devices",
        args.cuda_devices,
        "--name",
        episode_name,
        "--set",
        f"recipe_initial_data_file={case_dir / 'initial_data.gridinit'}",
        "--set",
        f"wormhole_support_strength={support}",
        "--set",
        "wormhole_phi_perturbation_amplitude=0.0",
        "--set",
        f"stop_time={args.stop_time}",
        "--set",
        f"plot_interval={args.plot_interval}",
        "reproduce",
    ]
    if args.dry_run:
        cmd.insert(-1, "--dry-run")
    if args.consume_plotfiles:
        cmd.insert(-1, "--consume-plotfiles")
        cmd.insert(-1, "--consumer-delete")
        cmd.insert(-1, "--consumer-radii")
        cmd.insert(-1, "12")
        cmd.insert(-1, "16")
        cmd.insert(-1, "20")
        cmd.insert(-1, "24")

    env = dict(os.environ)
    path_entries = [
        str(LOCAL_OPENMPI / "bin"),
        str(GRTRESNA_ENV / "bin"),
        env.get("PATH", ""),
    ]
    lib_entries = [
        str(LOCAL_OPENMPI / "lib"),
        str(GRTRESNA_ENV / "lib"),
        env.get("LD_LIBRARY_PATH", ""),
    ]
    env["PATH"] = os.pathsep.join(entry for entry in path_entries if entry)
    env["LD_LIBRARY_PATH"] = os.pathsep.join(entry for entry in lib_entries if entry)
    env["PYTHONPATH"] = str(WRAPPER_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(cmd, cwd=args.repo_root, text=True, env=env)
    episode_dir = args.runs_dir / episode_name
    return {
        "case": case_dir.name,
        "case_dir": str(case_dir),
        "episode": str(episode_dir),
        "support_strength": support,
        "returncode": completed.returncode,
        "dry_run": args.dry_run,
        "manifest": _load_manifest(case_dir),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run rotating-vs-non-rotating support scans for GRTresna gridinit cases."
    )
    parser.add_argument("cases", nargs="+", type=Path, help="Case dirs or parent dirs.")
    parser.add_argument(
        "--supports",
        type=_parse_float_list,
        default=[1.0, 0.8, 0.6, 0.5, 0.4, 0.3],
        help="Comma/space separated support_strength values.",
    )
    parser.add_argument("--runs-dir", type=Path, default=Path("runs/rotating_wormhole_support_scan"))
    parser.add_argument("--repo-root", type=Path, default=Path.cwd())
    parser.add_argument("--mpi-ranks", type=int, default=8)
    parser.add_argument("--cuda-devices", default="0,1,2,3,4,5,6,7")
    parser.add_argument("--stop-time", type=float, default=30.0)
    parser.add_argument("--plot-interval", type=int, default=1)
    parser.add_argument("--consume-plotfiles", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("runs/rotating_wormhole_support_scan/summary.jsonl"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.repo_root = args.repo_root.expanduser().resolve()
    args.runs_dir = (args.repo_root / args.runs_dir).resolve() if not args.runs_dir.is_absolute() else args.runs_dir
    args.summary = (args.repo_root / args.summary).resolve() if not args.summary.is_absolute() else args.summary
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    cases = _case_dirs([path.expanduser().resolve() for path in args.cases])
    if not cases:
        raise SystemExit("No case directories with initial_data.gridinit were found.")

    status = 0
    with args.summary.open("w", encoding="utf-8") as out:
        for case_dir in cases:
            for support in args.supports:
                record = run_case(args, case_dir, support)
                out.write(json.dumps(record, sort_keys=True) + "\n")
                out.flush()
                status = status or int(record["returncode"])
    print(f"Wrote scan summary: {args.summary}")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
