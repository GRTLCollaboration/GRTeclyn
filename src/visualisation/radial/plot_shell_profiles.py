from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _parse_shell_profiles(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"Empty shell profile file: {path}")
    header = lines[0].lstrip("#").split()
    if header[0] != "time":
        raise ValueError(f"Unexpected shell profile header in {path}")
    cols = {name: [] for name in header[1:]}
    times: list[float] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(header):
            continue
        times.append(float(parts[0]))
        for name, value in zip(header[1:], parts[1:]):
            cols[name].append(float(value))
    return np.asarray(times), {k: np.asarray(v) for k, v in cols.items()}


def plot_shell_profiles(path: Path, out_dir: Path, name: str = "shell_profiles") -> Path:
    t, cols = _parse_shell_profiles(path)
    out_dir.mkdir(parents=True, exist_ok=True)

    fields = sorted({key.split("_mean_")[0] for key in cols if "_mean_" in key})
    fig, axes = plt.subplots(len(fields), 1, figsize=(10, 3 * len(fields)), sharex=True)
    if len(fields) == 1:
        axes = [axes]

    for ax, field in zip(axes, fields):
        mean_keys = sorted(k for k in cols if k.startswith(f"{field}_mean_"))
        for key in mean_keys:
            radius = key.split("_mean_")[1]
            min_key = f"{field}_min_{radius}"
            max_key = f"{field}_max_{radius}"
            ax.plot(t, cols[key], label=f"{field} mean @ {radius}")
            if min_key in cols and max_key in cols:
                ax.fill_between(t, cols[min_key], cols[max_key], alpha=0.15)
        ax.set_ylabel(field)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)

    axes[-1].set_xlabel("time")
    fig.tight_layout()
    out_path = out_dir / f"{name}.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot shell_profiles.dat time series.")
    parser.add_argument("shell_profiles", type=Path, help="Path to shell_profiles.dat")
    parser.add_argument("--out", type=Path, required=True, help="Output directory")
    parser.add_argument("--name", default="shell_profiles", help="Output basename")
    args = parser.parse_args()
    out = plot_shell_profiles(args.shell_profiles, args.out, name=args.name)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
