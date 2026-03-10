#!/usr/bin/env python3
"""
Plot streaming collapse diagnostics written by WormholeCollapse.

Expected input: collapse_diagnostics.dat (SmallDataIO ASCII), typically located at:
  <run_output>/data/collapse_diagnostics.dat

Columns (current format):
  time  min_lapse  min_chi  max_abs_K  min_lapse_x  min_lapse_y  min_lapse_z

This script produces a publication-style multi-panel figure (one panel per value).
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np


def _default_run_dir() -> Path:
    # Mirror other visualisation scripts: <repo>/src/visualisation/<module>/... -> <repo>/data by default
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent.parent  # .../GRTeclyn/src
    return (project_root.parent / "data").resolve()


def _resolve_input_path(data_dir: Path, explicit_input: str | None) -> Path:
    if explicit_input:
        p = Path(explicit_input).expanduser().resolve()
        if not p.exists():
            raise SystemExit(f"Input not found: {p}")
        return p

    # Common layouts:
    # - <run>/data/collapse_diagnostics.dat  (recommended)
    # - <run>/collapse_diagnostics.dat       (fallback)
    candidates = [
        data_dir / "data" / "collapse_diagnostics.dat",
        data_dir / "collapse_diagnostics.dat",
    ]
    for c in candidates:
        if c.exists():
            return c.resolve()
    raise SystemExit(
        "Could not find collapse_diagnostics.dat. Tried:\n"
        + "\n".join([f"  - {c}" for c in candidates])
        + "\nProvide it explicitly as a positional argument, or set --data to your run directory."
    )


def load_collapse_diagnostics(path: Path) -> Dict[str, np.ndarray]:
    # SmallDataIO files may or may not include a header; we accept both.
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            rows.append(vals)

    if not rows:
        raise SystemExit(f"No data rows found in {path}")

    arr = np.asarray(rows, dtype=float)
    if arr.shape[1] not in (4, 7):
        raise SystemExit(
            f"Unexpected number of columns in {path}: got {arr.shape[1]}, expected 4 or 7"
        )

    t = arr[:, 0]
    out: Dict[str, np.ndarray] = {
        "t": t,
        "min_lapse": arr[:, 1],
        "min_chi": arr[:, 2],
        "max_abs_K": arr[:, 3],
    }
    if arr.shape[1] == 7:
        out.update(
            {
                "min_lapse_x": arr[:, 4],
                "min_lapse_y": arr[:, 5],
                "min_lapse_z": arr[:, 6],
            }
        )
    else:
        out.update(
            {
                "min_lapse_x": np.full_like(t, np.nan),
                "min_lapse_y": np.full_like(t, np.nan),
                "min_lapse_z": np.full_like(t, np.nan),
            }
        )

    # sort by time
    idx = np.argsort(out["t"])
    for k in list(out.keys()):
        out[k] = out[k][idx]
    return out


def _apply_scientific_style() -> None:
    # Match src/visualisation/hamiltonian/__main__.py style.
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "serif"],
            "mathtext.fontset": "stix",
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.linewidth": 1.0,
            "grid.alpha": 0.5,
        }
    )


def plot_collapse_diagnostics(data: Dict[str, np.ndarray], out_path: Path) -> None:
    _apply_scientific_style()

    t = data["t"]

    fig, axes = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    axes = np.asarray(axes)

    # Left column: scalar collapse indicators (log scale)
    left_specs: Tuple[Tuple[str, str, str], ...] = (
        ("min_lapse", r"$\min(\alpha)$", r"Minimum lapse: $\alpha$"),
        ("min_chi", r"$\min(\chi)$", r"Minimum conformal factor: $\chi$"),
        ("max_abs_K", r"$\max(|K|)$", r"Maximum curvature: $|K|$"),
    )

    for i, (key, ylabel, title) in enumerate(left_specs):
        ax = axes[i, 0]
        y = np.asarray(data[key], dtype=float)
        # Avoid semilogy issues if zeros appear (shouldn't, but safe).
        y_plot = np.where(y > 0, y, np.nan)
        ax.semilogy(t, y_plot, color="blue", linewidth=1.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", alpha=0.6)
        ax.tick_params(axis="both", which="major", direction="in")

    # Right column: location of min lapse (linear)
    right_specs: Tuple[Tuple[str, str, str], ...] = (
        ("min_lapse_x", r"$x_{\min\alpha}$", r"Min-lapse location: $x$"),
        ("min_lapse_y", r"$y_{\min\alpha}$", r"Min-lapse location: $y$"),
        ("min_lapse_z", r"$z_{\min\alpha}$", r"Min-lapse location: $z$"),
    )
    for i, (key, ylabel, title) in enumerate(right_specs):
        ax = axes[i, 1]
        y = np.asarray(data[key], dtype=float)
        ax.plot(t, y, color="green", linewidth=1.5)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, which="both", ls="--", alpha=0.6)
        ax.tick_params(axis="both", which="major", direction="in")

    axes[-1, 0].set_xlabel(r"$t$")
    axes[-1, 1].set_xlabel(r"$t$")

    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


def main() -> None:
    script_dir = Path(__file__).resolve().parent
    default_data = _default_run_dir()

    parser = argparse.ArgumentParser(
        description="Plot collapse_diagnostics.dat (min lapse/chi/max|K| and min-lapse location)."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=None,
        help="Path to collapse_diagnostics.dat (optional; otherwise inferred from --data).",
    )
    parser.add_argument(
        "--data",
        default=str(default_data),
        help="Run output directory (expects data/collapse_diagnostics.dat inside).",
    )
    parser.add_argument(
        "--out",
        default=str(script_dir),
        help="Output directory for collapse_diagnostics_plot.png",
    )
    parser.add_argument(
        "--name",
        default="collapse_diagnostics_plot.png",
        help="Output filename (default: collapse_diagnostics_plot.png)",
    )
    args = parser.parse_args()

    data_dir = Path(args.data).expanduser().resolve()
    in_path = _resolve_input_path(data_dir, args.input)
    data = load_collapse_diagnostics(in_path)

    out_dir = Path(args.out).expanduser().resolve()
    out_path = out_dir / str(args.name)
    plot_collapse_diagnostics(data, out_path=out_path)


if __name__ == "__main__":
    main()