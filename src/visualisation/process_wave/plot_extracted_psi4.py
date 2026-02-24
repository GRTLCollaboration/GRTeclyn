#!/usr/bin/env python3
"""
Plot extracted Psi4 time-series produced by consume_plotfiles.py.

Input format (whitespace separated):
  time  Re(R=r1) Im(R=r1)  Re(R=r2) Im(R=r2)  ...

This avoids reading plotfiles entirely and is the intended "later plot" step.
"""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter, welch

matplotlib.use("Agg")


def _default_data_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent.parent.parent
    return (project_root.parent / "data").resolve()


def _smooth_psd(psd: np.ndarray, window: int, polyorder: int) -> np.ndarray:
    psd = np.asarray(psd, dtype=float)
    out = psd.copy()
    m = np.isfinite(psd) & (psd > 0)
    if np.sum(m) < 7:
        return out
    y = np.log10(psd[m])
    n = y.size
    w = int(window)
    if w < 5:
        return out
    if w % 2 == 0:
        w += 1
    if w > n:
        w = n if (n % 2 == 1) else (n - 1)
    p = int(polyorder)
    if p < 1:
        return out
    if p >= w:
        p = max(1, w - 2)
    y_s = savgol_filter(y, window_length=w, polyorder=p, mode="interp")
    out[m] = 10 ** y_s
    return out


def _parse_header_radii(header_line: str) -> List[float]:
    # Matches Re(R=14) / Im(R=14.0) style
    rs = re.findall(r"R=([0-9.+-eE]+)", header_line)
    radii: List[float] = []
    for x in rs[::2]:  # Re/Im pairs repeat the same R; take every other
        try:
            radii.append(float(x))
        except ValueError:
            continue
    # Deduplicate preserving order
    out: List[float] = []
    for r in radii:
        if not out or abs(out[-1] - r) > 0:
            out.append(r)
    return out


def load_extracted(path: Path) -> Tuple[np.ndarray, List[float], Dict[float, np.ndarray]]:
    header_radii: List[float] = []
    times: List[float] = []
    cols: List[List[float]] = []

    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if not header_radii and "R=" in s:
                    header_radii = _parse_header_radii(s)
                continue
            parts = s.split()
            try:
                t = float(parts[0])
            except ValueError:
                continue
            nums = []
            ok = True
            for p in parts[1:]:
                try:
                    nums.append(float(p))
                except ValueError:
                    ok = False
                    break
            if not ok:
                continue
            times.append(t)
            cols.append(nums)

    if not times:
        raise ValueError(f"No data rows found in {path}")

    t_arr = np.asarray(times, dtype=float)
    data = np.asarray(cols, dtype=float)

    if data.shape[1] % 2 != 0:
        raise ValueError(f"Expected Re/Im pairs; got {data.shape[1]} columns in {path}")

    n_r = data.shape[1] // 2
    if header_radii and len(header_radii) != n_r:
        # Header mismatch; ignore header and infer radii indices
        header_radii = []

    radii = header_radii if header_radii else [float(i) for i in range(n_r)]

    out: Dict[float, np.ndarray] = {}
    for i, r in enumerate(radii):
        re_col = data[:, 2 * i]
        im_col = data[:, 2 * i + 1]
        out[r] = re_col + 1j * im_col

    # sort by time
    idx = np.argsort(t_arr)
    t_arr = t_arr[idx]
    for r in list(out.keys()):
        out[r] = out[r][idx]

    return t_arr, radii, out


def main() -> None:
    default_data = _default_data_dir() / "small_data" / "psi4_mode_l2m0.dat"
    script_dir = Path(__file__).resolve().parent

    parser = argparse.ArgumentParser(description="Plot extracted Psi4 time-series (.dat) with PSD.")
    parser.add_argument("input", nargs="?", default=str(default_data), help=f"Input .dat (default: {default_data})")
    parser.add_argument("--out", default=str(script_dir), help="Output directory (default: this folder)")
    parser.add_argument("--radii", type=float, nargs="*", default=None, help="Subset of radii to plot (defaults to all)")
    parser.add_argument("--time-axis", choices=["simulation", "retarded"], default="simulation")
    parser.add_argument("--t-min", type=float, default=None)
    parser.add_argument("--t-max", type=float, default=None)
    parser.add_argument("--psd-smooth-window", type=int, default=21)
    parser.add_argument("--psd-smooth-polyorder", type=int, default=5)
    parser.add_argument("--psd-hide-raw", action="store_true")
    args = parser.parse_args()

    in_path = Path(args.input).expanduser().resolve()
    if not in_path.exists():
        raise SystemExit(f"Input not found: {in_path}")

    t, radii_all, series = load_extracted(in_path)

    radii = radii_all
    if args.radii:
        # Keep only requested radii that exist
        radii = [r for r in args.radii if r in series]
        if not radii:
            raise SystemExit(f"No requested radii found in file. Available: {radii_all}")

    # Publication-style LaTeX (mathtext) typography
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
            "mathtext.fontset": "cm",
            "axes.labelsize": 14,
            "axes.titlesize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "axes.linewidth": 1.2,
            "grid.alpha": 0.5,
            "legend.fontsize": 12,
            "figure.figsize": (10, 8),
        }
    )

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    n_r = len(radii)
    psi4_colors = plt.cm.Blues(np.linspace(0.5, 0.85, n_r))[::-1]
    psd_colors = plt.cm.Greens(np.linspace(0.5, 0.85, n_r))[::-1]

    # time step for PSD
    dt = np.median(np.diff(t)) if t.size >= 2 else np.nan
    fs = (1.0 / dt) if (np.isfinite(dt) and dt > 0) else None

    plotted_any = False
    for i, R in enumerate(radii):
        psi4 = series[R]
        x = t if args.time_axis == "simulation" else (t - float(R))
        y = np.real(psi4)

        m = np.ones_like(x, dtype=bool)
        if args.t_min is not None:
            m &= x >= args.t_min
        if args.t_max is not None:
            m &= x <= args.t_max

        ax1.plot(x[m], y[m], color=psi4_colors[i], linewidth=1.5, label=f"R={R:g}")

        if fs is not None and psi4.size >= 8:
            freqs, psd = welch(np.real(psi4), fs, nperseg=min(128, max(8, psi4.size // 2)))
            psd_s = _smooth_psd(psd, window=args.psd_smooth_window, polyorder=args.psd_smooth_polyorder)
            if not args.psd_hide_raw:
                ax2.semilogy(freqs, psd, "o", color="red", markersize=3.0, markeredgewidth=0, alpha=0.35)
            ax2.semilogy(freqs, psd_s, "-", color=psd_colors[i], linewidth=1.5, label=f"R={R:g}")

        plotted_any = True

    ax1.set_xlabel(r"$t$" if args.time_axis == "simulation" else r"$t - R_{\mathrm{ext}}$")
    ax1.set_ylabel(r"$r\,\mathrm{Re}\!\left(\Psi_4^{2,0}\right)$")
    ax1.set_title(r"Gravitational-wave curvature $\Psi_4$ (mode $\ell=2$, $m=0$)")
    if plotted_any:
        ax1.legend(loc="upper right", frameon=True, framealpha=0.9)
    ax1.grid(True, which="both", ls="--", alpha=0.6)
    ax1.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    ax2.set_xlabel(r"$f\,(M^{-1})$")
    ax2.set_ylabel(r"$\mathrm{PSD}\left[\Psi_4\right]$")
    ax2.grid(True, which="major", ls="--", alpha=0.6)
    ax2.tick_params(axis="both", which="major", direction="in", top=True, right=True)

    out_dir = Path(args.out).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_".join([f"R{r:g}" for r in radii])
    out_path = out_dir / f"psi4_extracted_{suffix}.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved to {out_path}")


if __name__ == "__main__":
    main()

