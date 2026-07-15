#!/usr/bin/env python3
"""Extrapolate rΨ₄ → r→∞ with a 1/r polynomial fit across extraction radii.

Also reports relative constraint-norm diagnostics from constraint_norms.dat.

Usage:
  python psi4_extrapolate.py path/to/run/output [--radii 12 18 24]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def load_psi4(path: Path, radii: list[float]) -> tuple[np.ndarray, dict[float, np.ndarray]]:
    d = np.loadtxt(path, comments="#")
    t = d[:, 0]
    series = {}
    for i, R in enumerate(radii):
        series[R] = d[:, 1 + 2 * i] + 1j * d[:, 2 + 2 * i]
    return t, series


def extrapolate_peak(t: np.ndarray, series: dict[float, np.ndarray],
                     t_lo: float = 13.0, t_hi: float = 24.0) -> dict:
    """Fit |r Ψ₄|(R) = a0 + a1/R  (or a0 + a1/R + a2/R²) over the burst window."""
    radii = sorted(series)
    # Peak |Ψ4| in burst window at each radius
    mask = (t >= t_lo) & (t <= t_hi)
    peaks = []
    rpsi_peaks = []
    for R in radii:
        amp = np.abs(series[R][mask])
        pk = float(np.max(amp)) if amp.size else float("nan")
        peaks.append(pk)
        rpsi_peaks.append(pk * R)
    R = np.asarray(radii, dtype=float)
    y = np.asarray(rpsi_peaks, dtype=float)
    # Linear fit y = a0 + a1/R  →  extrapolant a0 = lim r→∞ rΨ₄
    A = np.column_stack([np.ones_like(R), 1.0 / R])
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    a0, a1 = coef
    # Quadratic if ≥3 radii
    a0_q = a0
    if len(R) >= 3:
        A2 = np.column_stack([np.ones_like(R), 1.0 / R, 1.0 / R**2])
        coef2, *_ = np.linalg.lstsq(A2, y, rcond=None)
        a0_q = float(coef2[0])
    return {
        "radii": radii,
        "peak_psi4": peaks,
        "peak_rpsi4": rpsi_peaks,
        "rpsi4_inf_linear": float(a0),
        "rpsi4_inf_quad": float(a0_q),
        "fit_a1": float(a1),
    }


def relative_constraints(path: Path) -> dict:
    d = np.loadtxt(path, comments="#")
    t, ham, mom = d[:, 0], d[:, 1], d[:, 2]
    ham0, mom0 = ham[0], mom[0]
    return {
        "t0": float(t[0]),
        "tend": float(t[-1]),
        "Ham_L2_max": float(np.max(ham)),
        "Mom_L2_max": float(np.max(mom)),
        "Ham_L2_t0": float(ham0),
        "Mom_L2_t0": float(mom0),
        "Ham_rel_max": float(np.max(ham) / ham0) if ham0 > 0 else float("nan"),
        "Mom_rel_max": float(np.max(mom) / mom0) if mom0 > 0 else float("nan"),
        "Ham_L2_end": float(ham[-1]),
        "Mom_L2_end": float(mom[-1]),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_output", type=Path)
    ap.add_argument("--radii", type=float, nargs="+", default=[12.0, 18.0, 24.0])
    args = ap.parse_args()
    psi4 = args.run_output / "small_data" / "psi4_mode_l2m0.dat"
    cons = args.run_output / "data" / "constraint_norms.dat"
    if psi4.exists():
        t, series = load_psi4(psi4, args.radii)
        ext = extrapolate_peak(t, series)
        print(f"Psi4 file: {psi4}")
        for R, p, rp in zip(ext["radii"], ext["peak_psi4"], ext["peak_rpsi4"]):
            print(f"  R={R:g}: peak|Ψ₄|={p:.4e}  peak|rΨ₄|={rp:.4e}")
        print(f"  lim r→∞ |rΨ₄| (linear  a0+a1/r) = {ext['rpsi4_inf_linear']:.4e}")
        print(f"  lim r→∞ |rΨ₄| (quad a0+a1/r+a2/r²) = {ext['rpsi4_inf_quad']:.4e}")
        print("  NOTE: extraction radii are near-zone (R < λ/2 for f≲0.1).")
    if cons.exists():
        c = relative_constraints(cons)
        print(f"Constraints: {cons}")
        print(f"  Ham L2: t0={c['Ham_L2_t0']:.3e} max={c['Ham_L2_max']:.3e} "
              f"rel_max={c['Ham_rel_max']:.2f}×")
        print(f"  Mom L2: t0={c['Mom_L2_t0']:.3e} max={c['Mom_L2_max']:.3e} "
              f"rel_max={c['Mom_rel_max']:.2f}×")
        print("  (No Θ/Z_i or L∞ columns in this file; relative = L2/L2(t0).)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
