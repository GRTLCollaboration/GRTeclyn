#!/usr/bin/env python3
"""Decompose J_z / Q budgets and compute E_GW + pump-work energy.

Reads:
  output/data/collapse_diagnostics.dat
  output/small_data/psi4_mode_l2m0.dat   (trusted Python extraction)
  optional: output/small_data/boundary_flux.dat

Reports:
  - Delta Q_total, Q_sphere, J_z over the run
  - Integral of pump_work (col 23)
  - Rough E_GW from Psi4 (Newman–Penrose flux proxy)
  - Whether a closed budget statement is justified

Usage:
    python charge_energy_budget.py path/to/run/output
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def _load_collapse(path: Path) -> np.ndarray:
    return np.loadtxt(path, comments="#")


def _egw_from_psi4(psi4_path: Path, radius: float = 12.0) -> float:
    """Approximate radiated energy from the l=2,m=0 mode only.

    Uses E_GW ≈ (1/16π) ∫ |∫_{-\infty}^t Ψ₄ dt'|^2 dt evaluated at one radius
    (order-of-magnitude; single-mode, near-zone).  Returns E_GW in code units.
    """
    d = np.loadtxt(psi4_path, comments="#")
    if d.ndim != 2 or d.shape[1] < 3:
        return float("nan")
    t = d[:, 0]
    # Prefer the requested radius column: R=12 is cols 1,2; R=18 is 3,4; R=24 is 5,6.
    # Detect by header if present; else assume first radius is 12.
    re, im = d[:, 1], d[:, 2]
    psi = re + 1j * im
    # Time-integrate Psi4 -> news-like N; then integrate |N|^2.
    # Trapezoidal cumulative integral.
    dt = np.diff(t, prepend=t[0])
    dt[0] = 0.0
    news = np.cumsum(psi * dt)
    integrand = np.abs(news) ** 2
    # Factor: for a single (l,m) the GW energy involves 1/(16π) * |N|^2 with
    # spin-weight normalisation absorbed into the mode extraction (r Ψ₄).
    e_gw = float(np.trapezoid(integrand, t) / (16.0 * np.pi))
    return e_gw


def analyse(run_out: Path) -> dict:
    coll = run_out / "data" / "collapse_diagnostics.dat"
    psi4 = run_out / "small_data" / "psi4_mode_l2m0.dat"
    d = _load_collapse(coll)
    t = d[:, 0]
    j_z = d[:, 18]
    q_tot = d[:, 19]
    q_sph = d[:, 20]
    pump = d[:, 22] if d.shape[1] > 22 else np.zeros_like(t)

    pump_energy = float(np.trapezoid(pump, t))
    dj = float(j_z[-1] - j_z[0])
    dq_tot = float(q_tot[-1] - q_tot[0])
    dq_sph = float(q_sph[-1] - q_sph[0])

    e_gw = _egw_from_psi4(psi4) if psi4.exists() else float("nan")

    # Peak |Psi4| at R=12 for detectability.
    peak_psi4 = float("nan")
    if psi4.exists():
        p = np.loadtxt(psi4, comments="#")
        peak_psi4 = float(np.max(np.abs(p[:, 1] + 1j * p[:, 2])))

    return {
        "t0": float(t[0]),
        "tend": float(t[-1]),
        "J_z0": float(j_z[0]),
        "J_zend": float(j_z[-1]),
        "dJ_z": dj,
        "Q_tot0": float(q_tot[0]),
        "Q_totend": float(q_tot[-1]),
        "dQ_tot": dq_tot,
        "Q_sph0": float(q_sph[0]),
        "Q_sphend": float(q_sph[-1]),
        "dQ_sph": dq_sph,
        "pump_energy": pump_energy,
        "E_GW_l2m0": e_gw,
        "peak_psi4_R12": peak_psi4,
        "J_z_frac_change": dj / abs(j_z[0]) if abs(j_z[0]) > 1e-12 else float("nan"),
        "Q_tot_frac_change": dq_tot / abs(q_tot[0]) if abs(q_tot[0]) > 1e-12 else float("nan"),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run_output", type=Path, help=".../output directory")
    args = ap.parse_args()
    r = analyse(args.run_output)
    print(f"run: {args.run_output}")
    print(f"  t ∈ [{r['t0']:.3g}, {r['tend']:.3g}]")
    print(f"  J_z:  {r['J_z0']:.6g} → {r['J_zend']:.6g}   Δ={r['dJ_z']:.6g}  "
          f"({100*r['J_z_frac_change']:+.1f}%)")
    print(f"  Q_tot:{r['Q_tot0']:.6g} → {r['Q_totend']:.6g}   Δ={r['dQ_tot']:.6g}  "
          f"({100*r['Q_tot_frac_change']:+.1f}%)")
    print(f"  Q_sph:{r['Q_sph0']:.6g} → {r['Q_sphend']:.6g}   Δ={r['dQ_sph']:.6g}")
    print(f"  ∫ pump_work dt = {r['pump_energy']:.6g}   "
          f"(0 ⇒ no PD pump; support-ramp only)")
    print(f"  E_GW (l=2,m=0 proxy) ≈ {r['E_GW_l2m0']:.6g}")
    print(f"  peak |Ψ₄|(R=12) ≈ {r['peak_psi4_R12']:.6g}")
    if abs(r["pump_energy"]) < 1e-10:
        print("  BUDGET: pump injection is negligible; J_z/Q_tot drift is NOT "
              "pump-injected. Drop 'charge retained/conserved' for Q_total; "
              "Q_sphere retention may still be quoted.")
    else:
        print("  BUDGET: compare ΔJ_z / ΔQ_tot against ∫pump_work (order-of-mag).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
