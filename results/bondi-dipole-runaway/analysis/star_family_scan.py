#!/usr/bin/env python3
"""Scan the dressed-star families of BOTH sectors at the ultraweak rung.

For each requested frequency omega the fixed-frequency solver returns the
gravitationally dressed sextic star (amplitude shooting at fixed omega);
``gravity_sign = +1`` is the canonical sector, ``-1`` the phantom sector whose
self-gravity is repulsive (negative ADM mass, lapse > 1 in the core).

Writes ``stars/star_family.csv`` -- the M(omega) curve behind the equal-|ADM|
pairing: the phantom star at omega = 0.56598 weighs exactly as much as the
canonical star at omega = 0.550.

Run with the wrapper venv (needs numpy/scipy):
    PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
        results/bondi-dipole-runaway/analysis/star_family_scan.py
"""

from __future__ import annotations

import csv
import pathlib
import sys

# Ultraweak rung: lambda^2/mu = 4.8 => omega_min = sqrt(1 - 3 lambda^2/16 mu) = 0.316.
MASS, LAM, MU = 1.0, 10240.0, 21845333.0
OMEGAS = [
    0.52, 0.53, 0.54, 0.55, 0.56, 0.56598, 0.575, 0.60,
    0.625, 0.65, 0.675, 0.70, 0.75, 0.80,
]


def main(argv: list[str]) -> int:
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )

    pack = pathlib.Path(argv[0]) if argv else pathlib.Path(__file__).resolve().parent.parent
    out = pack / "stars" / "star_family.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for omega in OMEGAS:
        for sector, sign in (("canonical", 1.0), ("phantom", -1.0)):
            try:
                star = cached_selfgrav_at_omega(MASS, LAM, MU, omega, sign)
            except Exception as exc:  # branch edge: no dressed star exists here
                rows.append({
                    "sector": sector, "omega_requested": omega, "status": f"no solution: {exc}",
                    "omega_solved": "", "phi_c": "", "adm_mass": "",
                    "alpha_central": "", "psi_central": "", "compactness": "",
                })
                print(f"[scan] {sector:9s} omega={omega:.5f}  NO SOLUTION ({exc})")
                continue
            rows.append({
                "sector": sector,
                "omega_requested": omega,
                "status": "ok",
                "omega_solved": f"{star.omega:.9f}",
                "phi_c": f"{star.phi_c:.9f}",
                "adm_mass": f"{star.adm_mass:.9f}",
                "alpha_central": f"{float(star.alpha[0]):.9f}",
                "psi_central": f"{float(star.psi[0]):.9f}",
                "compactness": f"{star.compactness:.9e}",
            })
            print(
                f"[scan] {sector:9s} omega={omega:.5f}  phi_c={star.phi_c:.6f}  "
                f"ADM={star.adm_mass:+.6f}  alpha(0)={float(star.alpha[0]):.5f}"
            )

    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[scan] wrote {out} ({len(rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
