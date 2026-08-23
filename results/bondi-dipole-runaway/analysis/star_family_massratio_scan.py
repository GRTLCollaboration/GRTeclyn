#!/usr/bin/env python3
"""Scan both star branches into the LIGHT end, for the mass-ratio ladder.

``star_family_scan.py`` covers omega 0.52-0.80, which is the range the matched
pair and the phase-3 cell live in.  The mass-ratio ladder (phase 3b) needs
lighter stars than that on both branches: a phantom light enough to sit at 0.60
and 0.40 of the matched mass, and a canonical light enough to be OUTWEIGHED by
the matched phantom.  This scan extends both branches to omega 1.0 and answers
the question the campaign plan asks it to record -- **whether the phantom
branch has a bound on how light its stars can be**.

It does.  |M-| falls monotonically to a MINIMUM of about 0.00791 near
omega = 0.94, then rises again, and at omega = 1.0 no bound star exists at all.
The lightest phantom star on this branch is therefore ~0.55 of the matched
mass, which is why the ladder has a 0.60 rung and no 0.40 rung.

Separate file, separate CSV: ``star_family.csv`` is the table the matched
pairing was chosen from and every published number is traceable to it, so it is
not rewritten to answer a later question.

Run with the wrapper venv (needs numpy/scipy):
    PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
        results/bondi-dipole-runaway/analysis/star_family_massratio_scan.py
"""

from __future__ import annotations

import csv
import pathlib
import sys

# Same rung as star_family_scan.py -- the tables must be comparable.
MASS, LAM, MU = 1.0, 10240.0, 21845333.0

# The mass every ratio in the ladder is quoted against: the canonical star at
# omega = 0.75, whose |ADM| the matched phantom at omega = 0.7603 reproduces.
MATCHED_MASS = 0.014350

OMEGAS = [
    0.75, 0.7603, 0.78, 0.80, 0.8040, 0.81, 0.82, 0.84, 0.86, 0.88,
    0.90, 0.92, 0.94, 0.96, 0.98, 1.00,
]


def main(argv: list[str]) -> int:
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )

    pack = (
        pathlib.Path(argv[0])
        if argv
        else pathlib.Path(__file__).resolve().parent.parent
    )
    out = pack / "stars" / "star_family_massratio.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for omega in OMEGAS:
        for sector, sign in (("canonical", 1.0), ("phantom", -1.0)):
            try:
                star = cached_selfgrav_at_omega(MASS, LAM, MU, omega, sign)
            except Exception as exc:  # branch edge: no dressed star exists here
                rows.append({
                    "sector": sector, "omega_requested": omega,
                    "status": f"no solution: {exc}", "omega_solved": "",
                    "phi_c": "", "adm_mass": "", "abs_mass_over_matched": "",
                })
                print(f"[massratio] {sector:9s} omega={omega:.5f}  NO SOLUTION ({exc})")
                continue
            rows.append({
                "sector": sector,
                "omega_requested": omega,
                "status": "ok",
                "omega_solved": f"{star.omega:.9f}",
                "phi_c": f"{star.phi_c:.9f}",
                "adm_mass": f"{star.adm_mass:.9f}",
                "abs_mass_over_matched": f"{abs(star.adm_mass) / MATCHED_MASS:.6f}",
            })
            print(
                f"[massratio] {sector:9s} omega={omega:.5f}  "
                f"M={star.adm_mass:+.6f}  |M|/matched={abs(star.adm_mass)/MATCHED_MASS:.4f}"
            )

    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"[massratio] wrote {out} ({len(rows)} rows)")

    phantom = [
        (r["omega_requested"], float(r["adm_mass"]))
        for r in rows
        if r["sector"] == "phantom" and r["status"] == "ok"
    ]
    if phantom:
        w_min, m_min = min(phantom, key=lambda t: abs(t[1]))
        print(
            f"[massratio] lightest phantom on this branch: |M| = {abs(m_min):.6f} "
            f"at omega = {w_min} ({abs(m_min)/MATCHED_MASS:.3f} of matched)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
