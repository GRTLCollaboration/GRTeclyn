#!/usr/bin/env python3
"""Point-mass Bondi reference trajectories, for comparison with the NR runs.

Bondi's setup in the bicomplex model: both sectors carry POSITIVE inertial and
passive gravitational mass (identical Klein-Gordon equations), while the ACTIVE
gravitational mass is signed (the sign lives only in the Einstein source).  So
each body's acceleration is set entirely by the OTHER body's active mass, and
is independent of its own mass (passive = inertial cancels):

    a_1 = M2_active * (x2 - x1) / |x2 - x1|^3
    a_2 = M1_active * (x1 - x2) / |x1 - x2|^3

With M1 = +M+ (canonical at +x) and M2 = -|M-| (phantom at -x) both
accelerations point +x -- the runaway.  Note the point-mass prediction for the
GAP: it grows when |M-| > |M+| (the canonical body is pushed harder than the
phantom is pulled), and is CONSTANT when the masses match.  Neither is what the
NR runs do, because at separation 8 the stars (rms radius ~5) overlap heavily.

Writes ``analysis/newtonian_reference.csv`` and prints the observed-vs-point-mass
comparison table used in FINDINGS.md.

Usage: python3 analysis/newtonian_reference.py [PACK_DIR]
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys

# Dressed-star ADM masses solved for this campaign (see stars/star_family.csv).
M_CANON = 0.06395
M_PHANT_UNEQUAL = -0.07696  # omega = 0.550
M_PHANT_EQUAL = -0.06395  # omega = 0.56598, |ADM| matched to the canonical star

# Finest grid of the convergence campaign (dx = 0.25).
CASES = [
    ("convA_pm_n256", M_CANON, M_PHANT_UNEQUAL),
    ("convA_pm_eqm_n256", M_CANON, M_PHANT_EQUAL),
]
X0_CANON, X0_PHANT = 4.0, -4.0  # grid centre subtracted (run coords: 36 / 28)
T_END, DT = 60.0, 0.01


def derivatives(state: list[float], m1: float, m2: float) -> list[float]:
    x1, v1, x2, v2 = state
    r = x2 - x1
    r3 = max(abs(r), 1.0e-6) ** 3
    return [v1, m2 * r / r3, v2, m1 * (-r) / r3]


def integrate(m1: float, m2: float) -> list[tuple[float, float, float]]:
    """RK4; returns [(t, x1, x2)] sampled every 1.0 in time."""
    state = [X0_CANON, 0.0, X0_PHANT, 0.0]
    out = [(0.0, state[0], state[2])]
    steps = int(round(T_END / DT))
    for step in range(1, steps + 1):
        k1 = derivatives(state, m1, m2)
        k2 = derivatives([s + 0.5 * DT * k for s, k in zip(state, k1)], m1, m2)
        k3 = derivatives([s + 0.5 * DT * k for s, k in zip(state, k2)], m1, m2)
        k4 = derivatives([s + DT * k for s, k in zip(state, k3)], m1, m2)
        state = [
            s + DT / 6.0 * (a + 2 * b + 2 * c + d)
            for s, a, b, c, d in zip(state, k1, k2, k3, k4)
        ]
        t = step * DT
        if abs(t - round(t)) < 0.5 * DT:
            out.append((round(t), state[0], state[2]))
        if abs(state[2] - state[0]) < 1.0e-3:
            break  # contact: the point-mass idealisation is over
    return out


def read_observed(pack: pathlib.Path, cell: str) -> dict[float, tuple[float, float]]:
    path = pack / "campaign" / cell / "sector_barycenters.dat"
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        parts = [float(v) for v in line.split()]
        rows.append((parts[0], parts[2], parts[7]))
    x0c, x0p = rows[0][1], rows[0][2]
    return {t: (xc - x0c, xp - x0p) for t, xc, xp in rows}


def nearest(series: dict[float, tuple[float, float]], t: float) -> tuple[float, float]:
    return series[min(series, key=lambda k: abs(k - t))]


def main(argv: list[str]) -> int:
    pack = pathlib.Path(argv[0]) if argv else pathlib.Path(__file__).resolve().parent.parent
    out_rows = []
    print(
        "cell          t   drift_canon (NR / point-mass)   drift_phantom (NR / point-mass)"
        "   separation (NR / point-mass)"
    )
    for cell, m1, m2 in CASES:
        model = integrate(m1, m2)
        observed = read_observed(pack, cell)
        for t, x1, x2 in model:
            d1_pm, d2_pm = x1 - X0_CANON, x2 - X0_PHANT
            d1_nr, d2_nr = nearest(observed, t)
            out_rows.append({
                "cell": cell,
                "t": t,
                "drift_canon_nr": f"{d1_nr:.4f}",
                "drift_canon_pointmass": f"{d1_pm:.4f}",
                "drift_phantom_nr": f"{d2_nr:.4f}",
                "drift_phantom_pointmass": f"{d2_pm:.4f}",
                "sep_nr": f"{8.0 + d1_nr - d2_nr:.4f}",
                "sep_pointmass": f"{x1 - x2:.4f}",
            })
            if t in (10, 20, 30, 40, 60):
                print(
                    f"{cell:12s} {t:3.0f}   {d1_nr:+7.3f} / {d1_pm:+7.3f}"
                    f"              {d2_nr:+7.3f} / {d2_pm:+7.3f}"
                    f"            {8.0 + d1_nr - d2_nr:6.3f} / {x1 - x2:6.3f}"
                )
        print()

    path = pack / "analysis" / "newtonian_reference.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"[newtonian] wrote {path} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
