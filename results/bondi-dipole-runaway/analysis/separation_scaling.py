#!/usr/bin/env python3
"""Separation-scaling test: does the runaway follow the inverse-square law?

Three cells with IDENTICAL stars and numerics, differing only in the initial
separation (8, 12, 16).  Each is compared against its OWN point-mass Bondi
integration (same masses, same geometry), so the figure of merit is the ratio

    NR drift / point-mass drift        (1.0 = pure inverse-square point masses)

Why not fit a power law to the displacements directly: the closest pair's gap
closes substantially during the run, so its acceleration GROWS with time.  A
naive exponent fit at fixed time then reports 2.4-4.6 depending on when it is
evaluated -- it is measuring runaway feedback, not the force law.  Comparing
each run to its own point-mass integration removes that, because the reference
closes its gap too.

The phantom sector is used: the canonical sector's barycentre is biased by
boundary loss (its shed matter is pushed toward +x, exits the domain, and the
missing leading edge drags the remaining centroid backwards -- visible as a
spurious NEGATIVE canonical drift in the widest-gap cell).  Use
``sector_dynamics.dat`` (core tracker) instead when a run has it.

Usage: python3 analysis/separation_scaling.py [PACK_DIR]
"""

from __future__ import annotations

import csv
import pathlib
import sys

# Dressed-star ADM masses at omega = 0.550 (stars/star_family.csv).
M_CANON, M_PHANT = 0.06395, -0.07696
# Finest grid of the convergence campaign (dx = 0.25), which is where the
# article's separation numbers come from.
CELLS = {8: "convA_pm_n256", 12: "convA_pm_sep12_n256", 16: "convA_pm_sep16_n256"}
SAMPLE_TIMES = (20, 30, 40, 50)
DT = 0.01


def point_mass(sep: float, t_end: float = 60.0) -> dict[float, tuple[float, float]]:
    """RK4 point-mass Bondi pair; returns {t: (drift_canon, drift_phantom)}."""
    state = [sep / 2.0, 0.0, -sep / 2.0, 0.0]

    def deriv(s: list[float]) -> list[float]:
        x1, v1, x2, v2 = s
        r = x2 - x1
        r3 = max(abs(r), 1.0e-6) ** 3
        return [v1, M_PHANT * r / r3, v2, M_CANON * (-r) / r3]

    out = {0.0: (0.0, 0.0)}
    for step in range(1, int(t_end / DT) + 1):
        k1 = deriv(state)
        k2 = deriv([q + 0.5 * DT * k for q, k in zip(state, k1)])
        k3 = deriv([q + 0.5 * DT * k for q, k in zip(state, k2)])
        k4 = deriv([q + DT * k for q, k in zip(state, k3)])
        state = [
            q + DT / 6.0 * (a + 2 * b + 2 * c + d)
            for q, a, b, c, d in zip(state, k1, k2, k3, k4)
        ]
        t = step * DT
        if abs(t - round(t)) < 0.5 * DT:
            out[float(round(t))] = (state[0] - sep / 2.0, state[2] + sep / 2.0)
    return out


def load_drift(path: pathlib.Path) -> list[tuple[float, float, float]]:
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        v = [float(x) for x in line.split()]
        rows.append((v[0], v[2], v[7]))
    x0c, x0p = rows[0][1], rows[0][2]
    return [(t, xc - x0c, xp - x0p) for t, xc, xp in rows]


def main(argv: list[str]) -> int:
    pack = pathlib.Path(argv[0]) if argv else pathlib.Path(__file__).resolve().parent.parent
    out_rows = []
    print("phantom drift, NR vs its own point-mass reference")
    print(" sep |" + "".join(f"   t={t:<3d}" for t in SAMPLE_TIMES))
    for sep, cell in CELLS.items():
        path = pack / "campaign" / cell / "sector_barycenters.dat"
        if not path.is_file():
            print(f" {sep:3d} | (missing {cell})")
            continue
        series = load_drift(path)
        reference = point_mass(float(sep))
        ratios = []
        for t in SAMPLE_TIMES:
            nr = min(series, key=lambda r: abs(r[0] - t))
            pm_c, pm_p = reference[float(t)]
            ratio = nr[2] / pm_p if pm_p else float("nan")
            ratios.append(ratio)
            out_rows.append({
                "separation": sep,
                "t": t,
                "drift_phantom_nr": f"{nr[2]:.6f}",
                "drift_phantom_pointmass": f"{pm_p:.6f}",
                "ratio_phantom": f"{ratio:.4f}",
                "drift_canon_nr": f"{nr[1]:.6f}",
                "drift_canon_pointmass": f"{pm_c:.6f}",
            })
        print(f" {sep:3d} |" + "".join(f"  {r:6.2f}" for r in ratios))

    path = pack / "analysis" / "separation_scaling.csv"
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"\n[scaling] wrote {path} ({len(out_rows)} rows)")
    print(
        "Read: ratio -> 1 means the drift IS the inverse-square point-mass result.\n"
        "The excess at sep 8 is the overlap enhancement (rms radius ~5 per star);\n"
        "it falls to ~10% by sep 12-16, and grows with time in every cell because\n"
        "the gap closes as the pair runs away."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
