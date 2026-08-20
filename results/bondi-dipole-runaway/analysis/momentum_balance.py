#!/usr/bin/env python3
"""Item B: is the runaway momentum-balanced in the MATTER sector?

Bondi's point-mass runaway conserves total momentum exactly: the canonical
body carries p = +Mv, the phantom body p = -|M|v (negative mass, positive
velocity), and the pair self-accelerates at zero net momentum.  The field
theory version has a third player -- the gravitational field -- so the matter
volume integral alone need not vanish.  This script measures how far from
zero it actually is, and whether matching the sector |ADM| masses (the _eqm
cell) tightens the balance.

Column conventions (extraction/sector_dynamics.py): ``px_canon`` and
``px_phantom`` are each the plain integral S_x sqrt(gamma) d3x of their own
complex field, BOTH evaluated with the canonical-sign formula.  The phantom
sector's stress tensor enters the constraints with the opposite sign, so its
physical (gravitating) momentum is -px_phantom and the total the metric feels
is ``px_total = px_canon - px_phantom``.

For scale, each sector's integral is also compared against the naive
"mass times core velocity" estimate M * dx_core/dt using the dressed ADM
masses -- the point-mass number the volume integral would match if the lump
moved rigidly and the field contribution were negligible.

Own script and own output file (``analysis/momentum_balance.csv``) on
purpose: the published trajectory/summary CSV contracts are not widened.

Usage: python3 analysis/momentum_balance.py [PACK_DIR]
"""

from __future__ import annotations

import csv
import pathlib
import sys

# Dressed-star ADM masses (stars/star_family.csv, star_radius.csv).
M_CANON = 0.06395
CELLS = {
    # cell dir -> phantom ADM mass (signed)
    "pair_pm_v2": -0.07696,      # standard pair: |M-| exceeds |M+| by 20%
    "pair_pm_eqm_v2": -0.06395,  # equal-|ADM| pair (phantom at omega=0.56598)
}
# The published science window ends at t=30 (bath contamination beyond ~60;
# tagger effects beyond ~47 in the mixed cell), so quote the balance there.
SAMPLE_TIMES = (10.0, 20.0, 30.0, 45.0, 60.0)


def _read_stream(path: pathlib.Path) -> dict[str, list[float]]:
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().lstrip("#").split()
        cols: dict[str, list[float]] = {name: [] for name in header}
        for line in fh:
            parts = line.split()
            if len(parts) != len(header):
                continue
            for name, tok in zip(header, parts):
                cols[name].append(float(tok))
    return cols


def _dxdt(ts: list[float], xs: list[float], i: int) -> float:
    lo, hi = max(0, i - 1), min(len(ts) - 1, i + 1)
    return (xs[hi] - xs[lo]) / (ts[hi] - ts[lo]) if ts[hi] > ts[lo] else 0.0


def _nearest(ts: list[float], t: float) -> int:
    return min(range(len(ts)), key=lambda i: abs(ts[i] - t))


def main(argv: list[str]) -> int:
    pack = (
        pathlib.Path(argv[0])
        if argv
        else pathlib.Path(__file__).resolve().parent.parent
    )
    out_rows: list[dict[str, object]] = []

    for cell, m_phantom in CELLS.items():
        stream = pack / "campaign" / cell / "sector_dynamics.dat"
        if not stream.is_file():
            print(f"[momB] {cell}: no sector_dynamics.dat -- skipped")
            continue
        cols = _read_stream(stream)
        ts = cols["time"]
        for i, t in enumerate(ts):
            v_c = _dxdt(ts, cols["core_x_canon"], i)
            v_p = _dxdt(ts, cols["core_x_phantom"], i)
            out_rows.append(
                {
                    "cell": cell,
                    "time": f"{t:.4f}",
                    "px_canon": f"{cols['px_canon'][i]:.6e}",
                    "px_phantom_raw": f"{cols['px_phantom'][i]:.6e}",
                    "px_total": f"{cols['px_total'][i]:.6e}",
                    "v_core_canon": f"{v_c:.6e}",
                    "v_core_phantom": f"{v_p:.6e}",
                    "mv_canon": f"{M_CANON * v_c:.6e}",
                    "mv_phantom": f"{m_phantom * v_p:.6e}",
                    "mv_total": f"{M_CANON * v_c + m_phantom * v_p:.6e}",
                    "coord_sep": f"{cols['coord_sep'][i]:.6f}",
                }
            )

        print(f"\n=== {cell} (M+={M_CANON:+.5f}, M-={m_phantom:+.5f}) ===")
        for t_q in SAMPLE_TIMES:
            i = _nearest(ts, t_q)
            px_c, px_t = cols["px_canon"][i], cols["px_total"][i]
            v_c = _dxdt(ts, cols["core_x_canon"], i)
            v_p = _dxdt(ts, cols["core_x_phantom"], i)
            mv_t = M_CANON * v_c + m_phantom * v_p
            scale = max(abs(px_c), abs(cols["px_phantom"][i]), 1e-300)
            print(
                f"  t={ts[i]:5.1f}  px_canon={px_c:+.3e}  "
                f"px_total={px_t:+.3e} ({px_t / scale:+7.1%} of sector scale)  "
                f"point-mass mv_total={mv_t:+.3e}  sep={cols['coord_sep'][i]:.2f}"
            )

    out = pack / "analysis" / "momentum_balance.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"\n[momB] wrote {out} ({len(out_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
