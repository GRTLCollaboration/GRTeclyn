#!/usr/bin/env python3
"""Score confinement runs in ABSOLUTE per-sector activity (Debug.md 19.1).

`confined_frac` is a ratio whose denominator (total_activity) grows during
every run -- unpumped through dispersal, pumped through injection -- so the
fraction moves when nothing is being lost. This scorer reads the same
confinement.dat and reports the quantity that answers the actual question,
"how much matter is still inside r_conf":

    abs_total = total_activity * confined_frac
    abs_canon = total_activity * canon_mass_frac * confined_frac_canon
    abs_phan  = total_activity * (1 - canon_mass_frac) * confined_frac_phantom

Caveats carried from Debug.md 18.1: activity is amplitude-weighted (not an
energy) and canon_mass_frac is an amplitude share, so the sector split is
approximate. Comparisons across runs of the same seed are safe; absolute
cross-seed comparisons are not.

Per run it prints the three-phase summary the ladder re-read established
(Debug.md 19.2): injection (first slice vs t=0), hold flatness, late decline
rate, end state -- plus the governor safety record from constraint_norms.dat
(an arm whose governor crosses 0.5 before t=27 bought its numbers with
constraint violation; baseline tp30 crosses at t~28.9).

Usage:
    score_confine_abs.py RUN_DIR [RUN_DIR ...] [--ref RUN_DIR] [--tsv]

    --ref   baseline for the gain columns (default: none)
    --tsv   machine-readable one-line-per-run output
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# confinement.dat columns (0-indexed), per the consumer header
C_TIME, C_TOTAL, C_FRAC = 0, 1, 4
C_FRAC_CANON, C_FRAC_PHAN, C_CMF, C_MINCHI = 13, 15, 16, 17

LATE_T0, LATE_T1 = 20.16, 27.36   # clean late window: past onset, pre-choke
HOLD_T0, HOLD_T1 = 1.44, 18.72    # hold phase


def read_confinement(run: Path) -> list[list[float]]:
    f = run / "small_data/confinement.dat"
    rows = []
    for line in f.read_text().splitlines():
        s = line.split()
        if not s or s[0].startswith("#"):
            continue
        rows.append([float(x) for x in s])
    return rows


def governor_cross(run: Path) -> tuple[str, str]:
    """(first t with governor < 0.5, min governor) from constraint_norms.dat."""
    f = run / "data/constraint_norms.dat"
    if not f.exists():
        return "n/a", "n/a"
    cross, gmin = None, 1.0
    for line in f.read_text().splitlines():
        s = line.split()
        if len(s) < 10:
            continue
        t, g = float(s[0]), float(s[9])
        gmin = min(gmin, g)
        if g < 0.5 and cross is None:
            cross = t
    return ("never" if cross is None else f"{cross:.2f}"), f"{gmin:.4f}"


def at(rows: list[list[float]], t: float) -> list[float] | None:
    best = min(rows, key=lambda r: abs(r[C_TIME] - t))
    return best if abs(best[C_TIME] - t) < 0.75 else None


def absolutes(r: list[float]) -> tuple[float, float, float]:
    tot = r[C_TOTAL] * r[C_FRAC]
    can = r[C_TOTAL] * r[C_CMF] * r[C_FRAC_CANON]
    pha = r[C_TOTAL] * (1.0 - r[C_CMF]) * r[C_FRAC_PHAN]
    return tot, can, pha


def rate(a0: float, a1: float, dt: float) -> float:
    """Exponential decline rate; positive = losing matter."""
    if a0 <= 0 or a1 <= 0 or dt <= 0:
        return math.nan
    return math.log(a0 / a1) / dt


def score(run: Path) -> dict | None:
    rows = read_confinement(run)
    if len(rows) < 2:
        return None
    r0 = rows[0]
    first = rows[1]
    t_end = rows[-1][C_TIME]
    out = {"name": run.name, "t_end": t_end, "n_rows": len(rows)}

    a0 = absolutes(r0)
    a1 = absolutes(first)
    out["inj_canon"] = a1[1] / a0[1] if a0[1] > 0 else math.nan
    out["inj_phan"] = a1[2] / a0[2] if a0[2] > 0 else math.nan

    h0, h1 = at(rows, HOLD_T0), at(rows, HOLD_T1)
    if h0 and h1:
        b0, b1 = absolutes(h0), absolutes(h1)
        out["hold_canon"] = rate(b0[1], b1[1], h1[C_TIME] - h0[C_TIME])
        out["hold_phan"] = rate(b0[2], b1[2], h1[C_TIME] - h0[C_TIME])

    l0, l1 = at(rows, LATE_T0), at(rows, LATE_T1)
    if l0 and l1:
        c0, c1 = absolutes(l0), absolutes(l1)
        out["late_canon"] = rate(c0[1], c1[1], l1[C_TIME] - l0[C_TIME])
        out["late_phan"] = rate(c0[2], c1[2], l1[C_TIME] - l0[C_TIME])

    e = absolutes(rows[-1])
    out["end_total"], out["end_canon"], out["end_phan"] = e
    out["min_chi"] = min(r[C_MINCHI] for r in rows)
    out["gov_cross"], out["gov_min"] = governor_cross(run)
    return out


def fnum(v, w=7, p=3):
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return " " * (w - 3) + "n/a"
    return f"{v:{w}.{p}f}"


def main() -> int:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tsv = "--tsv" in sys.argv
    ref_dir = None
    if "--ref" in sys.argv:
        ref_dir = Path(sys.argv[sys.argv.index("--ref") + 1])
        args = [a for a in args if Path(a) != ref_dir]
    if not args:
        print(__doc__)
        return 2

    ref = score(ref_dir) if ref_dir else None
    results = []
    for a in args:
        s = score(Path(a))
        if s is None:
            print(f"{a}: no confinement rows yet", file=sys.stderr)
            continue
        results.append(s)

    if tsv:
        keys = ["name", "t_end", "inj_canon", "inj_phan", "hold_canon",
                "hold_phan", "late_canon", "late_phan", "end_total",
                "end_canon", "end_phan", "min_chi", "gov_cross", "gov_min"]
        print("\t".join(keys))
        for s in results:
            print("\t".join(str(s.get(k, "n/a")) for k in keys))
        return 0

    hdr = (f"{'run':<14} {'t_end':>5} | {'inj_c':>6} {'inj_p':>6} | "
           f"{'late_c':>7} {'late_p':>7} | {'end_c':>6} {'end_p':>6}"
           f"{' | vs_ref_c vs_ref_p' if ref else ''} | {'gov<0.5':>7} {'min_chi':>8}")
    print(hdr)
    print("-" * len(hdr))
    for s in results:
        line = (f"{s['name']:<14} {s['t_end']:>5.1f} | "
                f"{fnum(s.get('inj_canon'), 6, 2)} {fnum(s.get('inj_phan'), 6, 2)} | "
                f"{fnum(s.get('late_canon'), 7, 4)} {fnum(s.get('late_phan'), 7, 4)} | "
                f"{fnum(s.get('end_canon'), 6, 1)} {fnum(s.get('end_phan'), 6, 1)}")
        if ref:
            rc = s.get("end_canon", math.nan) / ref["end_canon"] if ref.get("end_canon") else math.nan
            rp = s.get("end_phan", math.nan) / ref["end_phan"] if ref.get("end_phan") else math.nan
            line += f" | {fnum(rc, 8, 2)} {fnum(rp, 8, 2)}"
        line += f" | {s['gov_cross']:>7} {s['min_chi']:>8.2e}"
        print(line)
    print("\ninj_* = abs at first slice / t=0 (injection factor); late_* = decline"
          "\nrate per unit time over [20.2, 27.4] (positive = losing); end_* ="
          "\nabsolute confined activity at t_end. Baselines (tp30): inj_c 2.51,"
          "\ninj_p 1.05, late_p 0.0247, end_c 62.2*, end_p 22.9, gov cross 29.03."
          "\n(*canonical end contaminated past t~27.4 by the L2_Ham runaway.)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
