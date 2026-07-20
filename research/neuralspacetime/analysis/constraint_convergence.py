#!/usr/bin/env python3
"""Constraint-norm Richardson analysis for candidate-146 RC/RM/RF ladder.

Reads volume-weighted L2 Ham/Mom norms from the completed resolution-matrix
runs and:
  1. Emits a decimated three-resolution overlay table for pgfplots.
  2. Estimates the observed convergence order p(t) from the unequal-spacing
     Richardson relation already written in research.tex:
         (Q_c - Q_m) / (Q_m - Q_f) = (h_c^p - h_m^p) / (h_m^p - h_f^p)
  3. Prints summary maxima used in the paper.

Usage (from repo root):
  python research/neuralspacetime/analysis/constraint_convergence.py
"""

from __future__ import annotations

import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
PROMOTE = REPO / "runs" / "grtresna_promote"
ARTICLE_DATA = REPO / "research" / "neuralspacetime" / "article" / "data"

RUNS = {
    "RC": PROMOTE / "bcma_rc_L128_N192_t30_hq_eval000146",
    "RM": PROMOTE / "bcma_rm_L128_N256_t30_hq_eval000146",
    "RF": PROMOTE / "bcma_rf_L128_N384_t30_hq_eval000146",
}
H = {"RC": 128.0 / 192.0, "RM": 128.0 / 256.0, "RF": 128.0 / 384.0}


def load_constraints(path: Path) -> list[tuple[float, float, float]]:
    rows: list[tuple[float, float, float]] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            rows.append((float(parts[0]), float(parts[1]), float(parts[2])))
    return rows


def interpolate(
    series: list[tuple[float, float, float]], t_query: list[float]
) -> list[tuple[float, float, float]]:
    """Linear interpolate (ham, mom) onto a common time grid."""
    if not series:
        return []
    times = [r[0] for r in series]
    hams = [r[1] for r in series]
    moms = [r[2] for r in series]
    out: list[tuple[float, float, float]] = []
    j = 0
    for t in t_query:
        if t <= times[0]:
            out.append((t, hams[0], moms[0]))
            continue
        if t >= times[-1]:
            out.append((t, hams[-1], moms[-1]))
            continue
        while j + 1 < len(times) and times[j + 1] < t:
            j += 1
        t0, t1 = times[j], times[j + 1]
        a = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
        out.append(
            (
                t,
                (1 - a) * hams[j] + a * hams[j + 1],
                (1 - a) * moms[j] + a * moms[j + 1],
            )
        )
    return out


def richardson_order(qc: float, qm: float, qf: float, hc: float, hm: float, hf: float) -> float | None:
    """Solve (qc-qm)/(qm-qf) = (hc^p - hm^p)/(hm^p - hf^p) for p by bisection."""
    num = qc - qm
    den = qm - qf
    if abs(den) < 1e-30 or abs(num) < 1e-30:
        return None
    # Same-sign monotonicity required for a real positive order.
    if num * den <= 0:
        return None
    ratio = num / den

    def f(p: float) -> float:
        return (hc**p - hm**p) / (hm**p - hf**p) - ratio

    lo, hi = 0.05, 8.0
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        # Expand or give up.
        for hi_try in (10.0, 12.0):
            fhi = f(hi_try)
            if flo * fhi <= 0:
                hi = hi_try
                break
        else:
            return None
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if flo * fm <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def decimate(rows: list[tuple[float, float, float]], n_keep: int = 80) -> list[tuple[float, float, float]]:
    if len(rows) <= n_keep:
        return rows
    step = max(1, len(rows) // n_keep)
    kept = rows[::step]
    if kept[-1] != rows[-1]:
        kept.append(rows[-1])
    return kept


def main() -> None:
    series = {k: load_constraints(v / "data" / "constraint_norms.dat") for k, v in RUNS.items()}
    for k, rows in series.items():
        if not rows:
            raise SystemExit(f"missing constraint data for {k}: {RUNS[k]}")
        print(
            f"{k}: n={len(rows)}  t=[{rows[0][0]:.3f},{rows[-1][0]:.3f}]  "
            f"max||H||2={max(r[1] for r in rows):.4e}  "
            f"max||M||2={max(r[2] for r in rows):.4e}"
        )

    # Common time grid from RM (headline).
    t_full = [r[0] for r in series["RM"]]
    # Use ~0.1 cadence for order analysis.
    t_query = [t for t in t_full if abs(t * 10 - round(t * 10)) < 1e-9 or t in (t_full[0], t_full[-1])]
    if len(t_query) < 20:
        t_query = t_full[:: max(1, len(t_full) // 300)]

    interp = {k: interpolate(series[k], t_query) for k in RUNS}

    ARTICLE_DATA.mkdir(parents=True, exist_ok=True)
    overlay_path = ARTICLE_DATA / "constraints_resolution.txt"
    # Decimate RM times for typesetting; interpolate RC/RF onto those times.
    t_plot = [r[0] for r in decimate(series["RM"], n_keep=60)]
    plot = {k: interpolate(series[k], t_plot) for k in RUNS}
    with overlay_path.open("w") as fh:
        # First line is a pgfplots-readable header (no leading '#').
        fh.write("t ham_RC ham_RM ham_RF mom_RC mom_RM mom_RF\n")
        fh.write(
            "# Generated by research/neuralspacetime/analysis/constraint_convergence.py\n"
        )
        for i, t in enumerate(t_plot):
            fh.write(
                f"{t:.6f}  "
                f"{plot['RC'][i][1]:.6e}  {plot['RM'][i][1]:.6e}  {plot['RF'][i][1]:.6e}  "
                f"{plot['RC'][i][2]:.6e}  {plot['RM'][i][2]:.6e}  {plot['RF'][i][2]:.6e}\n"
            )
    print(f"wrote {overlay_path}")

    order_path = ARTICLE_DATA / "constraint_order.txt"
    orders_ham: list[float] = []
    orders_mom: list[float] = []
    with order_path.open("w") as fh:
        fh.write("t p_ham p_mom\n")
        fh.write("# Observed Richardson order from RC/RM/RF (unequal spacing)\n")
        for i, t in enumerate(t_query):
            hc, hm, hf = H["RC"], H["RM"], H["RF"]
            p_h = richardson_order(
                interp["RC"][i][1],
                interp["RM"][i][1],
                interp["RF"][i][1],
                hc,
                hm,
                hf,
            )
            p_m = richardson_order(
                interp["RC"][i][2],
                interp["RM"][i][2],
                interp["RF"][i][2],
                hc,
                hm,
                hf,
            )
            ph_s = f"{p_h:.4f}" if p_h is not None else "nan"
            pm_s = f"{p_m:.4f}" if p_m is not None else "nan"
            # Keep a moderate cadence for the figure.
            if i % max(1, len(t_query) // 80) == 0 or i == len(t_query) - 1:
                fh.write(f"{t:.6f}  {ph_s}  {pm_s}\n")
            if p_h is not None and 0.2 < p_h < 6.0:
                orders_ham.append(p_h)
            if p_m is not None and 0.2 < p_m < 6.0:
                orders_mom.append(p_m)
    print(f"wrote {order_path}")

    def stats(vals: list[float]) -> str:
        if not vals:
            return "n/a"
        vals_s = sorted(vals)
        med = vals_s[len(vals_s) // 2]
        return f"n={len(vals)} median={med:.2f} mean={sum(vals)/len(vals):.2f} [{min(vals):.2f},{max(vals):.2f}]"

    # Early window (driven / ignition phase) and free-evolution window.
    early_ham = [
        richardson_order(
            interpolate(series["RC"], [t])[0][1],
            interpolate(series["RM"], [t])[0][1],
            interpolate(series["RF"], [t])[0][1],
            H["RC"],
            H["RM"],
            H["RF"],
        )
        for t in t_query
        if t <= 4.0
    ]
    early_ham = [p for p in early_ham if p is not None and 0.2 < p < 6.0]
    late_ham = [
        richardson_order(
            interpolate(series["RC"], [t])[0][1],
            interpolate(series["RM"], [t])[0][1],
            interpolate(series["RF"], [t])[0][1],
            H["RC"],
            H["RM"],
            H["RF"],
        )
        for t in t_query
        if 4.0 < t <= 30.0
    ]
    late_ham = [p for p in late_ham if p is not None and 0.2 < p < 6.0]

    print("--- observed order summary ---")
    print(f"||H||2 all:   {stats(orders_ham)}")
    print(f"||M||2 all:   {stats(orders_mom)}")
    print(f"||H||2 t<=4:  {stats(early_ham)}")
    print(f"||H||2 t>4:   {stats(late_ham)}")

    # Peak values table for the paper.
    print("--- maxima ---")
    for k in ("RC", "RM", "RF"):
        rows = series[k]
        print(
            f"  {k}: max||H||2={max(r[1] for r in rows):.4e}  "
            f"max||M||2={max(r[2] for r in rows):.4e}  "
            f"final||H||2={rows[-1][1]:.4e}"
        )


if __name__ == "__main__":
    main()
