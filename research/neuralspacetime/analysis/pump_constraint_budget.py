#!/usr/bin/env python3
"""Duhamel upper bound on controller-attributed constraint violation.

Reads ``constraint_norms.dat`` columns appended by the always-on-pump
diagnostics:

    time, L2_Ham, L2_Mom, min_rho_req, max_rho_req, integral_neg_rho,
    [L2_Ham_rel, L2_Mom_rel, pump_force_L2, governor]

and forms the integral bounds

    ||H_pump||_2 ≲ 16 π ∫_0^t ||α f_⊥||_2 dt'
    ||M_pump||_2 ≲  8 π ∫_0^t ||α f_i||_2 dt'   (proxy: same pump_force_L2)

These are upper bounds, not identities: L² norms do not add linearly through
the constraint-propagation system. Exact cancellation requires the mode-1/2
controller reservoir (see ControllerReservoirMatter).

Usage:
    python pump_constraint_budget.py path/to/run/data/constraint_norms.dat
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np


def load_constraint_norms(path: Path) -> dict[str, np.ndarray]:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    ncols = data.shape[1]
    out = {
        "time": data[:, 0],
        "L2_Ham": data[:, 1],
        "L2_Mom": data[:, 2],
    }
    if ncols >= 10:
        out["L2_Ham_rel"] = data[:, 6]
        out["L2_Mom_rel"] = data[:, 7]
        out["pump_force_L2"] = data[:, 8]
        out["governor"] = data[:, 9]
    elif ncols >= 9:
        # tolerate missing governor
        out["L2_Ham_rel"] = data[:, 6]
        out["L2_Mom_rel"] = data[:, 7]
        out["pump_force_L2"] = data[:, 8]
        out["governor"] = np.ones_like(data[:, 0])
    else:
        raise SystemExit(
            f"{path}: need >= 9 columns (got {ncols}); re-run with the "
            "extended constraint diagnostics."
        )
    return out


def cumulative_trapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Inclusive cumulative trapezoid; length matches ``x`` (leading zero)."""
    out = np.zeros_like(x)
    for i in range(1, len(x)):
        out[i] = out[i - 1] + 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1])
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("constraint_norms", type=Path)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Optional TSV of t, L2_Ham, ham_bound, L2_Mom, mom_bound, governor",
    )
    args = parser.parse_args()

    cols = load_constraint_norms(args.constraint_norms)
    t = cols["time"]
    force = cols["pump_force_L2"]
    integ = cumulative_trapz(force, t)
    ham_bound = 16.0 * math.pi * integ
    mom_bound = 8.0 * math.pi * integ

    gov = cols["governor"]
    print(f"file: {args.constraint_norms}")
    print(f"rows: {len(t)}  t in [{t[0]:.3f}, {t[-1]:.3f}]")
    print(f"pump_force_L2: max={force.max():.3e}  mean={force.mean():.3e}")
    print(f"governor: min={gov.min():.3f}  max={gov.max():.3f}")
    print(
        f"at t={t[-1]:.2f}: L2_Ham={cols['L2_Ham'][-1]:.3e}  "
        f"ham_bound={ham_bound[-1]:.3e}  "
        f"L2_Mom={cols['L2_Mom'][-1]:.3e}  mom_bound={mom_bound[-1]:.3e}"
    )
    # Peak Ham during early window vs bound
    early = t <= 4.0
    if np.any(early):
        i_peak = int(np.argmax(cols["L2_Ham"][early]))
        print(
            f"early-window peak L2_Ham={cols['L2_Ham'][early][i_peak]:.3e} "
            f"at t={t[early][i_peak]:.2f}; "
            f"bound there={ham_bound[early][i_peak]:.3e}"
        )

    if args.output is not None:
        table = np.column_stack(
            [
                t,
                cols["L2_Ham"],
                ham_bound,
                cols["L2_Mom"],
                mom_bound,
                gov,
                force,
            ]
        )
        header = "time L2_Ham ham_bound L2_Mom mom_bound governor pump_force_L2"
        np.savetxt(args.output, table, header=header)
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
