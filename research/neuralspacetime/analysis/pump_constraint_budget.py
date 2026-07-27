#!/usr/bin/env python3
"""Duhamel upper bound on controller-attributed constraint violation.

Reads ``constraint_norms.dat`` columns appended by the always-on-pump
diagnostics:

    time, L2_Ham, L2_Mom, min_rho_req, max_rho_req, integral_neg_rho,
    [L2_Ham_rel, L2_Mom_rel, pump_force_L2, governor, pump_fi_L2]

and forms the integral bounds

    ||H_pump||_2 ≲ 16 π ∫_0^t ||f_⊥||_2 dt'   from pump_force_L2
    ||M_pump||_2 ≲  8 π ∫_0^t ||f_i||_2 dt'   from pump_fi_L2

Both force norms are measured with the same law and the same centre-relative
coordinates as the evolution RHS, and carry no lapse factor: the pump adds
S_A to ∂_t Π_A, a coordinate-time rate, so ∂_t ρ|_pump = f_⊥ exactly.  (Runs
predating the pump_fi_L2 column fall back to pump_force_L2 for the momentum
bound and say so.)

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


def load_constraint_norms(path: Path) -> dict[str, object]:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    ncols = data.shape[1]
    out = {
        "time": data[:, 0],
        "L2_Ham": data[:, 1],
        "L2_Mom": data[:, 2],
    }
    if ncols < 9:
        raise SystemExit(
            f"{path}: need >= 9 columns (got {ncols}); re-run with the "
            "extended constraint diagnostics."
        )
    out["L2_Ham_rel"] = data[:, 6]
    out["L2_Mom_rel"] = data[:, 7]
    out["pump_force_L2"] = data[:, 8]
    # tolerate runs written before the governor / pump_fi_L2 columns existed
    out["governor"] = data[:, 9] if ncols >= 10 else np.ones_like(data[:, 0])
    if ncols >= 11:
        out["pump_fi_L2"] = data[:, 10]
        out["fi_measured"] = True
    else:
        out["pump_fi_L2"] = out["pump_force_L2"]
        out["fi_measured"] = False
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
    force_i = cols["pump_fi_L2"]
    ham_bound = 16.0 * math.pi * cumulative_trapz(force, t)
    mom_bound = 8.0 * math.pi * cumulative_trapz(force_i, t)

    gov = cols["governor"]
    print(f"file: {args.constraint_norms}")
    print(f"rows: {len(t)}  t in [{t[0]:.3f}, {t[-1]:.3f}]")
    print(f"pump_force_L2: max={force.max():.3e}  mean={force.mean():.3e}")
    if cols["fi_measured"]:
        print(f"pump_fi_L2:    max={force_i.max():.3e}  mean={force_i.mean():.3e}")
    else:
        print("pump_fi_L2:    NOT MEASURED in this run -- momentum bound uses "
              "pump_force_L2 as a stand-in and is NOT a bound on ||M_pump||.")
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
                force_i,
            ]
        )
        header = (
            "time L2_Ham ham_bound L2_Mom mom_bound governor "
            "pump_force_L2 pump_fi_L2"
        )
        np.savetxt(args.output, table, header=header)
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
