#!/usr/bin/env python3
"""Compare the exotic-matter ("negative energy") budget of a discovered candidate
against an analytic Alcubierre warp bubble, on ONE common measure::

    E_exotic = | int_{rho<0} rho sqrt(det gamma) d^3x |        (geometric units)

with ``rho = T_{mu nu} n^mu n^nu`` the Eulerian energy density and ``T = G/8pi``.

Both sides go through ``grteclyn_wrapper.metrics.probes.warpfactory``
(:func:`stress_energy` + :func:`exotic_energy_budget`), so the comparison answers
"how much *less* exotic matter does our self-consistent candidate need than a
textbook warp drive?" apples-to-apples.

The Alcubierre side is built analytically (:func:`alcubierre_metric`).  The
candidate side is read from a finished eval directory:

  * ``score.json``            -> min_rho_required, energy_conditions.min_nec
  * ``data/constraint_norms.dat`` -> the integral_negative_rho profile
    (cols: t, ham_l2, mom_l2, min_rho, max_rho, integral_negative_rho), reported
    at t=0 (assembled), its peak, and the run-max.  NOTE the run-max can be
    inflated by late-time numerical drift, so the t=0..warp-peak range is the
    honest "exotic matter that supports the warp".

Examples
--------
    # Alcubierre control only (also a sanity check vs the documented 31.5% / -0.47):
    uv run python scripts/validation/exotic_energy_compare.py --v-s 2.0 --radius 4.0

    # Compare against eval 177:
    uv run python scripts/validation/exotic_energy_compare.py --v-s 2.0 --radius 4.0 \\
        --candidate-eval runs/grtresna_promote/l128n256t30_ftl_cmaes_v17_robust_qd_eval000177
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from grteclyn_wrapper.metrics.probes import warpfactory as wf


def alcubierre_budget(*, v_s: float, radius: float, sigma: float, n: int,
                      half_width: float, dt: float, crop: int) -> dict:
    """Exotic-energy budget + min NEC of an Alcubierre bubble (analytic)."""
    g, spacing = wf.alcubierre_metric(
        velocity=v_s, bubble_radius=radius, sigma=sigma,
        half_width=half_width, n_space=n, dt=dt,
    )
    budget = wf.exotic_energy_budget(g, spacing, crop=crop)
    ec = wf.evaluate_four_metric(g, spacing, crop=crop)
    return {
        "exotic_energy": budget.exotic_energy,
        "min_rho": budget.min_rho,
        "negative_fraction": budget.negative_fraction,
        "min_nec": ec.nec_min,
        "n_points": budget.n_points,
    }


def _constraint_profile(eval_dir: Path) -> dict | None:
    dat = eval_dir / "data" / "constraint_norms.dat"
    if not dat.is_file():
        return None
    a = np.loadtxt(dat)
    if a.ndim != 2 or a.shape[1] < 6:
        return None
    t, int_neg = a[:, 0], a[:, 5]
    return {
        "int_neg_t0": float(int_neg[0]),
        "int_neg_max": float(int_neg.max()),
        "t_at_max": float(t[int_neg.argmax()]),
        "int_neg_final": float(int_neg[-1]),
    }


def candidate_budget(eval_dir: Path) -> dict:
    """Read a candidate's exotic-matter diagnostics from a finished eval dir.

    These are GRTeclyn's evolved constraint diagnostics (integral of negative
    rho_required = (R+K^2-K_ij K^ij)/16pi), the *same* physical quantity as the
    Eulerian rho integrated for Alcubierre, computed during the evolution.
    """
    score = json.loads((eval_dir / "score.json").read_text(encoding="utf-8"))
    m = score.get("metrics", {})
    constraints = m.get("constraints", {})
    ec = m.get("energy_conditions", {})
    out = {
        "min_rho_required": constraints.get("min_rho_required"),
        "integral_negative_rho_runmax": constraints.get("integral_negative_rho"),
        "min_nec": ec.get("min_nec"),
    }
    profile = _constraint_profile(eval_dir)
    if profile:
        out.update(profile)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v-s", type=float, default=2.0, help="Alcubierre bubble speed")
    ap.add_argument("--radius", type=float, default=4.0, help="Alcubierre bubble radius R")
    ap.add_argument("--sigma", type=float, default=2.0, help="Alcubierre wall steepness")
    ap.add_argument("--n", type=int, default=97, help="grid points per axis")
    ap.add_argument("--half-width", type=float, default=10.0)
    ap.add_argument("--dt", type=float, default=0.1, help="time spacing for d_t terms")
    ap.add_argument("--crop", type=int, default=3, help="boundary cells dropped per axis")
    ap.add_argument("--candidate-eval", type=Path, default=None,
                    help="finished eval dir with score.json + data/constraint_norms.dat")
    args = ap.parse_args()

    alc = alcubierre_budget(
        v_s=args.v_s, radius=args.radius, sigma=args.sigma,
        n=args.n, half_width=args.half_width, dt=args.dt, crop=args.crop,
    )
    print("=" * 72)
    print(f"Alcubierre control  v_s={args.v_s}  R={args.radius}  sigma={args.sigma}  "
          f"{args.n}^3")
    print("=" * 72)
    print(f"  total exotic energy |E_neg| : {alc['exotic_energy']:.3f}")
    print(f"  min energy density (min rho): {alc['min_rho']:+.4f}")
    print(f"  min NEC                     : {alc['min_nec']:+.4f}")
    print(f"  negative-energy cell frac   : {alc['negative_fraction']:.2%}")

    if args.candidate_eval is None:
        return 0

    cand = candidate_budget(args.candidate_eval.expanduser().resolve())
    print()
    print("=" * 72)
    print(f"Candidate  {args.candidate_eval.name}")
    print("=" * 72)
    e_t0 = cand.get("int_neg_t0")
    e_pk = cand.get("int_neg_max")
    print(f"  exotic energy @ t=0 (built) : {e_t0}")
    print(f"  exotic energy run-max (t={cand.get('t_at_max')}) : {e_pk}  "
          f"(may include late-time numerical drift)")
    print(f"  min rho_required            : {cand.get('min_rho_required')}")
    print(f"  min NEC                     : {cand.get('min_nec')}")

    print()
    print("-" * 72)
    print("RATIO  (Alcubierre / candidate) -- 'how many times more exotic matter'")
    print("-" * 72)
    alc_E = alc["exotic_energy"]
    if e_t0:
        print(f"  total exotic energy (vs t=0 build) : {alc_E / e_t0:5.1f}x")
    if e_pk:
        print(f"  total exotic energy (vs run-max)   : {alc_E / e_pk:5.1f}x")
    if cand.get("min_nec"):
        print(f"  pointwise NEC intensity            : {alc['min_nec'] / cand['min_nec']:5.0f}x milder")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
