#!/usr/bin/env python3
"""Is R(omega) monotonic for the dressed sextic stars at the ultraweak rung?

The Discussion claims a cleaner constant-gap test could be built from "more
compact stars (higher omega)".  That is only true if the family radius really
decreases with frequency.  Sextic (solitonic) Q-balls are famously *not*
monotonic in R(omega): the radius diverges in both the thin-wall limit
(omega -> omega_min, the ball inflates to hold its surface tension) and the
thick-wall limit (omega -> m, the field spreads out as the binding vanishes),
so R(omega) has a minimum somewhere in between.  If the working frequency
omega = 0.55 sits on the *thin-wall* side of that minimum, "higher omega" is
the right advice; if it sits past it, the sentence is backwards.

Separate script and separate output file on purpose: ``star_family_scan.py``
writes the published ``stars/star_family.csv`` whose columns back Table/Fig.
"family", and that contract is not widened here.

Radii reported (all from the same solved profile):
  * ``r99_iso`` / ``r90_iso`` -- isotropic radius enclosing 99% / 90% of the
    integrated |phi0| r^2.  ``r99_iso`` is exactly the quantity behind the
    ``compactness`` column of ``star_family.csv``, reproduced here so the two
    files can be cross-checked.
  * ``r99_areal`` -- the physical (circumferential) radius psi^2 r at r99_iso.
    This is the one a compactness statement should quote.
  * ``compactness_areal`` = 2 |M_ADM| / r99_areal.  Absolute value because the
    phantom sector carries M < 0; the *size* comparison is the point.

Run with the wrapper venv (needs numpy/scipy):
    PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
        results/bondi-dipole-runaway/analysis/star_radius_scan.py
"""

from __future__ import annotations

import csv
import pathlib
import sys

# Same rung as the campaign: lambda^2/mu = 4.8 => omega_min = 0.316.
MASS, LAM, MU = 1.0, 10240.0, 21845333.0

# Denser than star_family_scan.py's grid: a minimum in R(omega) is easy to step
# over with 14 points.  0.52 is known to have no dressed solution in either
# sector, so the sweep starts above it and runs to the thick-wall end.
OMEGAS = [
    0.530, 0.535, 0.540, 0.545, 0.550, 0.555, 0.560, 0.56598, 0.570,
    0.575, 0.580, 0.590, 0.600, 0.615, 0.630, 0.650, 0.670, 0.690,
    0.710, 0.730, 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900,
]

# The two frequencies the paper actually uses, for the verdict line.
OMEGA_CANON, OMEGA_PHANTOM_EQM = 0.550, 0.56598


def _radii(star) -> dict[str, float]:
    """Enclosed-fraction radii of one solved profile, isotropic and areal."""
    import numpy as np

    r = np.asarray(star.r, dtype=float)
    integrand = np.abs(star.phi0) * r**2
    cumulative = np.cumsum(integrand) * np.gradient(r)
    total = float(cumulative[-1]) if cumulative.size else 1.0

    def enclosing(frac: float) -> float:
        idx = int(np.searchsorted(cumulative, frac * total))
        return float(r[min(idx, r.size - 1)])

    r99, r90 = enclosing(0.99), enclosing(0.90)
    psi99 = float(np.interp(r99, r, np.asarray(star.psi, dtype=float)))
    return {
        "r99_iso": r99,
        "r90_iso": r90,
        # Conformally flat data: areal radius = psi^2 * r_isotropic.
        "r99_areal": psi99**2 * r99,
    }


#: A step in R smaller than this fraction of R is treated as flat, not as a
#: trend reversal.  R_99 is read off a *grid point* of the ODE solution
#: (``r_full[searchsorted(...)]``), so it is quantised: neighbouring frequencies
#: can land on the same or an adjacent node and produce a step of ~1e-3 in
#: absolute terms with no physical content.  Without this guard the detector
#: reported a spurious minimum/maximum pair at omega = 0.575/0.580 from a
#: 0.01% wiggle, which would have moved the reported location of the real
#: radius minimum by 0.04 in omega.
TURN_REL_TOL = 2.0e-3


def _turning_points(
    xs: list[float], ys: list[float], rel_tol: float = TURN_REL_TOL
) -> list[tuple[float, str]]:
    """Interior sign changes of the discrete derivative of y(x).

    Steps below ``rel_tol`` (relative) are treated as flat and merged into the
    surrounding trend, so grid quantisation in R does not manufacture turns.
    """
    out: list[tuple[float, str]] = []
    trend = 0  # -1 falling, +1 rising, 0 not yet established
    for i in range(1, len(xs)):
        step = ys[i] - ys[i - 1]
        scale = max(abs(ys[i]), abs(ys[i - 1]), 1e-300)
        if abs(step) / scale < rel_tol:
            continue  # flat within quantisation -- carry the previous trend
        sign = 1 if step > 0 else -1
        if trend != 0 and sign != trend:
            out.append((xs[i - 1], "minimum" if sign > 0 else "maximum"))
        trend = sign
    return out


# r90 matters as much as r99 here and turns over at a different frequency: r99
# includes the exponential tail, which starts inflating (thick-wall limit) while
# the field core is still shrinking.  Reporting only r99 understates how much
# the *core* can be shrunk by raising omega (5% vs 32%).
CURVE_KEYS = ("r90_iso", "r99_iso", "r99_areal", "compactness_areal")


def _curves_from_rows(rows: list[dict]) -> dict[str, dict[str, list[float]]]:
    curves: dict[str, dict[str, list[float]]] = {
        s: {"omega": [], **{k: [] for k in CURVE_KEYS}}
        for s in ("canonical", "phantom")
    }
    for row in rows:
        if str(row.get("status", "")) != "ok":
            continue
        sec = str(row["sector"])
        curves[sec]["omega"].append(float(row["omega"]))
        for key in CURVE_KEYS:
            curves[sec][key].append(float(row[key]))
    return curves


def _solve_rows(argv_unused=None) -> list[dict[str, object]]:
    from grteclyn_wrapper.grtresna.profiles.boson_star_ode import (
        cached_selfgrav_at_omega,
    )

    rows: list[dict[str, object]] = []
    for omega in OMEGAS:
        for sector, sign in (("canonical", 1.0), ("phantom", -1.0)):
            try:
                star = cached_selfgrav_at_omega(MASS, LAM, MU, omega, sign)
            except Exception as exc:  # branch edge: no dressed star here
                rows.append(
                    {
                        "sector": sector,
                        "omega": omega,
                        "status": f"no solution: {exc}",
                        "phi_c": "",
                        "adm_mass": "",
                        "r99_iso": "",
                        "r90_iso": "",
                        "r99_areal": "",
                        "compactness_iso": "",
                        "compactness_areal": "",
                    }
                )
                print(f"[radius] {sector:9s} omega={omega:.5f}  NO SOLUTION")
                continue

            rad = _radii(star)
            m_abs = abs(float(star.adm_mass))
            comp_iso = 2.0 * m_abs / rad["r99_iso"] if rad["r99_iso"] > 0 else 0.0
            comp_areal = (
                2.0 * m_abs / rad["r99_areal"] if rad["r99_areal"] > 0 else 0.0
            )
            rows.append(
                {
                    "sector": sector,
                    "omega": omega,
                    "status": "ok",
                    "phi_c": f"{star.phi_c:.9f}",
                    "adm_mass": f"{star.adm_mass:.9f}",
                    "r99_iso": f"{rad['r99_iso']:.6f}",
                    "r90_iso": f"{rad['r90_iso']:.6f}",
                    "r99_areal": f"{rad['r99_areal']:.6f}",
                    "compactness_iso": f"{comp_iso:.6e}",
                    "compactness_areal": f"{comp_areal:.6e}",
                }
            )
            print(
                f"[radius] {sector:9s} omega={omega:.5f}  R99_iso={rad['r99_iso']:8.4f}"
                f"  R99_areal={rad['r99_areal']:8.4f}  |M|={m_abs:.6f}"
                f"  2|M|/R={comp_areal:.4e}"
            )
    return rows


def _report(curves: dict[str, dict[str, list[float]]]) -> None:
    """Print the verdict the Discussion sentence needs."""
    print("\n=== R(omega) monotonicity ===")
    for sector in ("canonical", "phantom"):
        om = curves[sector]["omega"]
        if len(om) < 3:
            print(f"{sector:9s}: too few solved points ({len(om)}) to judge")
            continue
        for key in CURVE_KEYS:
            ys = curves[sector][key]
            turns = _turning_points(om, ys)
            trend = "decreasing" if ys[-1] < ys[0] else "increasing"
            if not turns:
                print(
                    f"{sector:9s} {key:18s}: MONOTONIC ({trend}) over "
                    f"{om[0]:.3f}-{om[-1]:.3f}  "
                    f"[{ys[0]:.4g} -> {ys[-1]:.4g}, factor {ys[0] / ys[-1]:.1f}]"
                )
                continue
            where = ", ".join(f"{w:.4f} ({k})" for w, k in turns)
            print(f"{sector:9s} {key:18s}: NON-MONOTONIC -- turning at {where}")
            mins = [w for w, k in turns if k == "minimum"]
            if not mins or key == "compactness_areal":
                continue
            first_min = min(mins)
            for probe, label in (
                (OMEGA_CANON, "canonical run  omega=0.550   "),
                (OMEGA_PHANTOM_EQM, "phantom-eqm    omega=0.56598"),
            ):
                side = "SHRINKS" if probe < first_min else "GROWS"
                # How much is actually on the table before the turn?
                r_here = float(_interp(om, ys, probe))
                r_min = min(y for x, y in zip(om, ys) if x <= first_min)
                gain = 100.0 * (r_here - r_min) / r_here
                print(
                    f"{'':9s} {'':18s}  -> raising omega from {label} {side} it; "
                    f"best case R {r_here:.3f} -> {r_min:.3f} "
                    f"({gain:+.1f}%) before the minimum at omega={first_min:.3f}"
                )


def _interp(xs: list[float], ys: list[float], x: float) -> float:
    for i in range(1, len(xs)):
        if xs[i] >= x:
            f = (x - xs[i - 1]) / (xs[i] - xs[i - 1]) if xs[i] != xs[i - 1] else 0.0
            return ys[i - 1] + f * (ys[i] - ys[i - 1])
    return ys[-1]


def main(argv: list[str]) -> int:
    # ``--reuse`` re-derives the verdict from an existing CSV.  The 54 shooting
    # solves cost ~15 min; the verdict costs nothing, so iterating on the
    # analysis must not require re-solving the family.
    reuse = "--reuse" in argv
    positional = [a for a in argv if not a.startswith("-")]
    pack = (
        pathlib.Path(positional[0])
        if positional
        else pathlib.Path(__file__).resolve().parent.parent
    )
    out = pack / "stars" / "star_radius.csv"
    out.parent.mkdir(parents=True, exist_ok=True)

    if reuse and out.is_file():
        with out.open(newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        print(f"[radius] reusing {out} ({len(rows)} rows) -- no solves")
    else:
        rows = _solve_rows()
        with out.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"\n[radius] wrote {out} ({len(rows)} rows)")

    _report(_curves_from_rows(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
