#!/usr/bin/env python3
"""Did the elliptic solve actually reach the tolerance it was asked for?

WHY THIS EXISTS.  ``BONDI_NL_TOL`` is a *request*, not a guarantee.  The solver
leaves its nonlinear loop by three different doors and only one of them means
the initial data is as good as it was asked to be:

  converged  both relative errors fell below ``NL_exit_tolerance``.  Good.
  stalled    the per-iteration improvement fell below ``NL_stall_tolerance``,
             so the loop gave up.  It can stop ANYWHERE above the tolerance.
  cap        ``max_NL_iterations`` ran out (default 50).  Same problem, and
             the solver does not even print a line saying so.

Nothing downstream checks which door was used.  The Python layer parses the
final residuals and logs "GRTresna converged: ..." unconditionally -- that
message is a label, not a verdict -- and the evolution then runs happily on
initial data that never met its constraint target.  A run built that way looks
completely normal: the star sits still, the lapse is flat, the constraints are
"small".  They are simply not as small as the number in the launch script says,
which quietly turns a convergence ladder into a measurement of the error floor.

The defaults make this easy to hit.  ``GRTresnaConfig`` ships
``nl_exit_tolerance = 1.0`` -- one percent -- because it was tuned for
MAP-Elites throughput, where thousands of cheap solves matter more than any
single one being tight.  Anything that does not override it inherits a search
campaign's tolerance into a paper run.

WHAT IS CHECKED.  The residual table the solve writes
(``grtresna/Ham_and_Mom_errors.txt``) against the tolerance the solve was
actually given (``grtresna/params.txt``).  Both files survive archiving; the
per-rank ``pout`` logs do not, which is why the verdict is reconstructed from
the table rather than scraped from the log.

Exit status is the gate: 0 = converged, 1 = did not.

    grteclyn-wrapper/.venv/bin/python research/bondi_dipole/check_solve_exit.py \
        runs/<campaign>/<cell>
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

VERDICT_OK = "converged"


def _find_grtresna(target: Path) -> Path | None:
    """Accept a cell dir, an episode dir, or the grtresna dir itself."""
    if (target / "Ham_and_Mom_errors.txt").exists():
        return target
    direct = target / "grtresna"
    if (direct / "Ham_and_Mom_errors.txt").exists():
        return direct
    for child in sorted(target.glob("*/grtresna")):
        if (child / "Ham_and_Mom_errors.txt").exists():
            return child
    return None


def _tolerances(params: Path) -> tuple[float | None, float | None, int | None]:
    if not params.exists():
        return None, None, None
    text = params.read_text()

    def _num(key: str) -> str | None:
        hit = re.search(rf"^\s*{key}\s*=\s*(\S+)", text, re.M)
        return hit.group(1) if hit else None

    exit_tol = _num("NL_exit_tolerance")
    stall_tol = _num("NL_stall_tolerance")
    cap = _num("max_NL_iterations")
    return (
        float(exit_tol) if exit_tol else None,
        float(stall_tol) if stall_tol else None,
        int(cap) if cap else None,
    )


def _residuals(errors: Path) -> list[tuple[int, float, float]]:
    rows = []
    for line in errors.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit():
            try:
                rows.append((int(parts[0]), float(parts[1]), float(parts[2])))
            except ValueError:
                continue
    return rows


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("target", type=Path, help="cell dir, episode dir, or grtresna dir")
    ap.add_argument(
        "--tolerance",
        type=float,
        default=None,
        help="override the tolerance to judge against (default: the one the "
        "solve was given, read from its own params.txt)",
    )
    args = ap.parse_args(argv)

    work = _find_grtresna(args.target)
    if work is None:
        print(f"FAIL  no grtresna/Ham_and_Mom_errors.txt under {args.target}")
        return 1

    exit_tol, stall_tol, cap = _tolerances(work / "params.txt")
    tol = args.tolerance if args.tolerance is not None else exit_tol
    rows = _residuals(work / "Ham_and_Mom_errors.txt")
    if not rows:
        print(f"FAIL  {work}/Ham_and_Mom_errors.txt has no residual rows")
        return 1

    last_iter, ham, mom = rows[-1]
    print(f"solve      {work}")
    print(f"tolerance  exit={exit_tol}  stall={stall_tol}  cap={cap}")
    print(f"final      iteration {last_iter}  Ham {ham:.3e}%  Mom {mom:.3e}%")

    if tol is None:
        print("FAIL  no NL_exit_tolerance recorded -- cannot judge this solve")
        return 1

    met = ham < tol and mom < tol
    if met and cap is not None and last_iter >= cap:
        # Landing on the last permitted iteration is met-by-luck, not converged.
        print(
            f"FAIL  reached the tolerance only on the final permitted iteration "
            f"({last_iter} of {cap}) -- raise max_NL_iterations and re-solve so "
            f"the exit is the tolerance, not the cap"
        )
        return 1
    if not met:
        reason = "ran out of iterations" if (cap and last_iter >= cap) else "stalled"
        print(
            f"FAIL  {reason} at Ham {ham:.3e}% / Mom {mom:.3e}%, above the "
            f"requested {tol:g}% -- this cell's initial data is looser than its "
            f"launch script claims; do not put it in a convergence ladder"
        )
        return 1

    headroom = tol / max(ham, mom)
    print(
        f"PASS  converged at iteration {last_iter}"
        + (f" of {cap}" if cap else "")
        + f", {headroom:.0f}x inside the requested {tol:g}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
