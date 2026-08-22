#!/usr/bin/env python3
"""Read a cell's launch script and refuse the ones that cannot be compared.

WHY THIS EXISTS (2026-08-22).  Every other check in this campaign judges a cell
on its own: is its grid aligned, did its solve converge, is its matter intact.
A cell can pass all of them and still be worthless, because the thing that makes
it worthless is a difference from the cell it will be *read against*.  That is
what happened to the all-canonical null control.

The specific trap, and why a dry run does not catch it.  GRTresna switches to
the maximal-slicing (K = 0) path BY ITSELF whenever any lump carries negative
energy -- the CTTK ansatz K = sign*sqrt(24 pi G rho) is imaginary for rho < 0,
so it has no choice.  Switching it also relaxes psi, sets a psi floor and
changes the coefficient averaging.  A command copied from any phantom-bearing
cell therefore inherits that whole path INVISIBLY: it is nowhere in the launch
script, because nobody asked for it.  Delete the phantom star and all four
settings vanish just as silently, and the canonical cell is born on an
already-collapsing slice.  It looks fine.  Its solve converges.  Its grid is
aligned.  It is simply not a null for the run it is supposed to be a null for.

And `BONDI_DRYRUN=1` cannot show you: the solve parameter file does not exist
until the solve runs, and the dry run's own metadata does not record the choice.
The first moment the truth is written down is after the CPU time is spent.

So this reads the launch script itself -- the only artefact that exists before
anything is spent -- and answers the between-cell question directly.

    grteclyn-wrapper/.venv/bin/python research/bondi_dipole/preflight_cell.py \
        runs/bondi/staging/<cell>/launch.sh \
        --reference runs/bondi/runaway_paper/runaway_pair_d10_L64_N128_lev0/launch.sh

Exit status is the gate: 0 = safe to launch, 1 = do not spend time on it.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# Knobs a cell is *supposed* to vary against its reference.  Anything else that
# differs is reported, because the whole point is that surprises are the bug.
EXPECTED_TO_DIFFER = {
    "BONDI_GPU",
    "BONDI_RUNS_DIR",
    "BONDI_S0",
    "BONDI_S1",
    "BONDI_S0_OMEGA",
    "BONDI_S1_OMEGA",
    "BONDI_EXOTIC",
    "BONDI_OMEGA",
    "BONDI_SEP",
    "BONDI_GRTRESNA_MAXIMAL_SLICING",
    # Drawing pictures changes nothing a number depends on.
    "GRTECLYN_FRAMES",
    "GRTECLYN_FRAMES_CACHE_SLICES",
}

ASSIGN = re.compile(r"\b(BONDI_[A-Z0-9_]+|GRTECLYN_[A-Z0-9_]+)=(\"[^\"]*\"|'[^']*'|\S+)")


def read_env(script: Path) -> dict[str, str]:
    """Every BONDI_/GRTECLYN_ assignment in the file, last one winning."""
    env: dict[str, str] = {}
    for name, raw in ASSIGN.findall(script.read_text()):
        env[name] = raw.strip("\"'").rstrip("\\").strip()
    return env


def has_phantom(env: dict[str, str]) -> bool:
    """True when the cell contains at least one negative-energy star.

    Sector 1 is the phantom sector.  The pair script defaults to a mixed pair
    (S0=0, S1=1); the single-star script uses BONDI_EXOTIC.
    """
    if "BONDI_EXOTIC" in env:
        return env["BONDI_EXOTIC"] not in ("0", "")
    s0 = env.get("BONDI_S0", "0")
    s1 = env.get("BONDI_S1", "1")
    return "1" in (s0, s1)


def check(env: dict[str, str], where: str) -> list[str]:
    problems = []

    # 1. The between-cell trap this script exists for.
    if not has_phantom(env) and env.get("BONDI_GRTRESNA_MAXIMAL_SLICING", "0") != "1":
        problems.append(
            f"{where}: every star here is canonical, so nothing forces the flat "
            f"K = 0 start that every phantom-bearing cell gets for free. This "
            f"cell would be born on an already-collapsing slice and could not be "
            f"read against one. Set BONDI_GRTRESNA_MAXIMAL_SLICING=1."
        )

    # 2. The solve grid must land on the evolution grid (README rule 1).
    try:
        nfull = int(env.get("BONDI_NFULL", "128"))
        lfull = float(env.get("BONDI_LFULL", "64"))
        domain = float(env.get("BONDI_GRTRESNA_DOMAIN_L", "128"))
        want = nfull * (domain / lfull)
        got = float(env.get("BONDI_GRTRESNA_N", nfull))
        if abs(got - want) > 1e-9:
            problems.append(
                f"{where}: the solve grid does not land on the evolution grid. "
                f"With N={nfull} in a box of {lfull:g} solved in a box of "
                f"{domain:g}, the solve needs {want:g} cells, not {got:g}. The "
                f"metric would arrive displaced by a fraction of a cell."
            )
    except ValueError:
        problems.append(f"{where}: grid knobs are not numbers; cannot check alignment")

    if env.get("BONDI_GRTRESNA_MAXLEVEL", "3") != "0":
        problems.append(
            f"{where}: the solve is refined (BONDI_GRTRESNA_MAXLEVEL != 0), which "
            f"puts interpolation back into the transfer path"
        )

    # 3. The solve tolerance must be asked for (README rule 8).
    if "BONDI_NL_TOL" not in env:
        problems.append(
            f"{where}: no BONDI_NL_TOL, so the solve inherits the shipped default "
            f"of one percent -- a MAP-Elites throughput setting, not a paper one"
        )
    return problems


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("launch", type=Path, help="the cell's launch.sh")
    ap.add_argument(
        "--reference",
        type=Path,
        default=None,
        help="launch.sh of the cell this one will be read against; every knob "
        "that differs and is not expected to differ is reported",
    )
    args = ap.parse_args(argv)

    env = read_env(args.launch)
    if not env:
        print(f"FAIL  no BONDI_ assignments found in {args.launch}")
        return 1

    problems = check(env, args.launch.parent.name or str(args.launch))

    if args.reference:
        ref = read_env(args.reference)
        ref_name = args.reference.parent.name
        surprises = []
        for key in sorted(set(env) | set(ref)):
            if key in EXPECTED_TO_DIFFER:
                continue
            a, b = env.get(key), ref.get(key)
            if a != b:
                surprises.append(f"    {key}: this={a!r}  {ref_name}={b!r}")
        if surprises:
            problems.append(
                f"{args.launch.parent.name}: differs from {ref_name} in knobs that "
                f"were not meant to vary --\n" + "\n".join(surprises)
            )
        # The reference's slicing path must be reachable by this cell.
        if has_phantom(ref) != has_phantom(env):
            forced = env.get("BONDI_GRTRESNA_MAXIMAL_SLICING", "0") == "1"
            if not forced:
                problems.append(
                    f"{args.launch.parent.name}: this cell and {ref_name} do not "
                    f"contain the same kind of matter, so they do not get the same "
                    f"slicing by default -- force it or they are not comparable"
                )

    for p in problems:
        print(f"FAIL  {p}")
    if problems:
        return 1
    print(f"PASS  {args.launch.parent.name or args.launch}: safe to launch")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
