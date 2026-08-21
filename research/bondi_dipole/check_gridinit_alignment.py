#!/usr/bin/env python3
"""Is the metric in a gridinit file centred on the same place as the matter?

WHY THIS EXISTS (measured 2026-08-21).  The elliptic solve runs on its own
grid and the answer is copied cell-by-cell onto the evolution grid by
``grtresna/io/chombo.py::_target_span_slice``.  That copy is piecewise
constant, and where a solve cell is FINER than an evolution cell the last
source cell to touch a target cell wins -- always the one above it.  The metric
therefore lands displaced downward by a fraction of a cell, by the same amount
on all three axes because the arithmetic is identical per axis.

The matter does not move with it.  ``io/conversion.py`` repaints the lumps
analytically on the target grid, exact to machine precision.  So each star is
born sitting slightly off the centre of its own gravitational well.  A
canonical star falls back toward that centre; a phantom star, sitting on a hill
rather than in a well, is pushed the other way.  That sign flip is the
signature -- an external field would accelerate both sectors identically.

WHAT IS MEASURED.  A pair configuration is exactly symmetric under y -> C-y and
z -> C-z about the axis joining the stars, so both the matter and the metric
must centroid at exactly C on those axes.  The matter always does.  Any offset
in the metric is the bug, and it is read off directly rather than inferred.

The window is snapped to an exactly symmetric set of cell centres; get that
wrong by half a cell and the estimator invents an offset of its own (which is
how this was first mis-measured).  The assertion guards it.

THE FIX this checks for: set the solve cells to
``N_full * (solve_box / evolution_box)`` with no solve refinement, so the solve
spacing equals the evolution spacing and the transfer is a straight copy.
Offsets then come back at ~1e-3 or below instead of ~1e-1.

Usage:
    check_gridinit_alignment.py <initial_data.gridinit> [--centre 32.0]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

#: Components carrying each sector's scalar amplitude.  The bicomplex model
#: splits one star across a real and an imaginary slot and leaves the partner
#: slot identically zero, so squaring and summing the candidates is safe.
SECTOR_FIELDS: dict[str, tuple[str, ...]] = {
    "canonical": ("phi", "Pi", "phi2", "Pi2"),
    "phantom": ("phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2"),
}


def read_header(path: Path) -> tuple[int, list[str], int, float, np.ndarray]:
    """Return (payload offset, component names, n, dx, origin)."""
    with path.open("rb") as handle:
        raw = b""
        while not raw.endswith(b"END_HEADER\n"):
            chunk = handle.read(1)
            if not chunk:
                raise ValueError(f"no END_HEADER in {path}")
            raw += chunk
        offset = handle.tell()
    text = raw.decode()
    if not text.startswith("GRTECLYN_GRID_INIT_V2"):
        raise ValueError(f"not a GRTECLYN_GRID_INIT_V2 file: {path}")
    fields: dict[str, list[str]] = {}
    for line in text.splitlines()[1:]:
        if not line.strip() or line == "END_HEADER":
            continue
        key, _, rest = line.partition(" ")
        fields[key] = rest.split()
    names = fields["component_names"]
    n = int(fields["nx_ny_nz"][0])
    dx = float(fields["dx"][0])
    origin = np.array([float(v) for v in fields["origin"]])
    return offset, names, n, dx, origin


def symmetric_window(centre: float, half_cells: int, dx: float, origin: float,
                     n: int) -> slice:
    """Cells whose centres are placed symmetrically about ``centre``.

    Off by half a cell here and the estimator reports an offset for a perfectly
    symmetric field, so the result is asserted rather than trusted.
    """
    first_above = int(np.floor((centre - origin) / dx))
    window = slice(first_above - half_cells, first_above + half_cells)
    if window.start < 0 or window.stop > n:
        raise ValueError(f"window {half_cells} cells wide falls outside the grid")
    coords = origin + (np.arange(n) + 0.5) * dx
    mid = 0.5 * (coords[window][0] + coords[window][-1])
    if abs(mid - centre) > 1.0e-9:
        raise AssertionError(f"window midpoint {mid} is not {centre}")
    return window


def centroid_along(weight: np.ndarray, coords: np.ndarray, axis: int) -> float:
    """Centroid of a 3-D non-negative weight along one array axis."""
    other = tuple(a for a in range(3) if a != axis)
    profile = weight.sum(axis=other)
    total = profile.sum()
    if total <= 0:
        return float("nan")
    return float((profile * coords).sum() / total)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gridinit", type=Path)
    ap.add_argument("--centre", type=float, default=32.0,
                    help="transverse centre of the configuration (default 32)")
    ap.add_argument("--half-width", type=float, default=4.0,
                    help="window half-width in code units (default 4)")
    ap.add_argument("--tolerance", type=float, default=5.0e-3,
                    help="offset above which the transfer is called biased")
    args = ap.parse_args(argv)

    offset, names, n, dx, origin = read_header(args.gridinit)
    # Array axes are (z, y, x); the transverse axes are y (1) and z (0).
    data = np.memmap(args.gridinit, dtype="<f8", mode="r",
                     offset=offset).reshape(n, n, n, len(names))
    index = {name: i for i, name in enumerate(names)}
    coords = origin[0] + (np.arange(n) + 0.5) * dx

    print(f"{args.gridinit}")
    print(f"  grid {n}^3   spacing {dx:.6f}   origin {origin[0]:g}")

    half = int(round(args.half_width / dx))
    win = symmetric_window(args.centre, half, dx, origin[0], n)
    chi = np.abs(np.asarray(data[win, win, :, index["chi"]]) - 1.0)

    worst = 0.0
    for sector, candidates in SECTOR_FIELDS.items():
        present = [c for c in candidates if c in index]
        rho = np.zeros(chi.shape)
        for name in present:
            rho = rho + np.asarray(data[win, win, :, index[name]]) ** 2
        if rho.max() <= 0:
            continue
        # Restrict to this sector's own half of the box so the companion's
        # metric feature does not drag the centroid.
        peak_x = int(np.argmax(rho.sum(axis=(0, 1))))
        lo, hi = max(0, peak_x - 2 * half), min(n, peak_x + 2 * half)
        rho_s, chi_s = rho[:, :, lo:hi], chi[:, :, lo:hi]
        for axis, label in ((1, "y"), (0, "z")):
            c_mat = centroid_along(rho_s, coords[win], axis) - args.centre
            c_met = centroid_along(chi_s, coords[win], axis) - args.centre
            gap = c_met - c_mat
            worst = max(worst, abs(gap))
            print(f"  {sector:<10} {label}: matter {c_mat:+.3e}   "
                  f"metric {c_met:+.3e}   offset {gap:+.4f}")

    verdict = "ALIGNED" if worst <= args.tolerance else "BIASED"
    print(f"  -> worst offset {worst:.4f}  [{verdict}]  "
          f"(tolerance {args.tolerance:g})")
    return 0 if worst <= args.tolerance else 1


if __name__ == "__main__":
    sys.exit(main())
