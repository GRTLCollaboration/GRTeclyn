from __future__ import annotations

import math
import os

import numpy as np

# Time-resolved MATTER-CONFINEMENT stream.  One row per plotfile.
#
# Every other matter signal in the pipeline is a peak/range (single-point) or a
# geometry-health quantity (min_lapse, min_chi); none of them can see matter
# spreading out / flying away, because that is a *spatial* property.  A lump can
# disperse into a ragged halo while its peak density (even total energy) RISES
# under pump injection -- which is exactly how "everything is confined and good"
# was reported while the frames showed the matter blowing apart.
#
# This stream measures the spatial distribution of the matter directly, using the
# mass weight w = scalar_activity = sum_k sqrt(phi_k^2 + Pi_k^2) (model-aware:
# real, complex, and bicomplex matter all map onto it):
#   total          = sum(w * dV)                 -- total matter "activity"
#   peak           = max(w)                      -- peak (the old, blind, signal)
#   rms_radius     = sqrt(sum w |x - x_bary|^2 / sum w)  -- THE spread / "flew away"
#   confined_frac  = sum_{|x-x_bary| < R_conf} w / sum w -- fraction still in lumps
#   bary_{x,y,z}   = sum(w x) / sum(w)           -- true MATTER barycentre
# R_conf defaults to 4 * well_width (the lump scale).  Dispersal shows up as
# rms_radius growing and confined_frac collapsing, regardless of peak/total.
CONFINEMENT_TIMESERIES_HEADER = (
    "# time  total_activity  peak_activity  rms_radius  confined_frac  "
    "bary_x  bary_y  bary_z  r_conf"
)


def _matter_weight_on_grid(cg, ftype: str, available: set[str]) -> np.ndarray | None:
    """Build the matter weight w = sum_k sqrt(phi_k^2 + Pi_k^2) on a covering grid.

    Mirrors the ``scalar_activity`` derived field (fields.py): combine the base
    complex field (phi/Pi) with any per-lump components (phi_lump*/Pi_lump*),
    covering real, single-complex and bicomplex (canonical + phantom) matter.
    """

    def arr(name: str) -> np.ndarray | None:
        if name not in available:
            return None
        try:
            return np.asarray(cg[ftype, name], dtype=np.float64)
        except Exception:  # noqa: BLE001
            return None

    lump_pairs = [
        (f"phi_lump{k}", f"Pi_lump{k}")
        for k in range(5)
        if f"phi_lump{k}" in available and f"Pi_lump{k}" in available
    ]
    phi = arr("phi")
    pi = arr("Pi")

    # Single complex field: phi/Pi = Re, phi_lump0/Pi_lump0 = Im -> one |Phi| norm.
    if lump_pairs == [("phi_lump0", "Pi_lump0")] and phi is not None and pi is not None:
        phi2 = arr("phi_lump0")
        pi2 = arr("Pi_lump0")
        return np.sqrt(phi**2 + pi**2 + phi2**2 + pi2**2)

    if lump_pairs:
        total: np.ndarray | None = None
        for phi_name, pi_name in lump_pairs:
            a = arr(phi_name)
            b = arr(pi_name)
            if a is None or b is None:
                continue
            term = np.sqrt(a**2 + b**2)
            total = term if total is None else total + term
        if total is not None:
            return total

    if phi is not None and pi is not None:
        return np.sqrt(phi**2 + pi**2)
    return None


def _extract_confinement_line(
    p: str,
    *,
    t: float,
    well_width: float = 1.5,
    level: int = 0,
    verbose: bool = False,
) -> str | None:
    """Compute per-plotfile matter-confinement moments -> dat row.

    Uses a uniform covering grid at ``level`` (default the base grid) so the
    mass-weighted moments are AMR-consistent and cheap.  Returns None if no
    matter field is present (no row written).
    """
    try:
        import yt

        ds = yt.load(p)
        available = {name for (_ftype, name) in ds.field_list}
        ftype = "boxlib"
        # Resolve field type robustly (stream datasets in tests use "stream").
        if not any(f == "boxlib" for (f, _n) in ds.field_list):
            ftype = ds.field_list[0][0] if ds.field_list else "boxlib"

        lvl = max(0, min(int(level), int(getattr(ds, "max_level", 0) or 0)))
        dims = np.asarray(ds.domain_dimensions, dtype=int) * (ds.refine_by**lvl)
        cg = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge, dims=dims)

        w = _matter_weight_on_grid(cg, ftype, available)
        if w is None:
            if verbose:
                print(f"WARNING: confinement: no matter field in {os.path.basename(p)}")
            return None

        total = float(w.sum())
        peak = float(w.max())
        if not math.isfinite(total) or total <= 0.0:
            # No matter present: emit a row so the gap in time is explicit.
            return (
                f"{t:.16e}  {total:.16e}  {peak:.16e}  "
                f"{0.0:.16e}  {0.0:.16e}  {0.0:.16e}  {0.0:.16e}  {0.0:.16e}  "
                f"{4.0 * well_width:.16e}"
            )

        dx = np.asarray(ds.domain_width, dtype=np.float64) / dims
        le = np.asarray(ds.domain_left_edge, dtype=np.float64)
        xs = le[0] + (np.arange(dims[0]) + 0.5) * dx[0]
        ys = le[1] + (np.arange(dims[1]) + 0.5) * dx[1]
        zs = le[2] + (np.arange(dims[2]) + 0.5) * dx[2]
        X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")

        bx = float((w * X).sum() / total)
        by = float((w * Y).sum() / total)
        bz = float((w * Z).sum() / total)
        r2 = (X - bx) ** 2 + (Y - by) ** 2 + (Z - bz) ** 2
        rms_radius = float(math.sqrt(max(0.0, (w * r2).sum() / total)))

        r_conf = 4.0 * float(well_width)
        confined_frac = float(w[r2 < r_conf * r_conf].sum() / total)

        return (
            f"{t:.16e}  {total:.16e}  {peak:.16e}  {rms_radius:.16e}  "
            f"{confined_frac:.16e}  {bx:.16e}  {by:.16e}  {bz:.16e}  {r_conf:.16e}"
        )
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: confinement extraction failed for "
                f"{os.path.basename(p)}: {exc}"
            )
        return None
