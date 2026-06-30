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


def _format_confinement_row(
    *,
    t: float,
    total: float,
    peak: float,
    rms_radius: float,
    confined_frac: float,
    bx: float,
    by: float,
    bz: float,
    r_conf: float,
) -> str:
    return (
        f"{t:.16e}  {total:.16e}  {peak:.16e}  {rms_radius:.16e}  "
        f"{confined_frac:.16e}  {bx:.16e}  {by:.16e}  {bz:.16e}  {r_conf:.16e}"
    )


def _code_values(arr) -> np.ndarray:
    """Return a unyt/ndarray as plain float64 in code units when possible.

    Positions and cell volumes must be in *code units* so they are comparable to
    ``r_conf = 4 * well_width`` (also code units).  yt tags boxlib fields with cgs
    lengths by default, so fall back through code_length -> raw values.
    """
    for unit in ("code_length", "code_length**3"):
        try:
            return np.asarray(arr.to_value(unit), dtype=np.float64)
        except Exception:  # noqa: BLE001
            continue
    try:
        return np.asarray(arr.to_ndarray(), dtype=np.float64)
    except Exception:  # noqa: BLE001
        return np.asarray(arr, dtype=np.float64)


def _extract_confinement_line(
    p: str,
    *,
    t: float,
    well_width: float = 1.5,
    level: int | None = None,
    verbose: bool = False,
) -> str | None:
    """Compute per-plotfile matter-confinement moments -> dat row.

    AMR-AWARE by default (``level is None``): integrates over *all* cells weighted
    by their true cell volume, so each location is counted once at its FINEST
    available resolution.  This is what makes ``confined_frac`` actually sensitive
    to the AMR ``max_level`` -- the old level-0 covering grid was resolution-blind
    (ml=2 and ml=3 gave byte-identical numbers because only the base grid was
    sampled, never the refined soliton core).

    The mass weight is ``w = scalar_activity`` and the integrated mass element is
    ``w * dV``; on a uniform (single-level) grid dV is constant and cancels in the
    fraction, so this is backward-compatible with the old covering-grid result.

    Pass an integer ``level`` to force the legacy uniform covering grid at that
    level (used by the fallback path and by callers that want a fixed grid).
    Returns None if no matter field is present (no row written).
    """
    try:
        import yt

        ds = yt.load(p)
        available = {name for (_ftype, name) in ds.field_list}
        ftype = "boxlib"
        # Resolve field type robustly (stream datasets in tests use "stream").
        if not any(f == "boxlib" for (f, _n) in ds.field_list):
            ftype = ds.field_list[0][0] if ds.field_list else "boxlib"

        r_conf = 4.0 * float(well_width)
        max_level = int(getattr(ds, "max_level", 0) or 0)

        # ---- AMR-aware, cell-volume-weighted integral (default) -------------
        if level is None and max_level > 0:
            ad = ds.all_data()
            w = _matter_weight_on_grid(ad, ftype, available)
            if w is None:
                if verbose:
                    print(
                        f"WARNING: confinement: no matter field in "
                        f"{os.path.basename(p)}"
                    )
                return None
            w = np.asarray(w, dtype=np.float64).ravel()
            dV = _code_values(ad[("index", "cell_volume")]).ravel()
            mass = w * dV
            total = float(mass.sum())
            peak = float(w.max()) if w.size else 0.0
            if not math.isfinite(total) or total <= 0.0:
                return _format_confinement_row(
                    t=t, total=total, peak=peak, rms_radius=0.0,
                    confined_frac=0.0, bx=0.0, by=0.0, bz=0.0, r_conf=r_conf,
                )
            x = _code_values(ad[("index", "x")]).ravel()
            y = _code_values(ad[("index", "y")]).ravel()
            z = _code_values(ad[("index", "z")]).ravel()
            bx = float((mass * x).sum() / total)
            by = float((mass * y).sum() / total)
            bz = float((mass * z).sum() / total)
            r2 = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
            rms_radius = float(math.sqrt(max(0.0, (mass * r2).sum() / total)))
            confined_frac = float(mass[r2 < r_conf * r_conf].sum() / total)
            return _format_confinement_row(
                t=t, total=total, peak=peak, rms_radius=rms_radius,
                confined_frac=confined_frac, bx=bx, by=by, bz=bz, r_conf=r_conf,
            )

        # ---- Legacy uniform covering grid (forced level, or single-level) ---
        lvl = 0 if level is None else max(0, min(int(level), max_level))
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
            return _format_confinement_row(
                t=t, total=total, peak=peak, rms_radius=0.0,
                confined_frac=0.0, bx=0.0, by=0.0, bz=0.0, r_conf=r_conf,
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

        confined_frac = float(w[r2 < r_conf * r_conf].sum() / total)

        return _format_confinement_row(
            t=t, total=total, peak=peak, rms_radius=rms_radius,
            confined_frac=confined_frac, bx=bx, by=by, bz=bz, r_conf=r_conf,
        )
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: confinement extraction failed for "
                f"{os.path.basename(p)}: {exc}"
            )
        return None
