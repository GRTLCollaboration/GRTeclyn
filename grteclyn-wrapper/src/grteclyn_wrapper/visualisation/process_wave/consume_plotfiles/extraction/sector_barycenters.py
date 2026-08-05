"""Per-sector (canonical vs phantom) barycentre stream.

The Bondi-dipole runaway diagnostic: a mixed-sign pair translates while the
AGGREGATE matter barycentre in ``confinement.dat`` stays put (equal-magnitude
sectors cancel in the sum), so the runaway is only visible as the drift of each
sector's own centroid.  Plotfiles are purged after consumption, which means
this has to be measured live, in the streaming consumer -- it cannot be
recovered afterwards.

Deliberately a separate stream from ``confinement.py`` (own columns tuple, own
``sector_barycenters.dat``): the confinement contract is positional and shared
with live campaigns, and stays untouched.  Pure helpers (`_matter_sectors`,
`_moments`, `_code_values`) are reused by import, not copied.

Unlike the confinement stream, the sector split here takes the matter-model
tag from the caller: the field-sniff fallback cannot classify a run whose
canonical sector is identically zero (an all-phantom configuration looks like
independent scalars), and the Bondi control matrix contains exactly such runs.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .confinement import _code_values, _matter_sectors, _moments

SECTOR_BARYCENTERS_COLUMNS = (
    "time",
    "total_canon",
    "bary_x_canon",
    "bary_y_canon",
    "bary_z_canon",
    "rms_radius_canon",
    "total_phantom",
    "bary_x_phantom",
    "bary_y_phantom",
    "bary_z_phantom",
    "rms_radius_phantom",
)

SECTOR_BARYCENTERS_HEADER = "# " + "  ".join(SECTOR_BARYCENTERS_COLUMNS)

_SECTOR_KEYS = {"canon": "canonical", "phantom": "phantom"}


def _sector_values(
    w: np.ndarray,
    dV: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    r_conf: float,
) -> dict[str, float]:
    """Moments of one sector; nan barycentre when the sector holds no matter.

    ``_moments`` reports a (0,0,0) barycentre for an empty sector, which here
    would read as a real centroid sitting at the origin -- poison for a drift
    measurement.  An absent sector must be visibly absent.
    """
    total, bx, by, bz, rms, _frac = _moments(w, dV, x, y, z, r_conf)
    if not math.isfinite(total) or total <= 0.0:
        return {"total": total, "x": math.nan, "y": math.nan, "z": math.nan, "rms": math.nan}
    return {"total": total, "x": bx, "y": by, "z": bz, "rms": rms}


def _extract_sector_barycenters_line(
    p: str,
    *,
    t: float,
    matter_model: str = "",
    well_width: float = 1.5,
    verbose: bool = False,
) -> str | None:
    """Compute per-sector centroid moments for one plotfile -> dat row.

    AMR-aware: integrates every cell at its finest available resolution against
    the true cell volume (same convention as the confinement stream).  Returns
    None when the plotfile has no matter or no two-sector split is possible --
    a single-sector model produces no stream rather than a stream of nans.
    """
    try:
        import yt

        ds = yt.load(p)
        available = {name for (_ftype, name) in ds.field_list}
        ftype = "boxlib"
        if not any(f == "boxlib" for (f, _n) in ds.field_list):
            ftype = ds.field_list[0][0] if ds.field_list else "boxlib"

        src = ds.all_data()
        sectors = _matter_sectors(src, ftype, available, model=matter_model)
        if sectors is None or "canonical" not in sectors or "phantom" not in sectors:
            if verbose:
                print(
                    f"WARNING: sector_barycenters: no canonical/phantom split in "
                    f"{os.path.basename(p)} (model tag {matter_model!r})"
                )
            return None

        dV = _code_values(src[("index", "cell_volume")]).ravel()
        x = _code_values(src[("index", "x")]).ravel()
        y = _code_values(src[("index", "y")]).ravel()
        z = _code_values(src[("index", "z")]).ravel()
        r_conf = 4.0 * float(well_width)

        values: dict[str, float] = {"time": t}
        for suffix, sector_name in _SECTOR_KEYS.items():
            w = np.asarray(sectors[sector_name], dtype=np.float64).ravel()
            moments = _sector_values(w, dV, x, y, z, r_conf)
            values[f"total_{suffix}"] = moments["total"]
            values[f"bary_x_{suffix}"] = moments["x"]
            values[f"bary_y_{suffix}"] = moments["y"]
            values[f"bary_z_{suffix}"] = moments["z"]
            values[f"rms_radius_{suffix}"] = moments["rms"]

        return "  ".join(
            f"{float(values[name]):.16e}" for name in SECTOR_BARYCENTERS_COLUMNS
        )
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: sector_barycenters extraction failed for "
                f"{os.path.basename(p)}: {exc}"
            )
        return None
