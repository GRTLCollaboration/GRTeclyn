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
# mass weight w = scalar_activity (model-aware; see _matter_sectors):
#   total          = sum(w * dV)                 -- total matter "activity"
#   peak           = max(w)                      -- peak (the old, blind, signal)
#   rms_radius     = sqrt(sum w |x - x_bary|^2 / sum w)  -- THE spread / "flew away"
#   confined_frac  = sum_{|x-x_bary| < R_conf} w / sum w -- fraction still in lumps
#   bary_{x,y,z}   = sum(w x) / sum(w)           -- true MATTER barycentre
# R_conf defaults to 4 * well_width (the lump scale).  Dispersal shows up as
# rms_radius growing and confined_frac collapsing, regardless of peak/total.
#
# COORDINATE vs PROPER VOLUME (cols 9-11).  ``dV`` above is the COORDINATE cell
# volume.  The physical mass element is ``w sqrt(gamma) d^3x`` with
# ``sqrt(gamma) = chi^-3/2`` (gamma_ij = gammatilde_ij / chi, det gammatilde = 1).
# When a campaign compares a collapsing control run against pumped runs the two
# have central ``chi`` differing by orders of magnitude -- on the 2026-07-28 pump
# ladder tp0 reached chi=5.7e-4 (sqrt(gamma)=7.4e4) while tp24 sat at chi=0.58
# (sqrt(gamma)=2.2), so the control run's core was weighted ~3e4 times less per
# unit PROPER volume than its rivals'.  A coordinate-only confined_frac therefore
# overstates how much the control run dispersed.  Both are written; quote the
# proper one for any cross-run claim.
#
# SECTORS (cols 12-16).  For the bicomplex (canonical + phantom) model the two
# matter sectors disperse at different rates, and an aggregate weight hides it.
# Per-sector moments are written so the canonical/phantom split is measurable.

CONFINEMENT_COLUMNS = (
    "time",
    "total_activity",
    "peak_activity",
    "rms_radius",
    "confined_frac",
    "bary_x",
    "bary_y",
    "bary_z",
    "r_conf",
    # --- proper-volume (sqrt(gamma) = chi^-3/2) counterparts -----------------
    "total_activity_proper",
    "rms_radius_proper",
    "confined_frac_proper",
    # --- per-sector, coordinate volume (nan when the model has one sector) ---
    "rms_radius_canon",
    "confined_frac_canon",
    "rms_radius_phantom",
    "confined_frac_phantom",
    "canon_mass_frac",
    # --- geometry lever, so the coord/proper gap is auditable from the file --
    "min_chi",
)

CONFINEMENT_TIMESERIES_HEADER = "# " + "  ".join(CONFINEMENT_COLUMNS)

# Matter-model tags understood by _matter_sectors (params.txt recipe_matter_model).
_BICOMPLEX_TAGS = frozenset(
    {"grtresna_bicomplex_scalar", "bicomplex", "bicomplex_scalar"}
)
_COMPLEX_TAGS = frozenset(
    {"grtresna_complex_scalar", "complex", "complex_scalar", "boson_star"}
)


def _lump_pairs(available: set[str]) -> list[tuple[str, str]]:
    return [
        (f"phi_lump{k}", f"Pi_lump{k}")
        for k in range(5)
        if f"phi_lump{k}" in available and f"Pi_lump{k}" in available
    ]


def _matter_sectors(
    src, ftype: str, available: set[str], *, model: str = ""
) -> dict[str, np.ndarray] | None:
    """Split the matter into physical sectors and return their activity weights.

    Returns ``{"total": w, ...}`` with optional ``"canonical"`` / ``"phantom"``
    entries, or None when the plotfile carries no matter field.

    THE GROUPING MATTERS.  From ``Examples/RadialRecipe/StateVariables.hpp``:

      independent real scalars : lump k lives in ``phi_lumpk``/``Pi_lumpk``;
                                 ``phi``/``Pi`` are unused (identically zero).
      single complex scalar    : ``phi``/``Pi`` = Re, ``phi_lump0``/``Pi_lump0``
                                 = Im of ONE field.
      bicomplex (canonical +   : ``phi``/``Pi`` = canonical Phi+ Re,
      phantom)                   ``phi_lump0``/``Pi_lump0`` = Phi+ Im,
                                 ``phi_lump1``/``Pi_lump1`` = phantom Phi- Re,
                                 ``phi_lump2``/``Pi_lump2`` = phantom Phi- Im.

    The bicomplex and 3-independent-real-scalar layouts have IDENTICAL field
    lists, so the model tag is authoritative and the sniff below is the
    fallback: ``phi``/``Pi`` are exactly zero for independent real scalars, and
    non-zero whenever they carry a complex field's real part.

    Before 2026-07-28 this function fell through to the independent-scalar
    branch for bicomplex runs, so ``phi``/``Pi`` -- the canonical field's real
    part and its momentum, roughly a quarter of the matter -- were dropped from
    every dispersion number, and Re/Im of the same complex field were summed as
    ``|Re| + |Im|`` instead of one modulus (which makes a stationary U(1) Q-ball
    read as breathing at its phase frequency omega).
    """

    def arr(name: str) -> np.ndarray | None:
        if name not in available:
            return None
        try:
            return np.asarray(src[ftype, name], dtype=np.float64)
        except Exception:  # noqa: BLE001
            return None

    def modulus(*names: str) -> np.ndarray | None:
        acc: np.ndarray | None = None
        for name in names:
            a = arr(name)
            if a is None:
                return None
            acc = a * a if acc is None else acc + a * a
        return None if acc is None else np.sqrt(acc)

    pairs = _lump_pairs(available)
    phi = arr("phi")
    pi = arr("Pi")
    has_base = phi is not None and pi is not None
    base_active = bool(
        has_base and (np.abs(phi).max() > 0.0 or np.abs(pi).max() > 0.0)
    )

    tag = model.strip().lower()
    if tag in _BICOMPLEX_TAGS:
        layout = "bicomplex"
    elif tag in _COMPLEX_TAGS:
        layout = "complex"
    elif tag:
        layout = "independent"
    elif "phi2" in available and "Pi2" in available:
        layout = "complex_phi2"
    elif base_active and len(pairs) == 3:
        layout = "bicomplex"
    elif base_active and len(pairs) == 1:
        layout = "complex"
    elif pairs:
        layout = "independent"
    elif has_base:
        layout = "complex_bare"
    else:
        return None

    if layout == "bicomplex":
        canon = modulus("phi", "Pi", "phi_lump0", "Pi_lump0")
        phantom = modulus("phi_lump1", "Pi_lump1", "phi_lump2", "Pi_lump2")
        if canon is not None and phantom is not None:
            return {"canonical": canon, "phantom": phantom, "total": canon + phantom}
        layout = "independent" if pairs else "complex_bare"

    if layout == "complex_phi2":
        w = modulus("phi", "Pi", "phi2", "Pi2")
        if w is not None:
            return {"total": w}
        layout = "independent" if pairs else "complex_bare"

    if layout == "complex":
        w = modulus("phi", "Pi", "phi_lump0", "Pi_lump0")
        if w is not None:
            return {"total": w}
        layout = "independent" if pairs else "complex_bare"

    if layout == "independent" and pairs:
        # Genuinely independent real scalars: each lump is its own field, so the
        # moduli add.  phi/Pi are unused slots here, but include them when they
        # are non-zero rather than silently discarding matter.
        total: np.ndarray | None = None
        for phi_name, pi_name in pairs:
            term = modulus(phi_name, pi_name)
            if term is None:
                continue
            total = term if total is None else total + term
        if total is not None:
            if base_active:
                base = modulus("phi", "Pi")
                if base is not None:
                    total = total + base
            return {"total": total}

    w = modulus("phi", "Pi")
    return None if w is None else {"total": w}


def _format_confinement_row(values: dict[str, float]) -> str:
    return "  ".join(f"{float(values[name]):.16e}" for name in CONFINEMENT_COLUMNS)


def _empty_row(*, t: float, total: float, peak: float, r_conf: float) -> str:
    values = {name: 0.0 for name in CONFINEMENT_COLUMNS}
    values.update(time=t, total_activity=total, peak_activity=peak, r_conf=r_conf)
    for name in (
        "rms_radius_canon",
        "confined_frac_canon",
        "rms_radius_phantom",
        "confined_frac_phantom",
        "canon_mass_frac",
    ):
        values[name] = math.nan
    values["min_chi"] = math.nan
    return _format_confinement_row(values)


def _moments(
    w: np.ndarray,
    vol: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    r_conf: float,
) -> tuple[float, float, float, float, float, float]:
    """(total, bary_x, bary_y, bary_z, rms_radius, confined_frac)."""
    mass = w * vol
    total = float(mass.sum())
    if not math.isfinite(total) or total <= 0.0:
        return (total, 0.0, 0.0, 0.0, math.nan, math.nan)
    bx = float((mass * x).sum() / total)
    by = float((mass * y).sum() / total)
    bz = float((mass * z).sum() / total)
    r2 = (x - bx) ** 2 + (y - by) ** 2 + (z - bz) ** 2
    rms = float(math.sqrt(max(0.0, (mass * r2).sum() / total)))
    frac = float(mass[r2 < r_conf * r_conf].sum() / total)
    return (total, bx, by, bz, rms, frac)


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
    matter_model: str = "",
    verbose: bool = False,
) -> str | None:
    """Compute per-plotfile matter-confinement moments -> dat row.

    AMR-AWARE by default (``level is None``): integrates over *all* cells weighted
    by their true cell volume, so each location is counted once at its FINEST
    available resolution.  This is what makes ``confined_frac`` actually sensitive
    to the AMR ``max_level`` -- the old level-0 covering grid was resolution-blind
    (ml=2 and ml=3 gave byte-identical numbers because only the base grid was
    sampled, never the refined soliton core).

    The mass weight is the model-aware activity (see :func:`_matter_sectors`) and
    the integrated mass element is ``w * dV``; on a uniform (single-level) grid dV
    is constant and cancels in every fraction.

    Pass an integer ``level`` to force the legacy uniform covering grid at that
    level (used by the fallback path and by callers that want a fixed grid).
    Both paths integrate against the true cell volume -- before 2026-07-28 the
    covering-grid branch summed raw cell values with no ``dV``, so a run that
    momentarily de-refined emitted ``total_activity`` rows a factor ``1/dV``
    (8x at level 0, dx=0.5) larger than its neighbours in the same file.

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

        if level is None and max_level > 0:
            # ---- AMR-aware, cell-volume-weighted integral (default) ---------
            src = ds.all_data()
            sectors = _matter_sectors(src, ftype, available, model=matter_model)
            if sectors is None:
                if verbose:
                    print(
                        f"WARNING: confinement: no matter field in "
                        f"{os.path.basename(p)}"
                    )
                return None
            sectors = {k: np.asarray(v, dtype=np.float64).ravel() for k, v in sectors.items()}
            dV = _code_values(src[("index", "cell_volume")]).ravel()
            x = _code_values(src[("index", "x")]).ravel()
            y = _code_values(src[("index", "y")]).ravel()
            z = _code_values(src[("index", "z")]).ravel()
            chi = (
                np.asarray(src[ftype, "chi"], dtype=np.float64).ravel()
                if "chi" in available
                else None
            )
        else:
            # ---- Legacy uniform covering grid (forced level, or single-level)
            lvl = 0 if level is None else max(0, min(int(level), max_level))
            dims = np.asarray(ds.domain_dimensions, dtype=int) * (ds.refine_by**lvl)
            src = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge, dims=dims)

            sectors = _matter_sectors(src, ftype, available, model=matter_model)
            if sectors is None:
                if verbose:
                    print(
                        f"WARNING: confinement: no matter field in "
                        f"{os.path.basename(p)}"
                    )
                return None
            sectors = {k: np.asarray(v, dtype=np.float64).ravel() for k, v in sectors.items()}
            dx = np.asarray(ds.domain_width, dtype=np.float64) / dims
            le = np.asarray(ds.domain_left_edge, dtype=np.float64)
            xs = le[0] + (np.arange(dims[0]) + 0.5) * dx[0]
            ys = le[1] + (np.arange(dims[1]) + 0.5) * dx[1]
            zs = le[2] + (np.arange(dims[2]) + 0.5) * dx[2]
            X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
            x, y, z = X.ravel(), Y.ravel(), Z.ravel()
            # Same units as the AMR branch: a true volume integral, not a cell sum.
            dV = np.full(x.shape, float(dx[0] * dx[1] * dx[2]), dtype=np.float64)
            chi = (
                np.asarray(src[ftype, "chi"], dtype=np.float64).ravel()
                if "chi" in available
                else None
            )

        w = sectors["total"]
        peak = float(w.max()) if w.size else 0.0
        total_coord = float((w * dV).sum())
        if not math.isfinite(total_coord) or total_coord <= 0.0:
            return _empty_row(t=t, total=total_coord, peak=peak, r_conf=r_conf)

        # sqrt(gamma) = chi^-3/2 exactly (gamma_ij = gammatilde_ij / chi with
        # det gammatilde = 1), so this is the PROPER volume element.
        if chi is not None and chi.size == dV.size:
            chi_safe = np.maximum(chi, 1.0e-30)
            dV_proper = dV * np.power(chi_safe, -1.5)
            min_chi = float(chi.min())
        else:
            dV_proper = dV
            min_chi = math.nan

        total, bx, by, bz, rms, frac = _moments(w, dV, x, y, z, r_conf)
        total_pr, _, _, _, rms_pr, frac_pr = _moments(w, dV_proper, x, y, z, r_conf)

        values = {
            "time": t,
            "total_activity": total,
            "peak_activity": peak,
            "rms_radius": rms,
            "confined_frac": frac,
            "bary_x": bx,
            "bary_y": by,
            "bary_z": bz,
            "r_conf": r_conf,
            "total_activity_proper": total_pr,
            "rms_radius_proper": rms_pr,
            "confined_frac_proper": frac_pr,
            "rms_radius_canon": math.nan,
            "confined_frac_canon": math.nan,
            "rms_radius_phantom": math.nan,
            "confined_frac_phantom": math.nan,
            "canon_mass_frac": math.nan,
            "min_chi": min_chi,
        }

        if "canonical" in sectors and "phantom" in sectors:
            tot_c, _, _, _, rms_c, frac_c = _moments(
                sectors["canonical"], dV, x, y, z, r_conf
            )
            tot_p, _, _, _, rms_p, frac_p = _moments(
                sectors["phantom"], dV, x, y, z, r_conf
            )
            values["rms_radius_canon"] = rms_c
            values["confined_frac_canon"] = frac_c
            values["rms_radius_phantom"] = rms_p
            values["confined_frac_phantom"] = frac_p
            denom = tot_c + tot_p
            values["canon_mass_frac"] = tot_c / denom if denom > 0.0 else math.nan

        return _format_confinement_row(values)
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: confinement extraction failed for "
                f"{os.path.basename(p)}: {exc}"
            )
        return None
