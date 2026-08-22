"""Per-sector CORE dynamics: halo-free position, matter momentum, gauge check.

The scrutiny stream for the Bondi-dipole runaway.  ``sector_barycenters.dat``
answers "did the sectors move?" with a first moment over the WHOLE sector --
core plus shed radiation.  That is the right cheap diagnostic, but by late
times the halo carries a large share of the weight (in the reference pair run
the tracked activity falls ~35% and the rms radius doubles by t=60), so the
barycentre stops being a clean measure of where the star is.  Worse, a
barycentre is a COORDINATE quantity: a referee is entitled to ask whether the
slice moved rather than the matter.

This stream answers the three questions the barycentre cannot:

1. **Where is the core?**  Peak location and a core-weighted centroid
   (weight w^4, restricted to cells above a fraction of the peak) that the
   halo cannot pull.  Differentiating it gives velocity/acceleration directly
   instead of via an integrated displacement.
2. **Is momentum balanced?**  The matter momentum of each sector,
   ``P_i = integral S_i sqrt(gamma) d3x`` with
   ``S_i = -(Pi_1 d_i phi_1 + Pi_2 d_i phi_2)``, and the GRAVITATING total
   ``P_canon - P_phantom`` (the phantom's stress-energy enters Einstein with a
   flipped sign).  Bondi's runaway is momentum-conserving: both bodies
   accelerate while the total stays ~0.  Nothing else in the campaign measures
   that -- it is the signature that separates a physical runaway from a
   numerical kick.
   NB this is the MATTER momentum (a volume integral), not the ADM momentum
   (a surface integral); they differ by the gravitational field momentum.
3. **Is the motion gauge?**  The shift at each core and the domain maximum of
   |shift|, plus the PROPER separation between the cores
   (``integral sqrt(gamma_xx) dx`` along the pair axis), which does not care
   how the coordinates are laid out.

One module, one ``sector_dynamics.dat``, one opt-in flag: the existing
``confinement`` / ``sector_barycenters`` column contracts are positional and
shared with live campaigns, so they stay untouched.  The three quantities share
a single stream because they share the expensive part -- one plotfile load and
one covering grid -- and answer one question together ("is the drift real?").

Cost: this is the only consumer stream that builds a covering grid, so it is
off by default and enabled per campaign via ``GRTECLYN_SECTOR_DYNAMICS=1``
(~1.6 s per plotfile at level 0, ~4.6 s at level 1, on a 128^3 domain).

SCOPE.  The core tracker assumes ONE lump per sector, which is exactly the
mixed-pair geometry it exists for (pm / mp / equal-mass cells).  In a
same-sector control (pp, mm) both lumps live in one sector, and the reported
"core" is then a single centroid sitting between them -- still a valid
position, but not the position of either lump.  Read those cells with
``sector_barycenters.dat`` instead.
"""

from __future__ import annotations

import math
import os

import numpy as np

from .confinement import _code_values

SECTOR_DYNAMICS_COLUMNS = (
    "time",
    # Core position (halo-free) and peak amplitude, per sector.
    "core_x_canon",
    "core_y_canon",
    "core_z_canon",
    "peak_canon",
    "core_x_phantom",
    "core_y_phantom",
    "core_z_phantom",
    "peak_phantom",
    # Separation between cores: coordinate, and proper along the pair axis.
    "coord_sep",
    "proper_sep",
    # Matter momentum along x: per sector, and the gravitating total
    # (canonical - phantom).  Bondi: the total stays ~0 while both move.
    "px_canon",
    "px_phantom",
    "px_total",
    # Gauge: shift sampled at each core, and the domain maximum of |shift|.
    "shift_x_canon",
    "shift_x_phantom",
    "max_abs_shift",
)

SECTOR_DYNAMICS_HEADER = "# " + "  ".join(SECTOR_DYNAMICS_COLUMNS)

# Cells below this fraction of the sector peak are halo, not core.
_CORE_THRESHOLD = 0.25
# Centroid weight exponent; >1 suppresses whatever halo survives the threshold.
_CORE_POWER = 4.0

# Bicomplex layout (Examples/RadialRecipe/StateVariables.hpp): canonical Phi+ is
# (phi, Pi) Re + (phi_lump0, Pi_lump0) Im; phantom Phi- is (phi_lump1, Pi_lump1)
# Re + (phi_lump2, Pi_lump2) Im.
_SECTOR_FIELDS = {
    "canon": (("phi", "Pi"), ("phi_lump0", "Pi_lump0")),
    "phantom": (("phi_lump1", "Pi_lump1"), ("phi_lump2", "Pi_lump2")),
}

_NAN_CORE = {
    "x": math.nan, "y": math.nan, "z": math.nan, "peak": math.nan, "px": math.nan,
    "ix": None, "iy": None, "iz": None,
}


def _core_moments(
    w: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray
) -> dict:
    """Peak amplitude and core-weighted centroid of one sector.

    An absent sector (identically zero field, e.g. the phantom slot in a
    single-canonical run) must read as absent, not as a core sitting at the
    origin -- that would be poison for a drift measurement.
    """
    peak = float(np.max(w)) if w.size else 0.0
    if not math.isfinite(peak) or peak <= 0.0:
        return dict(_NAN_CORE)

    mask = w >= _CORE_THRESHOLD * peak
    weight = np.where(mask, w, 0.0) ** _CORE_POWER
    norm = float(weight.sum())
    if norm <= 0.0:
        return dict(_NAN_CORE)

    cx = float((weight * x).sum() / norm)
    cy = float((weight * y).sum() / norm)
    cz = float((weight * z).sum() / norm)
    peak_idx = np.unravel_index(int(np.argmax(w)), w.shape)
    return {
        "x": cx, "y": cy, "z": cz, "peak": peak,
        "ix": int(peak_idx[0]), "iy": int(peak_idx[1]), "iz": int(peak_idx[2]),
    }


def _momentum_x(
    phi1: np.ndarray, pi1: np.ndarray, phi2: np.ndarray, pi2: np.ndarray,
    sqrt_gamma: np.ndarray, dx: float,
) -> float:
    """integral S_x sqrt(gamma) d3x for one complex field.

    ``S_i = -(Pi_1 d_i phi_1 + Pi_2 d_i phi_2)`` -- the same expression the
    solver and the evolution use for the momentum source (GRTresna
    ``ComplexScalarField::compute_emtensor``, GRTeclyn
    ``BiComplexScalarField``), so this integral is the matter momentum those
    codes actually gravitate with.
    """
    d1 = np.gradient(phi1, dx, axis=0)
    d2 = np.gradient(phi2, dx, axis=0)
    s_x = -(pi1 * d1 + pi2 * d2)
    return float(np.sum(s_x * sqrt_gamma) * dx**3)


def _row_through(
    volume: np.ndarray, axis_y: np.ndarray, axis_z: np.ndarray,
    y_t: float, z_t: float,
) -> np.ndarray:
    """The x-row of ``volume`` through transverse position ``(y_t, z_t)``.

    Bilinear in the two transverse axes.  Snapping to the nearest cell instead
    would quantise the row to the grid, which is exactly the error this module
    used to make.
    """
    ny, nz = axis_y.size, axis_z.size
    if ny < 2 or nz < 2:
        return volume[:, min(ny - 1, 0), min(nz - 1, 0)]
    dy = float(axis_y[1] - axis_y[0])
    dz = float(axis_z[1] - axis_z[0])
    fy = (float(y_t) - float(axis_y[0])) / dy
    fz = (float(z_t) - float(axis_z[0])) / dz
    j = int(np.clip(np.floor(fy), 0, ny - 2))
    k = int(np.clip(np.floor(fz), 0, nz - 2))
    wy = float(np.clip(fy - j, 0.0, 1.0))
    wz = float(np.clip(fz - k, 0.0, 1.0))
    return (
        (1.0 - wy) * (1.0 - wz) * volume[:, j, k]
        + wy * (1.0 - wz) * volume[:, j + 1, k]
        + (1.0 - wy) * wz * volume[:, j, k + 1]
        + wy * wz * volume[:, j + 1, k + 1]
    )


def _proper_separation(
    h11: np.ndarray, chi: np.ndarray, core_a: dict, core_b: dict,
    axes: list[np.ndarray],
) -> float:
    """integral sqrt(gamma_xx) dx along the pair axis between the two cores.

    ``gamma_ij = h_ij / chi`` in the GRTeclyn variables, so the proper length
    element along x is ``sqrt(h11 / chi) dx``.

    BOTH ENDPOINTS ARE SUB-CELL.  Until 2026-08-21 this integrated between the
    two cores' integer PEAK CELL indices and sampled the row at the canonical
    core's integer transverse indices, which quantised the answer to whole
    cells -- dx = 0.5 in the campaign, coarser than the 0.1-ish signal it was
    supposed to resolve, so the column was useless for measuring whether a pair
    holds its separation.  The cores' continuous centroids are already computed
    (``core["x"]`` and friends); this uses them.

    Values written before that date are quantised and must not be compared with
    values written after it.
    """
    for core in (core_a, core_b):
        if not math.isfinite(core["x"]):
            return math.nan
    axis_x, axis_y, axis_z = axes
    if axis_x.size < 2:
        return math.nan

    lo, hi = sorted((float(core_a["x"]), float(core_b["x"])))
    if hi <= lo:
        return 0.0
    if lo < float(axis_x[0]) or hi > float(axis_x[-1]):
        return math.nan

    # The two cores sit on a common axis by construction; averaging their
    # transverse positions is the honest row to walk when they disagree
    # slightly, and is exact when they do not.
    y_t = 0.5 * (float(core_a["y"]) + float(core_b["y"]))
    z_t = 0.5 * (float(core_a["z"]) + float(core_b["z"]))

    ratio = h11 / np.maximum(chi, 1e-12)
    row = np.sqrt(np.maximum(_row_through(ratio, axis_y, axis_z, y_t, z_t), 0.0))

    # Cumulative trapezoid, then read off at the two sub-cell endpoints: exact
    # for the piecewise-linear interpolant the trapezoid rule already assumes.
    seg = 0.5 * (row[1:] + row[:-1]) * np.diff(axis_x)
    cum = np.concatenate(([0.0], np.cumsum(seg)))
    return float(np.interp(hi, axis_x, cum) - np.interp(lo, axis_x, cum))


def _extract_sector_dynamics_line(
    p: str, *, t: float, level: int = 0, verbose: bool = False,
) -> str | None:
    """Core dynamics / momentum / gauge row for one plotfile.

    Returns None when the plotfile carries no bicomplex matter, so a run
    without two sectors produces no stream rather than a stream of nans.
    """
    try:
        import yt

        ds = yt.load(p)
        available = {name for (_ftype, name) in ds.field_list}
        ftype = "boxlib"
        if not any(f == "boxlib" for (f, _n) in ds.field_list):
            ftype = ds.field_list[0][0] if ds.field_list else "boxlib"

        needed = {"phi", "Pi", "phi_lump0", "Pi_lump0", "phi_lump1", "Pi_lump1",
                  "phi_lump2", "Pi_lump2", "chi", "h11"}
        if not needed.issubset(available):
            if verbose:
                missing = sorted(needed - available)
                print(
                    f"WARNING: sector_dynamics: {os.path.basename(p)} lacks "
                    f"{missing} -- not a bicomplex run"
                )
            return None

        dims = ds.domain_dimensions * (2 ** int(level))
        grid = ds.covering_grid(int(level), left_edge=ds.domain_left_edge, dims=dims)

        def field(name: str) -> np.ndarray:
            return np.asarray(_code_values(grid[ftype, name]), dtype=np.float64)

        left = np.asarray(_code_values(ds.domain_left_edge), dtype=np.float64)
        width = np.asarray(_code_values(ds.domain_width), dtype=np.float64)
        dx = float(width[0] / dims[0])
        axes = [left[i] + (np.arange(dims[i]) + 0.5) * (width[i] / dims[i]) for i in range(3)]
        x, y, z = np.meshgrid(axes[0], axes[1], axes[2], indexing="ij")

        chi = field("chi")
        h11 = field("h11")
        # gamma_ij = h_ij / chi with det h = 1  =>  sqrt(gamma) = chi^(-3/2).
        sqrt_gamma = np.power(np.maximum(chi, 1e-12), -1.5)

        values: dict[str, float] = {"time": t}
        cores: dict[str, dict] = {}
        for suffix, ((re_phi, re_pi), (im_phi, im_pi)) in _SECTOR_FIELDS.items():
            phi1, pi1 = field(re_phi), field(re_pi)
            phi2, pi2 = field(im_phi), field(im_pi)
            w = np.sqrt(phi1 * phi1 + pi1 * pi1 + phi2 * phi2 + pi2 * pi2)
            core = _core_moments(w, x, y, z)
            cores[suffix] = core
            values[f"core_x_{suffix}"] = core["x"]
            values[f"core_y_{suffix}"] = core["y"]
            values[f"core_z_{suffix}"] = core["z"]
            values[f"peak_{suffix}"] = core["peak"]
            values[f"px_{suffix}"] = (
                _momentum_x(phi1, pi1, phi2, pi2, sqrt_gamma, dx)
                if core["ix"] is not None else math.nan
            )

        canon, phantom = cores["canon"], cores["phantom"]
        values["coord_sep"] = canon["x"] - phantom["x"]
        values["proper_sep"] = _proper_separation(h11, chi, canon, phantom, axes)
        # The phantom's stress-energy enters Einstein with a flipped sign, so
        # the GRAVITATING total subtracts its matter momentum.
        values["px_total"] = values["px_canon"] - values["px_phantom"]

        if {"shift1", "shift2", "shift3"}.issubset(available):
            shift1, shift2, shift3 = field("shift1"), field("shift2"), field("shift3")
            values["max_abs_shift"] = float(
                np.max(np.sqrt(shift1**2 + shift2**2 + shift3**2))
            )
            for suffix, core in cores.items():
                values[f"shift_x_{suffix}"] = (
                    float(shift1[core["ix"], core["iy"], core["iz"]])
                    if core["ix"] is not None else math.nan
                )
        else:
            values["max_abs_shift"] = math.nan
            values["shift_x_canon"] = math.nan
            values["shift_x_phantom"] = math.nan

        return "  ".join(
            f"{float(values[name]):.16e}" for name in SECTOR_DYNAMICS_COLUMNS
        )
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: sector_dynamics extraction failed for "
                f"{os.path.basename(p)}: {exc}"
            )
        return None
