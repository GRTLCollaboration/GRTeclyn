from __future__ import annotations

from typing import Sequence

import numpy as np
import yt

from ..config import _field_frame_config, _frames_auto_zlim_enabled
from ..fields import _canonical_field_name, _field_key, _register_derived_fields
from .center import _resolve_frame_physics_center

#: Marks a ``frame_zlims`` dict as measured over the whole series rather than
#: locked from one plotfile.  Not a field name, so it can never collide.
_SCANNED_KEY = "__scanned_series__"


def _auto_zlim_from_array(values: np.ndarray, field_name: str) -> tuple[float, float] | None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return None

    if field_name == "K":
        max_abs = float(np.nanpercentile(np.abs(finite), 99.5))
        max_abs = max(max_abs, 5.0e-6)
        return (-max_abs, max_abs)

    if field_name == "Weyl4_Mag":
        hi = float(np.nanpercentile(finite, 99.5))
        hi = max(hi, 1.0e-8)
        return (0.0, hi)

    if field_name == "chi":
        lo, hi = np.nanpercentile(finite, [0.5, 99.5])
        lo, hi = float(lo), float(hi)
        span = hi - lo
        min_span = 0.04
        if span < min_span:
            mid = 0.5 * (lo + hi)
            lo, hi = mid - 0.5 * min_span, mid + 0.5 * min_span
            span = min_span
        pad = max(0.05 * span, 0.005)
        return (max(0.05, lo - pad), min(1.05, hi + pad))

    if field_name in {"shift1", "shift2", "shift3"}:
        max_abs = float(np.nanpercentile(np.abs(finite), 99.5))
        max_abs = max(max_abs, 1.0e-4)
        return (-max_abs, max_abs)

    if field_name in {"Weyl4_Re", "Weyl4_Im", "phi", "Pi", "phi_lump0", "Pi_lump0", "chi_minus_1"}:
        max_abs = float(np.nanpercentile(np.abs(finite), 99.5))
        max_abs = max(max_abs, 1.0e-8)
        return (-max_abs, max_abs)

    if field_name in {"lump_activity", "scalar_activity"}:
        hi = float(np.nanpercentile(finite, 99.5))
        hi = max(hi, 1.0e-6)
        return (0.0, hi)

    lo, hi = np.nanpercentile(finite, [1.0, 99.0])
    lo = float(lo)
    hi = float(hi)
    span = hi - lo
    min_span = 1.0e-4 if field_name in {"lapse"} else 1.0e-12
    if span < min_span:
        mid = 0.5 * (lo + hi)
        half_span = 0.5 * min_span
        return (mid - half_span, mid + half_span)
    pad = 0.05 * span if span > 0.0 else min_span
    return (lo - pad, hi + pad)


def _resolve_plot_zlim(
    field: str,
    win: np.ndarray,
    cfg: dict,
    *,
    auto_zlim: bool | None,
    frame_zlims: dict[str, list[float]] | None,
    use_global_zlim: bool,
) -> tuple[float, float]:
    if _frames_auto_zlim_enabled(auto_zlim) or cfg.get("auto_zlim"):
        auto = _auto_zlim_from_array(win, field)
        if auto is not None:
            return auto

    # A scanned series range was measured over every frame of THIS run, so it
    # beats both the per-frame rescale (which makes the colourbar move) and the
    # hand-set presets (which are a guess made before the run).  The sentinel is
    # carried in the dict itself so the scan needs no new plumbing through the
    # worker.  See zlim_scan.scan_series_zlims.
    if frame_zlims and frame_zlims.get(_SCANNED_KEY):
        stored = frame_zlims.get(field)
        if stored is not None:
            return (float(stored[0]), float(stored[1]))

    if cfg.get("per_frame_zlim"):
        auto = _auto_zlim_from_array(win, field)
        if auto is not None:
            return auto

    preset = cfg["zlim"]
    if preset[0] is not None:
        # Fixed preset wins unconditionally: a per-frame fallback (previously used
        # for lump_activity/scalar_activity when the signal looked weak) makes the
        # colorbar limits change frame-to-frame, so the same color means a
        # different value in every frame and the movie "bounces".  For stable
        # movies the scale must be held fixed; opt into per-frame scaling via
        # GRTECLYN_FRAMES_AUTO_ZLIM or a per-field ``auto_zlim`` flag instead.
        return preset

    if use_global_zlim and frame_zlims is not None:
        stored = frame_zlims.get(field)
        if stored is not None:
            return (float(stored[0]), float(stored[1]))

    auto = _auto_zlim_from_array(win, field)
    return auto if auto is not None else (0.0, 1.0)


def _extract_native_slice_window(
    ds,
    field: str,
    axis: str,
    coord: float | None,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    corner: bool,
) -> np.ndarray:
    physics_center = _resolve_frame_physics_center(ds, axis, coord, zoom, center_xyz, corner)
    _register_derived_fields(ds, field)
    plot_field = _field_key(ds, field)

    lvl = int(ds.index.max_level)
    ds.force_periodicity()
    full_dims = [int(round(n * (2 ** lvl))) for n in ds.domain_dimensions]
    cg = ds.covering_grid(level=lvl, left_edge=ds.domain_left_edge, dims=full_dims)
    data = np.array(cg[plot_field])
    le = np.array(ds.domain_left_edge.d, dtype=float)
    re = np.array(ds.domain_right_edge.d, dtype=float)
    cell = (re - le) / np.array(full_dims, dtype=float)

    ai = {"x": 0, "y": 1, "z": 2}[axis]
    sidx = int(
        np.clip(round((float(physics_center[ai]) - le[ai]) / cell[ai] - 0.5), 0, full_dims[ai] - 1)
    )
    slab = np.take(data, sidx, axis=ai)
    h_ax, v_ax = {0: (1, 2), 1: (2, 0), 2: (0, 1)}[ai]
    remaining = [a for a in (0, 1, 2) if a != ai]
    arr = slab if remaining == [v_ax, h_ax] else slab.T
    extent = [le[h_ax], re[h_ax], le[v_ax], re[v_ax]]

    if zoom is not None:
        hh = float(zoom) / 2.0
        hcoords = np.linspace(extent[0], extent[1], arr.shape[1])
        vcoords = np.linspace(extent[2], extent[3], arr.shape[0])
        hmask = (hcoords >= physics_center[h_ax] - hh) & (hcoords <= physics_center[h_ax] + hh)
        vmask = (vcoords >= physics_center[v_ax] - hh) & (vcoords <= physics_center[v_ax] + hh)
        if hmask.any() and vmask.any():
            return arr[np.ix_(vmask, hmask)]
    return arr


def _lock_frame_zlims_from_plotfile(
    plot_path: str, args_dict: dict, *, include_per_frame_fields: bool = False
) -> dict[str, list[float]]:
    """Colour limits this one plotfile would choose, per field.

    ``include_per_frame_fields`` keeps the fields flagged ``per_frame_zlim``.
    Locking those from a single plotfile is wrong -- that is why they are
    skipped by default -- but the series scanner needs exactly this per-plotfile
    answer so it can take the envelope over the whole run.  See
    ``zlim_scan.scan_series_zlims``.
    """
    ds = yt.load(plot_path)
    frame_fields = [_canonical_field_name(f) for f in args_dict.get("frames_fields", [])]
    axis = args_dict.get("frames_axis", "z")
    zoom = args_dict.get("frames_zoom")
    center_xyz = args_dict.get("frames_center")
    corner = bool(args_dict.get("frames_corner"))
    zlims: dict[str, list[float]] = {}

    # AMR HQ plotfiles: covering_grid at max_level allocates a full uniform
    # 2048³ array (~200+ GB).  Use the same FRB SlicePlot path as rendering.
    use_frb = int(ds.index.max_level) > 0
    if use_frb:
        from .center import _frame_buff_size, _resolve_frame_physics_center
        from ..fields import _field_key, _register_derived_fields

        physics_center = _resolve_frame_physics_center(
            ds, axis, args_dict.get("frames_coord"), zoom, center_xyz, corner
        )
        plot_center = ds.arr(physics_center, "code_length")
        buff = _frame_buff_size(ds, zoom)

    for fld in frame_fields:
        try:
            if not include_per_frame_fields and _field_frame_config(fld).get("per_frame_zlim"):
                continue
            if use_frb:
                _register_derived_fields(ds, fld)
                plot_field = _field_key(ds, fld)
                available = list(getattr(ds, "field_list", [])) + list(
                    getattr(ds, "derived_field_list", [])
                )
                if plot_field not in available:
                    if args_dict.get("verbose", False):
                        print(f"WARNING: zlim-lock skipped missing field {fld!r}")
                    continue
                slc = yt.SlicePlot(
                    ds, axis, plot_field, center=plot_center, buff_size=(buff, buff)
                )
                if zoom is not None:
                    slc.set_width((float(zoom), "code_length"))
                slc.render()
                win = np.asarray(slc.frb[plot_field])
            else:
                win = _extract_native_slice_window(
                    ds,
                    fld,
                    axis,
                    args_dict.get("frames_coord"),
                    zoom,
                    center_xyz,
                    corner,
                )
            auto = _auto_zlim_from_array(win, fld)
            if auto is not None:
                zlims[fld] = [auto[0], auto[1]]
                if args_dict.get("verbose", False):
                    print(f"[zlim-lock] {fld}: {auto[0]:.6g} .. {auto[1]:.6g}")
        except Exception as exc:
            if args_dict.get("verbose", False):
                print(f"WARNING: zlim-lock skipped {fld!r}: {exc}")
            continue
    return zlims
