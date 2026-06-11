#!/usr/bin/env python3
"""
Consume AMReX plotfiles: extract Psi4 mode(s) -> append small-data -> (optional) delete.

This solves the "plot_interval=1 creates gigantic data" problem by keeping only a
tiny time-series on disk:
  - input: plotfile directories (e.g. WormholePlt0001234/)
  - output: ASCII .dat time-series (constraints-style)
  - optional: delete plotfile directories after successful extraction

Requirements:
  - plotfiles must contain ('boxlib','Weyl4_Re') and ('boxlib','Weyl4_Im')
    (in params: amr.derive_plot_vars = Weyl4)
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import yt
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 -- registers 3D projection
from scipy.integrate import cumulative_trapezoid

from grteclyn_wrapper.core.config import VISUALISATION_DIR, default_sim_data_dir

# Frame-rendering cost knobs.  These diagnostic frames (inspected, not
# publication figures) dominate the per-batch CPU post-processing that leaves
# the GPUs idle, so they are tuned for throughput: 110 dpi keeps an 8x7 figure
# legible (~880px) at ~2x fewer pixels than the old 150, and 2 FRB samples per
# simulation cell (down from 4) quarters the yt fixed-resolution-buffer fill
# while still over-resolving the native grid.
_FRAME_DPI = 110
_FRAME_SAMPLES_PER_CELL = 2
_FRAME_BUFF_CAP = 1024


def _default_data_dir() -> str:
    return str(default_sim_data_dir())


def _default_frames_out_dir() -> str:
    return str(VISUALISATION_DIR / "visualize")


def _frames_auto_zlim_enabled(explicit: bool | None = None) -> bool:
    """Per-frame percentile scaling (colorbar jumps). Default: fixed limits for movies."""
    if explicit is not None:
        return explicit
    return os.environ.get("GRTECLYN_FRAMES_AUTO_ZLIM", "").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }


# Fixed color limits keep movie colorbars stable across time steps.
# Override any field via GRTECLYN_FRAMES_ZLIM_<FIELD>=lo,hi (e.g. shift1=-0.05,0.05).
_FIELD_FRAME_CONFIGS: Dict[str, dict] = {
    # chi develops a deep, localised conformal well (min_chi can reach ~0.4 while
    # the far field stays ~1.0).  A narrow fixed window clamps the whole well to
    # the colormap floor and renders the frame blank, so both chi fields use
    # per-frame percentile scaling; the wide preset is only a fallback.
    "chi": {"zlim": (0.3, 1.05), "cmap": "magma", "label": r"Conformal Factor $\chi$", "auto_zlim": True},
    "chi_minus_1": {"zlim": (-0.6, 0.6), "cmap": "RdBu", "label": r"$\chi - 1$", "auto_zlim": True},
    "phi": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\phi$"},
    "Pi": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\Pi$"},
    "phi_lump_sum": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"$\sum_k\phi_k$"},
    "Pi_lump_sum": {"zlim": (-0.01, 0.01), "cmap": "RdBu", "label": r"$\sum_k\Pi_k$"},
    "scalar_activity": {"zlim": (0.0, 0.20), "cmap": "viridis", "label": r"$\sum_k\sqrt{\phi_k^2+\Pi_k^2}$"},
    "lump_activity": {"zlim": (0.0, 0.20), "cmap": "viridis", "label": r"$\sum_k\sqrt{\phi_k^2+\Pi_k^2}$"},
    "local_speed": {"zlim": (0.90, 1.30), "cmap": "magma", "label": r"Local Coordinate Speed"},
    "shift1": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"shift1"},
    "shift2": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"shift2"},
    "shift3": {"zlim": (-0.05, 0.05), "cmap": "RdBu", "label": r"shift3"},
    "rho_req": {"zlim": (-3.0e-3, 3.0e-3), "cmap": "RdBu", "label": r"$\rho_{\mathrm{req}}$"},
    "K": {"zlim": (-5.0e-4, 5.0e-4), "cmap": "RdBu", "label": r"Trace of Extrinsic Curvature $K$"},
    "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": r"Z4 Constraint $\Theta$"},
    "lapse": {"zlim": (0.995, 1.005), "cmap": "viridis", "label": r"Lapse $\alpha$"},
    "Ham": {"zlim": (-0.1, 0.1), "cmap": "RdBu", "label": r"Hamiltonian Constraint $\mathcal{H}$"},
    "A11": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{11}$"},
    "A12": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{12}$"},
    "A22": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{22}$"},
    "A33": {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{33}$"},
    "GW_Plus": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_+$"},
    "GW_Cross": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_\times$"},
    "Weyl4_Re": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Re}(\Psi_4)$"},
    "Weyl4_Im": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Im}(\Psi_4)$"},
    "Weyl4_Mag": {"zlim": (0.0, 1.0e-4), "cmap": "inferno", "label": r"$|\Psi_4|$"},
}


def _field_frame_config(field: str) -> dict:
    cfg = dict(_FIELD_FRAME_CONFIGS.get(field, {"zlim": (None, None), "cmap": "viridis", "label": field}))
    env_key = f"GRTECLYN_FRAMES_ZLIM_{field.upper()}"
    override = os.environ.get(env_key, "").strip()
    if override:
        parts = [p.strip() for p in override.split(",")]
        if len(parts) == 2:
            cfg["zlim"] = (float(parts[0]), float(parts[1]))
    return cfg


def _iter_plotfile_dirs(data_dir: str) -> List[str]:
    """Return sorted plotfile directories under data_dir."""
    out: List[str] = []
    if not os.path.isdir(data_dir):
        return out
    prefixes = (
        "WormholePlt",
        "SupportedWormholePlt",
        "RotatingWormholePlt",
        "RadialRecipePlt",
        "plt",
    )
    for name in os.listdir(data_dir):
        if not any(name.startswith(prefix) for prefix in prefixes):
            continue
        p = os.path.join(data_dir, name)
        if os.path.isdir(p):
            out.append(p)
    out.sort()
    return out


def _parse_plot_index(plot_dir_basename: str) -> int | None:
    """
    Parse trailing integer index from plotfile directory basename.
    Examples: WormholePlt00010 -> 10, plt000123 -> 123
    """
    m = re.search(r"(\d+)$", plot_dir_basename)
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None

def _canonical_field_name(name: str) -> str:
    # Accept a few convenient aliases.
    aliases = {
        "Weyl": "Weyl4_Re",
        "Weyl4": "Weyl4_Re",
        "Weyl_Re": "Weyl4_Re",
        "Weyl_Im": "Weyl4_Im",
        "Weyl_Mag": "Weyl4_Mag",
    }
    return aliases.get(name, name)


def _field_type_for(ds, *field_names: str) -> str:
    """Prefer BoxLib fields, but support yt stream datasets for isolated tests."""
    field_list = list(getattr(ds, "field_list", []))
    derived_field_list = list(getattr(ds, "derived_field_list", []))
    available = field_list + derived_field_list
    for name in field_names:
        if ("boxlib", name) in available:
            return "boxlib"
    for name in field_names:
        if ("stream", name) in available:
            return "stream"
    for name in field_names:
        for ftype, fname in available:
            if fname == name:
                return ftype
    return "boxlib"


def _field_key(ds, field_name: str) -> tuple[str, str]:
    return (_field_type_for(ds, field_name), field_name)


def _should_auto_reset(plot_dirs: List[str], state: Dict[str, bool]) -> bool:
    """
    Heuristic: if the output folder contains a "fresh" plot index (0) but the
    saved state references plotfiles that do not exist anymore, assume the user
    restarted a run in the same directory and reset outputs.
    """
    if not plot_dirs or not state:
        return False
    basenames = [os.path.basename(p) for p in plot_dirs]
    cur_set = set(basenames)
    state_set = set(state.keys())

    cur_idxs = [i for i in (_parse_plot_index(b) for b in basenames) if i is not None]
    if not cur_idxs:
        return False
    min_cur = min(cur_idxs)

    # Restart-like scenario: we see index 0 again, but state refers to old plotfiles.
    if min_cur == 0 and not state_set.issubset(cur_set):
        return True

    # Another restart-like scenario: indices restarted and current max is below
    # what we previously processed.
    state_idxs = [i for i in (_parse_plot_index(b) for b in state_set) if i is not None]
    if min_cur == 0 and state_idxs and max(cur_idxs) < max(state_idxs):
        return True

    return False


def _truncate_if_exists(path: Path) -> None:
    if path.exists() and path.is_file() and path.stat().st_size > 0:
        # Truncate content but keep the file path stable
        path.write_text("", encoding="utf-8")


def _is_plotfile_ready(plot_dir: str, stable_seconds: float) -> bool:
    """Best-effort check that the plotfile is not being written right now."""
    header = os.path.join(plot_dir, "Header")
    if not os.path.isfile(header):
        return False
    try:
        mtime = os.path.getmtime(header)
    except OSError:
        return False
    return (time.time() - mtime) >= stable_seconds


def get_extraction_points(radius: float, n_points: int) -> Tuple[np.ndarray, ...]:
    """Generate sphere points and weights (theta, phi grid)."""
    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")

    X = radius * np.sin(THETA) * np.cos(PHI)
    Y = radius * np.sin(THETA) * np.sin(PHI)
    Z = radius * np.cos(THETA)

    dtheta = np.pi / (n_points - 1)
    dphi = 2.0 * np.pi / n_points
    W = np.sin(THETA) * dtheta * dphi
    return X, Y, Z, THETA, PHI, W


def spin_weighted_sph_harm_s2_l2_m0(theta: np.ndarray) -> np.ndarray:
    """
    Exact (real) spin-weighted spherical harmonic:  _-2Y_{2,0}
    _-2Y20(theta,phi) = sqrt(15/(32*pi)) * sin^2(theta)
    """
    return np.sqrt(15.0 / (32.0 * np.pi)) * (np.sin(theta) ** 2)

def _register_derived_fields(ds, field: str) -> None:
    base_ftype = _field_type_for(ds, "phi", "Pi", "chi")
    if (base_ftype, field) in list(getattr(ds, "field_list", [])) + list(getattr(ds, "derived_field_list", [])):
        return

    available_names = {
        name for (_ftype, name) in list(getattr(ds, "field_list", []))
        + list(getattr(ds, "derived_field_list", []))
    }
    lump_pairs = [
        (f"phi_lump{k}", f"Pi_lump{k}")
        for k in range(5)
        if f"phi_lump{k}" in available_names and f"Pi_lump{k}" in available_names
    ]

    def _sum_lump_component(data, component: str):
        total = None
        for k in range(5):
            name = f"{component}_lump{k}"
            if name not in available_names:
                continue
            term = data[base_ftype, name]
            total = term if total is None else total + term
        return total

    # GW proxy fields
    if field == "GW_Plus":
        def _gw_plus(field, data):
            return data[base_ftype, "A11"] - data[base_ftype, "A22"]
        ds.add_field((base_ftype, "GW_Plus"), function=_gw_plus, sampling_type="cell", units="")
    elif field == "GW_Cross":
        def _gw_cross(field, data):
            return 2.0 * data[base_ftype, "A12"]
        ds.add_field((base_ftype, "GW_Cross"), function=_gw_cross, sampling_type="cell", units="")
    elif field == "Weyl4_Mag":
        def _weyl4_mag(field, data):
            re_v = data[base_ftype, "Weyl4_Re"]
            im_v = data[base_ftype, "Weyl4_Im"]
            return np.sqrt(re_v**2 + im_v**2)
        ds.add_field((base_ftype, "Weyl4_Mag"), function=_weyl4_mag, sampling_type="cell", units="")
    elif field == "chi_minus_1":
        def _chi_minus_1(field, data):
            return data[base_ftype, "chi"] - 1.0
        ds.add_field((base_ftype, "chi_minus_1"), function=_chi_minus_1, sampling_type="cell", units="")
    elif field == "phi_lump_sum":
        def _phi_lump_sum(field, data):
            total = _sum_lump_component(data, "phi")
            if total is not None:
                return total
            return data[base_ftype, "phi"]
        ds.add_field((base_ftype, "phi_lump_sum"), function=_phi_lump_sum, sampling_type="cell", units="")
    elif field == "Pi_lump_sum":
        def _pi_lump_sum(field, data):
            total = _sum_lump_component(data, "Pi")
            if total is not None:
                return total
            return data[base_ftype, "Pi"]
        ds.add_field((base_ftype, "Pi_lump_sum"), function=_pi_lump_sum, sampling_type="cell", units="")
    elif field == "scalar_activity":
        def _scalar_activity(field, data):
            if lump_pairs:
                total = None
                for phi_name, pi_name in lump_pairs:
                    term = np.sqrt(
                        data[base_ftype, phi_name] ** 2 + data[base_ftype, pi_name] ** 2
                    )
                    total = term if total is None else total + term
                return total
            phi = data[base_ftype, "phi"]
            pi = data[base_ftype, "Pi"]
            return np.sqrt(phi**2 + pi**2)
        ds.add_field((base_ftype, "scalar_activity"), function=_scalar_activity, sampling_type="cell", units="")
    elif field == "lump_activity":
        has_combined_scalar = "phi" in available_names and "Pi" in available_names
        ref_field = next(
            (name for name in ("rho_req", "chi", "K", "Weyl4_Re") if name in available_names),
            None,
        )

        def _lump_activity(field, data):
            if lump_pairs:
                total = None
                for phi_name, pi_name in lump_pairs:
                    term = np.sqrt(
                        data[base_ftype, phi_name] ** 2 + data[base_ftype, pi_name] ** 2
                    )
                    total = term if total is None else total + term
                return total
            if has_combined_scalar:
                phi = data[base_ftype, "phi"]
                pi = data[base_ftype, "Pi"]
                return np.sqrt(phi**2 + pi**2)
            if ref_field is None:
                raise RuntimeError(
                    "lump_activity requires phi_lump*/Pi_lump*, phi/Pi, or a reference field"
                )
            return np.zeros_like(data[base_ftype, ref_field])

        ds.add_field((base_ftype, "lump_activity"), function=_lump_activity, sampling_type="cell", units="")
    elif field == "local_speed":
        def _local_speed(field, data):
            eps = 1.0e-12
            chi = data[base_ftype, "chi"]
            lapse = data[base_ftype, "lapse"]
            c1 = np.abs(data[base_ftype, "shift1"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h11"], eps))
            c2 = np.abs(data[base_ftype, "shift2"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h22"], eps))
            c3 = np.abs(data[base_ftype, "shift3"]) + lapse * np.sqrt(chi / np.maximum(data[base_ftype, "h33"], eps))
            return np.maximum(np.maximum(c1, c2), c3)
        ds.add_field((base_ftype, "local_speed"), function=_local_speed, sampling_type="cell", units="")

    elif field == "Weyl4_Re":
        # Ensure base fields are available if asked for explicitly?
        # Usually they are just read from disk.
        pass


def _resolve_frame_physics_center(
    ds,
    axis: str,
    coord: float | None,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    corner: bool,
) -> list[float]:
    """Match visualize/__main__.py center logic; auto throat corner when unset."""
    mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
    mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)
    physics_center = [mid_x, mid_y, 0.0]

    if corner and zoom is not None:
        slice_plane_val = 0.0 if coord is None else float(coord)
        if center_xyz is not None:
            origin = np.array(center_xyz, dtype=float)
        else:
            origin = np.array(
                [mid_x - float(zoom) / 2.0, mid_y - float(zoom) / 2.0, slice_plane_val],
                dtype=float,
            )
        w = float(zoom)
        if axis == "z":
            physics_center = [origin[0] + w / 2.0, origin[1] + w / 2.0, slice_plane_val]
        elif axis == "y":
            physics_center = [origin[0] + w / 2.0, slice_plane_val, origin[2] + w / 2.0]
        elif axis == "x":
            physics_center = [slice_plane_val, origin[1] + w / 2.0, origin[2] + w / 2.0]
    elif center_xyz is not None:
        physics_center = [float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2])]

    if coord is not None:
        if axis == "z":
            physics_center[2] = float(coord)
        elif axis == "y":
            physics_center[1] = float(coord)
        elif axis == "x":
            physics_center[0] = float(coord)

    return physics_center


def _frame_buff_size(ds, zoom: float | None) -> int:
    """Pixels across the FRB: ``_FRAME_SAMPLES_PER_CELL`` samples per simulation
    cell in the zoom window (capped at ``_FRAME_BUFF_CAP``)."""
    if zoom is None:
        zoom = float(ds.domain_width[0])
    le = np.array(ds.domain_left_edge.d, dtype=float)
    re = np.array(ds.domain_right_edge.d, dtype=float)
    dims = np.array(ds.domain_dimensions, dtype=float) * (2 ** int(ds.index.max_level))
    dx = (re - le) / np.maximum(dims, 1.0)
    cells = max(int(round(float(zoom) / min(dx[0], dx[1]))), 64)
    return int(min(max(cells * _FRAME_SAMPLES_PER_CELL, 256), _FRAME_BUFF_CAP))


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

    if field_name in {"Weyl4_Re", "Weyl4_Im", "phi", "Pi", "chi_minus_1"}:
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
    min_span = 1.0e-4 if field_name in {"chi", "lapse"} else 1.0e-12
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

    preset = cfg["zlim"]
    if preset[0] is not None:
        auto = _auto_zlim_from_array(win, field)
        if auto is not None and field in {"lump_activity", "scalar_activity"}:
            signal_span = auto[1] - auto[0]
            preset_span = preset[1] - preset[0]
            if signal_span > 0.0 and signal_span < 0.25 * preset_span:
                return auto
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


def _lock_frame_zlims_from_plotfile(plot_path: str, args_dict: dict) -> dict[str, list[float]]:
    ds = yt.load(plot_path)
    frame_fields = [_canonical_field_name(f) for f in args_dict.get("frames_fields", [])]
    zlims: dict[str, list[float]] = {}
    for fld in frame_fields:
        win = _extract_native_slice_window(
            ds,
            fld,
            args_dict.get("frames_axis", "z"),
            args_dict.get("frames_coord"),
            args_dict.get("frames_zoom"),
            args_dict.get("frames_center"),
            bool(args_dict.get("frames_corner")),
        )
        auto = _auto_zlim_from_array(win, fld)
        if auto is not None:
            zlims[fld] = [auto[0], auto[1]]
            if args_dict.get("verbose", False):
                print(f"[zlim-lock] {fld}: {auto[0]:.6g} .. {auto[1]:.6g}")
    return zlims


def _render_native_slice_frame(
    ds,
    field: str,
    axis: str,
    coord: float | None,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    corner: bool,
    frames_out_dir: str,
    frame_idx: int,
    verbose: bool,
    auto_zlim: bool | None = None,
    frame_zlims: dict[str, list[float]] | None = None,
    use_global_zlim: bool = True,
) -> str:
    """Render from a uniform covering grid (native cell resolution, no FRB blocks)."""
    def _clean_zero(value: float, tol: float = 1.0e-10) -> float:
        value = float(value)
        return 0.0 if abs(value) < tol else value

    cfg = _field_frame_config(field)
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
        win = arr[np.ix_(vmask, hmask)] if hmask.any() and vmask.any() else arr
    else:
        win = arr

    zlim = _resolve_plot_zlim(
        field,
        win,
        cfg,
        auto_zlim=auto_zlim,
        frame_zlims=frame_zlims,
        use_global_zlim=use_global_zlim,
    )

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
    })

    xlabel_name, ylabel_name = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}[axis]
    coord_val = _clean_zero(physics_center[ai])
    title_text = r"%s $\quad t=%.2f \quad %s=%g$" % (
        cfg["label"],
        float(ds.current_time),
        axis,
        coord_val,
    )

    plot_arr = win
    plot_extent = list(extent)
    if zoom is not None and win.shape != arr.shape:
        hh = float(zoom) / 2.0
        hcoords = np.linspace(extent[0], extent[1], arr.shape[1])
        vcoords = np.linspace(extent[2], extent[3], arr.shape[0])
        hmask = (hcoords >= physics_center[h_ax] - hh) & (hcoords <= physics_center[h_ax] + hh)
        vmask = (vcoords >= physics_center[v_ax] - hh) & (vcoords <= physics_center[v_ax] + hh)
        if hmask.any() and vmask.any():
            plot_extent = [
                float(hcoords[hmask][0]),
                float(hcoords[hmask][-1]),
                float(vcoords[vmask][0]),
                float(vcoords[vmask][-1]),
            ]

    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(
        plot_arr,
        origin="lower",
        extent=plot_extent,
        aspect="equal",
        cmap=cfg["cmap"],
        vmin=zlim[0],
        vmax=zlim[1],
        interpolation="nearest",
    )
    ax.set_xlabel(r"$%s$" % xlabel_name)
    ax.set_ylabel(r"$%s$" % ylabel_name)
    ax.set_title(title_text, pad=8)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(cfg["label"])

    x_ticks = ax.get_xticks()
    left_x_native = float(x_ticks[0]) if len(x_ticks) else None

    def _fmt_x(val, _pos):
        return f"{float(val):g}"

    def _fmt_y(val, pos):
        if corner and abs(float(val)) < 1.0e-12:
            return ""
        display = f"{float(val):g}"
        if pos == 0 and left_x_native is not None and display == f"{left_x_native:g}":
            return ""
        return display

    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_x))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_y))

    output_dir = os.path.join(frames_out_dir, f"{field}_{axis}")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    frame_name = f"frame_{axis}_{frame_idx:04d}.png"
    out_path = os.path.join(frames_dir, frame_name)
    fig.savefig(out_path, dpi=_FRAME_DPI, bbox_inches="tight")
    plt.close(fig)

    if verbose:
        print(f"[frame-native] {field} {full_dims} -> {out_path}")

    return out_path


def _render_slice_frame(
    ds,
    field: str,
    axis: str,
    coord: float | None,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    corner: bool,
    frames_out_dir: str,
    frame_idx: int,
    verbose: bool,
    auto_zlim: bool | None = None,
    frame_zlims: dict[str, list[float]] | None = None,
    use_global_zlim: bool = True,
) -> str:
    """
    Render a SlicePlot frame and save it under:
      <frames_out_dir>/<field>_<axis>/frames/frame_<axis>_<idx>.png
    Returns the saved path.
    """
    def _clean_zero(value: float, tol: float = 1.0e-10) -> float:
        value = float(value)
        return 0.0 if abs(value) < tol else value

    def _format_tick(value: float, _pos=None) -> str:
        value = _clean_zero(value)
        return f"{value:g}"

    if int(ds.index.max_level) == 0:
        return _render_native_slice_frame(
            ds,
            field,
            axis,
            coord,
            zoom,
            center_xyz,
            corner,
            frames_out_dir,
            frame_idx,
            verbose,
            auto_zlim=auto_zlim,
            frame_zlims=frame_zlims,
            use_global_zlim=use_global_zlim,
        )

    cfg = _field_frame_config(field)
    physics_center = _resolve_frame_physics_center(ds, axis, coord, zoom, center_xyz, corner)
    plot_center = ds.arr(physics_center, "code_length")

    _register_derived_fields(ds, field)
    plot_field = _field_key(ds, field)

    # Apply scientific style
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
    })

    buff = _frame_buff_size(ds, zoom)
    slc = yt.SlicePlot(ds, axis, plot_field, center=plot_center, buff_size=(buff, buff))
    # Use dataset-native coordinates (e.g. [0,40]) on axes for symmetry-reduced domains.
    slc.set_origin("native")
    slc.set_axes_unit("code_length")
    if zoom is not None:
        slc.set_width((float(zoom), "code_length"))
    slc.set_log(plot_field, False)

    # Force axis labels to be purely LaTeX (overriding yt's defaults)
    # The axes mapping depends on the slice axis.
    # For z-slice: x is x-axis, y is y-axis.
    # For x-slice: y is x-axis, z is y-axis.
    # For y-slice: z is x-axis, x is y-axis.
    # (yt standard right-handed cyclic permutation)
    # Actually, let's verify yt's specific mapping for SlicePlot.
    # z-axis -> x (horizontal), y (vertical)
    # y-axis -> z (horizontal), x (vertical) -- wait, standard yt is usually:
    #   x-axis normal -> y horizontal, z vertical
    #   y-axis normal -> z horizontal, x vertical
    #   z-axis normal -> x horizontal, y vertical
    # Let's assume standard behavior first.

    try:
        zlim = _resolve_plot_zlim(
            field,
            np.asarray(slc.frb[plot_field]),
            cfg,
            auto_zlim=auto_zlim,
            frame_zlims=frame_zlims,
            use_global_zlim=use_global_zlim,
        )
    except Exception:
        zlim = cfg["zlim"]
    if zlim[0] is not None:
        slc.set_zlim(plot_field, zlim[0], zlim[1])
    slc.set_cmap(plot_field, cfg["cmap"])

    coord_val = _clean_zero(physics_center[{"x": 0, "y": 1, "z": 2}[axis]])

    # Format title with LaTeX
    title_text = r"%s $\quad t=%.2f \quad %s=%g$" % (cfg["label"], float(ds.current_time), axis, coord_val)

    # Force LaTeX labels for axes
    axis_map = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
    xlabel_name, ylabel_name = axis_map[axis]
    slc.set_xlabel(r"$%s$" % xlabel_name)
    slc.set_ylabel(r"$%s$" % ylabel_name)
    slc.set_colorbar_label(plot_field, cfg['label'])

    # Render now so we can tweak tick labels on the matplotlib axes.
    # NOTE: yt may re-render inside slc.save(); to make tick tweaks stick we save
    # via matplotlib after rendering.
    slc.render()
    plot = slc.plots[plot_field]
    ax = plot.axes
    ax.set_title(title_text, pad=8)

    for artist in ax.collections:
        if hasattr(artist, "set_linewidth"):
            artist.set_linewidth(0.0)
        if hasattr(artist, "set_edgecolor"):
            artist.set_edgecolor("none")
        if hasattr(artist, "set_antialiased"):
            artist.set_antialiased(False)

    # For centered full-domain RadialRecipe runs, display coordinates relative
    # to the physical object center instead of raw AMReX coordinates.
    if not corner:
        display_offsets = {
            "x": (physics_center[1], physics_center[2]),
            "y": (physics_center[2], physics_center[0]),
            "z": (physics_center[0], physics_center[1]),
        }[axis]
        x_offset, y_offset = display_offsets
        x_ticks = ax.get_xticks()
        left_x_native = float(x_ticks[0]) if len(x_ticks) else None

        def _fmt_x(val, pos):
            return _format_tick(val - x_offset, pos)

        def _fmt_y(val, pos):
            display = _format_tick(val - y_offset, pos)
            # Bottom-left corner: x- and y-axis min ticks overlap (e.g. "-48 -48").
            if pos == 0 and left_x_native is not None:
                left_x_display = _format_tick(left_x_native - x_offset)
                if display == left_x_display:
                    return ""
            return display

        ax.xaxis.set_major_formatter(FuncFormatter(_fmt_x))
        ax.yaxis.set_major_formatter(FuncFormatter(_fmt_y))
        ax.set_xlabel(r"$%s-%s_0$" % (xlabel_name, xlabel_name))
        ax.set_ylabel(r"$%s-%s_0$" % (ylabel_name, ylabel_name))

    # Remove the duplicated "0" tick label at the symmetry origin corner.
    # Keep the x-axis "0" label, hide ONLY the y-axis tick label at y=0.
    # Using a formatter is robust against yt/matplotlib regenerating tick labels.
    if corner:
        def _fmt_y(val, _pos):
            if abs(float(val)) < 1.0e-12:
                return ""
            return f"{val:g}"

        ax.yaxis.set_major_formatter(FuncFormatter(_fmt_y))
        try:
            plot.figure.canvas.draw()
        except Exception:
            pass

    output_dir = os.path.join(frames_out_dir, f"{field}_{axis}")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    frame_name = f"frame_{axis}_{frame_idx:04d}.png"
    out_path = os.path.join(frames_dir, frame_name)
    plot.figure.savefig(out_path, dpi=_FRAME_DPI, bbox_inches="tight")
    plt.close(plot.figure)

    if verbose:
        print(f"[frame] {field} -> {out_path}")

    return out_path


def _render_projection_frame(
    ds,
    field: str,
    axis: str,
    method: str,
    zoom: float | None,
    center_xyz: Sequence[float] | None,
    frames_out_dir: str,
    frame_idx: int,
    verbose: bool,
) -> str:
    """Render a line-of-sight projection frame for 3D matter placement."""
    cfg = _field_frame_config(field)

    mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
    mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)
    physics_center = [mid_x, mid_y, 0.0]
    if center_xyz is not None:
        physics_center = [float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2])]

    _register_derived_fields(ds, field)
    plot_field = _field_key(ds, field)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
    })

    plot_center = ds.arr(physics_center, "code_length")
    proj = yt.ProjectionPlot(
        ds,
        axis,
        plot_field,
        center=plot_center,
        method=method,
    )
    proj.set_origin("native")
    proj.set_axes_unit("code_length")
    if zoom is not None:
        width = float(zoom)
        z_left = float(ds.domain_left_edge[2])
        z_right = float(ds.domain_right_edge[2])
        # RadialRecipe/GRTresna evolves a half-z domain with reflective z=0.
        # For x/y projections, a symmetric width around z0=0 asks yt to render
        # z<0, which appears as a white cut-off band. Use the stored z extent
        # instead, while preserving the requested transverse zoom.
        if axis in {"x", "y"} and physics_center[2] <= z_left + 1.0e-12:
            z_width = min(width, z_right - z_left)
            if axis == "x":
                proj.set_width((width, z_width))
            else:
                proj.set_width((z_width, width))
        else:
            proj.set_width((width, "code_length"))
    proj.set_log(plot_field, False)
    if cfg["zlim"][0] is not None:
        proj.set_zlim(plot_field, cfg["zlim"][0], cfg["zlim"][1])
    proj.set_cmap(plot_field, cfg["cmap"])
    proj.set_colorbar_label(plot_field, cfg["label"])

    axis_map = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
    xlabel_name, ylabel_name = axis_map[axis]
    proj.set_xlabel(r"$%s$" % xlabel_name)
    proj.set_ylabel(r"$%s$" % ylabel_name)
    title_text = r"%s $\quad t=%.2f \quad \mathrm{%s\ projection\ along}\ %s$" % (
        cfg["label"],
        float(ds.current_time),
        method,
        axis,
    )

    proj.render()
    plot = proj.plots[plot_field]
    ax = plot.axes
    ax.set_title(title_text, pad=8)

    display_offsets = {
        "x": (physics_center[1], physics_center[2]),
        "y": (physics_center[2], physics_center[0]),
        "z": (physics_center[0], physics_center[1]),
    }[axis]
    x_offset, y_offset = display_offsets
    x_ticks = ax.get_xticks()
    left_x_native = float(x_ticks[0]) if len(x_ticks) else None

    def _proj_fmt_x(val, pos):
        return f"{float(val - x_offset):g}"

    def _proj_fmt_y(val, pos):
        display = f"{float(val - y_offset):g}"
        if pos == 0 and left_x_native is not None:
            left_x_display = f"{float(left_x_native - x_offset):g}"
            if display == left_x_display:
                return ""
        return display

    ax.xaxis.set_major_formatter(FuncFormatter(_proj_fmt_x))
    ax.yaxis.set_major_formatter(FuncFormatter(_proj_fmt_y))
    ax.set_xlabel(r"$%s-%s_0$" % (xlabel_name, xlabel_name))
    ax.set_ylabel(r"$%s-%s_0$" % (ylabel_name, ylabel_name))

    output_dir = os.path.join(frames_out_dir, f"{field}_proj_{axis}")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    frame_name = f"frame_proj_{axis}_{frame_idx:04d}.png"
    out_path = os.path.join(frames_dir, frame_name)
    plot.figure.savefig(out_path, dpi=_FRAME_DPI, bbox_inches="tight")
    plt.close(plot.figure)

    if verbose:
        print(f"[projection] {field} -> {out_path}")

    return out_path

def _cleanup_existing_frames(frames_out_dir: str, fields: Sequence[str], axis: str, verbose: bool) -> None:
    """
    Remove existing PNG frames (and movie mp4) for the requested field+axis outputs.

    This mirrors the behavior of visualize/__main__.py which clears frames before
    creating a new animation run, but is scoped to the fields requested via
    --frames-fields (so it won't touch other outputs).
    """
    base = Path(frames_out_dir)
    for fld in fields:
        out_dir = base / f"{fld}_{axis}"
        frames_dir = out_dir / "frames"
        if frames_dir.is_dir():
            for p in frames_dir.glob("*.png"):
                try:
                    p.unlink()
                except FileNotFoundError:
                    pass
            if verbose:
                print(f"[clean] cleared frames in {frames_dir}")
        # Also remove a previously stitched movie for this field/axis if present.
        movie = out_dir / f"movie_{fld}_{axis}.mp4"
        if movie.exists():
            try:
                movie.unlink()
                if verbose:
                    print(f"[clean] removed {movie}")
            except FileNotFoundError:
                pass


def _cleanup_projection_frames(frames_out_dir: str, fields: Sequence[str], axes: Sequence[str], verbose: bool) -> None:
    """Remove existing projection PNG frames for requested field+axis outputs."""
    base = Path(frames_out_dir)
    for fld in fields:
        for axis in axes:
            out_dir = base / f"{fld}_proj_{axis}"
            frames_dir = out_dir / "frames"
            if frames_dir.is_dir():
                for p in frames_dir.glob("*.png"):
                    try:
                        p.unlink()
                    except FileNotFoundError:
                        pass
                if verbose:
                    print(f"[clean] cleared projection frames in {frames_dir}")
            movie = out_dir / f"movie_{fld}_proj_{axis}.mp4"
            if movie.exists():
                try:
                    movie.unlink()
                    if verbose:
                        print(f"[clean] removed {movie}")
                except FileNotFoundError:
                    pass


def _extract_mode_amps_l2m0(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
) -> List[complex]:
    """
    Extract (l=2,m=0) mode amplitudes for multiple radii in ONE batched yt query.

    This is much faster than querying radii one-by-one because AMR indexing and
    disk reads are amortized.
    """
    if not radii:
        return []

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

    # Precompute unit-sphere grid (theta,phi) and weights (independent of radius).
    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    sinT = np.sin(THETA)
    X1 = sinT * np.cos(PHI)
    Y1 = sinT * np.sin(PHI)
    Z1 = np.cos(THETA)
    dtheta = np.pi / (n_points - 1)
    dphi = 2.0 * np.pi / n_points
    W = sinT * dtheta * dphi
    Wf = W.ravel()
    Y20 = spin_weighted_sph_harm_s2_l2_m0(THETA)
    shape = X1.shape

    # Build one big list of points for all radii.
    starts: List[int] = []
    ends: List[int] = []
    sx_all: List[np.ndarray] = []
    sy_all: List[np.ndarray] = []
    sz_all: List[np.ndarray] = []
    cursor = 0
    for r in radii:
        starts.append(cursor)
        n = shape[0] * shape[1]
        cursor += n
        ends.append(cursor)
        sx_all.append((float(r) * X1).ravel() + center[0])
        sy_all.append((float(r) * Y1).ravel() + center[1])
        sz_all.append((float(r) * Z1).ravel() + center[2])

    sx = np.concatenate(sx_all)
    sy = np.concatenate(sy_all)
    sz = np.concatenate(sz_all)

    in_domain = (
        (sx >= left[0])
        & (sx <= right[0])
        & (sy >= left[1])
        & (sy <= right[1])
        & (sz >= left[2])
        & (sz <= right[2])
    )
    idxs = np.where(in_domain)[0]
    if idxs.size < 4:
        raise RuntimeError("Too few sphere points inside domain (all radii).")

    weyl = np.full(sx.shape, np.nan + 1j * np.nan, dtype=np.complex128)
    pts = np.column_stack((sx[idxs], sy[idxs], sz[idxs]))

    def _coerce_scalar_samples(values) -> np.ndarray:
        arr = np.asarray(values)
        if arr.dtype != object and arr.ndim == 1:
            return arr.astype(float, copy=False)

        out = np.empty(len(values), dtype=float)
        for jj, v in enumerate(values):
            vv = np.asarray(v, dtype=float).reshape(-1)
            out[jj] = vv[0] if vv.size else np.nan
        return out

    try:
        vals = ds.find_field_values_at_points(
            [("boxlib", "Weyl4_Re"), ("boxlib", "Weyl4_Im")], pts
        )
        re_vals = _coerce_scalar_samples(vals[0])
        im_vals = _coerce_scalar_samples(vals[1])
        weyl[idxs] = re_vals + 1j * im_vals
    except Exception:
        # Fallback: per-point sampling (slow, but should still work)
        for ii, i in enumerate(idxs):
            pt = ds.point([pts[ii, 0], pts[ii, 1], pts[ii, 2]])
            re = np.asarray(pt[("boxlib", "Weyl4_Re")])
            im = np.asarray(pt[("boxlib", "Weyl4_Im")])
            if re.size and im.size:
                weyl[i] = float(re.flat[0]) + 1j * float(im.flat[0])

    amps: List[complex] = []
    for r, s, e in zip(radii, starts, ends):
        weyl_r = weyl[s:e]
        mask = np.isfinite(weyl_r.real) & np.isfinite(weyl_r.imag)
        if np.sum(mask) < 4:
            raise RuntimeError(f"Too few valid Weyl4 samples at R={r}.")

        omega_tot = np.sum(Wf[mask])
        if omega_tot <= 0.0:
            raise RuntimeError(f"Invalid integration weights at R={r} (omega_tot <= 0).")

        W_renorm = (Wf / omega_tot * 4.0 * np.pi).reshape(shape)
        psi4 = np.where(mask, weyl_r, 0.0).reshape(shape)
        amp = np.sum(psi4 * np.conj(Y20) * W_renorm) * float(r)
        amps.append(complex(amp))

    return amps


def _shell_stats_header(radii: Sequence[float], fields: Sequence[str]) -> str:
    cols = ["time"]
    for radius in radii:
        tag = f"R{radius:g}"
        for field in fields:
            cols.extend(
                [
                    f"{field}_mean_{tag}",
                    f"{field}_min_{tag}",
                    f"{field}_max_{tag}",
                ]
            )
    return "# " + "  ".join(cols)


def _extract_shell_field_stats(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
    fields: Sequence[str],
) -> Dict[str, Tuple[float, float, float]]:
    """Sample mean/min/max of fields on spherical shells at given radii."""
    if not radii or not fields:
        return {}

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    sinT = np.sin(THETA)
    X1 = sinT * np.cos(PHI)
    Y1 = sinT * np.sin(PHI)
    Z1 = np.cos(THETA)

    yt_fields = []
    for field in fields:
        key = ("boxlib", field)
        if key not in ds.field_list:
            raise RuntimeError(f"Plotfile missing field {field!r} for shell extraction.")
        yt_fields.append(key)

    out: Dict[str, Tuple[float, float, float]] = {}
    for radius in radii:
        try:
            sx = (float(radius) * X1).ravel() + center[0]
            sy = (float(radius) * Y1).ravel() + center[1]
            sz = (float(radius) * Z1).ravel() + center[2]
            in_domain = (
                (sx >= left[0])
                & (sx <= right[0])
                & (sy >= left[1])
                & (sy <= right[1])
                & (sz >= left[2])
                & (sz <= right[2])
            )
            idxs = np.where(in_domain)[0]
            if idxs.size < 4:
                continue

            pts = np.column_stack((sx[idxs], sy[idxs], sz[idxs]))
            vals = ds.find_field_values_at_points(yt_fields, pts)
            tag = f"R{radius:g}"
            for field, samples in zip(fields, vals):
                arr = np.asarray(samples, dtype=float).reshape(-1)
                arr = arr[np.isfinite(arr)]
                if arr.size < 4:
                    continue
                out[f"{field}_mean_{tag}"] = (
                    float(np.mean(arr)),
                    float(np.min(arr)),
                    float(np.max(arr)),
                )
        except Exception:
            continue
    if not out:
        raise RuntimeError("Shell extraction produced no valid samples at any radius.")
    return out


def _format_shell_stats_line(
    t: float,
    stats: Dict[str, Tuple[float, float, float]],
    radii: Sequence[float],
    fields: Sequence[str],
) -> str:
    parts = [f"{t:.16e}"]
    for radius in radii:
        tag = f"R{radius:g}"
        for field in fields:
            key = f"{field}_mean_{tag}"
            if key not in stats:
                parts.extend(["nan", "nan", "nan"])
                continue
            mean_v, min_v, max_v = stats[key]
            parts.extend([f"{mean_v:.16e}", f"{min_v:.16e}", f"{max_v:.16e}"])
    return "  ".join(parts)


def _extract_areal_radius_min(
    ds,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    chi_floor: float = 1.0e-8,
) -> Tuple[float, float]:
    """Extract the minimum areal radius R_areal = r / sqrt(chi) along the x-axis.

    Returns (R_areal_min, r_at_min).
    """
    right = float(ds.domain_right_edge[0])
    c = np.asarray(center, dtype=float)
    ray = ds.ray(c, [right, c[1], c[2]])

    r_arr = np.asarray(ray[("index", "x")], dtype=float) - c[0]
    chi_arr = np.asarray(ray[("boxlib", "chi")], dtype=float)

    order = np.argsort(r_arr)
    r_arr = r_arr[order]
    chi_arr = chi_arr[order]

    chi_arr = np.maximum(chi_arr, chi_floor)
    R_areal = r_arr / np.sqrt(chi_arr)

    skip = r_arr > 1.0e-12
    if not np.any(skip):
        return (0.0, 0.0)
    R_areal_valid = R_areal[skip]
    r_valid = r_arr[skip]
    i_min = np.argmin(R_areal_valid)
    return (float(R_areal_valid[i_min]), float(r_valid[i_min]))


def _extract_embedding_profile(
    ds,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    chi_floor: float = 1.0e-8,
    r_max: float | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute the isometric embedding profile (R_areal, z_embed) from chi along x-axis.

    The spatial metric is dl^2 = (1/chi) dr^2.  The embedding surface satisfies
    dR^2 + dz^2 = dl^2,  where R_areal = r / sqrt(chi).

    In isotropic coordinates, the throat (minimum R_areal) is at r = b0/2, NOT
    at the grid origin r=0 (which represents the compactified second universe's
    asymptotic infinity).  We locate the throat, then integrate z_embed outward
    in both directions to produce a symmetric two-funnel embedding with the
    throat at z_embed=0.

    The raw AMR ray data is interpolated onto a uniform fine grid to eliminate
    step artifacts at refinement-level boundaries.

    Returns arrays (R_areal, z_embed) suitable for surface-of-revolution plotting.
    """
    right = float(ds.domain_right_edge[0])
    c = np.asarray(center, dtype=float)
    ray = ds.ray(c, [right, c[1], c[2]])

    r_raw = np.asarray(ray[("index", "x")], dtype=float) - c[0]
    chi_raw = np.asarray(ray[("boxlib", "chi")], dtype=float)

    order = np.argsort(r_raw)
    r_raw = r_raw[order]
    chi_raw = chi_raw[order]

    chi_raw = np.maximum(chi_raw, chi_floor)

    valid = r_raw > 1.0e-12
    r_raw = r_raw[valid]
    chi_raw = chi_raw[valid]

    if len(r_raw) < 3:
        return np.array([]), np.array([])

    r_end = float(r_max) if r_max is not None else float(r_raw[-1])
    n_pts = max(1024, len(r_raw) * 4)
    r_arr = np.linspace(float(r_raw[0]), r_end, n_pts)
    chi_arr = np.interp(r_arr, r_raw, chi_raw)
    chi_arr = np.maximum(chi_arr, chi_floor)

    from scipy.ndimage import gaussian_filter1d
    dr = r_arr[1] - r_arr[0]
    dx_fine = float(np.min(np.diff(r_raw))) if len(r_raw) > 1 else dr
    sigma_pts = max(2.0, 3.0 * dx_fine / dr)
    chi_arr = gaussian_filter1d(chi_arr, sigma=sigma_pts, mode="nearest")
    chi_arr = np.maximum(chi_arr, chi_floor)

    R_areal = r_arr / np.sqrt(chi_arr)
    R_areal = gaussian_filter1d(R_areal, sigma=sigma_pts, mode="nearest")

    dR_dr = np.gradient(R_areal, r_arr)
    dl_dr_sq = 1.0 / chi_arr
    dz_dr_sq = dl_dr_sq - dR_dr**2
    dz_dr = np.sqrt(np.maximum(dz_dr_sq, 0.0))

    # Find the throat (minimum R_areal) and build the outer funnel.
    # The Ellis-Bronnikov geometry is symmetric under r <-> b0^2/(4r),
    # so both funnels have identical embedding profiles when parameterized
    # by R_areal.  The inner funnel (r -> 0, second universe) is poorly
    # resolved numerically (chi hits the floor), so we use the well-resolved
    # outer funnel and mirror it to produce the classic symmetric embedding.
    i_throat = int(np.argmin(R_areal))

    R_outer = R_areal[i_throat:]
    dz_outer = dz_dr[i_throat:]
    z_outer = np.zeros(len(R_outer))
    if len(R_outer) > 1:
        z_outer[1:] = cumulative_trapezoid(dz_outer, r_arr[i_throat:])

    R_full = np.concatenate([R_outer[::-1], R_outer[1:]])
    z_full = np.concatenate([-z_outer[::-1], z_outer[1:]])

    return R_full, z_full


def _render_embedding_frame(
    ds,
    frames_out_dir: str,
    frame_idx: int,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    r_max: float | None = None,
    verbose: bool = False,
) -> str:
    """Render a 3D embedding-diagram frame (surface of revolution) and save as PNG."""
    R_areal, z_embed = _extract_embedding_profile(ds, center=center, r_max=r_max)
    if len(R_areal) < 3:
        raise RuntimeError("Too few points for embedding diagram.")

    if r_max is None:
        r_max = float(ds.domain_right_edge[0]) - center[0]

    lim_xy = float(r_max) * 1.15
    lim_z = float(r_max) * 0.8

    # Clip embedding data to the z-axis display range. Near the grid origin
    # (compactified second universe) chi -> 0 produces enormous dl/dr,
    # making z_embed blow up while R_areal remains large. Matplotlib 3D does
    # not clip surfaces to axis limits, so unclipped data obscures the throat.
    mask = np.abs(z_embed) <= lim_z * 1.05
    if mask.sum() >= 3:
        R_areal = R_areal[mask]
        z_embed = z_embed[mask]
        lim_xy = min(lim_xy, float(np.max(R_areal)) * 1.15)

    phi = np.linspace(0, 2 * np.pi, 80)
    PHI, _ = np.meshgrid(phi, np.arange(len(R_areal)))
    X_3d = R_areal[:, None] * np.cos(PHI)
    Y_3d = R_areal[:, None] * np.sin(PHI)
    Z_3d = np.broadcast_to(z_embed[:, None], X_3d.shape)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif", "Times New Roman", "serif"],
        "mathtext.fontset": "cm",
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "axes.linewidth": 1.2,
    })

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")

    norm = plt.Normalize(vmin=-lim_z, vmax=lim_z)
    colors = plt.cm.viridis(norm(Z_3d))
    ax.plot_surface(X_3d, Y_3d, Z_3d, facecolors=colors, shade=True, alpha=0.85,
                    rstride=2, cstride=2, linewidth=0)

    ax.set_xlim(-lim_xy, lim_xy)
    ax.set_ylim(-lim_xy, lim_xy)
    ax.set_zlim(-lim_z, lim_z)

    ax.view_init(elev=25, azim=-60)
    t_sim = float(ds.current_time)
    ax.set_title(r"$\mathrm{Embedding\;Diagram}\quad t=%.3f$" % t_sim, fontsize=16)
    ax.set_xlabel(r"$x$", fontsize=14)
    ax.set_ylabel(r"$y$", fontsize=14)
    ax.set_zlabel(r"$z_{\mathrm{embed}}$", fontsize=14)

    from matplotlib.ticker import FuncFormatter
    math_fmt = FuncFormatter(lambda v, _: r"$%g$" % v)
    ax.xaxis.set_major_formatter(math_fmt)
    ax.yaxis.set_major_formatter(math_fmt)
    ax.zaxis.set_major_formatter(math_fmt)

    output_dir = os.path.join(frames_out_dir, "embedding")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    frame_name = f"frame_{frame_idx:04d}.png"
    out_path = os.path.join(frames_dir, frame_name)
    fig.savefig(out_path, dpi=_FRAME_DPI, bbox_inches="tight", pad_inches=0.3)
    plt.close(fig)

    if verbose:
        print(f"[embedding] t={t_sim:.4f} -> {out_path}")

    return out_path


def _cleanup_embedding_frames(frames_out_dir: str, verbose: bool) -> None:
    """Remove existing embedding PNG frames and movie."""
    base = Path(frames_out_dir) / "embedding" / "frames"
    if base.is_dir():
        for p in base.glob("*.png"):
            try:
                p.unlink()
            except FileNotFoundError:
                pass
        if verbose:
            print(f"[clean] cleared embedding frames in {base}")
    movie = Path(frames_out_dir) / "embedding" / "movie_embedding.mp4"
    if movie.exists():
        try:
            movie.unlink()
        except FileNotFoundError:
            pass


def _load_state(state_path: Path) -> Dict[str, bool]:
    if not state_path.exists():
        return {}
    try:
        return json.loads(state_path.read_text())
    except Exception:
        return {}


def _save_state(state_path: Path, state: Dict[str, bool]) -> None:
    state_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = state_path.with_suffix(state_path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True))
    tmp.replace(state_path)


def _append_line(path: Path, header: str, line: str) -> None:
    is_new = not path.exists()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        if is_new:
            f.write(header.rstrip() + "\n")
        f.write(line.rstrip() + "\n")


def _process_single_plotfile(p: str, args_dict: dict, protected: set, fallback_frame_idx: int) -> dict:
    import yt
    
    try:
        yt.funcs.mylog.setLevel(30 if args_dict.get("verbose") else 40)
    except Exception:
        pass

    t0 = time.time()
    result = {
        "p": p,
        "key": os.path.basename(p),
        "t": 0.0,
        "psi4_line": None,
        "areal_line": None,
        "shell_line": None,
        "boundary_flux_line": None,
        "success": False,
        "deleted": False,
        "status_str": "",
        "dt_s": 0.0,
        "error": None,
    }

    try:
        ds = yt.load(p)
        t = float(ds.current_time)
        result["t"] = t
        key = result["key"]

        if args_dict.get("boundary_flux", True):
            try:
                from grteclyn_wrapper.metrics.probes.boundary import extract_scalar_boundary_flux

                flux = extract_scalar_boundary_flux(p)
                if flux is not None:
                    net, psi4_amp = flux
                    psi4_val = psi4_amp if psi4_amp is not None else 0.0
                    result["boundary_flux_line"] = f"{t:.16e}  {net:.16e}  {psi4_val:.16e}"
            except Exception as exc:
                if args_dict.get("verbose", False):
                    print(f"WARNING: boundary flux extraction failed for {key}: {exc}")

        if args_dict.get("psi4"):
            if ("boxlib", "Weyl4_Re") not in ds.field_list or ("boxlib", "Weyl4_Im") not in ds.field_list:
                raise RuntimeError("Plotfile missing Weyl4_Re/Im. Set: amr.derive_plot_vars = Weyl4 and re-run.")
            amps = _extract_mode_amps_l2m0(
                ds,
                radii=[float(r) for r in args_dict["radii"]],
                n_points=int(args_dict["n_points"]),
                center=args_dict["center"],
            )
            result["psi4_line"] = f"{t:.16e}  " + "  ".join([f"{a.real:.16e}  {a.imag:.16e}" for a in amps])

        shell_fields = list(args_dict.get("shell_fields") or [])
        if shell_fields:
            try:
                stats = _extract_shell_field_stats(
                    ds,
                    radii=[float(r) for r in args_dict["radii"]],
                    n_points=int(args_dict["n_points"]),
                    center=args_dict["center"],
                    fields=shell_fields,
                )
                result["shell_line"] = _format_shell_stats_line(
                    t,
                    stats,
                    [float(r) for r in args_dict["radii"]],
                    shell_fields,
                )
            except Exception as exc:
                if args_dict.get("verbose", False):
                    print(f"WARNING: shell extraction failed for {key}: {exc}")

        frame_fields = [_canonical_field_name(f) for f in args_dict.get("frames_fields", [])]
        if frame_fields:
            idx = _parse_plot_index(key)
            frame_idx = idx if idx is not None else fallback_frame_idx
            for fld in frame_fields:
                _render_slice_frame(
                    ds,
                    field=fld,
                    axis=args_dict["frames_axis"],
                    coord=args_dict.get("frames_coord"),
                    zoom=args_dict.get("frames_zoom"),
                    center_xyz=args_dict.get("frames_center"),
                    corner=args_dict.get("frames_corner"),
                    frames_out_dir=args_dict["frames_out"],
                    frame_idx=int(frame_idx),
                    verbose=args_dict.get("verbose", False),
                    auto_zlim=args_dict.get("frames_auto_zlim"),
                    frame_zlims=args_dict.get("frame_zlims"),
                    use_global_zlim=args_dict.get("frames_global_zlim", True),
                )

        projection_fields = [_canonical_field_name(f) for f in args_dict.get("projection_fields", [])]
        projection_axes = list(args_dict.get("projection_axes", []) or [])
        if projection_fields and projection_axes:
            idx = _parse_plot_index(key)
            frame_idx = idx if idx is not None else fallback_frame_idx
            for fld in projection_fields:
                for axis in projection_axes:
                    _render_projection_frame(
                        ds,
                        field=fld,
                        axis=axis,
                        method=args_dict.get("projection_method", "mip"),
                        zoom=args_dict.get("frames_zoom"),
                        center_xyz=args_dict.get("frames_center"),
                        frames_out_dir=args_dict["frames_out"],
                        frame_idx=int(frame_idx),
                        verbose=args_dict.get("verbose", False),
                    )

        if args_dict.get("areal_radius"):
            if ("boxlib", "chi") in ds.field_list:
                try:
                    R_min, r_min = _extract_areal_radius_min(ds, center=args_dict["center"])
                    result["areal_line"] = f"{t:.16e}  {R_min:.16e}  {r_min:.16e}"
                except Exception as exc:
                    if args_dict.get("verbose", False):
                        print(f"WARNING: areal extraction failed for {key}: {exc}")
            elif args_dict.get("verbose", False):
                print(f"WARNING: plotfile {key} missing chi field; skipping areal radius.")

        if args_dict.get("embedding"):
            if ("boxlib", "chi") in ds.field_list:
                e_idx = _parse_plot_index(key)
                embed_frame_idx = e_idx if e_idx is not None else fallback_frame_idx
                _render_embedding_frame(
                    ds,
                    frames_out_dir=args_dict["frames_out"],
                    frame_idx=int(embed_frame_idx),
                    center=args_dict["center"],
                    r_max=args_dict.get("embedding_rmax"),
                    verbose=args_dict.get("verbose", False),
                )
            else:
                if args_dict.get("verbose", False):
                    print(f"WARNING: plotfile {key} missing chi field; skipping embedding.")

        result["success"] = bool(
            result["psi4_line"]
            or result["areal_line"]
            or result["shell_line"]
            or frame_fields
            or projection_fields
        )

        if args_dict.get("delete") and (p not in protected):
            shutil.rmtree(p)
            result["deleted"] = True

    except Exception as e:
        result["error"] = str(e)

    result["dt_s"] = time.time() - t0
    status = "deleted" if result["deleted"] else ("kept" if args_dict.get("delete") else "kept(no-delete)")
    result["status_str"] = status
    
    return result


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract Psi4 (l=2,m=0) from plotfiles and optionally delete them."
    )
    parser.add_argument("--data", default=_default_data_dir(), help="Directory containing WormholePlt*/plt*")
    parser.add_argument("--out", default=None, help="Output directory (default: <data>/small_data)")
    parser.add_argument("--radii", type=float, nargs="+", default=[14.0, 30.0], help="Extraction radii")
    parser.add_argument("--n-points", type=int, default=64, help="Angular resolution N (N×N points)")
    parser.add_argument(
        "--center",
        type=float,
        nargs=3,
        default=[0.0, 0.0, 0.0],
        help="Extraction center (x y z) in code units",
    )
    parser.add_argument("--stable-seconds", type=float, default=5.0, help="Require Header mtime older than this")
    parser.add_argument("--poll-seconds", type=float, default=2.0, help="Polling interval when --watch")
    parser.add_argument("--watch", action="store_true", help="Keep running and process new plotfiles")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-plotfile details (name, time, delete/keep, timings).",
    )
    parser.add_argument(
        "--psi4",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Enable/disable Psi4 mode extraction to .dat (default: enabled).",
    )
    parser.add_argument(
        "--shell-fields",
        nargs="+",
        default=[],
        help="Extract mean/min/max on spherical shells for these fields (e.g. chi lapse K).",
    )
    parser.add_argument(
        "--frames-fields",
        nargs="+",
        default=[],
        help="Render SlicePlot frames for these fields (e.g. chi K Weyl4_Re).",
    )
    parser.add_argument("--frames-axis", default="z", choices=["x", "y", "z"], help="Slice axis for frames.")
    parser.add_argument("--frames-coord", type=float, default=None, help="Slice coordinate for frames (axis coordinate).")
    parser.add_argument("--frames-zoom", type=float, default=None, help="Zoom width in code units for frames.")
    parser.add_argument("--frames-center", type=float, nargs=3, default=None, help="Center (x y z) for frames.")
    parser.add_argument("--frames-corner", action="store_true", help="Corner mode for symmetry-reduced domains (frames).")
    parser.add_argument(
        "--frames-auto-zlim",
        action="store_true",
        help="Per-frame colorbar from slice percentiles (faint fields only; causes movie flicker).",
    )
    parser.add_argument(
        "--frames-global-zlim",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Lock per-field colorbar limits from the first plotfile (stable movies, visible faint fields).",
    )
    parser.add_argument("--frames-out", default=_default_frames_out_dir(), help="Frames output base dir (default: grteclyn_wrapper/visualisation/visualize).")
    parser.add_argument(
        "--projection-fields",
        nargs="+",
        default=[],
        help="Render 3D line-of-sight projection frames for these fields.",
    )
    parser.add_argument(
        "--projection-axes",
        nargs="+",
        choices=["x", "y", "z"],
        default=[],
        help="Axes to project along for --projection-fields.",
    )
    parser.add_argument(
        "--projection-method",
        choices=["mip", "integrate", "sum"],
        default="mip",
        help="Projection method. mip is best for locating compact blobs.",
    )
    parser.add_argument(
        "--areal-radius",
        action="store_true",
        help="Extract minimum areal radius R_areal = r/sqrt(chi) along x-axis to areal_radius.dat.",
    )
    parser.add_argument(
        "--embedding",
        action="store_true",
        help="Render 3D embedding-diagram frames (surface of revolution from chi profile).",
    )
    parser.add_argument(
        "--embedding-rmax",
        type=float,
        default=None,
        help="Maximum coordinate radius for embedding diagram (default: full domain).",
    )
    parser.add_argument(
        "--delete",
        action="store_true",
        help="Delete plotfile directory after successful extraction",
    )
    parser.add_argument(
        "--keep-last",
        type=int,
        default=2,
        help="Never delete the newest N plotfiles (safety, default: 2)",
    )
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=1,
        help="Number of parallel worker processes to use",
    )
    parser.add_argument(
        "--keep-existing-frames",
        action="store_true",
        help="Do not delete existing frames at startup.",
    )
    args = parser.parse_args()

    # Reduce yt logging overhead/spam (can be noisy in watch mode).
    try:
        # WARNING by default; INFO is extremely verbose.
        yt.funcs.mylog.setLevel(30 if args.verbose else 40)  # WARNING/ERROR
    except Exception:
        pass

    data_dir = os.path.abspath(args.data)
    out_dir = Path(args.out) if args.out else Path(data_dir) / "small_data"
    state_path = out_dir / "consume_state.json"
    out_path = out_dir / "psi4_mode_l2m0.dat"
    areal_out_path = out_dir / "areal_radius.dat"
    shell_out_path = out_dir / "shell_profiles.dat"
    boundary_flux_out_path = out_dir / "boundary_flux.dat"
    header = "# time  " + "  ".join([f"Re(R={R:g})  Im(R={R:g})" for R in args.radii])
    areal_header = "# time  R_areal_min  r_at_R_areal_min"
    shell_header = _shell_stats_header(args.radii, args.shell_fields)

    state = _load_state(state_path)

    # Auto-reset on simulation restart (same output dir reused):
    plot_dirs_now = _iter_plotfile_dirs(data_dir)
    if _should_auto_reset(plot_dirs_now, state):
        print("Detected a likely simulation restart in the same output directory.")
        print(f"Resetting: {out_path} and {state_path}")
        _truncate_if_exists(out_path)
        _truncate_if_exists(areal_out_path)
        if args.shell_fields:
            _truncate_if_exists(shell_out_path)
        _save_state(state_path, {})
        state = {}

    # If rendering frames, clear existing frames for the requested fields/axis at startup.
    frame_fields_startup = [_canonical_field_name(f) for f in args.frames_fields]
    if frame_fields_startup and not args.keep_existing_frames:
        _cleanup_existing_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            fields=frame_fields_startup,
            axis=args.frames_axis,
            verbose=args.verbose,
        )

    projection_fields_startup = [_canonical_field_name(f) for f in args.projection_fields]
    if projection_fields_startup and args.projection_axes and not args.keep_existing_frames:
        _cleanup_projection_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            fields=projection_fields_startup,
            axes=args.projection_axes,
            verbose=args.verbose,
        )

    if args.embedding and not args.keep_existing_frames:
        _cleanup_embedding_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            verbose=args.verbose,
        )

    def process_once() -> int:
        plot_dirs = _iter_plotfile_dirs(data_dir)
        if not plot_dirs:
            return 0
        processed_count = 0

        # Never delete newest keep_last plotfiles. With keep_last=0, protect none.
        keep_last = max(0, int(args.keep_last))
        protected = set(plot_dirs[-keep_last:]) if keep_last > 0 else set()

        to_process = []
        for p in plot_dirs:
            key = os.path.basename(p)
            if state.get(key):
                continue
            if not _is_plotfile_ready(p, stable_seconds=float(args.stable_seconds)):
                continue
            to_process.append(p)

        if not to_process:
            return 0

        print(
            f"Processing {len(to_process)} plotfile(s) "
            f"(~1–5 min each on large AMR; use -j 1 for reliability)...",
            flush=True,
        )

        args_dict = vars(args)
        args_dict["frames_out"] = os.path.abspath(args.frames_out)
        frame_zlims = state.get("frame_zlims", {})
        use_global_zlim = bool(
            args.frames_global_zlim and not _frames_auto_zlim_enabled(args.frames_auto_zlim)
        )
        if frame_fields_startup and use_global_zlim and not frame_zlims and to_process:
            first = min(
                to_process,
                key=lambda p: _parse_plot_index(os.path.basename(p)) if _parse_plot_index(os.path.basename(p)) is not None else 10**12,
            )
            frame_zlims = _lock_frame_zlims_from_plotfile(first, args_dict)
            state["frame_zlims"] = frame_zlims
            _save_state(state_path, state)
        args_dict["frame_zlims"] = frame_zlims
        args_dict["frames_global_zlim"] = use_global_zlim

        if args.jobs > 1:
            from concurrent.futures import ProcessPoolExecutor, as_completed

            with ProcessPoolExecutor(max_workers=args.jobs) as executor:
                futures = {
                    executor.submit(
                        _process_single_plotfile, p, args_dict, protected, processed_count + i
                    ): p
                    for i, p in enumerate(to_process)
                }
                done = 0
                for f in as_completed(futures):
                    p = futures[f]
                    done += 1
                    key = os.path.basename(p)
                    print(f"[{done}/{len(to_process)}] finished {key}", flush=True)
                    try:
                        res = f.result()
                        if res["success"]:
                            if res["psi4_line"]:
                                _append_line(out_path, header=header, line=res["psi4_line"])
                            if res["areal_line"]:
                                _append_line(areal_out_path, header=areal_header, line=res["areal_line"])
                            if res["shell_line"]:
                                _append_line(shell_out_path, header=shell_header, line=res["shell_line"])
                            if res.get("boundary_flux_line"):
                                _append_line(
                                    boundary_flux_out_path,
                                    header="time  net_outward_flux  psi4_boundary_amp",
                                    line=res["boundary_flux_line"],
                                )
                            
                            state[res["key"]] = True
                            _save_state(state_path, state)
                            processed_count += 1
                            
                            if args.verbose:
                                print(f"[ok] {res['key']}  t={res['t']:.6g}  {res['status_str']}  ({res['dt_s']:.2f}s)")
                        else:
                            print(f"WARNING: failed to process {p}: {res.get('error', 'Unknown error')}")
                    except Exception as e:
                        print(f"WARNING: worker failed for {p}: {e}")
        else:
            for i, p in enumerate(to_process):
                key = os.path.basename(p)
                print(f"[{i + 1}/{len(to_process)}] loading {key} ...", flush=True)
                res = _process_single_plotfile(p, args_dict, protected, processed_count + i)
                if res["success"]:
                    if res["psi4_line"]:
                        _append_line(out_path, header=header, line=res["psi4_line"])
                    if res["areal_line"]:
                        _append_line(areal_out_path, header=areal_header, line=res["areal_line"])
                    if res["shell_line"]:
                        _append_line(shell_out_path, header=shell_header, line=res["shell_line"])
                    
                    state[res["key"]] = True
                    _save_state(state_path, state)
                    processed_count += 1
                    
                    if args.verbose:
                        print(f"[ok] {res['key']}  t={res['t']:.6g}  {res['status_str']}  ({res['dt_s']:.2f}s)")
                else:
                    print(f"WARNING: failed to process {p}: {res.get('error', 'Unknown error')}")

        # Cleanup pass: delete plotfiles that were processed earlier but were
        # inside keep-last at the time. Once they are no longer protected, we
        # can safely delete them.
        if args.delete:
            for p in plot_dirs:
                if p in protected:
                    continue
                key = os.path.basename(p)
                if not state.get(key):
                    continue
                if not _is_plotfile_ready(p, stable_seconds=float(args.stable_seconds)):
                    continue
                try:
                    shutil.rmtree(p)
                    if args.verbose:
                        print(f"[gc] deleted previously-processed {key}")
                except Exception as e:
                    print(f"WARNING: failed to delete {p}: {e}")

        return processed_count

    if args.watch:
        print(f"Watching {data_dir} for plotfiles. Writing {out_path}")
        while True:
            n = process_once()
            if n:
                if not args.verbose:
                    print(f"Processed {n} plotfile(s).")
            time.sleep(float(args.poll_seconds))
    else:
        n = process_once()
        print(f"Processed {n} plotfile(s). Output: {out_path}")


if __name__ == "__main__":
    main()

