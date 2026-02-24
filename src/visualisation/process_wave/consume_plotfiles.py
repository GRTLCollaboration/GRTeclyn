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


def _default_data_dir() -> str:
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(script_dir)))
    return os.path.abspath(os.path.join(project_root, "..", "data"))

def _default_frames_out_dir() -> str:
    # .../src/visualisation/process_wave -> .../src/visualisation/visualize
    script_dir = os.path.dirname(os.path.abspath(__file__))
    visualisation_dir = os.path.dirname(script_dir)
    return os.path.join(visualisation_dir, "visualize")


def _iter_plotfile_dirs(data_dir: str) -> List[str]:
    """Return sorted plotfile directories under data_dir."""
    out: List[str] = []
    if not os.path.isdir(data_dir):
        return out
    for name in os.listdir(data_dir):
        if not (name.startswith("WormholePlt") or name.startswith("plt")):
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
    }
    return aliases.get(name, name)


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
    # GW proxy fields
    if field == "GW_Plus":
        def _gw_plus(_field, data):
            return data["boxlib", "A11"] - data["boxlib", "A22"]
        ds.add_field(("boxlib", "GW_Plus"), function=_gw_plus, sampling_type="cell", units="")
    elif field == "GW_Cross":
        def _gw_cross(_field, data):
            return 2.0 * data["boxlib", "A12"]
        ds.add_field(("boxlib", "GW_Cross"), function=_gw_cross, sampling_type="cell", units="")
    elif field == "Weyl4_Mag":
        def _weyl4_mag(_field, data):
            re_v = data["boxlib", "Weyl4_Re"]
            im_v = data["boxlib", "Weyl4_Im"]
            return np.sqrt(re_v**2 + im_v**2)
        ds.add_field(("boxlib", "Weyl4_Mag"), function=_weyl4_mag, sampling_type="cell", units="")

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
) -> str:
    """
    Render a SlicePlot frame and save it under:
      <frames_out_dir>/<field>_<axis>/frames/frame_<axis>_<idx>.png
    Returns the saved path.
    """
    configs = {
        "chi":   {"zlim": (0.0, 1.0),  "cmap": "magma", "label": r"Conformal Factor $\chi$"},
        "K":     {"zlim": (-0.1, 0.1), "cmap": "RdBu",  "label": r"Trace of Extrinsic Curvature $K$"},
        "Theta": {"zlim": (-0.005, 0.005), "cmap": "RdBu", "label": r"Z4 Constraint $\Theta$"},
        "lapse": {"zlim": (0.0, 1.0),  "cmap": "viridis", "label": r"Lapse $\alpha$"},
        "Ham":   {"zlim": (-0.1, 0.1), "cmap": "RdBu", "label": r"Hamiltonian Constraint $\mathcal{H}$"},
        "A11":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{11}$"},
        "A12":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{12}$"},
        "A22":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{22}$"},
        "A33":   {"zlim": (-0.05, 0.05), "cmap": "PuOr", "label": r"Extrinsic Curvature $\tilde{A}_{33}$"},
        "GW_Plus": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_+$"},
        "GW_Cross": {"zlim": (-5.0e-6, 5.0e-6), "cmap": "PiYG", "label": r"GW Strain $h_\times$"},
        "Weyl4_Re": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Re}(\Psi_4)$"},
        "Weyl4_Im": {"zlim": (-1.0e-4, 1.0e-4), "cmap": "PiYG", "label": r"$\mathrm{Im}(\Psi_4)$"},
        "Weyl4_Mag": {"zlim": (0.0, 1.0e-4),    "cmap": "inferno", "label": r"$|\Psi_4|$"},
    }
    cfg = configs.get(field, {"zlim": (None, None), "cmap": "viridis", "label": field})

    # Center logic matches visualize/__main__.py
    mid_x = float((ds.domain_right_edge[0] + ds.domain_left_edge[0]) / 2.0)
    mid_y = float((ds.domain_right_edge[1] + ds.domain_left_edge[1]) / 2.0)
    physics_center = [mid_x, mid_y, 0.0]

    # Corner mode: treat `center_xyz` as the symmetry origin (corner) and set
    # plot center to origin + zoom/2 in the slice plane.
    if corner and zoom is not None:
        slice_plane_val = 0.0 if coord is None else float(coord)
        origin = np.array(center_xyz, dtype=float) if center_xyz is not None else np.array(ds.domain_left_edge, dtype=float)
        w = float(zoom)
        if axis == "z":
            physics_center = [origin[0] + w / 2.0, origin[1] + w / 2.0, slice_plane_val]
        elif axis == "y":
            physics_center = [origin[0] + w / 2.0, slice_plane_val, origin[2] + w / 2.0]
        elif axis == "x":
            physics_center = [slice_plane_val, origin[1] + w / 2.0, origin[2] + w / 2.0]
    elif center_xyz is not None:
        physics_center = [float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2])]

    # Apply coord override
    if coord is not None:
        if axis == "z":
            physics_center[2] = float(coord)
        elif axis == "y":
            physics_center[1] = float(coord)
        elif axis == "x":
            physics_center[0] = float(coord)

    plot_center = ds.arr(physics_center, "code_length")

    _register_derived_fields(ds, field)

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

    slc = yt.SlicePlot(ds, axis, ("boxlib", field), center=plot_center)
    # Use dataset-native coordinates (e.g. [0,40]) on axes for symmetry-reduced domains.
    slc.set_origin("native")
    slc.set_axes_unit("code_length")
    if zoom is not None:
        slc.set_width((float(zoom), "code_length"))
    slc.set_log(("boxlib", field), False)

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

    if cfg["zlim"][0] is not None:
        slc.set_zlim(("boxlib", field), cfg["zlim"][0], cfg["zlim"][1])
    slc.set_cmap(("boxlib", field), cfg["cmap"])

    coord_val = physics_center[{"x": 0, "y": 1, "z": 2}[axis]]

    # Format title with LaTeX
    slc.annotate_title(r"%s $\quad t=%.2f \quad %s=%.1f$" % (cfg['label'], float(ds.current_time), axis, coord_val))

    # Force LaTeX labels for axes
    axis_map = {"x": ("y", "z"), "y": ("z", "x"), "z": ("x", "y")}
    xlabel_name, ylabel_name = axis_map[axis]
    slc.set_xlabel(r"$%s$" % xlabel_name)
    slc.set_ylabel(r"$%s$" % ylabel_name)
    slc.set_colorbar_label(("boxlib", field), cfg['label'])

    # Remove the duplicated "0" tick label at the symmetry origin corner.
    if corner:
        slc.render()
        plot = slc.plots[("boxlib", field)]
        ax = plot.axes
        ymin = float(ax.get_ylim()[0])
        for tick, lbl in zip(ax.get_yticks(), ax.get_yticklabels()):
            if abs(float(tick) - ymin) < 1.0e-9:
                lbl.set_visible(False)

    output_dir = os.path.join(frames_out_dir, f"{field}_{axis}")
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    frame_name = f"frame_{axis}_{frame_idx:04d}.png"
    out_path = os.path.join(frames_dir, frame_name)
    slc.save(out_path)

    if verbose:
        print(f"[frame] {field} -> {out_path}")

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

    try:
        vals = ds.find_field_values_at_points(
            [("boxlib", "Weyl4_Re"), ("boxlib", "Weyl4_Im")], pts
        )
        re_vals = np.asarray(vals[0], dtype=float)
        im_vals = np.asarray(vals[1], dtype=float)
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


def _load_state(state_path: Path) -> Dict[str, bool]:
    if not state_path.exists():
        return {}
    try:
        return json.loads(state_path.read_text())
    except Exception:
        return {}


def _save_state(state_path: Path, state: Dict[str, bool]) -> None:
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
    parser.add_argument("--frames-out", default=_default_frames_out_dir(), help="Frames output base dir (default: src/visualisation/visualize).")
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
    header = "# time  " + "  ".join([f"Re(R={R:g})  Im(R={R:g})" for R in args.radii])

    state = _load_state(state_path)

    # Auto-reset on simulation restart (same output dir reused):
    # If we detect plotfiles starting again from 0 but our saved state refers to
    # plotfiles that don't exist, truncate outputs and start fresh.
    plot_dirs_now = _iter_plotfile_dirs(data_dir)
    if _should_auto_reset(plot_dirs_now, state):
        print("Detected a likely simulation restart in the same output directory.")
        print(f"Resetting: {out_path} and {state_path}")
        _truncate_if_exists(out_path)
        _save_state(state_path, {})
        state = {}

    # If rendering frames, clear existing frames for the requested fields/axis at startup.
    frame_fields_startup = [_canonical_field_name(f) for f in args.frames_fields]
    if frame_fields_startup:
        _cleanup_existing_frames(
            frames_out_dir=os.path.abspath(args.frames_out),
            fields=frame_fields_startup,
            axis=args.frames_axis,
            verbose=args.verbose,
        )

    def process_once() -> int:
        plot_dirs = _iter_plotfile_dirs(data_dir)
        if not plot_dirs:
            return 0
        processed = 0

        # Never delete newest keep_last plotfiles
        protected = set(plot_dirs[-max(0, int(args.keep_last)) :])

        for p in plot_dirs:
            key = os.path.basename(p)
            if state.get(key):
                continue
            if not _is_plotfile_ready(p, stable_seconds=float(args.stable_seconds)):
                continue

            try:
                t0 = time.time()
                ds = yt.load(p)
                t = float(ds.current_time)

                # --- Psi4 extraction (optional) ---
                if args.psi4:
                    # Validate fields once per plotfile
                    if ("boxlib", "Weyl4_Re") not in ds.field_list or ("boxlib", "Weyl4_Im") not in ds.field_list:
                        raise RuntimeError(
                            "Plotfile missing Weyl4_Re/Im. Set: amr.derive_plot_vars = Weyl4 and re-run."
                        )
                    amps = _extract_mode_amps_l2m0(
                        ds,
                        radii=[float(r) for r in args.radii],
                        n_points=int(args.n_points),
                        center=args.center,
                    )
                    line = f"{t:.16e}  " + "  ".join([f"{a.real:.16e}  {a.imag:.16e}" for a in amps])
                    _append_line(out_path, header=header, line=line)

                # --- Frame rendering (optional) ---
                frame_fields = [_canonical_field_name(f) for f in args.frames_fields]
                if frame_fields:
                    idx = _parse_plot_index(key)
                    frame_idx = idx if idx is not None else processed
                    for fld in frame_fields:
                        _render_slice_frame(
                            ds,
                            field=fld,
                            axis=args.frames_axis,
                            coord=args.frames_coord,
                            zoom=args.frames_zoom,
                            center_xyz=args.frames_center,
                            corner=args.frames_corner,
                            frames_out_dir=os.path.abspath(args.frames_out),
                            frame_idx=int(frame_idx),
                            verbose=args.verbose,
                        )

                state[key] = True
                _save_state(state_path, state)
                processed += 1

                deleted = False
                if args.delete and (p not in protected):
                    shutil.rmtree(p)
                    deleted = True

                if args.verbose:
                    dt_s = time.time() - t0
                    status = "deleted" if deleted else ("kept" if args.delete else "kept(no-delete)")
                    print(f"[ok] {key}  t={t:.6g}  {status}  ({dt_s:.2f}s)")
            except Exception as e:
                # Do not mark as processed; do not delete.
                print(f"WARNING: failed to process {p}: {e}")
                continue

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

        return processed

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

