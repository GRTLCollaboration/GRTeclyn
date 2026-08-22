from __future__ import annotations

import os
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np
import yt
from matplotlib.ticker import FuncFormatter

from ..config import _FRAME_DPI, _field_frame_config
from ..fields import _field_key, _register_derived_fields
from .center import _frame_buff_size, _resolve_frame_physics_center
from .slice_cache import cache_slice
from .zlim import _resolve_plot_zlim


def draw_slice_png(
    plot_arr,
    plot_extent: Sequence[float],
    *,
    field: str,
    cfg: dict,
    axis: str,
    coord_val: float,
    time: float,
    zlim: Sequence[float],
    frames_out_dir: str,
    frame_idx: int,
    corner: bool = False,
    verbose: bool = False,
    note: str = "",
) -> str:
    """Draw one slice PNG from an array that is already windowed.

    Split out of the live render path so the same picture can be drawn later
    from a cached slice, against a colour scale measured over the whole run.
    That is the only way to get a colourbar that does not move when the
    plotfiles are deleted as they are consumed.  Both callers must draw
    identically, so there is exactly one copy of this code.
    """
    plot_extent = list(plot_extent)
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
    title_text = r"%s $\quad t=%.2f \quad %s=%g$" % (
        cfg["label"],
        float(time),
        axis,
        coord_val,
    )

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
        print(f"[frame-native] {field}{note} -> {out_path}")

    return out_path


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
    cache_slices: bool = False,
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

    coord_val = _clean_zero(physics_center[ai])

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

    if cache_slices:
        # Free: the slice is already in memory and is ~1e4 times smaller than
        # the plotfile it came from, so caching it costs nothing worth saving.
        cache_slice(
            frames_out_dir, field, axis, frame_idx, plot_arr, plot_extent,
            time=float(ds.current_time), coord_val=coord_val,
        )

    return draw_slice_png(
        plot_arr,
        plot_extent,
        field=field,
        cfg=cfg,
        axis=axis,
        coord_val=coord_val,
        time=float(ds.current_time),
        zlim=zlim,
        frames_out_dir=frames_out_dir,
        frame_idx=frame_idx,
        corner=corner,
        verbose=verbose,
        note=f" {full_dims}",
    )


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
    cache_slices: bool = False,
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
            cache_slices=cache_slices,
        )

    # AMR (level > 0) draws through yt's fixed-resolution buffer, which this
    # cache does not cover; those movies keep the old per-frame behaviour.
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
