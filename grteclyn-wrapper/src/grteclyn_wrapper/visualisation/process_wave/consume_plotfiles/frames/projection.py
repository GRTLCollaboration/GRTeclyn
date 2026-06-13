from __future__ import annotations

import os

import matplotlib.pyplot as plt
import yt
from matplotlib.ticker import FuncFormatter

from ..config import _FRAME_DPI, _field_frame_config
from ..fields import _field_key, _register_derived_fields


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
