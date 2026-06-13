from __future__ import annotations

import os
from typing import Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

from ..config import _FRAME_DPI


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
