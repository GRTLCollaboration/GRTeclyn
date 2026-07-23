"""Multi-panel midplane figures for a rendered geometry-atlas elite."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from numpy.typing import NDArray

from ...search.geometry_atlas.render import RenderConfig, RenderedGeometry, render_genome
from .embedding import intrinsic_embedding, plot_intrinsic_embedding
from .elites import EliteRecord


_PUBLICATION_STYLE = {
    # Matplotlib's bundled STIX renderer gives LaTeX-style text without requiring
    # a LaTeX installation.
    "text.usetex": False,
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "DejaVu Serif", "serif"],
    "mathtext.fontset": "stix",
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "axes.linewidth": 0.8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
}
_DPI = 300
_IMAGE_OPTIONS = {"origin": "lower", "interpolation": "none", "rasterized": True}


def render_elite(
    elite: EliteRecord,
    *,
    n: int | None = None,
    L: float | None = None,
) -> tuple[RenderedGeometry, RenderConfig]:
    """Re-render an elite genome onto a Cartesian grid for plotting."""
    cfg_n = 64 if n is None else int(n)
    if L is None:
        # Prefer the campaign box encoded on the genome; fall back to 32.
        half = float(elite.genome.config.box_half_width)
        cfg_L = 2.0 * half if half > 0.0 else 32.0
    else:
        cfg_L = float(L)
    render_cfg = RenderConfig(n=cfg_n, L=cfg_L)
    return render_genome(elite.genome, render_cfg), render_cfg


def _midplane_index(n: int) -> int:
    return n // 2


def _extent(cfg: RenderConfig) -> tuple[float, float, float, float]:
    half = 0.5 * cfg.L
    return (-half, half, -half, half)


def _gamma_deviation(gamma: NDArray[np.float64]) -> NDArray[np.float64]:
    """Frobenius ||γ − I|| at each grid point."""
    eye = np.eye(3).reshape((1,) * (gamma.ndim - 2) + (3, 3))
    return np.linalg.norm(gamma - eye, axis=(-2, -1))


def _slice_xy(field: NDArray, k: int) -> NDArray:
    """Return the z=midplane slice. Field is (nz, ny, nx[, ...])."""
    return np.asarray(field[k])


def _slice_coordinate(k: int, cfg: RenderConfig) -> float:
    """Physical coordinate of a cell-centred slice."""
    return -0.5 * cfg.L + (k + 0.5) * cfg.dx


def _full_sequential_norm(field: NDArray) -> Normalize:
    """Colour normalization spanning every finite data value."""
    finite = np.asarray(field)[np.isfinite(field)]
    vmin = float(finite.min())
    vmax = float(finite.max())
    if np.isclose(vmin, vmax):
        pad = max(abs(vmin) * 1.0e-12, 1.0e-12)
        vmin -= pad
        vmax += pad
    return Normalize(vmin=vmin, vmax=vmax)


def _full_centered_norm(field: NDArray, *, center: float = 0.0) -> TwoSlopeNorm:
    """Symmetric normalization containing the full range around ``center``."""
    finite = np.asarray(field)[np.isfinite(field)]
    span = float(np.max(np.abs(finite - center)))
    span = max(span, 1.0e-12)
    return TwoSlopeNorm(vcenter=center, vmin=center - span, vmax=center + span)


def _style_axes(axes: NDArray, cfg: RenderConfig) -> None:
    half = 0.5 * cfg.L
    for ax in np.asarray(axes).ravel():
        ax.set_xlim(-half, half)
        ax.set_ylim(-half, half)
        ax.set_aspect("equal")
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")


def _add_colorbar(fig: plt.Figure, image, ax: plt.Axes):
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.035)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 3))
    colorbar.formatter = formatter
    colorbar.update_ticks()
    colorbar.ax.tick_params(labelsize=8)
    return colorbar


def plot_elite_panels(
    elite: EliteRecord,
    rendered: RenderedGeometry,
    cfg: RenderConfig,
    *,
    out_path: Path,
    quiver_stride: int = 4,
) -> Path:
    """Write a publication-quality 2×3 central-slice figure for one elite.

    Every colour limit is derived from the complete field; no percentile
    clipping or image interpolation is applied.
    """
    k = _midplane_index(cfg.n)
    z = _slice_coordinate(k, cfg)
    extent = _extent(cfg)
    alpha = _slice_xy(rendered.alpha, k)
    beta = _slice_xy(rendered.beta, k)  # (ny, nx, 3)
    beta_mag = np.linalg.norm(beta, axis=-1)
    gamma_dev = _slice_xy(_gamma_deviation(rendered.gamma), k)
    rho = _slice_xy(rendered.grid.rho, k)

    modes = ", ".join(elite.active_modes())
    title = (
        rf"Elite #{elite.rank}: $f_{{\rm geo}}={elite.f_geo:.3f}$, "
        rf"cell $({elite.cell[0]},{elite.cell[1]})$, {elite.family_label}, "
        rf"$f_{{\rm shift}}={elite.shift_fraction:.2f}$"
        "\n"
        rf"{modes}; $\int\rho_-\,d^3x={elite.integral_negative_rho:.2f}$; "
        rf"central slice $z={z:+.2f}$"
    )

    with plt.rc_context(_PUBLICATION_STYLE):
        fig, axes = plt.subplots(
            2, 3, figsize=(13.8, 8.8), constrained_layout=True
        )
        fig.suptitle(title, fontsize=13)

        # Field layout is (nz, ny, nx); imshow(rows=y, cols=x) maps it directly.
        ax = axes[0, 0]
        im = ax.imshow(
            alpha,
            extent=extent,
            cmap="RdBu_r",
            norm=_full_centered_norm(alpha, center=1.0),
            **_IMAGE_OPTIONS,
        )
        ax.set_title(r"Lapse $\alpha$")
        _add_colorbar(fig, im, ax)

        ax = axes[0, 1]
        beta_norm = _full_sequential_norm(beta_mag)
        im = ax.imshow(
            beta_mag,
            extent=extent,
            cmap="inferno",
            norm=beta_norm,
            **_IMAGE_OPTIONS,
        )
        ax.set_title(r"Shift magnitude $\|\boldsymbol{\beta}\|_2$")
        _add_colorbar(fig, im, ax)

        ax = axes[0, 2]
        ax.imshow(
            beta_mag,
            extent=extent,
            cmap="inferno",
            norm=beta_norm,
            alpha=0.92,
            **_IMAGE_OPTIONS,
        )
        s = max(1, int(quiver_stride))
        half = 0.5 * cfg.L
        xs = np.linspace(-half, half, cfg.n, endpoint=False) + 0.5 * cfg.dx
        ys = np.linspace(-half, half, cfg.n, endpoint=False) + 0.5 * cfg.dx
        X, Y = np.meshgrid(xs[::s], ys[::s], indexing="xy")
        U = beta[::s, ::s, 0]
        V = beta[::s, ::s, 1]
        max_vector = float(np.max(np.hypot(U, V)))
        scale = max_vector / (0.75 * s * cfg.dx) if max_vector > 0.0 else 1.0
        quiver = ax.quiver(
            X,
            Y,
            U,
            V,
            color="white",
            edgecolor="black",
            linewidth=0.25,
            angles="xy",
            scale_units="xy",
            scale=scale,
            width=0.004,
            headwidth=3.5,
        )
        if max_vector > 0.0:
            ax.quiverkey(
                quiver,
                0.93,
                1.035,
                max_vector,
                rf"$\|\boldsymbol{{\beta}}_{{xy}}\|_{{\max}}={max_vector:.2g}$",
                coordinates="axes",
                labelpos="W",
                fontproperties={"size": 8},
            )
        ax.set_title(r"Shift field $(\beta^x,\beta^y)$")

        ax = axes[1, 0]
        im = ax.imshow(
            gamma_dev,
            extent=extent,
            cmap="cividis",
            norm=_full_sequential_norm(gamma_dev),
            **_IMAGE_OPTIONS,
        )
        ax.set_title(r"Spatial deformation $\|\gamma-I\|_{\rm F}$")
        _add_colorbar(fig, im, ax)

        ax = axes[1, 1]
        im = ax.imshow(
            rho,
            extent=extent,
            cmap="RdBu_r",
            norm=_full_centered_norm(rho),
            **_IMAGE_OPTIONS,
        )
        ax.set_title(r"Required energy density $\rho$ (exotic: $\rho<0$)")
        _add_colorbar(fig, im, ax)

        _style_axes(axes.ravel()[:5], cfg)
        axes[1, 2].remove()
        ax = fig.add_subplot(2, 3, 6, projection="3d")
        axes[1, 2] = ax
        embedding = intrinsic_embedding(rendered, cfg)
        mappable = plot_intrinsic_embedding(
            ax,
            embedding,
            title=(
                rf"Intrinsic spatial geometry "
                rf"(${embedding.plane}$ plane; embedding-style)"
            ),
        )
        _add_colorbar(fig, mappable, ax)

        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=_DPI, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
    return out_path


def plot_gallery(
    elites: list[EliteRecord],
    *,
    out_path: Path,
    n: int | None = None,
    L: float | None = None,
) -> Path:
    """Comparison gallery with slices and an intrinsic-geometry surface."""
    if not elites:
        raise ValueError("no elites to plot")

    rendered_rows = []
    for elite in elites:
        rendered, cfg = render_elite(elite, n=n, L=L)
        k = _midplane_index(cfg.n)
        rendered_rows.append(
            (
                elite,
                cfg,
                np.linalg.norm(_slice_xy(rendered.beta, k), axis=-1),
                _slice_xy(_gamma_deviation(rendered.gamma), k),
                _slice_xy(rendered.grid.rho, k),
                intrinsic_embedding(rendered, cfg),
            )
        )

    cmaps = ("inferno", "cividis", "RdBu_r")
    column_titles = (
        r"Shift magnitude $\|\boldsymbol{\beta}\|_2$",
        r"Spatial deformation $\|\gamma-I\|_{\rm F}$",
        r"Required energy density $\rho$",
    )

    with plt.rc_context(_PUBLICATION_STYLE):
        fig = plt.figure(
            figsize=(16.5, 2.8 * len(elites)),
            constrained_layout=True,
        )
        grid = fig.add_gridspec(
            len(elites),
            4,
            width_ratios=(1.0, 1.0, 1.0, 1.18),
        )
        axes = np.empty((len(elites), 4), dtype=object)
        for row_index in range(len(elites)):
            for column in range(3):
                axes[row_index, column] = fig.add_subplot(grid[row_index, column])
            axes[row_index, 3] = fig.add_subplot(
                grid[row_index, 3],
                projection="3d",
            )

        z = _slice_coordinate(_midplane_index(rendered_rows[0][1].n), rendered_rows[0][1])
        fig.suptitle(
            rf"Geometry-atlas top elites: central slices at $z={z:+.2f}$ "
            r"and intrinsic spatial geometry",
            fontsize=14,
        )

        for row_index, (
            elite,
            cfg,
            beta_mag,
            gamma_dev,
            rho,
            embedding,
        ) in enumerate(rendered_rows):
            fields = (beta_mag, gamma_dev, rho)
            norms = (
                _full_sequential_norm(beta_mag),
                _full_sequential_norm(gamma_dev),
                _full_centered_norm(rho),
            )
            for column, (field, cmap, norm) in enumerate(zip(fields, cmaps, norms)):
                ax = axes[row_index, column]
                image = ax.imshow(
                    field,
                    extent=_extent(cfg),
                    cmap=cmap,
                    norm=norm,
                    **_IMAGE_OPTIONS,
                )
                if row_index == 0:
                    ax.set_title(column_titles[column], pad=8)
                ax.set_aspect("equal")
                ax.xaxis.set_major_locator(MaxNLocator(5))
                ax.yaxis.set_major_locator(MaxNLocator(5))
                if column == 0:
                    ax.set_ylabel(
                        rf"#{elite.rank}  $f_{{\rm geo}}={elite.f_geo:.3f}$"
                        "\n"
                        rf"{elite.family_label}  $(f_{{\rm shift}}={elite.shift_fraction:.2f})$"
                        "\n"
                        r"$y$",
                        labelpad=10,
                    )
                else:
                    ax.tick_params(labelleft=False)
                if row_index == len(elites) - 1:
                    ax.set_xlabel(r"$x$")
                _add_colorbar(fig, image, ax)

            embedding_ax = axes[row_index, 3]
            plot_intrinsic_embedding(
                embedding_ax,
                embedding,
                title="Intrinsic area distortion" if row_index == 0 else "",
                show_note=False,
            )
            embedding_ax.set_xlabel("")
            embedding_ax.set_ylabel("")
            embedding_ax.set_zlabel("")
            embedding_ax.set_xticklabels([])
            embedding_ax.set_yticklabels([])
            embedding_ax.set_zticklabels([])
            embedding_ax.tick_params(length=0, pad=0)
            embedding_ax.text2D(
                0.5,
                -0.02,
                (
                    rf"${embedding.plane}$ plane; "
                    rf"$\max|\ln(dA/dA_0)|={embedding.max_abs_distortion:.2g}$"
                ),
                transform=embedding_ax.transAxes,
                fontsize=7,
                color="0.2",
                ha="center",
            )

        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=_DPI, bbox_inches="tight", pad_inches=0.08)
        plt.close(fig)
    return out_path


__all__ = ["plot_elite_panels", "plot_gallery", "render_elite"]
