"""Embedding-style views of an elite's intrinsic spatial geometry.

An arbitrary non-axisymmetric 3-metric need not admit a global isometric
embedding in Euclidean 3-space.  The surface here therefore displays the
logarithmic proper-area distortion of the most deformed central coordinate
plane.  It is a data-faithful height map, not a wormhole cartoon.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import MaxNLocator
from numpy.typing import NDArray

from ...search.geometry_atlas.render import RenderConfig, RenderedGeometry


@dataclass(frozen=True)
class IntrinsicEmbedding:
    """Proper-area distortion sampled on one central coordinate plane."""

    u: NDArray[np.float64]
    v: NDArray[np.float64]
    height: NDArray[np.float64]
    plane: str
    u_label: str
    v_label: str
    normal_label: str
    normal_coordinate: float

    @property
    def max_abs_distortion(self) -> float:
        return float(np.max(np.abs(self.height)))


def _plane_embedding(
    plane_gamma: NDArray[np.float64],
    *,
    metric_axes: tuple[int, int],
    coordinates: NDArray[np.float64],
    plane: str,
    labels: tuple[str, str, str],
    normal_coordinate: float,
) -> IntrinsicEmbedding:
    """Build ``log(dA_proper / dA_flat)`` from an induced 2-metric."""
    a, b = metric_axes
    determinant = (
        plane_gamma[..., a, a] * plane_gamma[..., b, b]
        - plane_gamma[..., a, b] * plane_gamma[..., b, a]
    )
    area_ratio = np.sqrt(np.clip(determinant, np.finfo(float).tiny, None))
    height = np.log(area_ratio)
    u, v = np.meshgrid(coordinates, coordinates, indexing="xy")
    return IntrinsicEmbedding(
        u=u,
        v=v,
        height=height,
        plane=plane,
        u_label=labels[0],
        v_label=labels[1],
        normal_label=labels[2],
        normal_coordinate=normal_coordinate,
    )


def intrinsic_embedding(
    rendered: RenderedGeometry,
    cfg: RenderConfig,
) -> IntrinsicEmbedding:
    """Return the central plane with the largest RMS proper-area distortion."""
    k = cfg.n // 2
    coordinates = (
        np.linspace(-0.5 * cfg.L, 0.5 * cfg.L, cfg.n, endpoint=False)
        + 0.5 * cfg.dx
    )
    normal_coordinate = float(coordinates[k])
    gamma = rendered.gamma
    candidates = (
        _plane_embedding(
            gamma[k, :, :, :, :],
            metric_axes=(0, 1),
            coordinates=coordinates,
            plane="xy",
            labels=("x", "y", "z"),
            normal_coordinate=normal_coordinate,
        ),
        _plane_embedding(
            gamma[:, k, :, :, :],
            metric_axes=(0, 2),
            coordinates=coordinates,
            plane="xz",
            labels=("x", "z", "y"),
            normal_coordinate=normal_coordinate,
        ),
        _plane_embedding(
            gamma[:, :, k, :, :],
            metric_axes=(1, 2),
            coordinates=coordinates,
            plane="yz",
            labels=("y", "z", "x"),
            normal_coordinate=normal_coordinate,
        ),
    )
    return max(
        candidates,
        key=lambda surface: float(np.sqrt(np.mean(surface.height * surface.height))),
    )


def plot_intrinsic_embedding(
    ax,
    embedding: IntrinsicEmbedding,
    *,
    title: str | None = None,
    show_note: bool = True,
) -> ScalarMappable:
    """Plot an intrinsic proper-area distortion as a 3D surface."""
    span = max(embedding.max_abs_distortion, 1.0e-12)
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-span, vmax=span)
    surface = ax.plot_surface(
        embedding.u,
        embedding.v,
        embedding.height,
        cmap="RdBu_r",
        norm=norm,
        linewidth=0.0,
        antialiased=True,
        rasterized=True,
    )
    ax.contour(
        embedding.u,
        embedding.v,
        embedding.height,
        levels=8,
        cmap="Greys",
        linewidths=0.45,
        alpha=0.55,
    )
    ax.set_xlabel(rf"${embedding.u_label}$", labelpad=4)
    ax.set_ylabel(rf"${embedding.v_label}$", labelpad=4)
    ax.set_zlabel(
        r"$\ln(dA_{\rm proper}/dA_{\rm flat})$",
        labelpad=5,
    )
    ax.zaxis.set_major_locator(MaxNLocator(4))
    ax.view_init(elev=28.0, azim=-55.0)
    ax.set_box_aspect((1.0, 1.0, 0.48))
    ax.set_title(
        (
            rf"Intrinsic area distortion: ${embedding.plane}$ plane "
            rf"at ${embedding.normal_label}={embedding.normal_coordinate:+.2f}$"
        )
        if title is None
        else title,
        pad=10,
    )
    if show_note:
        ax.text2D(
            0.02,
            0.02,
            "height map; not a global isometric embedding",
            transform=ax.transAxes,
            fontsize=7,
            color="0.25",
        )
    return ScalarMappable(norm=norm, cmap=surface.cmap)


__all__ = [
    "IntrinsicEmbedding",
    "intrinsic_embedding",
    "plot_intrinsic_embedding",
]
