"""Render a stationary geometry genome onto a GRTeclyn gridinit."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from ...grtresna.io.gridinit import write_gridinit
from ...initial_data.adm_stationary import (
    StationaryGrid,
    constraint_residuals_from_adm,
    pack_ccz4_grid,
)
from .genome import GeometryGenome, decode_fields_at_points, minkowski_deviation


@dataclass(frozen=True)
class RenderConfig:
    """Cartesian grid used to materialize a genome."""

    n: int = 32
    L: float = 64.0  # full box length; domain is [0, L]^3 with center at L/2

    @property
    def dx(self) -> float:
        return self.L / self.n

    @property
    def origin(self) -> tuple[float, float, float]:
        return (0.0, 0.0, 0.0)

    @property
    def center(self) -> tuple[float, float, float]:
        half = 0.5 * self.L
        return (half, half, half)


@dataclass
class RenderedGeometry:
    """Rendered ADM / CCZ4 state plus diagnostics."""

    grid: StationaryGrid
    alpha: NDArray[np.float64]
    beta: NDArray[np.float64]
    gamma: NDArray[np.float64]
    kij_extra: NDArray[np.float64]
    diagnostics: dict[str, float]


def _shift_fraction(
    alpha: NDArray[np.float64],
    beta: NDArray[np.float64],
    gamma: NDArray[np.float64],
) -> float:
    """Morphology descriptor: fraction of the Minkowski deviation carried by shift.

    Separates geometry *families* for the atlas: ``1`` ≈ pure shift channel
    (Alcubierre-like warp tubes), ``0`` ≈ pure lapse/curvature deformation
    (conformal lenses). Complements the exotic-energy cost axis.
    """
    dev_shift = float(np.sum(beta * beta))
    dev_lapse = float(np.sum((alpha - 1.0) ** 2))
    eye = np.eye(3).reshape((1,) * (gamma.ndim - 2) + (3, 3))
    dev_gamma = float(np.sum((gamma - eye) ** 2))
    total = dev_shift + dev_lapse + dev_gamma
    if total <= 1.0e-30:
        return 0.0
    return dev_shift / total


def _mesh(cfg: RenderConfig) -> NDArray[np.float64]:
    """Return cell-center coordinates relative to the box center, shape (nz,ny,nx,3)."""
    dx = cfg.dx
    cx, cy, cz = cfg.center
    x = (np.arange(cfg.n, dtype=np.float64) + 0.5) * dx - cx
    y = (np.arange(cfg.n, dtype=np.float64) + 0.5) * dx - cy
    z = (np.arange(cfg.n, dtype=np.float64) + 0.5) * dx - cz
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    return np.stack([xx, yy, zz], axis=-1)


def render_genome(
    genome: GeometryGenome,
    cfg: RenderConfig | None = None,
    *,
    include_einstein_source: bool = True,
) -> RenderedGeometry:
    """Decode genome → ADM fields → CCZ4 gridinit payload."""
    cfg = cfg or RenderConfig()
    points = _mesh(cfg)
    alpha, beta, gamma, kij_extra = decode_fields_at_points(genome, points)

    # Signature / SPD sanity.
    eig = np.linalg.eigvalsh(gamma)
    if float(np.min(eig)) <= 0.0:
        raise ValueError("Rendered spatial metric is not positive-definite")
    if float(np.min(alpha)) <= 0.0:
        raise ValueError("Rendered lapse is non-positive")

    dx = (cfg.dx, cfg.dx, cfg.dx)
    packed = pack_ccz4_grid(
        alpha,
        beta,
        gamma,
        dx,
        cfg.origin,
        include_einstein_source=include_einstein_source,
        kij_extra=kij_extra,
    )

    residuals = constraint_residuals_from_adm(
        alpha,
        beta,
        gamma,
        dx,
        center=cfg.center,
        trim_cells=max(2, cfg.n // 16),
        puncture_cells=0.0,
        origin=cfg.origin,
    )
    cell_vol = cfg.dx**3
    neg_mask = packed.rho < 0.0
    integral_neg = float(-np.sum(packed.rho[neg_mask]) * cell_vol)
    integral_abs = float(np.sum(np.abs(packed.rho)) * cell_vol)

    diagnostics = {
        **minkowski_deviation(alpha, beta, gamma),
        "shift_fraction": _shift_fraction(alpha, beta, gamma),
        **residuals.as_dict(),
        "min_rho": float(np.min(packed.rho)),
        "max_rho": float(np.max(packed.rho)),
        "integral_negative_rho": integral_neg,
        "integral_abs_rho": integral_abs,
        "min_gamma_eig": float(np.min(eig)),
        "max_abs_j": float(np.max(np.abs(packed.j))),
        "boundary_max_abs_alpha_m1": float(
            max(
                np.max(np.abs(alpha[0] - 1.0)),
                np.max(np.abs(alpha[-1] - 1.0)),
                np.max(np.abs(alpha[:, 0] - 1.0)),
                np.max(np.abs(alpha[:, -1] - 1.0)),
                np.max(np.abs(alpha[:, :, 0] - 1.0)),
                np.max(np.abs(alpha[:, :, -1] - 1.0)),
            )
        ),
        "boundary_max_abs_beta": float(
            max(
                np.max(np.abs(beta[0])),
                np.max(np.abs(beta[-1])),
                np.max(np.abs(beta[:, 0])),
                np.max(np.abs(beta[:, -1])),
                np.max(np.abs(beta[:, :, 0])),
                np.max(np.abs(beta[:, :, -1])),
            )
        ),
        "max_abs_kij_extra": float(np.max(np.abs(kij_extra))),
    }
    return RenderedGeometry(
        grid=packed,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        kij_extra=kij_extra,
        diagnostics=diagnostics,
    )


def write_rendered_gridinit(
    rendered: RenderedGeometry,
    path: str | Path,
) -> Path:
    """Write a rendered geometry to a ``.gridinit`` file."""
    return write_gridinit(
        rendered.grid.data,
        rendered.grid.comp_names,
        rendered.grid.dx_xyz,
        rendered.grid.origin,
        path,
    )


def render_and_write(
    genome: GeometryGenome,
    path: str | Path,
    cfg: RenderConfig | None = None,
) -> tuple[RenderedGeometry, Path]:
    rendered = render_genome(genome, cfg)
    out = write_rendered_gridinit(rendered, path)
    return rendered, out


__all__ = [
    "RenderConfig",
    "RenderedGeometry",
    "render_and_write",
    "render_genome",
    "write_rendered_gridinit",
]
