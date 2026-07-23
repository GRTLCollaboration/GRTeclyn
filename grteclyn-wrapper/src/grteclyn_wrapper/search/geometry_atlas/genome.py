"""Compact-support RBF genome for broad stationary 4-metrics.

The genome parameterizes smooth, asymptotically flat deformations of Minkowski
with no spherical or axial symmetry assumption.  Each control center carries:

* a lapse amplitude,
* a 3-vector shift amplitude,
* six independent components of a symmetric log-metric matrix,
* six independent components of a symmetric extrinsic-curvature add-on.

A trailing *analytic block* appends coherent, physically-motivated topologies
(warp-bubble shift tube, compression tunnel, conformal lens, wormhole throat)
that the isotropic RBF bumps cannot form on their own — see :mod:`ansatz`.
Independent ``K_ij`` modes do not change the frozen 4-metric ``f_geo`` (which
depends only on ``alpha, beta, gamma``); they diversify the Stage-2 handoff.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from .ansatz import (
    ANALYTIC_PARAMS,
    analytic_bounds,
    analytic_defaults,
    ansatz_slice,
    decode_analytic,
    sample_analytic,
)


# Per-center layout:
# alpha, beta_x, beta_y, beta_z, S00, S01, S02, S11, S12, S22, K00..K22
PARAMS_PER_CENTER = 16
# Legacy trailing block size (Alcubierre-only genomes); kept for migration.
ALC_PARAMS = 4
# Slice of the Alcubierre sub-block inside the analytic tail.
_ALC_SLICE = ansatz_slice("alcubierre")
_SYMMETRIC_IDX = ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))


@dataclass(frozen=True)
class GeometryGenomeConfig:
    """Bounds and layout for the stationary geometry genome."""

    n_centers: int = 7
    support_radius: float = 12.0
    rbf_width: float = 4.0
    alpha_amp: float = 0.15
    beta_amp: float = 0.25
    log_metric_amp: float = 0.12
    kij_amp: float = 0.05
    mutation_sigma: float = 0.15
    # Box half-width used when placing centers (must match render L/2).
    box_half_width: float = 32.0
    # Alcubierre control bounds.
    alc_velocity_max: float = 2.0
    alc_radius_min: float = 2.0
    alc_radius_max: float = 8.0
    alc_sigma_min: float = 0.5
    alc_sigma_max: float = 3.0
    # Shared radial/wall bounds for the non-shift topologies.
    feature_radius_min: float = 2.0
    feature_radius_max: float = 10.0
    feature_sigma_min: float = 0.5
    feature_sigma_max: float = 4.0
    # Deformation amplitudes for the non-shift topologies.
    tunnel_amp: float = 0.6
    lens_amp: float = 0.6
    throat_amp: float = 0.6
    # Per-topology enable flags.  When False the block is pinned to its neutral
    # default (no deformation, not searchable).
    enable_alcubierre: bool = True
    enable_tunnel: bool = True
    enable_lens: bool = True
    enable_throat: bool = True


@dataclass(frozen=True)
class GeometryGenome:
    """Immutable genome: flat coefficient vector + fixed control centers."""

    coeffs: NDArray[np.float64]
    centers: NDArray[np.float64]  # (n_centers, 3)
    config: GeometryGenomeConfig

    def __post_init__(self) -> None:
        object.__setattr__(self, "coeffs", np.asarray(self.coeffs, dtype=np.float64))
        object.__setattr__(self, "centers", np.asarray(self.centers, dtype=np.float64))
        expected = self.config.n_centers * PARAMS_PER_CENTER + ANALYTIC_PARAMS
        if self.coeffs.shape != (expected,):
            raise ValueError(
                f"coeffs shape {self.coeffs.shape} != ({expected},) "
                f"for n_centers={self.config.n_centers}"
            )
        if self.centers.shape != (self.config.n_centers, 3):
            raise ValueError(
                f"centers shape {self.centers.shape} != "
                f"({self.config.n_centers}, 3)"
            )

    @property
    def ndim(self) -> int:
        return int(self.coeffs.size)

    @property
    def rbf_coeffs(self) -> NDArray[np.float64]:
        return self.coeffs[: self.config.n_centers * PARAMS_PER_CENTER]

    @property
    def analytic_coeffs(self) -> NDArray[np.float64]:
        """The whole trailing analytic block (all ansätze)."""
        return self.coeffs[self.config.n_centers * PARAMS_PER_CENTER :]

    @property
    def alc_coeffs(self) -> NDArray[np.float64]:
        """The Alcubierre sub-block within the analytic tail."""
        return self.analytic_coeffs[_ALC_SLICE]

    def clipped(self) -> GeometryGenome:
        """Return a copy with coefficients clipped to configured amplitudes."""
        cfg = self.config
        lo, hi = _bounds_vectors(cfg)
        return GeometryGenome(
            coeffs=np.clip(self.coeffs, lo, hi),
            centers=self.centers.copy(),
            config=cfg,
        )

    def to_dict(self) -> dict:
        return {
            "coeffs": self.coeffs.tolist(),
            "centers": self.centers.tolist(),
            "config": {
                "n_centers": self.config.n_centers,
                "support_radius": self.config.support_radius,
                "rbf_width": self.config.rbf_width,
                "alpha_amp": self.config.alpha_amp,
                "beta_amp": self.config.beta_amp,
                "log_metric_amp": self.config.log_metric_amp,
                "kij_amp": self.config.kij_amp,
                "mutation_sigma": self.config.mutation_sigma,
                "box_half_width": self.config.box_half_width,
                "alc_velocity_max": self.config.alc_velocity_max,
                "alc_radius_min": self.config.alc_radius_min,
                "alc_radius_max": self.config.alc_radius_max,
                "alc_sigma_min": self.config.alc_sigma_min,
                "alc_sigma_max": self.config.alc_sigma_max,
                "feature_radius_min": self.config.feature_radius_min,
                "feature_radius_max": self.config.feature_radius_max,
                "feature_sigma_min": self.config.feature_sigma_min,
                "feature_sigma_max": self.config.feature_sigma_max,
                "tunnel_amp": self.config.tunnel_amp,
                "lens_amp": self.config.lens_amp,
                "throat_amp": self.config.throat_amp,
                "enable_alcubierre": self.config.enable_alcubierre,
                "enable_tunnel": self.config.enable_tunnel,
                "enable_lens": self.config.enable_lens,
                "enable_throat": self.config.enable_throat,
            },
        }

    @classmethod
    def from_dict(cls, payload: dict) -> GeometryGenome:
        cfg_raw = dict(payload["config"])
        # Backward-compatible defaults for older elite JSON.
        cfg_raw.setdefault("kij_amp", 0.05)
        cfg_raw.setdefault("alc_velocity_max", 0.9)
        cfg_raw.setdefault("alc_radius_min", 2.0)
        cfg_raw.setdefault("alc_radius_max", 8.0)
        cfg_raw.setdefault("alc_sigma_min", 0.5)
        cfg_raw.setdefault("alc_sigma_max", 3.0)
        cfg_raw.setdefault("feature_radius_min", 2.0)
        cfg_raw.setdefault("feature_radius_max", 10.0)
        cfg_raw.setdefault("feature_sigma_min", 0.5)
        cfg_raw.setdefault("feature_sigma_max", 4.0)
        cfg_raw.setdefault("tunnel_amp", 0.6)
        cfg_raw.setdefault("lens_amp", 0.6)
        cfg_raw.setdefault("throat_amp", 0.6)
        cfg_raw.setdefault("enable_alcubierre", True)
        cfg_raw.setdefault("enable_tunnel", True)
        cfg_raw.setdefault("enable_lens", True)
        cfg_raw.setdefault("enable_throat", True)
        cfg = GeometryGenomeConfig(**cfg_raw)
        coeffs = np.asarray(payload["coeffs"], dtype=np.float64)
        coeffs = _migrate_coeffs(coeffs, cfg)
        return cls(
            coeffs=coeffs,
            centers=np.asarray(payload["centers"], dtype=np.float64),
            config=cfg,
        )


def _migrate_coeffs(coeffs: NDArray[np.float64], cfg: GeometryGenomeConfig) -> NDArray:
    """Upgrade legacy genomes (10-param centers, Alcubierre-only tail) in place."""
    expected = cfg.n_centers * PARAMS_PER_CENTER + ANALYTIC_PARAMS
    if coeffs.shape == (expected,):
        return coeffs

    defaults = analytic_defaults(cfg)

    # Legacy Alcubierre-only tail (4 trailing params) with current centers.
    legacy_alc = cfg.n_centers * PARAMS_PER_CENTER + ALC_PARAMS
    if coeffs.shape == (legacy_alc,):
        out = np.zeros(expected, dtype=np.float64)
        out[: cfg.n_centers * PARAMS_PER_CENTER] = coeffs[: cfg.n_centers * PARAMS_PER_CENTER]
        tail = defaults.copy()
        tail[_ALC_SLICE] = coeffs[-ALC_PARAMS:]
        out[cfg.n_centers * PARAMS_PER_CENTER :] = tail
        return out

    # Very legacy: 10 params per center, no analytic tail.
    legacy_10 = cfg.n_centers * 10
    if coeffs.shape == (legacy_10,):
        out = np.zeros(expected, dtype=np.float64)
        for c in range(cfg.n_centers):
            out[c * PARAMS_PER_CENTER : c * PARAMS_PER_CENTER + 10] = coeffs[
                c * 10 : c * 10 + 10
            ]
        out[cfg.n_centers * PARAMS_PER_CENTER :] = defaults
        return out

    raise ValueError(
        f"Cannot migrate coeffs of shape {coeffs.shape} to expected {expected}"
    )


def _bounds_vectors(cfg: GeometryGenomeConfig) -> tuple[NDArray, NDArray]:
    n = cfg.n_centers * PARAMS_PER_CENTER + ANALYTIC_PARAMS
    lo = np.zeros(n, dtype=np.float64)
    hi = np.zeros_like(lo)
    for c in range(cfg.n_centers):
        base = c * PARAMS_PER_CENTER
        lo[base] = -cfg.alpha_amp
        hi[base] = cfg.alpha_amp
        lo[base + 1 : base + 4] = -cfg.beta_amp
        hi[base + 1 : base + 4] = cfg.beta_amp
        lo[base + 4 : base + 10] = -cfg.log_metric_amp
        hi[base + 4 : base + 10] = cfg.log_metric_amp
        lo[base + 10 : base + 16] = -cfg.kij_amp
        hi[base + 10 : base + 16] = cfg.kij_amp
    a0 = cfg.n_centers * PARAMS_PER_CENTER
    a_lo, a_hi = analytic_bounds(cfg)
    lo[a0:] = a_lo
    hi[a0:] = a_hi
    return lo, hi


def fibonacci_centers(n: int, radius: float) -> NDArray[np.float64]:
    """Place ``n`` centers: one at origin (if n>=1) plus Fibonacci sphere."""
    if n <= 0:
        raise ValueError("n_centers must be positive")
    if n == 1:
        return np.zeros((1, 3), dtype=np.float64)
    pts = [np.zeros(3, dtype=np.float64)]
    n_sphere = n - 1
    golden = (1.0 + np.sqrt(5.0)) / 2.0
    for i in range(n_sphere):
        z = 1.0 - 2.0 * (i + 0.5) / n_sphere
        r_xy = np.sqrt(max(0.0, 1.0 - z * z))
        theta = 2.0 * np.pi * i / golden
        pts.append(
            np.array(
                [radius * r_xy * np.cos(theta), radius * r_xy * np.sin(theta), radius * z],
                dtype=np.float64,
            )
        )
    return np.stack(pts, axis=0)


def zero_genome(cfg: GeometryGenomeConfig | None = None) -> GeometryGenome:
    """Exact Minkowski genome (all coefficients zero / neutral Alcubierre)."""
    cfg = cfg or GeometryGenomeConfig()
    centers = fibonacci_centers(cfg.n_centers, 0.45 * cfg.support_radius)
    coeffs = np.zeros(cfg.n_centers * PARAMS_PER_CENTER + ANALYTIC_PARAMS, dtype=np.float64)
    coeffs[cfg.n_centers * PARAMS_PER_CENTER :] = analytic_defaults(cfg)
    return GeometryGenome(coeffs=coeffs, centers=centers, config=cfg)


def sample_genome(
    rng: np.random.Generator,
    cfg: GeometryGenomeConfig | None = None,
) -> GeometryGenome:
    """Sample a bounded random genome."""
    cfg = cfg or GeometryGenomeConfig()
    lo, hi = _bounds_vectors(cfg)
    # Prefer milder amplitudes so most samples stay near Minkowski.
    scale = 0.5
    coeffs = rng.uniform(lo * scale, hi * scale)
    # Analytic tail: each topology neutral by default, occasionally activated.
    a0 = cfg.n_centers * PARAMS_PER_CENTER
    coeffs[a0:] = sample_analytic(rng, cfg)
    centers = fibonacci_centers(cfg.n_centers, 0.45 * cfg.support_radius)
    # Small random jitter of centers inside the support ball.
    jitter = rng.normal(0.0, 0.15 * cfg.support_radius, size=centers.shape)
    centers = centers + jitter
    norms = np.linalg.norm(centers, axis=1, keepdims=True)
    max_r = 0.7 * cfg.support_radius
    too_far = norms[:, 0] > max_r
    centers[too_far] *= (max_r / norms[too_far])
    return GeometryGenome(coeffs=coeffs, centers=centers, config=cfg).clipped()


def mutate_genome(
    parent: GeometryGenome,
    rng: np.random.Generator,
    *,
    sigma: float | None = None,
) -> GeometryGenome:
    """Gaussian mutation of coefficients with clipping to bounds."""
    cfg = parent.config
    sigma = cfg.mutation_sigma if sigma is None else float(sigma)
    lo, hi = _bounds_vectors(cfg)
    spans = hi - lo
    # Disabled ansätze have zero span, so their blocks are left untouched.
    noise = rng.normal(0.0, sigma, size=parent.coeffs.shape) * spans
    child = GeometryGenome(
        coeffs=parent.coeffs + noise,
        centers=parent.centers.copy(),
        config=cfg,
    )
    return child.clipped()


def unpack_center(
    coeffs: NDArray[np.float64], center_idx: int
) -> tuple[float, NDArray, NDArray, NDArray]:
    """Return (alpha_amp, beta[3], S_sym[3,3], K_sym[3,3]) for one center."""
    base = center_idx * PARAMS_PER_CENTER
    alpha = float(coeffs[base])
    beta = coeffs[base + 1 : base + 4].astype(np.float64)
    s_flat = coeffs[base + 4 : base + 10]
    k_flat = coeffs[base + 10 : base + 16]
    S = np.zeros((3, 3), dtype=np.float64)
    K = np.zeros((3, 3), dtype=np.float64)
    for k, (i, j) in enumerate(_SYMMETRIC_IDX):
        S[i, j] = S[j, i] = float(s_flat[k])
        K[i, j] = K[j, i] = float(k_flat[k])
    return alpha, beta, S, K


def unpack_alcubierre(genome: GeometryGenome) -> tuple[float, float, float, int]:
    """Return (velocity, bubble_radius, sigma, axis) from the trailing block."""
    alc = genome.alc_coeffs
    if not genome.config.enable_alcubierre:
        return 0.0, float(alc[1]), float(alc[2]), 0
    velocity = float(max(alc[0], 0.0))
    radius = float(alc[1])
    sigma = float(alc[2])
    axis = int(np.clip(np.floor(alc[3]), 0, 2))
    return velocity, radius, sigma, axis


def genome_bounds(cfg: GeometryGenomeConfig) -> tuple[NDArray, NDArray]:
    """Public accessor for coefficient bounds."""
    return _bounds_vectors(cfg)


def decode_fields_at_points(
    genome: GeometryGenome,
    points: NDArray[np.float64],
) -> tuple[NDArray, NDArray, NDArray, NDArray]:
    """Evaluate (alpha, beta, gamma, kij_extra) at Cartesian points (..., 3).

    The compact envelope uses a C^2 bump that vanishes (with derivatives) at
    ``support_radius``, guaranteeing asymptotic flatness by construction.  The
    analytic topologies (shift tube, tunnel, lens, throat) are added on top.
    """
    cfg = genome.config
    shape = points.shape[:-1]
    alpha = np.ones(shape, dtype=np.float64)
    beta = np.zeros(shape + (3,), dtype=np.float64)
    S = np.zeros(shape + (3, 3), dtype=np.float64)
    kij = np.zeros(shape + (3, 3), dtype=np.float64)

    # Radial distance from domain center (assumed at origin for the genome).
    r = np.linalg.norm(points, axis=-1)
    env = compact_envelope(r, cfg.support_radius)
    w2 = cfg.rbf_width * cfg.rbf_width

    for c in range(cfg.n_centers):
        a_c, b_c, S_c, K_c = unpack_center(genome.rbf_coeffs, c)
        d = points - genome.centers[c]
        r2 = np.sum(d * d, axis=-1)
        phi = np.exp(-0.5 * r2 / w2) * env
        alpha = alpha + a_c * phi
        beta = beta + b_c.reshape((1,) * len(shape) + (3,)) * phi[..., None]
        S = S + S_c.reshape((1,) * len(shape) + (3, 3)) * phi[..., None, None]
        kij = kij + K_c.reshape((1,) * len(shape) + (3, 3)) * phi[..., None, None]

    # Analytic topologies contribute to lapse, shift, and the log-metric.
    d_alpha, d_beta, d_S = decode_analytic(genome.analytic_coeffs, points, r, env, cfg)
    alpha = alpha + d_alpha
    beta = beta + d_beta
    S = S + d_S

    # Floor the lapse away from zero / negative.
    alpha = np.clip(alpha, 1.0e-3, None)
    gamma = _batch_expm_spd(S)
    return alpha, beta, gamma, kij


def compact_envelope(r: NDArray[np.float64], support_radius: float) -> NDArray[np.float64]:
    """Smooth compact support: 1 at r=0, 0 for r >= support_radius."""
    x = np.clip(r / max(support_radius, 1.0e-12), 0.0, 1.0)
    # Polynomial bump with vanishing value and first derivative at the boundary.
    return (1.0 - x * x) ** 3


def _batch_expm_spd(S: NDArray[np.float64]) -> NDArray[np.float64]:
    """Matrix exponential of a symmetric field via eigendecomposition."""
    # S is (... ,3,3).  Eigen-decompose the last two axes.
    # For numerical safety clip eigenvalues into a mild range.
    vals, vecs = np.linalg.eigh(S)
    vals = np.clip(vals, -2.0, 2.0)
    return np.einsum("...ik,...k,...jk->...ij", vecs, np.exp(vals), vecs)


def minkowski_deviation(
    alpha: NDArray, beta: NDArray, gamma: NDArray
) -> dict[str, float]:
    """Scalar diagnostics of departure from Minkowski."""
    return {
        "max_abs_alpha_m1": float(np.max(np.abs(alpha - 1.0))),
        "max_abs_beta": float(np.max(np.abs(beta))),
        "max_abs_gamma_m_I": float(
            np.max(np.abs(gamma - np.eye(3).reshape((1,) * (gamma.ndim - 2) + (3, 3))))
        ),
    }


__all__ = [
    "ALC_PARAMS",
    "PARAMS_PER_CENTER",
    "GeometryGenome",
    "GeometryGenomeConfig",
    "compact_envelope",
    "decode_fields_at_points",
    "fibonacci_centers",
    "genome_bounds",
    "minkowski_deviation",
    "mutate_genome",
    "sample_genome",
    "unpack_alcubierre",
    "unpack_center",
    "zero_genome",
]
