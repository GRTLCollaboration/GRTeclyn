"""Analytic deformation ansätze for the geometry genome.

The genome carries a trailing *analytic block*: a fixed, ordered set of named
ansätze that each paint a coherent, physically-motivated deformation on top of
the compact-support RBF field.  Unlike the RBF bumps (isotropic Gaussians that
cannot organise into a tunnel or throat), these express whole *topologies*:

* :class:`AlcubierreAnsatz` — a shift tube (warp bubble).  Shortcut via ``beta``.
* :class:`TunnelAnsatz` — anisotropic spatial contraction along an axis: the
  pure-curvature analog of a warp corridor.  Shortcut via ``gamma``.
* :class:`LensAnsatz` — a spherical conformal well (``psi^4`` factor).  Shortcut
  via isotropic ``gamma`` contraction.
* :class:`ThroatAnsatz` — a radial contraction shell (Morris-Thorne-flavoured
  throat).  Shortcut via anisotropic radial ``gamma`` contraction.

Each ansatz contributes ``(dalpha, dbeta, dS)`` where ``S`` is the symmetric
log-metric (``gamma = expm(S)``), guaranteeing positive-definiteness and — via
the shared compact envelope — asymptotic flatness by construction.

Adding a new topology means appending one :class:`AnalyticAnsatz` to
:data:`ANSATZE`; nothing else in the genome/render/score stack changes
(open/closed).
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ...initial_data.motif import alcubierre_shape_function

if TYPE_CHECKING:  # pragma: no cover - typing only, avoids circular import
    from .genome import GeometryGenomeConfig

_EPS = 1.0e-8
_EYE3 = np.eye(3)


def _axis(value: float) -> int:
    """Map a continuous axis parameter in ``[0, 3)`` to ``{0, 1, 2}``."""
    return int(np.clip(np.floor(value), 0, 2))


class AnalyticAnsatz(ABC):
    """One named analytic deformation with a fixed parameter block."""

    name: str
    n_params: int

    @abstractmethod
    def bounds(self, cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
        """Return ``(lo, hi)`` coefficient bounds, each length ``n_params``."""

    @abstractmethod
    def default(self, cfg: "GeometryGenomeConfig") -> NDArray:
        """Return the neutral (no-deformation) parameter block."""

    @abstractmethod
    def enabled(self, cfg: "GeometryGenomeConfig") -> bool:
        """Whether this ansatz is searchable for the given config."""

    @abstractmethod
    def decode(
        self,
        params: NDArray[np.float64],
        points: NDArray[np.float64],
        r: NDArray[np.float64],
        env: NDArray[np.float64],
        cfg: "GeometryGenomeConfig",
    ) -> tuple[float | NDArray, NDArray, float | NDArray]:
        """Return ``(dalpha, dbeta, dS)`` contributions at ``points``."""

    def sample(
        self,
        rng: np.random.Generator,
        cfg: "GeometryGenomeConfig",
        *,
        probability: float,
    ) -> NDArray[np.float64]:
        """Sample a parameter block: neutral default, or (rarely) within bounds."""
        if not self.enabled(cfg) or rng.random() >= probability:
            return self.default(cfg)
        lo, hi = self.bounds(cfg)
        return rng.uniform(lo, hi)


class AlcubierreAnsatz(AnalyticAnsatz):
    """Warp-bubble shift tube: ``beta`` along an axis, no curvature/lapse term."""

    name = "alcubierre"
    n_params = 4  # velocity, bubble_radius, sigma, axis

    def bounds(self, cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
        lo = np.array(
            [0.0, cfg.alc_radius_min, cfg.alc_sigma_min, 0.0], dtype=np.float64
        )
        hi = np.array(
            [cfg.alc_velocity_max, cfg.alc_radius_max, cfg.alc_sigma_max, 2.999],
            dtype=np.float64,
        )
        return lo, hi

    def default(self, cfg: "GeometryGenomeConfig") -> NDArray:
        return np.array(
            [
                0.0,
                0.5 * (cfg.alc_radius_min + cfg.alc_radius_max),
                0.5 * (cfg.alc_sigma_min + cfg.alc_sigma_max),
                0.0,
            ],
            dtype=np.float64,
        )

    def enabled(self, cfg: "GeometryGenomeConfig") -> bool:
        return bool(cfg.enable_alcubierre)

    def decode(self, params, points, r, env, cfg):
        velocity = float(max(params[0], 0.0))
        if velocity <= _EPS:
            return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), 0.0
        radius = float(params[1])
        sigma = float(params[2])
        axis = _axis(float(params[3]))
        # Shape function already decays at large r; do not re-apply the RBF
        # envelope (it would clip the bubble wall by ~50%).
        f = alcubierre_shape_function(r, bubble_radius=radius, sigma=sigma)
        dbeta = np.zeros(points.shape[:-1] + (3,), dtype=np.float64)
        dbeta[..., axis] = -velocity * f
        return 0.0, dbeta, 0.0


class TunnelAnsatz(AnalyticAnsatz):
    """Anisotropic contraction along an axis: a pure-curvature warp corridor."""

    name = "tunnel"
    n_params = 3  # strength, width, axis

    def bounds(self, cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
        lo = np.array(
            [-cfg.tunnel_amp, cfg.feature_radius_min, 0.0], dtype=np.float64
        )
        hi = np.array(
            [cfg.tunnel_amp, cfg.feature_radius_max, 2.999], dtype=np.float64
        )
        return lo, hi

    def default(self, cfg: "GeometryGenomeConfig") -> NDArray:
        return np.array(
            [0.0, 0.5 * (cfg.feature_radius_min + cfg.feature_radius_max), 0.0],
            dtype=np.float64,
        )

    def enabled(self, cfg: "GeometryGenomeConfig") -> bool:
        return bool(cfg.enable_tunnel)

    def decode(self, params, points, r, env, cfg):
        strength = float(params[0])
        if abs(strength) <= _EPS:
            return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), 0.0
        width = max(float(params[1]), 0.5)
        axis = _axis(float(params[2]))
        # Perpendicular distance from the axis line through the origin.
        perp2 = np.clip(r * r - points[..., axis] ** 2, 0.0, None)
        prof = np.exp(-0.5 * perp2 / (width * width)) * env
        dS = np.zeros(points.shape[:-1] + (3, 3), dtype=np.float64)
        # Contract the along-axis metric component: gamma_|| = exp(-strength*prof).
        dS[..., axis, axis] = -strength * prof
        return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), dS


class LensAnsatz(AnalyticAnsatz):
    """Spherical conformal well: isotropic ``gamma`` contraction in a shell."""

    name = "lens"
    n_params = 3  # strength, radius, sigma

    def bounds(self, cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
        lo = np.array([-cfg.lens_amp, 0.0, cfg.feature_sigma_min], dtype=np.float64)
        hi = np.array(
            [cfg.lens_amp, cfg.feature_radius_max, cfg.feature_sigma_max],
            dtype=np.float64,
        )
        return lo, hi

    def default(self, cfg: "GeometryGenomeConfig") -> NDArray:
        return np.array(
            [0.0, 0.0, 0.5 * (cfg.feature_sigma_min + cfg.feature_sigma_max)],
            dtype=np.float64,
        )

    def enabled(self, cfg: "GeometryGenomeConfig") -> bool:
        return bool(cfg.enable_lens)

    def decode(self, params, points, r, env, cfg):
        strength = float(params[0])
        if abs(strength) <= _EPS:
            return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), 0.0
        radius = float(params[1])
        sigma = max(float(params[2]), 0.25)
        prof = np.exp(-0.5 * ((r - radius) / sigma) ** 2) * env
        s_iso = (strength * prof)[..., None, None]
        eye = _EYE3.reshape((1,) * (points.ndim - 1) + (3, 3))
        dS = s_iso * eye
        return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), dS


class ThroatAnsatz(AnalyticAnsatz):
    """Radial contraction shell: a Morris-Thorne-flavoured wormhole throat."""

    name = "throat"
    n_params = 3  # strength, throat_radius, sigma

    def bounds(self, cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
        lo = np.array(
            [-cfg.throat_amp, cfg.feature_radius_min, cfg.feature_sigma_min],
            dtype=np.float64,
        )
        hi = np.array(
            [cfg.throat_amp, cfg.feature_radius_max, cfg.feature_sigma_max],
            dtype=np.float64,
        )
        return lo, hi

    def default(self, cfg: "GeometryGenomeConfig") -> NDArray:
        return np.array(
            [
                0.0,
                0.5 * (cfg.feature_radius_min + cfg.feature_radius_max),
                0.5 * (cfg.feature_sigma_min + cfg.feature_sigma_max),
            ],
            dtype=np.float64,
        )

    def enabled(self, cfg: "GeometryGenomeConfig") -> bool:
        return bool(cfg.enable_throat)

    def decode(self, params, points, r, env, cfg):
        strength = float(params[0])
        if abs(strength) <= _EPS:
            return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), 0.0
        radius = float(params[1])
        sigma = max(float(params[2]), 0.25)
        prof = np.exp(-0.5 * ((r - radius) / sigma) ** 2) * env
        rhat = points / np.clip(r[..., None], _EPS, None)
        outer = rhat[..., :, None] * rhat[..., None, :]
        # Contract the radial metric component across the shell.
        dS = (-strength * prof)[..., None, None] * outer
        return 0.0, np.zeros(points.shape[:-1] + (3,), dtype=np.float64), dS


# Ordered registry.  Alcubierre stays first so its 4-param block keeps a stable
# offset for legacy elite JSON (see genome._migrate_coeffs).
ANSATZE: tuple[AnalyticAnsatz, ...] = (
    AlcubierreAnsatz(),
    TunnelAnsatz(),
    LensAnsatz(),
    ThroatAnsatz(),
)

ANALYTIC_PARAMS: int = sum(a.n_params for a in ANSATZE)


def ansatz_offset(name: str) -> int:
    """Return the start index of ``name``'s block within the analytic tail."""
    off = 0
    for a in ANSATZE:
        if a.name == name:
            return off
        off += a.n_params
    raise KeyError(name)


def ansatz_slice(name: str) -> slice:
    """Return the slice of ``name``'s block within the analytic tail."""
    off = ansatz_offset(name)
    for a in ANSATZE:
        if a.name == name:
            return slice(off, off + a.n_params)
    raise KeyError(name)  # pragma: no cover - unreachable


def analytic_bounds(cfg: "GeometryGenomeConfig") -> tuple[NDArray, NDArray]:
    """Concatenated ``(lo, hi)`` for the whole analytic tail.

    Disabled ansätze are pinned (``lo == hi == default``) so sampling and
    mutation leave them untouched.
    """
    los: list[NDArray] = []
    his: list[NDArray] = []
    for a in ANSATZE:
        if a.enabled(cfg):
            lo, hi = a.bounds(cfg)
        else:
            lo = hi = a.default(cfg)
        los.append(lo)
        his.append(hi)
    return np.concatenate(los), np.concatenate(his)


def analytic_defaults(cfg: "GeometryGenomeConfig") -> NDArray[np.float64]:
    """Concatenated neutral parameter blocks (renders exact Minkowski)."""
    return np.concatenate([a.default(cfg) for a in ANSATZE])


def sample_analytic(
    rng: np.random.Generator,
    cfg: "GeometryGenomeConfig",
    *,
    probability: float = 0.4,
) -> NDArray[np.float64]:
    """Sample the analytic tail: each ansatz neutral, or occasionally active."""
    return np.concatenate(
        [a.sample(rng, cfg, probability=probability) for a in ANSATZE]
    )


def decode_analytic(
    analytic_coeffs: NDArray[np.float64],
    points: NDArray[np.float64],
    r: NDArray[np.float64],
    env: NDArray[np.float64],
    cfg: "GeometryGenomeConfig",
) -> tuple[float | NDArray, NDArray, float | NDArray]:
    """Sum ``(dalpha, dbeta, dS)`` over all enabled ansätze."""
    shape = points.shape[:-1]
    dalpha: float | NDArray = 0.0
    dbeta = np.zeros(shape + (3,), dtype=np.float64)
    dS: float | NDArray = 0.0
    off = 0
    for a in ANSATZE:
        block = analytic_coeffs[off : off + a.n_params]
        off += a.n_params
        if not a.enabled(cfg):
            continue
        da, db, ds = a.decode(block, points, r, env, cfg)
        dalpha = dalpha + da
        dbeta = dbeta + db
        dS = dS + ds
    return dalpha, dbeta, dS


__all__ = [
    "ANALYTIC_PARAMS",
    "ANSATZE",
    "AlcubierreAnsatz",
    "AnalyticAnsatz",
    "LensAnsatz",
    "ThroatAnsatz",
    "TunnelAnsatz",
    "analytic_bounds",
    "analytic_defaults",
    "ansatz_offset",
    "ansatz_slice",
    "decode_analytic",
    "sample_analytic",
]
