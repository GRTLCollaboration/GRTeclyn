"""Fit GRTresna scalar-lump ansätze from geometry-first motifs.

Matter and momentum fitting are separate stages so static-lens support and
transport motors can be tested independently.  Output is always a GRTresnaConfig
candidate -- never final GRTeclyn initial data.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ...initial_data.motif import GeometryMotif, MomentumTarget, SupportRegion
from ..profiles.boson_star import PROFILE_ODE_BOUND
from ..profiles.qball_couplings import QBallCouplings
from ..solver import GRTresnaConfig, apply_exotic_safe_solver

DEFAULT_AMP_SCALE = 0.15
DEFAULT_WIDTH = 3.0
MAX_LUMP_AMP = 0.35
MIN_LUMP_WIDTH = 1.5
MAX_LUMP_WIDTH = 5.0
MAX_VELOCITY = 0.6

# Q-ball compact preset for toroidal (warp-drive) exotic lumps.
# Uses the same couplings as the QD MAP-Elites campaign (m=1, λ=640, μ=85333,
# ω=0.8) so GRTresna sees a bound soliton, not a free Klein-Gordon wave packet.
WARP_QBALL_COUPLINGS = QBallCouplings.compact()
WARP_QBALL_AMP = WARP_QBALL_COUPLINGS.with_equilibrium_paint().phi_core


@dataclass(frozen=True)
class FittedMatter:
    """GRTresna lump basis fitted from a geometry motif."""

    lumps: tuple[dict[str, Any], ...]
    scalar_mass: float
    maximal_slicing: bool
    static_lens_only: bool
    momentum_target: MomentumTarget
    notes: tuple[str, ...]


def _lump_dict(
    *,
    amp: float,
    width: float,
    center: tuple[float, float, float],
    velocity: tuple[float, float, float] = (0.0, 0.0, 0.0),
    omega: float = 0.0,
    mode: int = 0,
    exotic: int = 0,
) -> dict[str, Any]:
    return {
        "amp": float(amp),
        "width": float(width),
        "center": tuple(float(v) for v in center),
        "velocity": tuple(float(v) for v in velocity),
        "omega": float(omega),
        "mode": int(mode),
        "exotic": int(exotic),
    }


def _amp_from_rho(peak_rho: float) -> float:
    amp = min(MAX_LUMP_AMP, max(0.05, abs(peak_rho) * 8.0))
    return float(amp)


def _region_to_lump(region: SupportRegion) -> dict[str, Any]:
    return _lump_dict(
        amp=_amp_from_rho(region.peak_rho),
        width=max(MIN_LUMP_WIDTH, min(MAX_LUMP_WIDTH, region.width)),
        center=region.center,
        exotic=1 if region.exotic else 0,
    )


def fit_matter_from_motif(
    motif: GeometryMotif,
    *,
    max_lumps: int = 3,
) -> FittedMatter:
    """Fit scalar lumps to geometry-first rho_req support.

    When fewer support regions than *max_lumps* are available, the dominant
    region is split into multiple smaller lumps distributed around its center
    (ring layout in the xy-plane).  This gives the solver more degrees of
    freedom to realize strong-curvature targets that a single Gaussian cannot
    support within the convergence basin.

    For **warp-drive (toroidal)** motifs — identified by
    ``momentum_target.template == "toroidal"`` — the ring is placed in the
    plane **perpendicular** to the transport axis (yz-plane for x-axis
    travel) at the cylindrical radius of the exotic support region, matching
    the toroidal negative-energy ring of an Alcubierre bubble.
    """
    notes: list[str] = []
    regions = list(motif.support_regions)
    if not regions:
        notes.append("no support regions found; placing one canonical lump at origin")
        regions = [
            SupportRegion(
                center=(0.0, 0.0, 0.0),
                width=DEFAULT_WIDTH,
                peak_rho=DEFAULT_AMP_SCALE / 8.0,
                exotic=False,
                radial_center=0.0,
            )
        ]

    exotic_regions = [region for region in regions if region.exotic]
    canonical_regions = [region for region in regions if not region.exotic]
    ordered = exotic_regions + canonical_regions
    selected = ordered[: max(1, max_lumps)]

    is_toroidal = motif.momentum_target.template == "toroidal"

    # If we have fewer regions than max_lumps, split the dominant region
    # into additional lumps arranged in a ring around its center.
    if len(selected) < max_lumps and len(selected) >= 1:
        n_extra = max_lumps - len(selected)
        dominant = selected[0]
        ring_radius = max(dominant.width * 0.8, 2.0)
        # Amplitude is split evenly; each sub-lump gets 1/n of the original
        sub_amp = _amp_from_rho(dominant.peak_rho) / max_lumps
        sub_width = max(MIN_LUMP_WIDTH, min(MAX_LUMP_WIDTH, dominant.width * 0.7))
        cx, cy, cz = dominant.center

        if is_toroidal:
            # Toroidal ring in the plane perpendicular to the transport axis.
            # For transport along x, the ring is in the yz-plane at x=cx,
            # at cylindrical radius = |dominant.radial_center| (the bubble
            # wall radius).  This matches the Alcubierre negative-energy ring.
            # The ring is centred on the transport axis (y=0, z=0), NOT on the
            # original support region's centre, so all lumps sit at the same
            # cylindrical radius from the axis.
            ring_radius = max(dominant.radial_center, ring_radius)
            for i in range(n_extra):
                angle = 2.0 * math.pi * (i + 1) / max_lumps
                offset_y = ring_radius * math.cos(angle)
                offset_z = ring_radius * math.sin(angle)
                selected.append(SupportRegion(
                    center=(cx, offset_y, offset_z),
                    width=sub_width,
                    peak_rho=dominant.peak_rho / max_lumps,
                    exotic=dominant.exotic,
                    radial_center=dominant.radial_center,
                ))
            notes.append(
                f"split dominant support into {max_lumps} toroidal lumps "
                f"(ring_radius={ring_radius:.2f}, plane=yz)"
            )
        else:
            for i in range(n_extra):
                angle = 2.0 * math.pi * (i + 1) / max_lumps
                offset_x = ring_radius * math.cos(angle)
                offset_y = ring_radius * math.sin(angle)
                selected.append(SupportRegion(
                    center=(cx + offset_x, cy + offset_y, cz),
                    width=sub_width,
                    peak_rho=dominant.peak_rho / max_lumps,
                    exotic=dominant.exotic,
                    radial_center=dominant.radial_center,
                ))
            notes.append(
                f"split dominant support region into {max_lumps} ring-distributed lumps "
                f"(ring_radius={ring_radius:.2f})"
            )

    lumps = tuple(_region_to_lump(region) for region in selected)
    has_exotic = any(lump["exotic"] for lump in lumps)

    # For toroidal (warp-drive) motifs, upgrade lumps to Q-ball solitons
    # using the compact preset.  Plain Gaussian lumps with scalar_mass=0 have
    # no binding and barely curve spacetime (Ham ~100%), so GRTresna cannot
    # converge.  Q-ball lumps with the compact preset (m=1, λ=640, μ=85333,
    # ω=0.8) are bound solitons that the QD MAP-Elites campaign converges
    # routinely.
    scalar_mass = 0.0
    scalar_lambda = 0.0
    scalar_mu = 0.0
    if is_toroidal:
        qball = WARP_QBALL_COUPLINGS
        scalar_mass = qball.mass
        scalar_lambda = qball.lam
        scalar_mu = qball.mu
        qball_amp = WARP_QBALL_AMP
        updated_lumps = []
        for lump in lumps:
            updated = dict(lump)
            updated["amp"] = qball_amp
            updated["omega"] = qball.omega
            updated["profile"] = PROFILE_ODE_BOUND
            updated["qball_mass"] = qball.mass
            updated["qball_lam"] = qball.lam
            updated["qball_mu"] = qball.mu
            updated["qball_omega"] = qball.omega
            updated_lumps.append(updated)
        lumps = tuple(updated_lumps)
        notes.append(
            f"toroidal lumps upgraded to Q-ball solitons "
            f"(m={qball.mass}, λ={qball.lam}, μ={qball.mu}, ω={qball.omega}, "
            f"amp={qball_amp:.4f}, profile=ODE)"
        )

    if has_exotic:
        notes.append("exotic lumps present; maximal_slicing required")

    return FittedMatter(
        lumps=lumps,
        scalar_mass=scalar_mass,
        maximal_slicing=has_exotic or motif.exotic_needed,
        static_lens_only=motif.static_lens_only,
        momentum_target=motif.momentum_target,
        notes=tuple(notes),
    )


def _estimate_lump_pi(
    lump: Mapping[str, Any],
    *,
    sample_point: tuple[float, float, float] | None = None,
) -> float:
    """Python mirror of GRTresna MatterParams::lump_pi for unit tests."""
    amp = float(lump.get("amp", 0.0))
    width = float(lump.get("width", 5.0))
    center = tuple(float(v) for v in lump.get("center", (0.0, 0.0, 0.0)))
    velocity = tuple(float(v) for v in lump.get("velocity", (0.0, 0.0, 0.0)))
    omega = float(lump.get("omega", 0.0))
    exotic = int(lump.get("exotic", 0))
    if amp == 0.0:
        return 0.0
    if not any(abs(v) > 0.0 for v in velocity) and abs(omega) < 1.0e-12:
        return 0.0

    effective_amp = 0.25 * amp if exotic else amp
    loc = np.asarray(sample_point if sample_point is not None else center, dtype=float)
    dx = loc - np.asarray(center, dtype=float)
    r2 = float(np.dot(dx, dx))
    env = math.exp(-r2 / (2.0 * width * width))
    mode = int(lump.get("mode", 0))
    if mode == 1:
        angular = dx[0] / width
    elif mode == 2:
        angular = (dx[0] * dx[0] - dx[1] * dx[1]) / (width * width)
    else:
        angular = 1.0
    phi = effective_amp * angular * env

    eps = 1.0e-3 * width
    grad = np.zeros(3, dtype=float)
    for axis in range(3):
        lp = loc.copy()
        lm = loc.copy()
        lp[axis] += eps
        lm[axis] -= eps

        def phi_at(point: np.ndarray) -> float:
            d = point - np.asarray(center, dtype=float)
            r2_local = float(np.dot(d, d))
            env_local = math.exp(-r2_local / (2.0 * width * width))
            if mode == 1:
                ang = d[0] / width
            elif mode == 2:
                ang = (d[0] * d[0] - d[1] * d[1]) / (width * width)
            else:
                ang = 1.0
            return effective_amp * ang * env_local

        grad[axis] = (phi_at(lp) - phi_at(lm)) / (2.0 * eps)

    boost = -float(np.dot(velocity, grad))
    rot = -omega * (dx[0] * grad[1] - dx[1] * grad[0])
    return boost + rot


def _lump_phi_grad(
    lump: Mapping[str, Any],
    point: np.ndarray,
) -> tuple[float, np.ndarray]:
    amp = float(lump.get("amp", 0.0))
    width = float(lump.get("width", 5.0))
    center = np.asarray(lump.get("center", (0.0, 0.0, 0.0)), dtype=float)
    exotic = int(lump.get("exotic", 0))
    mode = int(lump.get("mode", 0))
    effective_amp = 0.25 * amp if exotic else amp
    if effective_amp == 0.0:
        return 0.0, np.zeros(3, dtype=float)

    eps = 1.0e-3 * width
    grad = np.zeros(3, dtype=float)

    def phi_at(pt: np.ndarray) -> float:
        d = pt - center
        r2_local = float(np.dot(d, d))
        env_local = math.exp(-r2_local / (2.0 * width * width))
        if mode == 1:
            ang = d[0] / width
        elif mode == 2:
            ang = (d[0] * d[0] - d[1] * d[1]) / (width * width)
        else:
            ang = 1.0
        return effective_amp * ang * env_local

    for axis in range(3):
        lp = point.copy()
        lm = point.copy()
        lp[axis] += eps
        lm[axis] -= eps
        grad[axis] = (phi_at(lp) - phi_at(lm)) / (2.0 * eps)
    return phi_at(point), grad


def estimate_momentum_source(
    lumps: Sequence[Mapping[str, Any]],
    *,
    sample_point: tuple[float, float, float] | None = None,
) -> tuple[float, float, float]:
    """Estimate S_i ~ -Pi grad_i phi at a sample point from fitted lumps."""
    point = np.asarray(sample_point if sample_point is not None else (0.0, 0.0, 0.0), dtype=float)
    grad_total = np.zeros(3, dtype=float)
    pi_total = 0.0
    for lump in lumps:
        pi = _estimate_lump_pi(lump, sample_point=tuple(point))
        _phi, grad = _lump_phi_grad(lump, point)
        pi_total += pi
        grad_total += grad

    source = -pi_total * grad_total
    return (float(source[0]), float(source[1]), float(source[2]))


def fit_momentum_from_motif(
    motif: GeometryMotif,
    fitted_matter: FittedMatter,
) -> FittedMatter:
    """Attach velocity/omega so GRTresna can solve both H and M_i."""
    target = motif.momentum_target
    notes = list(fitted_matter.notes)

    if motif.static_lens_only or not target.credible:
        notes.append("momentum fit skipped: static_lens_only")
        return FittedMatter(
            lumps=fitted_matter.lumps,
            scalar_mass=fitted_matter.scalar_mass,
            maximal_slicing=fitted_matter.maximal_slicing,
            static_lens_only=True,
            momentum_target=target,
            notes=tuple(notes),
        )

    direction = np.asarray(target.direction, dtype=float)
    norm = float(np.linalg.norm(direction))
    if norm < 1.0e-12:
        notes.append("momentum fit skipped: zero target direction")
        return FittedMatter(
            lumps=fitted_matter.lumps,
            scalar_mass=fitted_matter.scalar_mass,
            maximal_slicing=fitted_matter.maximal_slicing,
            static_lens_only=True,
            momentum_target=target,
            notes=tuple(notes),
        )
    direction = direction / norm
    velocity_mag = min(MAX_VELOCITY, 0.25 + 0.35 * target.strength)

    updated_lumps: list[dict[str, Any]] = []
    for idx, lump in enumerate(fitted_matter.lumps):
        center = tuple(float(v) for v in lump["center"])
        velocity = (
            float(direction[0] * velocity_mag),
            float(direction[1] * velocity_mag),
            float(direction[2] * velocity_mag),
        )
        if target.template == "counter_moving_pair" and idx % 2 == 1:
            velocity = tuple(-v for v in velocity)
        updated = dict(lump)
        updated["velocity"] = velocity
        if target.template == "toroidal":
            # For Q-ball lumps, keep the Q-ball omega (set in fit_matter).
            # Only set omega from strength if it wasn't already set.
            if "qball_omega" not in lump:
                updated["omega"] = 0.2 * target.strength
            updated["mode"] = max(int(lump.get("mode", 0)), 1)
        updated_lumps.append(updated)

    notes.append(
        f"momentum fit applied: template={target.template}, |v|={velocity_mag:.3f}"
    )
    return FittedMatter(
        lumps=tuple(updated_lumps),
        scalar_mass=fitted_matter.scalar_mass,
        maximal_slicing=fitted_matter.maximal_slicing,
        static_lens_only=False,
        momentum_target=target,
        notes=tuple(notes),
    )


def build_grtresna_config_from_fitted(
    fitted: FittedMatter,
    base: GRTresnaConfig | None = None,
) -> GRTresnaConfig:
    """Build a GRTresnaConfig without touching matter-first search helpers."""
    import dataclasses

    cfg = dataclasses.replace(base) if base is not None else GRTresnaConfig()
    cfg.scalar_mass = fitted.scalar_mass
    # Propagate Q-ball self-interaction couplings for toroidal (warp) lumps.
    # FittedMatter only carries scalar_mass directly; lambda/mu come from the
    # lump dicts (set by the toroidal Q-ball upgrade in fit_matter_from_motif).
    if fitted.lumps and "qball_lam" in fitted.lumps[0]:
        cfg.scalar_lambda = float(fitted.lumps[0]["qball_lam"])
        cfg.scalar_mu = float(fitted.lumps[0]["qball_mu"])
        cfg.bs_omega = float(fitted.lumps[0].get("qball_omega", 0.0))
    cfg.lumps = [dict(lump) for lump in fitted.lumps]
    if fitted.maximal_slicing:
        cfg.maximal_slicing = True
        apply_exotic_safe_solver(cfg)
    return cfg


def fitted_matter_to_dict(fitted: FittedMatter) -> dict[str, Any]:
    payload = asdict(fitted)
    payload["momentum_target"] = asdict(fitted.momentum_target)
    return payload


def fitted_matter_from_dict(payload: Mapping[str, Any]) -> FittedMatter:
    momentum = payload["momentum_target"]
    return FittedMatter(
        lumps=tuple(payload.get("lumps", ())),
        scalar_mass=float(payload.get("scalar_mass", 0.0)),
        maximal_slicing=bool(payload.get("maximal_slicing", False)),
        static_lens_only=bool(payload.get("static_lens_only", True)),
        momentum_target=MomentumTarget(**momentum),
        notes=tuple(payload.get("notes", ())),
    )


def write_fitted_matter_json(fitted: FittedMatter, path: str | Path) -> None:
    path = Path(path)
    path.write_text(
        json.dumps(fitted_matter_to_dict(fitted), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_fitted_matter_json(path: str | Path) -> FittedMatter:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return fitted_matter_from_dict(payload)


def write_momentum_target_json(target: MomentumTarget, path: str | Path) -> None:
    path = Path(path)
    path.write_text(
        json.dumps(asdict(target), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
