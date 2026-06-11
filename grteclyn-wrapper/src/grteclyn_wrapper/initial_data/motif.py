"""Extract compact geometry-first motifs from RadialRecipe search episodes.

Geometry-first scouts produce static lens geometries with algebraically
reconstructed matter.  This module distils those episodes into a small target
descriptor that the GRTresna projection pipeline can fit without knowing how the
motif was found.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..initial_data.constrained_recipe import RecipeBasis, constrained_phi
from ..metrics.diagnostics.constraints import read_constraint_metrics
from ..metrics.probes.ftl.analytic import compute_ftl_metrics, load_overrides_from_episode

MOTOR_BETA_EPS = 0.05
RHO_SUPPORT_FRAC = 0.15
DEFAULT_TRANSPORT_AXIS = (1.0, 0.0, 0.0)


@dataclass(frozen=True)
class SupportRegion:
    """Localized energy-support target derived from rho_req."""

    center: tuple[float, float, float]
    width: float
    peak_rho: float
    exotic: bool
    radial_center: float


@dataclass(frozen=True)
class MomentumTarget:
    """Desired scalar momentum motor for GRTresna lump fitting."""

    direction: tuple[float, float, float]
    support_center: tuple[float, float, float]
    strength: float
    template: str
    credible: bool
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class GeometryMotif:
    """Compact target distilled from a geometry-first episode."""

    episode_path: str
    overrides: dict[str, Any]
    transport_axis: tuple[float, float, float]
    polarity: float
    f_shortcut: float
    f_op: float | None
    f_null: float
    f_portal: float
    f_throat: float
    f_asymmetry: float
    beta_max: float
    beta_mean: float
    exotic_needed: bool
    min_rho_required: float | None
    integral_negative_rho: float | None
    static_lens_only: bool
    support_regions: tuple[SupportRegion, ...]
    momentum_target: MomentumTarget
    notes: tuple[str, ...] = field(default_factory=tuple)


def _unit(vec: Sequence[float]) -> tuple[float, float, float]:
    arr = np.asarray(vec, dtype=float)
    norm = float(np.linalg.norm(arr))
    if norm < 1.0e-12:
        return (0.0, 0.0, 0.0)
    arr = arr / norm
    return (float(arr[0]), float(arr[1]), float(arr[2]))


def _beta_profile(overrides: Mapping[str, Any], *, L: float) -> tuple[np.ndarray, np.ndarray]:
    from ..metrics.probes.ftl.analytic import _axis_profiles

    x, _r, _chi, _alpha, beta_x = _axis_profiles(overrides, L=L, n_points=512)
    return x, beta_x


def _support_regions_from_rho(
    r: np.ndarray,
    rho_required: np.ndarray,
    *,
    basis_radius_max: float,
) -> list[SupportRegion]:
    peak = float(np.max(np.abs(rho_required)))
    if peak < 1.0e-10:
        return []

    threshold = RHO_SUPPORT_FRAC * peak
    mask = np.abs(rho_required) >= threshold
    if not mask.any():
        return []

    regions: list[SupportRegion] = []
    indices = np.where(mask)[0]
    start = int(indices[0])
    prev = start
    for idx in indices[1:]:
        if idx != prev + 1:
            regions.append(_region_from_segment(r, rho_required, start, prev, basis_radius_max))
            start = int(idx)
        prev = int(idx)
    regions.append(_region_from_segment(r, rho_required, start, prev, basis_radius_max))
    return regions


def _region_from_segment(
    r: np.ndarray,
    rho_required: np.ndarray,
    start: int,
    end: int,
    basis_radius_max: float,
) -> SupportRegion:
    segment_r = r[start : end + 1]
    segment_rho = rho_required[start : end + 1]
    peak_idx = int(np.argmax(np.abs(segment_rho)))
    radial_center = float(segment_r[peak_idx])
    peak_rho = float(segment_rho[peak_idx])
    width = max(1.0, float(segment_r[-1] - segment_r[0]) * 0.5 + 0.5)
    width = min(width, 0.5 * basis_radius_max)
    return SupportRegion(
        center=(radial_center, 0.0, 0.0),
        width=width,
        peak_rho=peak_rho,
        exotic=peak_rho < 0.0,
        radial_center=radial_center,
    )


def _infer_momentum_target(
    overrides: Mapping[str, Any],
    *,
    L: float,
    beta_max: float,
    beta_mean: float,
    support_regions: Sequence[SupportRegion],
) -> MomentumTarget:
    notes: list[str] = []
    if beta_max < MOTOR_BETA_EPS:
        notes.append("shift amplitude below motor threshold")
        return MomentumTarget(
            direction=(0.0, 0.0, 0.0),
            support_center=(0.0, 0.0, 0.0),
            strength=0.0,
            template="none",
            credible=False,
            notes=tuple(notes),
        )

    x, beta_x = _beta_profile(overrides, L=L)
    weights = np.abs(beta_x)
    dipole = float(np.trapezoid(beta_x * x, x))
    peak_idx = int(np.argmax(weights))
    peak_sign = float(np.sign(beta_x[peak_idx])) if weights[peak_idx] > 0.0 else 0.0
    signed = dipole if abs(dipole) >= MOTOR_BETA_EPS else peak_sign

    if abs(signed) < 1.0e-3:
        notes.append("no coherent signed shift dipole")
        return MomentumTarget(
            direction=(0.0, 0.0, 0.0),
            support_center=(0.0, 0.0, 0.0),
            strength=0.0,
            template="none",
            credible=False,
            notes=tuple(notes),
        )

    direction = _unit((math.copysign(1.0, signed), 0.0, 0.0))
    support_center = support_regions[0].center if support_regions else (0.0, 0.0, 0.0)
    strength = min(1.0, beta_max / 0.5)
    notes.append("axial boost inferred from shift dipole")
    return MomentumTarget(
        direction=direction,
        support_center=support_center,
        strength=strength,
        template="axial_boost",
        credible=True,
        notes=tuple(notes),
    )


def extract_motif_from_episode(
    episode_dir: str | Path,
    *,
    phantom: bool = True,
    ftl_L: float | None = None,
) -> GeometryMotif:
    """Load a geometry-first episode and extract a projection target."""
    episode_path = Path(episode_dir).expanduser().resolve()
    overrides = load_overrides_from_episode(episode_path)
    if overrides is None:
        raise FileNotFoundError(
            f"No overrides found in {episode_path} (metadata.json or params.txt)"
        )

    basis_radius_max = float(overrides.get("recipe_basis_radius_max", 8.0))
    integration_L = ftl_L if ftl_L is not None else basis_radius_max
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=basis_radius_max,
    )
    chi_coeffs = [
        float(overrides.get(f"recipe_chi_coeff_{n}", 0.0))
        for n in range(num_bases)
    ]
    constrained = constrained_phi(
        basis=basis,
        chi_asymptotic=float(overrides.get("recipe_chi_asymptotic", 1.0)),
        chi_coeffs=chi_coeffs,
        K_const=float(overrides.get("recipe_K_asymptotic", 0.0)),
        phantom=phantom,
    )
    support_regions = tuple(
        _support_regions_from_rho(
            constrained.r,
            constrained.rho_required,
            basis_radius_max=basis_radius_max,
        )
    )

    ftl = compute_ftl_metrics(overrides, L=integration_L)
    x, beta_x = _beta_profile(overrides, L=integration_L)
    beta_max = float(np.max(np.abs(beta_x)))
    beta_mean = float(np.mean(beta_x[np.abs(beta_x) >= 0.1 * max(beta_max, 1.0e-8)])) if beta_max >= 1.0e-8 else 0.0

    constraint_path = episode_path / "data" / "constraint_norms.dat"
    constraints = read_constraint_metrics(constraint_path)
    min_rho_required = constraints.min_rho_required if constraints else None
    integral_negative_rho = constraints.integral_negative_rho if constraints else None
    exotic_needed = (
        (min_rho_required is not None and min_rho_required < 0.0)
        or any(region.exotic for region in support_regions)
        or phantom
    )

    momentum_target = _infer_momentum_target(
        overrides,
        L=integration_L,
        beta_max=beta_max,
        beta_mean=beta_mean,
        support_regions=support_regions,
    )
    static_lens_only = not momentum_target.credible

    score_path = episode_path / "score.json"
    f_op: float | None = None
    if score_path.exists():
        score_payload = json.loads(score_path.read_text(encoding="utf-8"))
        general = score_payload.get("metrics", {}).get("general_ftl") or {}
        if isinstance(general, dict):
            f_op = general.get("f_op")

    notes = list(ftl.notes)
    if static_lens_only:
        notes.append("static_lens_only: geometry-first reconstruction has no credible motor")

    return GeometryMotif(
        episode_path=str(episode_path),
        overrides=dict(overrides),
        transport_axis=DEFAULT_TRANSPORT_AXIS,
        polarity=ftl.f_asymmetry,
        f_shortcut=ftl.f_shortcut,
        f_op=f_op,
        f_null=ftl.f_null,
        f_portal=ftl.f_portal,
        f_throat=ftl.f_throat,
        f_asymmetry=ftl.f_asymmetry,
        beta_max=beta_max,
        beta_mean=beta_mean,
        exotic_needed=exotic_needed,
        min_rho_required=min_rho_required,
        integral_negative_rho=integral_negative_rho,
        static_lens_only=static_lens_only,
        support_regions=support_regions,
        momentum_target=momentum_target,
        notes=tuple(notes),
    )


def motif_to_dict(motif: GeometryMotif) -> dict[str, Any]:
    return asdict(motif)


def motif_from_dict(payload: Mapping[str, Any]) -> GeometryMotif:
    support_regions = tuple(
        SupportRegion(**region) for region in payload.get("support_regions", ())
    )
    momentum = MomentumTarget(**payload["momentum_target"])
    data = dict(payload)
    data["support_regions"] = support_regions
    data["momentum_target"] = momentum
    data["notes"] = tuple(data.get("notes", ()))
    return GeometryMotif(**data)


def write_motif_json(motif: GeometryMotif, path: str | Path) -> None:
    path = Path(path)
    path.write_text(
        json.dumps(motif_to_dict(motif), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_motif_json(path: str | Path) -> GeometryMotif:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return motif_from_dict(payload)
