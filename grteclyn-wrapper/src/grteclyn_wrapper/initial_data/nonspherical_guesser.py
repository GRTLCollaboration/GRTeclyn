"""First-stage non-spherical metric guesser validation.

This module deliberately stops short of a full 3D constraint solve.  It
extends the radial conformal factor with axisymmetric angular modes and
validates the proposed metric along 1D radial rays:

    chi(r, theta) = chi_0(r) + sum_a A_a exp(-(r-r_a)^2/(2 w_a^2)) P_l(cos theta)

The checks here answer a narrow question: does the proposed non-spherical
geometry remain a valid positive conformal metric along sampled rays, and how
large is the angular deformation?  Hamiltonian/momentum constraints for these
non-spherical data still require the future 3D C++/elliptic layer.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .constrained_recipe import RecipeBasis


ACCEPTED_RAY_SANE = "accepted_ray_sane"
REJECTED_NEGATIVE_CHI = "rejected_negative_chi"
REJECTED_HIGH_ANGULAR_CONTRAST = "rejected_high_angular_contrast"


@dataclass(frozen=True)
class AngularMode:
    """One axisymmetric angular deformation."""

    ell: int
    amplitude: float
    radial_center: float
    radial_width: float


@dataclass(frozen=True)
class NonSphericalCandidate:
    """A radial base profile plus angular perturbations."""

    candidate_id: str
    description: str
    basis: RecipeBasis
    chi_asymptotic: float
    chi_coeffs: list[float]
    modes: list[AngularMode]


@dataclass(frozen=True)
class RayValidationRecord:
    candidate_id: str
    description: str
    chi_min: float
    chi_max: float
    negative_chi_fraction: float
    max_angular_contrast: float
    max_stretch: float
    min_stretch: float
    classification: str
    reason: str


@dataclass(frozen=True)
class RayValidationSummary:
    total: int
    by_classification: dict[str, int]
    seed: int
    timestamp: str
    output_csv: str | None = None
    output_json: str | None = None


def legendre_p(ell: int, mu: NDArray[np.float64]) -> NDArray[np.float64]:
    """Low-order Legendre polynomials for axisymmetric angular modes."""
    if ell == 0:
        return np.ones_like(mu)
    if ell == 1:
        return mu
    if ell == 2:
        return 0.5 * (3.0 * mu * mu - 1.0)
    if ell == 3:
        return 0.5 * (5.0 * mu**3 - 3.0 * mu)
    raise ValueError(f"Unsupported ell={ell}; supported: 0, 1, 2, 3")


def angular_bump(
    r: NDArray[np.float64],
    theta: float,
    mode: AngularMode,
) -> NDArray[np.float64]:
    mu = np.array(np.cos(theta), dtype=np.float64)
    radial = np.exp(-((r - mode.radial_center) ** 2) / (2.0 * mode.radial_width**2))
    return mode.amplitude * radial * legendre_p(mode.ell, mu)


def nonspherical_to_overrides(candidate: NonSphericalCandidate) -> dict[str, Any]:
    """Map a non-spherical candidate to RadialRecipe params overrides."""
    overrides: dict[str, Any] = {
        "recipe_num_bases": candidate.basis.num_bases,
        "recipe_basis_width": candidate.basis.basis_width,
        "recipe_basis_radius_max": candidate.basis.basis_radius_max,
        "recipe_chi_asymptotic": candidate.chi_asymptotic,
        "recipe_alpha_asymptotic": 1.0,
        "recipe_beta_asymptotic": 0.0,
        "recipe_K_asymptotic": 0.0,
        "recipe_phi_asymptotic": 0.0,
        "recipe_Pi_asymptotic": 0.0,
        "recipe_num_chi_angular_modes": len(candidate.modes),
    }
    for n in range(candidate.basis.num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = candidate.chi_coeffs[n]
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = 0.0
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0
        overrides[f"recipe_alpha_coeff_{n}"] = 0.0
        overrides[f"recipe_beta_coeff_{n}"] = 0.0
    for idx, mode in enumerate(candidate.modes):
        overrides[f"recipe_chi_mode_ell_{idx}"] = mode.ell
        overrides[f"recipe_chi_mode_amp_{idx}"] = mode.amplitude
        overrides[f"recipe_chi_mode_rc_{idx}"] = mode.radial_center
        overrides[f"recipe_chi_mode_rw_{idx}"] = mode.radial_width
    return overrides


def evaluate_chi_ray(
    candidate: NonSphericalCandidate,
    r: NDArray[np.float64],
    theta: float,
) -> NDArray[np.float64]:
    """Evaluate chi(r, theta) along one radial ray."""
    chi = candidate.basis.evaluate(r, candidate.chi_asymptotic, candidate.chi_coeffs)
    for mode in candidate.modes:
        chi = chi + angular_bump(r, theta, mode)
    return chi


def evaluate_ray_bundle(
    candidate: NonSphericalCandidate,
    *,
    n_points: int = 2048,
    thetas: tuple[float, ...] = (0.0, np.pi / 4.0, np.pi / 2.0, 3.0 * np.pi / 4.0, np.pi),
) -> tuple[NDArray[np.float64], NDArray[np.float64], tuple[float, ...]]:
    r_max = 2.0 * candidate.basis.basis_radius_max
    r = np.linspace(r_max / n_points, r_max, n_points)
    bundle = np.vstack([evaluate_chi_ray(candidate, r, theta) for theta in thetas])
    return r, bundle, thetas


def classify_ray_bundle(
    *,
    chi_min: float,
    negative_chi_fraction: float,
    max_angular_contrast: float,
    contrast_threshold: float,
) -> tuple[str, str]:
    if chi_min <= 0.0 or negative_chi_fraction > 0.0:
        return (
            REJECTED_NEGATIVE_CHI,
            f"chi_min={chi_min:.4g}; neg_frac={negative_chi_fraction:.2%}",
        )
    if max_angular_contrast > contrast_threshold:
        return (
            REJECTED_HIGH_ANGULAR_CONTRAST,
            f"max angular contrast={max_angular_contrast:.4g}>{contrast_threshold}",
        )
    return ACCEPTED_RAY_SANE, "chi positive on all sampled rays"


def validate_candidate_rays(
    candidate: NonSphericalCandidate,
    *,
    contrast_threshold: float = 0.5,
) -> RayValidationRecord:
    """Run one-step metric validation along sampled 1D rays."""
    _, bundle, _ = evaluate_ray_bundle(candidate)
    chi_min = float(np.min(bundle))
    chi_max = float(np.max(bundle))
    negative_chi_fraction = float(np.mean(bundle <= 0.0))
    per_radius_contrast = np.max(bundle, axis=0) - np.min(bundle, axis=0)
    max_angular_contrast = float(np.max(per_radius_contrast))
    chi_safe = np.clip(bundle, 1.0e-12, None)
    stretch = 1.0 / np.sqrt(chi_safe)
    min_stretch = float(np.min(stretch))
    max_stretch = float(np.max(stretch))
    classification, reason = classify_ray_bundle(
        chi_min=chi_min,
        negative_chi_fraction=negative_chi_fraction,
        max_angular_contrast=max_angular_contrast,
        contrast_threshold=contrast_threshold,
    )
    return RayValidationRecord(
        candidate_id=candidate.candidate_id,
        description=candidate.description,
        chi_min=chi_min,
        chi_max=chi_max,
        negative_chi_fraction=negative_chi_fraction,
        max_angular_contrast=max_angular_contrast,
        max_stretch=max_stretch,
        min_stretch=min_stretch,
        classification=classification,
        reason=reason,
    )


def generate_nonspherical_candidates(seed: int = 42) -> list[NonSphericalCandidate]:
    """Deterministic non-spherical proposals for first-step validation."""
    rng = np.random.default_rng(seed)
    basis = RecipeBasis(num_bases=6, basis_width=1.0, basis_radius_max=8.0)
    base = [-0.12, -0.04, 0.0, 0.0, 0.0, 0.0]

    candidates = [
        NonSphericalCandidate(
            candidate_id="dipole_lopsided_000",
            description="mild l=1 lopsided radial compression",
            basis=basis,
            chi_asymptotic=1.0,
            chi_coeffs=base,
            modes=[AngularMode(ell=1, amplitude=0.08, radial_center=3.0, radial_width=1.2)],
        ),
        NonSphericalCandidate(
            candidate_id="quadrupole_bubble_001",
            description="mild l=2 prolate/oblate bubble-wall deformation",
            basis=basis,
            chi_asymptotic=1.0,
            chi_coeffs=base,
            modes=[AngularMode(ell=2, amplitude=-0.10, radial_center=4.0, radial_width=1.0)],
        ),
        NonSphericalCandidate(
            candidate_id="mixed_modes_002",
            description="combined l=1 and l=2 deformation",
            basis=basis,
            chi_asymptotic=1.0,
            chi_coeffs=base,
            modes=[
                AngularMode(ell=1, amplitude=0.06, radial_center=2.5, radial_width=1.0),
                AngularMode(ell=2, amplitude=-0.07, radial_center=5.0, radial_width=1.4),
            ],
        ),
        NonSphericalCandidate(
            candidate_id="strong_quadrupole_bad_003",
            description="bad control: angular mode drives chi negative on some rays",
            basis=basis,
            chi_asymptotic=1.0,
            chi_coeffs=[-0.55, -0.20, 0.0, 0.0, 0.0, 0.0],
            modes=[AngularMode(ell=2, amplitude=-0.70, radial_center=0.0, radial_width=1.5)],
        ),
    ]

    # Add a few randomized but bounded angular candidates.
    for idx in range(4, 8):
        ell = int(rng.choice([1, 2, 3]))
        amplitude = float(rng.uniform(-0.12, 0.12))
        center = float(rng.uniform(2.0, 6.0))
        width = float(rng.uniform(0.8, 1.8))
        candidates.append(
            NonSphericalCandidate(
                candidate_id=f"random_angular_{idx:03d}",
                description=f"random axisymmetric ell={ell} angular deformation",
                basis=basis,
                chi_asymptotic=1.0,
                chi_coeffs=base,
                modes=[AngularMode(ell=ell, amplitude=amplitude, radial_center=center, radial_width=width)],
            )
        )

    return candidates


def summarize(records: list[RayValidationRecord], *, seed: int, output_csv: Path | None, output_json: Path | None) -> RayValidationSummary:
    counts: dict[str, int] = {}
    for record in records:
        counts[record.classification] = counts.get(record.classification, 0) + 1
    return RayValidationSummary(
        total=len(records),
        by_classification=counts,
        seed=seed,
        timestamp=datetime.now(timezone.utc).isoformat(),
        output_csv=str(output_csv) if output_csv else None,
        output_json=str(output_json) if output_json else None,
    )


def write_csv(records: list[RayValidationRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(records[0]).keys()))
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))


def write_json(summary: RayValidationSummary, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")


def run_nonspherical_validation(
    *,
    seed: int = 42,
    output_dir: Path | None = None,
    contrast_threshold: float = 0.5,
) -> RayValidationSummary:
    candidates = generate_nonspherical_candidates(seed=seed)
    records = [
        validate_candidate_rays(candidate, contrast_threshold=contrast_threshold)
        for candidate in candidates
    ]

    csv_path = json_path = None
    if output_dir is not None:
        output_dir = output_dir.expanduser().resolve()
        csv_path = output_dir / "nonspherical_ray_validation.csv"
        json_path = output_dir / "nonspherical_ray_validation_summary.json"
        write_csv(records, csv_path)

    summary = summarize(records, seed=seed, output_csv=csv_path, output_json=json_path)
    if json_path is not None:
        write_json(summary, json_path)

    print("=" * 60)
    print("NON-SPHERICAL 1D RAY VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Total candidates: {summary.total}")
    for label, count in sorted(summary.by_classification.items()):
        print(f"  {label:36s} {count:4d}")
    if csv_path:
        print(f"CSV: {csv_path}")
    if json_path:
        print(f"JSON: {json_path}")
    print("=" * 60)
    return summary

