"""Batch validation harness for the constraint-satisfying metric guesser.

Generates synthetic candidate geometries (not known exact metrics), applies
``constrained_overrides`` once, runs ``preflight_check``, and classifies
each outcome.  No optimizer feedback loop is used.

Run via CLI::

    uv run python -m grteclyn_wrapper validate --output-dir ./validation_out

Or programmatically::

    from grteclyn_wrapper.initial_data.validate_guesser import run_validation
    summary = run_validation(output_dir=Path("validation_out"))
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, Literal

import numpy as np

from .constrained_recipe import (
    RecipeBasis,
    constrained_overrides,
    constrained_phi,
    fit_gaussian_basis,
)
from .preflight import PreflightResult, preflight_check

FieldType = Literal["normal", "phantom"]

# Classification labels (plan-defined).
ACCEPTED_NORMAL = "accepted_normal"
ACCEPTED_PHANTOM = "accepted_phantom"
REJECTED_NEGATIVE_CHI = "rejected_negative_chi"
REJECTED_UNSATISFIABLE_SOURCE = "rejected_unsatisfiable_source"
REJECTED_HIGH_CONSTRAINT = "rejected_high_constraint"
REJECTED_FIT_ERROR = "rejected_fit_error"
REJECTED_EXOTIC_ENERGY = "rejected_exotic_energy"

# Thresholds for classification (aligned with preflight defaults where applicable).
DEFAULT_FIT_ERROR_THRESHOLD = 0.5
DEFAULT_UNSATISFIABLE_SOURCE_THRESHOLD = 0.8
DEFAULT_EXOTIC_ENERGY_INTEGRAL_THRESHOLD = 1.0e-3


@dataclass(frozen=True)
class CandidateSpec:
    """One synthetic geometry proposal before guessing."""

    candidate_id: str
    family: str
    field_type: FieldType
    num_bases: int
    basis_width: float
    basis_radius_max: float
    chi_asymptotic: float
    chi_coeffs: list[float]
    K_asymptotic: float
    K_coeffs: list[float]
    phi_asymptotic: float
    phi_coeffs: list[float]
    Pi_asymptotic: float
    Pi_coeffs: list[float]
    apply_constrained: bool = True
    notes: str = ""


@dataclass
class ValidationRecord:
    """Full validation result for one candidate."""

    candidate_id: str
    family: str
    field_type: FieldType
    num_bases: int
    basis_width: float
    basis_radius_max: float
    apply_constrained: bool
    notes: str

    chi_min: float
    chi_max: float
    negative_chi_fraction: float

    phi_fit_residual: float
    integrand_negative_fraction: float
    rho_required_min: float
    rho_required_max: float
    integral_negative_rho: float

    raw_hamiltonian_l2: float
    raw_momentum_l2: float
    raw_preflight_passed: bool
    raw_preflight_reason: str

    constrained_hamiltonian_l2: float
    constrained_momentum_l2: float
    constrained_preflight_passed: bool
    constrained_preflight_reason: str

    hamiltonian_improvement_ratio: float | None
    classification: str
    classification_reason: str


@dataclass
class ValidationSummary:
    """Aggregate counts from a validation batch."""

    total: int
    by_classification: dict[str, int]
    by_family: dict[str, int]
    by_family_and_classification: dict[str, dict[str, int]]
    seed: int
    timestamp: str
    output_csv: str | None = None
    output_json: str | None = None


def alcubierre_shape(
    r: np.ndarray,
    *,
    radius: float,
    sigma: float,
) -> np.ndarray:
    """Alcubierre top-hat shape function (radial proxy only)."""
    normalizer = 2.0 * np.tanh(sigma * radius)
    return (
        np.tanh(sigma * (r + radius)) - np.tanh(sigma * (r - radius))
    ) / normalizer


def _shell_profile(
    r: np.ndarray,
    *,
    center: float,
    width: float,
    amplitude: float,
) -> np.ndarray:
    return amplitude * np.exp(-((r - center) ** 2) / (2.0 * width**2))


def _profile_to_coeffs(
    profile: np.ndarray,
    r: np.ndarray,
    basis: RecipeBasis,
    asymptotic: float,
) -> list[float]:
    coeffs, _ = fit_gaussian_basis(r, profile, basis, asymptotic=asymptotic)
    return [float(c) for c in coeffs]


def _blank_recipe_coeffs(num_bases: int) -> dict[str, Any]:
    overrides: dict[str, Any] = {}
    for n in range(num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = 0.0
        overrides[f"recipe_K_coeff_{n}"] = 0.0
        overrides[f"recipe_phi_coeff_{n}"] = 0.0
        overrides[f"recipe_Pi_coeff_{n}"] = 0.0
    return overrides


def lookup_candidate(candidate_id: str, *, seed: int = 42) -> CandidateSpec:
    """Return a deterministic validation candidate by ID."""
    for spec in generate_candidates(seed):
        if spec.candidate_id == candidate_id:
            return spec
    raise KeyError(
        f"Unknown candidate_id={candidate_id!r} (validation seed={seed})."
    )


def spec_to_overrides(spec: CandidateSpec) -> dict[str, Any]:
    """Convert a candidate spec to RadialRecipe params overrides."""
    overrides: dict[str, Any] = {
        "recipe_num_bases": spec.num_bases,
        "recipe_basis_width": spec.basis_width,
        "recipe_basis_radius_max": spec.basis_radius_max,
        "recipe_chi_asymptotic": spec.chi_asymptotic,
        "recipe_K_asymptotic": spec.K_asymptotic,
        "recipe_phi_asymptotic": spec.phi_asymptotic,
        "recipe_Pi_asymptotic": spec.Pi_asymptotic,
    }
    overrides.update(_blank_recipe_coeffs(spec.num_bases))
    for n in range(spec.num_bases):
        overrides[f"recipe_chi_coeff_{n}"] = spec.chi_coeffs[n]
        overrides[f"recipe_K_coeff_{n}"] = spec.K_coeffs[n]
        overrides[f"recipe_phi_coeff_{n}"] = spec.phi_coeffs[n]
        overrides[f"recipe_Pi_coeff_{n}"] = spec.Pi_coeffs[n]
    return overrides


def generate_candidates(seed: int = 42) -> list[CandidateSpec]:
    """Build a fixed deterministic set of synthetic candidate geometries."""
    rng = np.random.default_rng(seed)
    candidates: list[CandidateSpec] = []
    idx = 0

    def next_id(family: str) -> str:
        nonlocal idx
        cid = f"{family}_{idx:03d}"
        idx += 1
        return cid

    # --- Random smooth Gaussian profiles ---
    for i in range(12):
        n_bases = int(rng.integers(4, 9))
        width = float(rng.uniform(0.4, 2.0))
        rmax = float(rng.uniform(6.0, 10.0))
        chi_coeffs = (rng.standard_normal(n_bases) * 0.12).tolist()
        k_const = float(rng.uniform(-0.2, 0.2))
        field: FieldType = "phantom" if i % 3 == 0 else "normal"
        phi_rand = (rng.standard_normal(n_bases) * 0.08).tolist()
        k_rand = (rng.standard_normal(n_bases) * 0.05).tolist()
        pi_rand = (rng.standard_normal(n_bases) * 0.03).tolist()
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("random"),
                family="random",
                field_type=field,
                num_bases=n_bases,
                basis_width=width,
                basis_radius_max=rmax,
                chi_asymptotic=1.0,
                chi_coeffs=chi_coeffs,
                K_asymptotic=k_const,
                K_coeffs=k_rand,
                phi_asymptotic=0.0,
                phi_coeffs=phi_rand,
                Pi_asymptotic=0.0,
                Pi_coeffs=pi_rand,
            )
        )

    # --- Shell profiles ---
    fit_r = np.linspace(0.01, 20.0, 4096)
    for i, (center, width, amp) in enumerate(
        [(2.0, 0.5, -0.08), (4.0, 0.8, -0.12), (6.0, 1.0, 0.06)]
    ):
        basis = RecipeBasis(num_bases=10, basis_width=0.6, basis_radius_max=10.0)
        profile = 1.0 + amp * _shell_profile(fit_r, center=center, width=width, amplitude=1.0)
        chi_coeffs = _profile_to_coeffs(profile, fit_r, basis, asymptotic=1.0)
        field = "phantom" if amp > 0 else "normal"
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("shell"),
                family="shell",
                field_type=field,
                num_bases=10,
                basis_width=0.6,
                basis_radius_max=10.0,
                chi_asymptotic=1.0,
                chi_coeffs=chi_coeffs,
                K_asymptotic=0.0,
                K_coeffs=[0.0] * 10,
                phi_asymptotic=0.0,
                phi_coeffs=[0.0] * 10,
                Pi_asymptotic=0.0,
                Pi_coeffs=[0.0] * 10,
                notes=f"center={center}, width={width}, amp={amp}",
            )
        )

    # --- Bubble-wall (Alcubierre-inspired radial proxy) ---
    for i, (sigma, amp, field) in enumerate(
        [(1.5, 0.08, "normal"), (2.0, 0.10, "phantom"), (3.0, 0.15, "normal")]
    ):
        basis = RecipeBasis(num_bases=12, basis_width=0.6, basis_radius_max=10.0)
        shape = alcubierre_shape(fit_r, radius=4.0, sigma=sigma)
        wall = shape * (1.0 - shape)
        sign = 1.0 if field == "phantom" else -1.0
        profile = 1.0 + sign * amp * wall
        chi_coeffs = _profile_to_coeffs(profile, fit_r, basis, asymptotic=1.0)
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("bubble_wall"),
                family="bubble_wall",
                field_type=field,
                num_bases=12,
                basis_width=0.6,
                basis_radius_max=10.0,
                chi_asymptotic=1.0,
                chi_coeffs=chi_coeffs,
                K_asymptotic=0.0,
                K_coeffs=[0.0] * 12,
                phi_asymptotic=0.0,
                phi_coeffs=[0.0] * 12,
                Pi_asymptotic=0.0,
                Pi_coeffs=[0.0] * 12,
                notes=f"sigma={sigma}, amp={amp}",
            )
        )

    # --- Multi-shell profiles ---
    for i, amps in enumerate([(-0.06, 0.04), (-0.10, 0.08), (-0.05, -0.03)]):
        basis = RecipeBasis(num_bases=14, basis_width=0.5, basis_radius_max=12.0)
        profile = 1.0
        profile += amps[0] * _shell_profile(fit_r, center=3.0, width=0.6, amplitude=1.0)
        profile += amps[1] * _shell_profile(fit_r, center=7.0, width=1.0, amplitude=1.0)
        chi_coeffs = _profile_to_coeffs(profile, fit_r, basis, asymptotic=1.0)
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("multi_shell"),
                family="multi_shell",
                field_type="phantom" if sum(amps) > 0 else "normal",
                num_bases=14,
                basis_width=0.5,
                basis_radius_max=12.0,
                chi_asymptotic=1.0,
                chi_coeffs=chi_coeffs,
                K_asymptotic=0.0,
                K_coeffs=[0.0] * 14,
                phi_asymptotic=0.0,
                phi_coeffs=[0.0] * 14,
                Pi_asymptotic=0.0,
                Pi_coeffs=[0.0] * 14,
                notes=f"amps={amps}",
            )
        )

    # --- Near-pathological profiles ---
    pathological_cases = [
        ("narrow_basis", 4, 0.15, 4.0, [-0.25, -0.10, 0.05, 0.0]),
        ("steep_gradient", 8, 0.35, 8.0, [-0.45, 0.30, -0.20, 0.10, 0.0, 0.0, 0.0, 0.0]),
        ("depressed_asymptotic", 6, 0.8, 8.0, [-0.08, -0.04, 0.0, 0.0, 0.0, 0.0]),
    ]
    for name, n_b, w, rmax, chi_c in pathological_cases:
        chi_asym = 0.8 if name == "depressed_asymptotic" else 1.0
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("pathological"),
                family="pathological",
                field_type="phantom" if name == "depressed_asymptotic" else "normal",
                num_bases=n_b,
                basis_width=w,
                basis_radius_max=rmax,
                chi_asymptotic=chi_asym,
                chi_coeffs=chi_c + [0.0] * (n_b - len(chi_c)),
                K_asymptotic=0.0,
                K_coeffs=[0.0] * n_b,
                phi_asymptotic=0.0,
                phi_coeffs=[0.0] * n_b,
                Pi_asymptotic=0.0,
                Pi_coeffs=[0.0] * n_b,
                notes=name,
            )
        )

    # --- Bad controls (should be rejected) ---
    bad_controls = [
        ("negative_chi", [-1.2, 0.3, -0.2, 0.0], [0.0] * 4, [0.0] * 4, [0.0] * 4),
        ("spatial_K", [-0.15, -0.05, 0.0, 0.0], [0.4, -0.2, 0.1, 0.0], [0.0] * 4, [0.0] * 4),
        ("random_pi", [-0.10, -0.05, 0.0, 0.0], [0.0] * 4, [0.05, 0.0, 0.0, 0.0], [0.08, -0.04, 0.02, 0.0]),
    ]
    for name, chi_c, k_c, phi_c, pi_c in bad_controls:
        candidates.append(
            CandidateSpec(
                candidate_id=next_id("bad_control"),
                family="bad_control",
                field_type="normal",
                num_bases=4,
                basis_width=1.0,
                basis_radius_max=8.0,
                chi_asymptotic=1.0,
                chi_coeffs=chi_c,
                K_asymptotic=0.0,
                K_coeffs=k_c,
                phi_asymptotic=0.0,
                phi_coeffs=phi_c,
                Pi_asymptotic=0.0,
                Pi_coeffs=pi_c,
                apply_constrained=name != "negative_chi",
                notes=name,
            )
        )

    return candidates


def _chi_diagnostics(
    spec: CandidateSpec,
    *,
    n_points: int = 2048,
) -> tuple[float, float, float]:
    basis = RecipeBasis(
        num_bases=spec.num_bases,
        basis_width=spec.basis_width,
        basis_radius_max=spec.basis_radius_max,
    )
    r_max = 2.0 * spec.basis_radius_max
    r = np.linspace(r_max / n_points, r_max, n_points)
    chi = basis.evaluate(r, spec.chi_asymptotic, spec.chi_coeffs)
    neg_frac = float(np.mean(chi < 1.0e-10))
    return float(np.min(chi)), float(np.max(chi)), neg_frac


def _integral_negative_rho(rho: np.ndarray, r: np.ndarray) -> float:
    dr = r[1] - r[0]
    shell = 4.0 * math.pi * r * r * dr
    neg = np.where(rho < 0.0, -rho * shell, 0.0)
    return float(np.sum(neg))


def classify_candidate(
    *,
    spec: CandidateSpec,
    chi_min: float,
    negative_chi_fraction: float,
    phi_fit_residual: float,
    integrand_negative_fraction: float,
    integral_negative_rho: float,
    rho_required_min: float,
    constrained_preflight: PreflightResult,
    fit_error_threshold: float = DEFAULT_FIT_ERROR_THRESHOLD,
    unsatisfiable_threshold: float = DEFAULT_UNSATISFIABLE_SOURCE_THRESHOLD,
    exotic_energy_threshold: float = DEFAULT_EXOTIC_ENERGY_INTEGRAL_THRESHOLD,
) -> tuple[str, str]:
    """Return (classification_label, reason)."""
    if negative_chi_fraction > 0.1 or chi_min < 0.0:
        return REJECTED_NEGATIVE_CHI, f"neg_chi={negative_chi_fraction:.2%}, chi_min={chi_min:.4g}"

    if integrand_negative_fraction >= unsatisfiable_threshold:
        return (
            REJECTED_UNSATISFIABLE_SOURCE,
            f"integrand_negative={integrand_negative_fraction:.2%}",
        )

    if phi_fit_residual > fit_error_threshold:
        return REJECTED_FIT_ERROR, f"phi_fit_residual={phi_fit_residual:.4g}"

    if (
        spec.field_type == "normal"
        and (rho_required_min < 0.0 or integral_negative_rho > exotic_energy_threshold)
    ):
        return (
            REJECTED_EXOTIC_ENERGY,
            f"rho_min={rho_required_min:.4g}, int_neg_rho={integral_negative_rho:.4g}",
        )

    if not constrained_preflight.passed:
        return REJECTED_HIGH_CONSTRAINT, constrained_preflight.reason

    if spec.field_type == "phantom":
        return ACCEPTED_PHANTOM, "preflight ok; phantom field required"

    return ACCEPTED_NORMAL, "preflight ok; normal field"


def validate_candidate(
    spec: CandidateSpec,
    *,
    fit_error_threshold: float = DEFAULT_FIT_ERROR_THRESHOLD,
    unsatisfiable_threshold: float = DEFAULT_UNSATISFIABLE_SOURCE_THRESHOLD,
    exotic_energy_threshold: float = DEFAULT_EXOTIC_ENERGY_INTEGRAL_THRESHOLD,
) -> ValidationRecord:
    """Validate one candidate through raw and constrained pipelines."""
    raw_overrides = spec_to_overrides(spec)
    phantom = spec.field_type == "phantom"
    raw_pf = preflight_check(raw_overrides, phantom=phantom)

    chi_min, chi_max, neg_chi_frac = _chi_diagnostics(spec)

    basis = RecipeBasis(
        num_bases=spec.num_bases,
        basis_width=spec.basis_width,
        basis_radius_max=spec.basis_radius_max,
    )

    if spec.apply_constrained:
        constrained_result = constrained_phi(
            basis,
            spec.chi_asymptotic,
            spec.chi_coeffs,
            K_const=spec.K_asymptotic,
            phantom=phantom,
        )
        constrained_overrides_dict = spec_to_overrides(spec)
        constrained_overrides(constrained_overrides_dict, phantom=phantom)
        constrained_pf = preflight_check(constrained_overrides_dict, phantom=phantom)

        phi_fit = constrained_result.fit_residual
        integrand_neg = float(constrained_result.integrand_negative_mask.mean())
        rho_min = float(constrained_result.rho_required.min())
        rho_max = float(constrained_result.rho_required.max())
        int_neg_rho = _integral_negative_rho(
            constrained_result.rho_required,
            constrained_result.r,
        )
    else:
        constrained_pf = raw_pf
        phi_fit = float("nan")
        integrand_neg = float("nan")
        rho_min = float("nan")
        rho_max = float("nan")
        int_neg_rho = float("nan")

    if raw_pf.hamiltonian_l2 > 0:
        ham_ratio = constrained_pf.hamiltonian_l2 / raw_pf.hamiltonian_l2
    else:
        ham_ratio = None

    classification, class_reason = classify_candidate(
        spec=spec,
        chi_min=chi_min,
        negative_chi_fraction=neg_chi_frac,
        phi_fit_residual=phi_fit if math.isfinite(phi_fit) else 999.0,
        integrand_negative_fraction=integrand_neg if math.isfinite(integrand_neg) else 1.0,
        integral_negative_rho=int_neg_rho if math.isfinite(int_neg_rho) else 999.0,
        rho_required_min=rho_min if math.isfinite(rho_min) else -999.0,
        constrained_preflight=constrained_pf,
        fit_error_threshold=fit_error_threshold,
        unsatisfiable_threshold=unsatisfiable_threshold,
        exotic_energy_threshold=exotic_energy_threshold,
    )

    return ValidationRecord(
        candidate_id=spec.candidate_id,
        family=spec.family,
        field_type=spec.field_type,
        num_bases=spec.num_bases,
        basis_width=spec.basis_width,
        basis_radius_max=spec.basis_radius_max,
        apply_constrained=spec.apply_constrained,
        notes=spec.notes,
        chi_min=chi_min,
        chi_max=chi_max,
        negative_chi_fraction=neg_chi_frac,
        phi_fit_residual=phi_fit,
        integrand_negative_fraction=integrand_neg,
        rho_required_min=rho_min,
        rho_required_max=rho_max,
        integral_negative_rho=int_neg_rho,
        raw_hamiltonian_l2=raw_pf.hamiltonian_l2,
        raw_momentum_l2=raw_pf.momentum_l2,
        raw_preflight_passed=raw_pf.passed,
        raw_preflight_reason=raw_pf.reason,
        constrained_hamiltonian_l2=constrained_pf.hamiltonian_l2,
        constrained_momentum_l2=constrained_pf.momentum_l2,
        constrained_preflight_passed=constrained_pf.passed,
        constrained_preflight_reason=constrained_pf.reason,
        hamiltonian_improvement_ratio=ham_ratio,
        classification=classification,
        classification_reason=class_reason,
    )


CSV_FIELDS = [
    "candidate_id",
    "family",
    "field_type",
    "num_bases",
    "basis_width",
    "basis_radius_max",
    "apply_constrained",
    "notes",
    "chi_min",
    "chi_max",
    "negative_chi_fraction",
    "phi_fit_residual",
    "integrand_negative_fraction",
    "rho_required_min",
    "rho_required_max",
    "integral_negative_rho",
    "raw_hamiltonian_l2",
    "raw_momentum_l2",
    "raw_preflight_passed",
    "raw_preflight_reason",
    "constrained_hamiltonian_l2",
    "constrained_momentum_l2",
    "constrained_preflight_passed",
    "constrained_preflight_reason",
    "hamiltonian_improvement_ratio",
    "classification",
    "classification_reason",
]


def _record_to_row(record: ValidationRecord) -> dict[str, Any]:
    row = asdict(record)
    ratio = row.get("hamiltonian_improvement_ratio")
    if ratio is not None and not math.isfinite(ratio):
        row["hamiltonian_improvement_ratio"] = ""
    for key, value in row.items():
        if isinstance(value, float) and not math.isfinite(value):
            row[key] = ""
    return row


def write_validation_csv(records: list[ValidationRecord], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for record in records:
            writer.writerow(_record_to_row(record))


def summarize_validation(
    records: list[ValidationRecord],
    *,
    seed: int,
    output_csv: Path | None = None,
    output_json: Path | None = None,
) -> ValidationSummary:
    by_class: dict[str, int] = {}
    by_family: dict[str, int] = {}
    by_family_class: dict[str, dict[str, int]] = {}

    for rec in records:
        by_class[rec.classification] = by_class.get(rec.classification, 0) + 1
        by_family[rec.family] = by_family.get(rec.family, 0) + 1
        fam_map = by_family_class.setdefault(rec.family, {})
        fam_map[rec.classification] = fam_map.get(rec.classification, 0) + 1

    return ValidationSummary(
        total=len(records),
        by_classification=by_class,
        by_family=by_family,
        by_family_and_classification=by_family_class,
        seed=seed,
        timestamp=datetime.now(timezone.utc).isoformat(),
        output_csv=str(output_csv) if output_csv else None,
        output_json=str(output_json) if output_json else None,
    )


def write_validation_summary(summary: ValidationSummary, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(summary), indent=2), encoding="utf-8")


def print_validation_summary(summary: ValidationSummary) -> None:
    print("=" * 60)
    print("METRIC GUESSER VALIDATION SUMMARY")
    print("=" * 60)
    print(f"Total candidates: {summary.total}")
    print(f"Seed: {summary.seed}")
    print()
    print("By classification:")
    for label, count in sorted(summary.by_classification.items()):
        print(f"  {label:32s} {count:4d}")
    print()
    print("By family:")
    for family, count in sorted(summary.by_family.items()):
        print(f"  {family:32s} {count:4d}")
    print()
    if summary.output_csv:
        print(f"CSV: {summary.output_csv}")
    if summary.output_json:
        print(f"JSON: {summary.output_json}")
    print("=" * 60)


def run_validation(
    *,
    seed: int = 42,
    output_dir: Path | None = None,
    fit_error_threshold: float = DEFAULT_FIT_ERROR_THRESHOLD,
    unsatisfiable_threshold: float = DEFAULT_UNSATISFIABLE_SOURCE_THRESHOLD,
    exotic_energy_threshold: float = DEFAULT_EXOTIC_ENERGY_INTEGRAL_THRESHOLD,
) -> ValidationSummary:
    """Run the full validation batch and optionally write outputs."""
    candidates = generate_candidates(seed=seed)
    records = [
        validate_candidate(
            spec,
            fit_error_threshold=fit_error_threshold,
            unsatisfiable_threshold=unsatisfiable_threshold,
            exotic_energy_threshold=exotic_energy_threshold,
        )
        for spec in candidates
    ]

    csv_path = json_path = None
    if output_dir is not None:
        output_dir = output_dir.expanduser().resolve()
        csv_path = output_dir / "guesser_validation.csv"
        json_path = output_dir / "guesser_validation_summary.json"
        write_validation_csv(records, csv_path)

    summary = summarize_validation(
        records,
        seed=seed,
        output_csv=csv_path,
        output_json=json_path,
    )
    if json_path is not None:
        write_validation_summary(summary, json_path)

    print_validation_summary(summary)
    return summary
