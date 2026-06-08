"""Compare GRTresna-projected geometry against a geometry-first motif target."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..grtresna.motif_fit import estimate_momentum_source, read_fitted_matter_json
from ..initial_data.motif import GeometryMotif, read_motif_json
from .ftl_solved_geometry import build_xz_slice_from_gridinit, compute_solved_geometry_ftl
from .ftl_metrics import calculate_expansion_asymmetry
from ..grtresna.io import read_gridinit


POLARITY_TOLERANCE = 0.35
LOCALIZATION_TOLERANCE = 0.5
MOMENTUM_ALIGNMENT_MIN = 0.2
F_OP_RETENTION_MIN = 0.25


@dataclass(frozen=True)
class PreservationReport:
    """How well a solved .gridinit preserves the scout motif."""

    passed: bool
    retention_score: float
    f_op_target: float
    f_op_solved: float
    f_op_retention: float
    polarity_target: float
    polarity_solved: float
    polarity_retention: float
    beta_max_solved: float
    shift_alignment: float
    momentum_alignment: float
    support_localized: bool
    mechanism_descriptor: float
    notes: tuple[str, ...]


def _polarity_from_slice(
    alpha: np.ndarray,
    beta: np.ndarray,
    gamma: np.ndarray,
    spacing: tuple[float, float],
) -> float:
    n = alpha.shape[0]
    x = (np.arange(n) - n / 2) * spacing[0]
    chi = 1.0 / np.clip(gamma[:, :, 0, 0], 1.0e-10, None)
    dbeta = np.gradient(beta[:, :, 0], x, axis=0)
    dln_sqrt_g = -0.5 * np.gradient(np.log(chi), x, axis=0)
    theta = dbeta + dln_sqrt_g
    return calculate_expansion_asymmetry(x, theta)


def _retention_ratio(target: float, solved: float) -> float:
    if target <= 1.0e-8:
        return 1.0 if solved <= 1.0e-8 else 0.0
    return float(np.clip(solved / target, 0.0, 1.5))


def compare_motif_preservation(
    motif: GeometryMotif,
    gridinit_path: str | Path,
    *,
    fitted_matter_path: str | Path | None = None,
    ftl_L: float | None = None,
) -> PreservationReport:
    """Score whether GRTresna projection preserved the scout motif."""
    notes: list[str] = []
    solved = compute_solved_geometry_ftl(gridinit_path, L=ftl_L)
    if solved is None:
        return PreservationReport(
            passed=False,
            retention_score=0.0,
            f_op_target=max(motif.f_op or 0.0, motif.f_shortcut),
            f_op_solved=0.0,
            f_op_retention=0.0,
            polarity_target=motif.polarity,
            polarity_solved=0.0,
            polarity_retention=0.0,
            beta_max_solved=0.0,
            shift_alignment=0.0,
            momentum_alignment=0.0,
            support_localized=False,
            mechanism_descriptor=0.0,
            notes=("failed to read solved geometry",),
        )

    f_op_target = max(motif.f_op or 0.0, motif.f_shortcut, motif.f_null)
    f_op_solved = solved.operational.f_op
    f_op_retention = _retention_ratio(f_op_target, f_op_solved)

    grid = read_gridinit(gridinit_path)
    alpha, beta, gamma, spacing = build_xz_slice_from_gridinit(grid, L=ftl_L)
    polarity_solved = _polarity_from_slice(alpha, beta, gamma, spacing)
    polarity_retention = 1.0 - min(
        1.0,
        abs(polarity_solved - motif.polarity) / max(motif.polarity, POLARITY_TOLERANCE),
    )

    beta_mag = np.sqrt(beta[:, :, 0] ** 2 + beta[:, :, 1] ** 2)
    beta_max_solved = float(np.max(beta_mag))
    shift_alignment = 0.0
    if motif.momentum_target.credible and motif.momentum_target.direction[0] != 0.0:
        mid = beta.shape[1] // 2
        signed_beta = float(np.mean(beta[:, mid, 0]))
        shift_alignment = float(
            np.clip(
                signed_beta / max(beta_max_solved, 1.0e-8) * np.sign(motif.momentum_target.direction[0]),
                -1.0,
                1.0,
            )
        )

    momentum_alignment = 0.0
    if fitted_matter_path is not None and Path(fitted_matter_path).exists():
        fitted = read_fitted_matter_json(fitted_matter_path)
        source = np.asarray(estimate_momentum_source(fitted.lumps), dtype=float)
        target = np.asarray(motif.momentum_target.direction, dtype=float)
        source_norm = float(np.linalg.norm(source))
        target_norm = float(np.linalg.norm(target))
        if source_norm > 1.0e-8 and target_norm > 1.0e-8:
            momentum_alignment = float(
                np.clip(np.dot(source, target) / (source_norm * target_norm), -1.0, 1.0)
            )

    chi = 1.0 / np.clip(gamma[:, :, 0, 0], 1.0e-10, None)
    chi_dev = float(np.max(np.abs(chi - 1.0)))
    support_localized = chi_dev >= LOCALIZATION_TOLERANCE

    retention_score = float(
        np.mean(
            [
                f_op_retention,
                polarity_retention,
                max(0.0, shift_alignment),
                max(0.0, momentum_alignment),
            ]
        )
    )

    passed = (
        f_op_retention >= F_OP_RETENTION_MIN
        and polarity_retention >= 0.5
        and support_localized
    )
    if motif.momentum_target.credible and momentum_alignment < MOMENTUM_ALIGNMENT_MIN:
        notes.append("momentum motor did not survive projection")
        passed = False
    if f_op_retention < F_OP_RETENTION_MIN:
        notes.append("operational FTL motif eroded after GRTresna projection")

    return PreservationReport(
        passed=passed,
        retention_score=retention_score,
        f_op_target=f_op_target,
        f_op_solved=f_op_solved,
        f_op_retention=f_op_retention,
        polarity_target=motif.polarity,
        polarity_solved=polarity_solved,
        polarity_retention=polarity_retention,
        beta_max_solved=beta_max_solved,
        shift_alignment=shift_alignment,
        momentum_alignment=momentum_alignment,
        support_localized=support_localized,
        mechanism_descriptor=solved.mechanisms.mechanism_descriptor,
        notes=tuple(notes),
    )


def write_preservation_report(report: PreservationReport, path: str | Path) -> None:
    Path(path).write_text(
        json.dumps(asdict(report), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_preservation_report(path: str | Path) -> PreservationReport:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    payload["notes"] = tuple(payload.get("notes", ()))
    return PreservationReport(**payload)
