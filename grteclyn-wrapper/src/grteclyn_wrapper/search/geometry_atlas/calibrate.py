"""Probe calibration anchors for the geometry atlas.

Scores known shortcuts with the *same* frozen null-geodesic probe the atlas
uses, so search results can be interpreted against analytic Alcubierre and
(optionally) candidate-146 initial data.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from ...metrics.probes.ftl.geodesic import build_metric_3d_from_gridinit
from ...projection.warp_gridinit import write_alcubierre_gridinit
from .score import probe_half_length, score_metric_g


# Default paper-champion IVP (t=0 solved slice).  Frozen-peak f_geo in the
# article is from an evolved snapshot (~29%); this file is the available
# Stage-0 handoff for apples-to-apples atlas-probe calibration.
DEFAULT_CAND146_GRIDINIT = (
    Path(__file__).resolve().parents[5]
    / "runs"
    / "grtresna_promote"
    / "_cache"
    / "bcma_rm_L128_N256_t30_hq_eval000146"
    / "initial_data.gridinit"
)


@dataclass(frozen=True)
class CalibrationCase:
    name: str
    f_geo: float
    t_flat: float
    t_min: float | None
    h_quality_ok: bool
    n_reached: int
    n_rays: int
    half_length: float | None
    path: str | None = None
    notes: tuple[str, ...] = ()
    expected_note: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "f_geo": self.f_geo,
            "t_flat": self.t_flat,
            "t_min": self.t_min,
            "h_quality_ok": self.h_quality_ok,
            "n_reached": self.n_reached,
            "n_rays": self.n_rays,
            "half_length": self.half_length,
            "path": self.path,
            "notes": list(self.notes),
            "expected_note": self.expected_note,
        }


def _score_gridinit(
    path: Path,
    *,
    name: str,
    n_rays: int,
    half_length: float | None,
    expected_note: str | None = None,
    subsample_n: int | None = None,
    subsample_half_width: float | None = None,
) -> CalibrationCase:
    g, origin, spacing = build_metric_3d_from_gridinit(
        path, n=subsample_n, half_width=subsample_half_width
    )
    geo, _ff = score_metric_g(
        g,
        origin,
        spacing,
        n_rays=n_rays,
        compute_ff=False,
        half_length=half_length,
        directions=("x", "y", "z"),
    )
    return CalibrationCase(
        name=name,
        f_geo=float(geo.f_geo),
        t_flat=float(geo.t_flat),
        t_min=float(geo.t_min) if geo.t_min is not None else None,
        h_quality_ok=bool(geo.h_quality_ok),
        n_reached=int(geo.n_reached),
        n_rays=int(geo.n_rays),
        half_length=half_length,
        path=str(path),
        notes=tuple(geo.notes),
        expected_note=expected_note,
    )


def calibrate_atlas_probe(
    out_dir: Path,
    *,
    n_rays: int = 3,
    alc_n: int = 48,
    alc_L: float = 16.0,
    alc_velocity: float = 0.9,
    alc_radius: float = 4.0,
    alc_sigma: float = 1.5,
    cand146_path: Path | None = DEFAULT_CAND146_GRIDINIT,
    localise: bool = True,
) -> dict[str, Any]:
    """Run Alcubierre (+ optional cand-146) calibration and write a report."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # write_alcubierre_gridinit uses box [0, 2L]^3 with centre at (L,L,L).
    box_length = 2.0 * float(alc_L)
    support = float(alc_radius) + 3.0 / max(float(alc_sigma), 1.0e-6)
    half = (
        probe_half_length(support_radius=support, box_length=box_length)
        if localise
        else None
    )

    alc_path = out_dir / "analytic_alcubierre.gridinit"
    write_alcubierre_gridinit(
        alc_path,
        n=int(alc_n),
        L=float(alc_L),
        velocity=float(alc_velocity),
        bubble_radius=float(alc_radius),
        sigma=float(alc_sigma),
    )

    cases: list[CalibrationCase] = [
        _score_gridinit(
            alc_path,
            name="alcubierre_analytic",
            n_rays=n_rays,
            half_length=half,
            expected_note=(
                "Positive control: soft Alcubierre should yield f_geo ≳ 0.15 "
                "with localised endpoints (paper frozen peak ~0.29 on denser grids)."
            ),
        ),
        _score_gridinit(
            alc_path,
            name="alcubierre_analytic_fullbox",
            n_rays=n_rays,
            half_length=None,
            expected_note="Same metric with historical 5–95% box span (diluted).",
        ),
    ]

    if cand146_path is not None and Path(cand146_path).exists():
        # Subsample a centred window so the heavy HQ grid stays CPU-cheap.
        cases.append(
            _score_gridinit(
                Path(cand146_path),
                name="candidate_146_initial_data",
                n_rays=n_rays,
                half_length=probe_half_length(support_radius=16.0, box_length=128.0)
                if localise
                else None,
                expected_note=(
                    "t=0 IVP slice of cand.146 (not the article's evolved "
                    "frozen-peak ~29.45%).  Useful as a matter-sourced geometry "
                    "anchor under the atlas probe."
                ),
                subsample_n=64,
                subsample_half_width=32.0,
            )
        )
    else:
        cases.append(
            CalibrationCase(
                name="candidate_146_initial_data",
                f_geo=float("nan"),
                t_flat=float("nan"),
                t_min=None,
                h_quality_ok=False,
                n_reached=0,
                n_rays=n_rays,
                half_length=None,
                path=str(cand146_path) if cand146_path is not None else None,
                notes=("gridinit not found; skipped",),
                expected_note="Optional cand.146 calibration skipped.",
            )
        )

    report = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "localise": localise,
        "alcubierre": {
            "n": alc_n,
            "L_half": alc_L,
            "box_length": box_length,
            "velocity": alc_velocity,
            "bubble_radius": alc_radius,
            "sigma": alc_sigma,
            "probe_half_length": half,
        },
        "cases": [c.to_dict() for c in cases],
        "verdict": _verdict(cases),
    }
    (out_dir / "calibration_report.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def _verdict(cases: list[CalibrationCase]) -> dict[str, Any]:
    by_name = {c.name: c for c in cases}
    alc = by_name.get("alcubierre_analytic")
    alc_full = by_name.get("alcubierre_analytic_fullbox")
    ok = alc is not None and alc.h_quality_ok and alc.f_geo >= 0.10
    return {
        "alcubierre_localised_ok": bool(ok),
        "alcubierre_localised_f_geo": None if alc is None else alc.f_geo,
        "alcubierre_fullbox_f_geo": None if alc_full is None else alc_full.f_geo,
        "localisation_gain": (
            None
            if alc is None or alc_full is None
            else float(alc.f_geo - alc_full.f_geo)
        ),
        "message": (
            "Probe detects Alcubierre shortcut under localised endpoints."
            if ok
            else "Probe failed Alcubierre positive control — fix endpoints/probe before hunting."
        ),
    }


__all__ = [
    "DEFAULT_CAND146_GRIDINIT",
    "CalibrationCase",
    "calibrate_atlas_probe",
]
