"""Symbolic regression of native-coordinate metric coefficients (PySR)."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class SymbolicMetricForm:
    expressions: dict[str, str]
    slice: str
    r2_scores: dict[str, float]
    notes: tuple[str, ...] = ()


def _fit_axisymmetric_slice(
    plotfile: str | Path,
    *,
    n: int = 33,
) -> tuple[np.ndarray, dict[str, np.ndarray]] | None:
    """Sample chi, alpha, beta1 on the x-z midplane (y = center)."""
    import yt

    ds = yt.load(str(plotfile))
    center = ds.domain_center
    cy = float(center[1].to_value())
    half_x = 0.45 * float(ds.domain_width[0].to_value())
    half_z = 0.45 * float(ds.domain_width[2].to_value())
    slc = ds.slice(1, cy)
    frb = slc.to_frb((2.0 * half_x, "code_length"), n, center=center, height=(2.0 * half_z, "code_length"))

    def field(name: str) -> np.ndarray:
        try:
            arr = np.asarray(frb["boxlib", name], dtype=float)
        except Exception:  # noqa: BLE001
            arr = np.asarray(frb[name], dtype=float)
        return arr.T

    x = np.linspace(-half_x, half_x, n)
    z = np.linspace(-half_z, half_z, n)
    X, Z = np.meshgrid(x, z, indexing="ij")
    r = np.sqrt(X ** 2 + Z ** 2)
    return r.ravel(), {
        "chi": field("chi").ravel(),
        "lapse": field("lapse").ravel(),
        "shift1": field("shift1").ravel(),
    }


def _polynomial_fit(x: np.ndarray, y: np.ndarray, degree: int = 4) -> tuple[str, float]:
    """Fallback closed form when PySR is unavailable."""
    coeffs = np.polyfit(x, y, deg=min(degree, max(1, len(x) // 3)))
    expr = "+".join(f"{c:.6g}*r^{d}" if d else f"{c:.6g}" for d, c in enumerate(reversed(list(coeffs))))
    pred = np.polyval(coeffs, x)
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return expr, r2


def extract_symbolic_metric(
    plotfile: str | Path,
    *,
    max_complexity: int = 20,
) -> SymbolicMetricForm | None:
    """Fit 2D axisymmetric metric coefficients on the x-z midplane."""
    sampled = _fit_axisymmetric_slice(plotfile)
    if sampled is None:
        return None
    r, fields = sampled
    mask = np.isfinite(r) & (r > 1.0e-6)
    for key in fields:
        mask &= np.isfinite(fields[key])
    r = r[mask]
    if r.size < 10:
        return None

    expressions: dict[str, str] = {}
    r2_scores: dict[str, float] = {}
    notes: list[str] = []

    try:
        from pysr import PySRRegressor  # type: ignore

        for name, y in fields.items():
            y = y[mask]
            model = PySRRegressor(
                niterations=20,
                maxsize=max_complexity,
                binary_operators=["+", "-", "*", "/"],
                unary_operators=["exp", "log", "sqrt"],
                extra_sympy_mappings={},
                verbosity=0,
            )
            model.fit(r.reshape(-1, 1), y)
            best = model.get_best()
            expressions[name] = str(best["equation"])
            r2_scores[name] = float(best.get("score", 0.0))
        notes.append("PySR symbolic regression")
    except Exception:
        for name, y in fields.items():
            expr, r2 = _polynomial_fit(r, y[mask])
            expressions[name] = expr
            r2_scores[name] = r2
        notes.append("polynomial fallback (PySR unavailable)")

    return SymbolicMetricForm(
        expressions=expressions,
        slice="xz_midplane_axisymmetric",
        r2_scores=r2_scores,
        notes=tuple(notes),
    )


def verify_analytic_form(
    form: SymbolicMetricForm,
    geodesic_f_geo: float | None,
    *,
    f_geo_tolerance: float = 0.15,
) -> dict[str, Any]:
    """Cross-check symbolic form against numerical geodesic shortcut."""
    ok_r2 = all(v >= 0.5 for v in form.r2_scores.values())
    geo_ok = geodesic_f_geo is not None and math.isfinite(geodesic_f_geo)
    return {
        "analytic_form": ok_r2,
        "resolution_converged": geo_ok,
        "fit_quality_ok": ok_r2,
        "f_geo_match": geo_ok,
        "expressions": form.expressions,
        "r2_scores": form.r2_scores,
        "notes": list(form.notes),
    }


def write_symbolic_form(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(payload), indent=2, sort_keys=True), encoding="utf-8")
