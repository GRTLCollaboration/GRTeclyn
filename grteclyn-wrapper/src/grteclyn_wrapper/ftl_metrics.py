"""1D FTL shortcut metrics from RadialRecipe initial-data profiles.

Computes warp null-ray travel time, portal proper-distance shortcut,
wormhole throat pinch, and anti-flatness from t=0 recipe coefficients
without running a 3D geodesic tracer.
"""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from numpy.typing import NDArray

from .constrained_recipe import RecipeBasis

_ASSIGNMENT_RE = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*=\s*(.+?)(?:\s+#.*)?$")


@dataclass(frozen=True)
class FtlMetrics:
    f_null: float
    f_portal: float
    f_throat: float
    s_nonflat: float
    f_shortcut: float
    t_curved: float | None
    t_flat: float | None
    d_proper: float | None
    d_coordinate: float | None
    r_areal_min: float | None
    r_asymp: float | None
    path_valid: bool
    asymptotic_ok: bool
    notes: tuple[str, ...]


def _parse_params_file(path: Path) -> dict[str, Any]:
    overrides: dict[str, Any] = {}
    if not path.exists():
        return overrides
    for line in path.read_text(encoding="utf-8").splitlines():
        match = _ASSIGNMENT_RE.match(line)
        if not match:
            continue
        key, raw = match.group(1), match.group(2).strip()
        if raw.startswith('"') and raw.endswith('"'):
            overrides[key] = raw.strip('"')
            continue
        if raw.lower() in {"true", "false"}:
            overrides[key] = raw.lower() == "true"
            continue
        try:
            if "." in raw or "e" in raw.lower():
                overrides[key] = float(raw)
            else:
                overrides[key] = int(raw)
        except ValueError:
            overrides[key] = raw
    return overrides


def load_overrides_from_episode(episode_dir: Path) -> dict[str, Any] | None:
    episode_dir = episode_dir.expanduser().resolve()
    meta_path = episode_dir / "metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        overrides = meta.get("overrides")
        if isinstance(overrides, dict) and overrides:
            return dict(overrides)
    params_path = episode_dir / "params.txt"
    if params_path.exists():
        parsed = _parse_params_file(params_path)
        return parsed if parsed else None
    return None


def _coeff_list(overrides: Mapping[str, Any], prefix: str, num_bases: int) -> list[float]:
    return [float(overrides.get(f"{prefix}_{n}", 0.0)) for n in range(num_bases)]


def _axis_profiles(
    overrides: Mapping[str, Any],
    *,
    L: float,
    n_points: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    num_bases = int(overrides.get("recipe_num_bases", 4))
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=float(overrides.get("recipe_basis_radius_max", 8.0)),
    )

    chi_coeffs = _coeff_list(overrides, "recipe_chi_coeff", num_bases)
    alpha_coeffs = _coeff_list(overrides, "recipe_alpha_coeff", num_bases)
    beta_coeffs = _coeff_list(overrides, "recipe_beta_coeff", num_bases)

    x = np.linspace(-L, L, n_points)
    r = np.abs(x)

    chi_asym = float(overrides.get("recipe_chi_asymptotic", 1.0))
    alpha_asym = float(overrides.get("recipe_alpha_asymptotic", 1.0))
    beta_asym = float(overrides.get("recipe_beta_asymptotic", 0.0))

    chi = basis.evaluate(r, chi_asym, chi_coeffs)
    alpha = basis.evaluate(r, alpha_asym, alpha_coeffs)
    beta_x = basis.evaluate(r, beta_asym, beta_coeffs)

    chi = np.clip(chi, 1.0e-10, None)
    alpha = np.clip(alpha, 1.0e-10, None)
    return x, r, chi, alpha, beta_x


def compute_ftl_metrics(
    overrides: Mapping[str, Any],
    *,
    L: float | None = None,
    n_points: int = 1024,
    velocity_floor: float = 1.0e-5,
    asymptotic_delta: float = 0.15,
    nonflat_tau: float = 0.08,
) -> FtlMetrics:
    """Compute FTL shortcut metrics from recipe overrides at t=0."""
    basis_radius_max = float(overrides.get("recipe_basis_radius_max", 8.0))
    integration_L = L if L is not None else basis_radius_max
    num_bases = int(overrides.get("recipe_num_bases", 4))
    chi_asym = float(overrides.get("recipe_chi_asymptotic", 1.0))
    chi_coeffs = _coeff_list(overrides, "recipe_chi_coeff", num_bases)
    basis = RecipeBasis(
        num_bases=num_bases,
        basis_width=float(overrides.get("recipe_basis_width", 1.0)),
        basis_radius_max=basis_radius_max,
    )

    x, r, chi, alpha, beta_x = _axis_profiles(overrides, L=integration_L, n_points=n_points)

    notes: list[str] = []

    chi_dev = float(np.max(np.abs(chi - 1.0)))
    beta_max = float(np.max(np.abs(beta_x)))
    s_nonflat = min(1.0, (chi_dev + beta_max) / max(nonflat_tau, 1.0e-12))

    chi_edge = float(chi[-1])
    alpha_edge = float(alpha[-1])
    beta_edge = float(beta_x[-1])
    asymptotic_ok = (
        abs(chi_edge - 1.0) <= asymptotic_delta
        and abs(alpha_edge - 1.0) <= asymptotic_delta
        and abs(beta_edge) <= asymptotic_delta
    )
    if not asymptotic_ok:
        notes.append("asymptotic flatness check failed at integration boundary")

    coordinate_velocity = alpha * np.sqrt(chi) - beta_x
    path_valid = bool(np.all(coordinate_velocity >= velocity_floor))
    if not path_valid:
        notes.append("null-ray path invalid: coordinate velocity below floor")

    t_flat = 2.0 * integration_L
    t_curved: float | None = None
    f_null = 0.0
    if path_valid and asymptotic_ok:
        integrand = 1.0 / np.clip(coordinate_velocity, velocity_floor, None)
        t_curved = float(np.trapezoid(integrand, x))
        f_null = max(0.0, (t_flat - t_curved) / t_flat)

    sqrt_chi = np.sqrt(chi)
    d_proper = float(np.trapezoid(1.0 / sqrt_chi, x))
    d_coordinate = t_flat
    f_portal = max(0.0, (d_coordinate - d_proper) / d_coordinate)

    r_pos = np.linspace(max(integration_L / n_points, 1.0e-6), integration_L, max(n_points // 2, 64))
    chi_pos = basis.evaluate(r_pos, chi_asym, chi_coeffs)
    chi_pos = np.clip(chi_pos, 1.0e-10, None)
    r_areal_pos = r_pos / np.sqrt(chi_pos)
    r_asymp = integration_L / max(float(np.sqrt(chi_pos[-1])), 1.0e-8)
    r_areal_min = float(np.min(r_areal_pos))
    f_throat = 0.0
    if r_asymp > 1.0e-8 and len(r_areal_pos) >= 5:
        for idx in range(1, len(r_areal_pos) - 1):
            if (
                r_areal_pos[idx] <= r_areal_pos[idx - 1]
                and r_areal_pos[idx] <= r_areal_pos[idx + 1]
                and r_areal_pos[idx] < 0.95 * r_asymp
            ):
                pinch = max(0.0, 1.0 - float(r_areal_pos[idx]) / r_asymp)
                f_throat = max(f_throat, pinch)

    raw_shortcut = max(f_null, f_portal, f_throat)
    f_shortcut = raw_shortcut * s_nonflat if asymptotic_ok else 0.0
    if raw_shortcut > 0.0 and s_nonflat < 1.0e-6:
        notes.append("shortcut suppressed: geometry is effectively flat")

    return FtlMetrics(
        f_null=f_null,
        f_portal=f_portal,
        f_throat=f_throat,
        s_nonflat=s_nonflat,
        f_shortcut=f_shortcut,
        t_curved=t_curved,
        t_flat=t_flat,
        d_proper=d_proper,
        d_coordinate=d_coordinate,
        r_areal_min=r_areal_min,
        r_asymp=r_asymp,
        path_valid=path_valid,
        asymptotic_ok=asymptotic_ok,
        notes=tuple(notes),
    )
