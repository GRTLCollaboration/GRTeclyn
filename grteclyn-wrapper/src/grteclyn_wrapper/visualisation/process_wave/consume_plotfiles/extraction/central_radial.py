from __future__ import annotations

import math
import os
from typing import Sequence

import numpy as np

from .central import _default_ball_radius, _resolve_field_name, _scalar_activity_from_values


CENTRAL_RADIAL_PROFILE_HEADER = "# r  rho_req  lapse  scalar_activity"


def _log_spaced_radii(r_min: float, r_max: float, n_bins: int = 12) -> list[float]:
    r_min = max(float(r_min), 0.0)
    r_max = max(float(r_max), r_min + 1.0e-6)
    if n_bins <= 1:
        return [r_min, r_max]
    edges = np.geomspace(max(r_min, 1.0e-6), r_max, num=n_bins)
    return [float(r) for r in edges]


def _shell_means(
    ds,
    center: Sequence[float],
    radius: float,
    n_points: int,
) -> tuple[float | None, float | None, float | None]:
    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center_arr = np.asarray(center, dtype=float)
    theta = np.linspace(0.0, math.pi, n_points)
    phi = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    sin_t = np.sin(theta_grid)
    sx = float(radius) * sin_t * np.cos(phi_grid) + center_arr[0]
    sy = float(radius) * sin_t * np.sin(phi_grid) + center_arr[1]
    sz = float(radius) * np.cos(theta_grid) + center_arr[2]
    in_domain = (
        (sx >= left[0])
        & (sx <= right[0])
        & (sy >= left[1])
        & (sy <= right[1])
        & (sz >= left[2])
        & (sz <= right[2])
    )
    idxs = np.where(in_domain)[0]
    if idxs.size < 4:
        return None, None, None

    pts = np.column_stack((sx.ravel()[idxs], sy.ravel()[idxs], sz.ravel()[idxs]))
    field_map = {
        "rho_req": _resolve_field_name(ds, "rho_req"),
        "lapse": _resolve_field_name(ds, "lapse"),
        "scalar_activity": _resolve_field_name(ds, "scalar_activity"),
        "phi_re": _resolve_field_name(ds, "phi_re"),
        "phi_im": _resolve_field_name(ds, "phi_im"),
        "pi_re": _resolve_field_name(ds, "pi_re"),
        "pi_im": _resolve_field_name(ds, "pi_im"),
    }
    yt_fields = [("boxlib", name) for name in field_map.values() if name]
    if not yt_fields:
        return None, None, None
    vals = ds.find_field_values_at_points(yt_fields, pts)
    names = [name for name in field_map.values() if name]
    values: dict[str, float | None] = {}
    for name, samples in zip(names, vals):
        arr = np.asarray(samples, dtype=float).reshape(-1)
        arr = arr[np.isfinite(arr)]
        values[name] = float(np.mean(arr)) if arr.size else None

    rho = values.get(field_map["rho_req"]) if field_map["rho_req"] else None
    lapse = values.get(field_map["lapse"]) if field_map["lapse"] else None
    activity = values.get(field_map["scalar_activity"]) if field_map["scalar_activity"] else None
    if activity is None:
        activity = _scalar_activity_from_values(
            {
                "phi_re": values.get(field_map["phi_re"]) if field_map["phi_re"] else None,
                "phi_im": values.get(field_map["phi_im"]) if field_map["phi_im"] else None,
                "pi_re": values.get(field_map["pi_re"]) if field_map["pi_re"] else None,
                "pi_im": values.get(field_map["pi_im"]) if field_map["pi_im"] else None,
            }
        )
    return rho, lapse, activity


def _extract_central_radial_profile_block(
    p: str,
    *,
    t: float,
    center: Sequence[float],
    r_max: float = 6.0,
    r_min: float | None = None,
    n_bins: int = 12,
    n_points: int = 16,
    verbose: bool = False,
) -> str | None:
    import yt

    try:
        ds = yt.load(p)
        r_min_eff = r_min if r_min is not None else _default_ball_radius(ds)
        radii = _log_spaced_radii(r_min_eff, r_max, n_bins=n_bins)
        lines = [f"# time={t:.16e}", CENTRAL_RADIAL_PROFILE_HEADER]
        wrote = False
        for radius in radii:
            rho, lapse, activity = _shell_means(ds, center, radius, n_points)
            if rho is None and lapse is None and activity is None:
                continue
            wrote = True
            lines.append(
                f"{radius:.16e}  "
                f"{(rho if rho is not None else float('nan')):.16e}  "
                f"{(lapse if lapse is not None else float('nan')):.16e}  "
                f"{(activity if activity is not None else float('nan')):.16e}"
            )
        if not wrote:
            return None
        return "\n".join(lines) + "\n"
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: central-radial-profile failed for {os.path.basename(p)}: {exc}"
            )
        return None
