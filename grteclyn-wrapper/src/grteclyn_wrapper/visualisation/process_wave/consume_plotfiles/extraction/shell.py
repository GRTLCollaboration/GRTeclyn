from __future__ import annotations

from typing import Dict, Sequence, Tuple

import numpy as np


def _shell_stats_header(radii: Sequence[float], fields: Sequence[str]) -> str:
    cols = ["time"]
    for radius in radii:
        tag = f"R{radius:g}"
        for field in fields:
            cols.extend(
                [
                    f"{field}_mean_{tag}",
                    f"{field}_min_{tag}",
                    f"{field}_max_{tag}",
                ]
            )
    return "# " + "  ".join(cols)


def _extract_shell_field_stats(
    ds,
    radii: Sequence[float],
    n_points: int,
    center: Sequence[float],
    fields: Sequence[str],
) -> Dict[str, Tuple[float, float, float]]:
    """Sample mean/min/max of fields on spherical shells at given radii."""
    if not radii or not fields:
        return {}

    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    center = np.asarray(center, dtype=float)

    theta = np.linspace(0.0, np.pi, n_points)
    phi = np.linspace(0.0, 2.0 * np.pi, n_points, endpoint=False)
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")
    sinT = np.sin(THETA)
    X1 = sinT * np.cos(PHI)
    Y1 = sinT * np.sin(PHI)
    Z1 = np.cos(THETA)

    yt_fields = []
    for field in fields:
        key = ("boxlib", field)
        if key not in ds.field_list:
            raise RuntimeError(f"Plotfile missing field {field!r} for shell extraction.")
        yt_fields.append(key)

    out: Dict[str, Tuple[float, float, float]] = {}
    for radius in radii:
        try:
            sx = (float(radius) * X1).ravel() + center[0]
            sy = (float(radius) * Y1).ravel() + center[1]
            sz = (float(radius) * Z1).ravel() + center[2]
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
                continue

            pts = np.column_stack((sx[idxs], sy[idxs], sz[idxs]))
            vals = ds.find_field_values_at_points(yt_fields, pts)
            tag = f"R{radius:g}"
            for field, samples in zip(fields, vals):
                arr = np.asarray(samples, dtype=float).reshape(-1)
                arr = arr[np.isfinite(arr)]
                if arr.size < 4:
                    continue
                out[f"{field}_mean_{tag}"] = (
                    float(np.mean(arr)),
                    float(np.min(arr)),
                    float(np.max(arr)),
                )
        except Exception:
            continue
    if not out:
        raise RuntimeError("Shell extraction produced no valid samples at any radius.")
    return out


def _format_shell_stats_line(
    t: float,
    stats: Dict[str, Tuple[float, float, float]],
    radii: Sequence[float],
    fields: Sequence[str],
) -> str:
    parts = [f"{t:.16e}"]
    for radius in radii:
        tag = f"R{radius:g}"
        for field in fields:
            key = f"{field}_mean_{tag}"
            if key not in stats:
                parts.extend(["nan", "nan", "nan"])
                continue
            mean_v, min_v, max_v = stats[key]
            parts.extend([f"{mean_v:.16e}", f"{min_v:.16e}", f"{max_v:.16e}"])
    return "  ".join(parts)
