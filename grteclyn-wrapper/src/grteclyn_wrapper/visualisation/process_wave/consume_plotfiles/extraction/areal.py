from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np


def _extract_areal_radius_min(
    ds,
    center: Sequence[float] = (0.0, 0.0, 0.0),
    chi_floor: float = 1.0e-8,
) -> Tuple[float, float]:
    """Extract the minimum areal radius R_areal = r / sqrt(chi) along the x-axis.

    Returns (R_areal_min, r_at_min).
    """
    right = float(ds.domain_right_edge[0])
    c = np.asarray(center, dtype=float)
    ray = ds.ray(c, [right, c[1], c[2]])

    r_arr = np.asarray(ray[("index", "x")], dtype=float) - c[0]
    chi_arr = np.asarray(ray[("boxlib", "chi")], dtype=float)

    order = np.argsort(r_arr)
    r_arr = r_arr[order]
    chi_arr = chi_arr[order]

    chi_arr = np.maximum(chi_arr, chi_floor)
    R_areal = r_arr / np.sqrt(chi_arr)

    skip = r_arr > 1.0e-12
    if not np.any(skip):
        return (0.0, 0.0)
    R_areal_valid = R_areal[skip]
    r_valid = r_arr[skip]
    i_min = np.argmin(R_areal_valid)
    return (float(R_areal_valid[i_min]), float(r_valid[i_min]))
