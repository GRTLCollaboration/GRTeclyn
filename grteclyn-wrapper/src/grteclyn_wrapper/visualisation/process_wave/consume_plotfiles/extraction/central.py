from __future__ import annotations

import math
import os

import numpy as np


CENTRAL_TIMESERIES_HEADER = (
    "# time  rho_req  lapse  scalar_activity  phi_re  phi_im"
)

_FIELD_CANDIDATES = {
    "rho_req": ("rho_req",),
    "lapse": ("lapse",),
    "phi_re": ("phi", "phi_lump_sum"),
    "phi_im": ("phi_lump0",),
    "pi_re": ("Pi", "Pi_lump_sum"),
    "pi_im": ("Pi_lump0",),
    "scalar_activity": ("scalar_activity",),
}


def _read_scalar(field_obj) -> float | None:
    try:
        value = float(np.asarray(field_obj, dtype=float).reshape(-1)[0])
    except Exception:  # noqa: BLE001
        return None
    if math.isfinite(value):
        return value
    return None


def _field_at_point(ds, point, name: str) -> float | None:
    for candidate in _FIELD_CANDIDATES.get(name, (name,)):
        key = ("boxlib", candidate)
        if key not in ds.field_list:
            continue
        value = _read_scalar(point[key])
        if value is not None:
            return value
    return None


def _scalar_activity(ds, point) -> float | None:
    direct = _field_at_point(ds, point, "scalar_activity")
    if direct is not None:
        return direct
    parts: list[float] = []
    for key in ("phi_re", "phi_im", "pi_re", "pi_im"):
        value = _field_at_point(ds, point, key)
        if value is not None:
            parts.append(value)
    if not parts:
        return None
    return float(math.sqrt(sum(value * value for value in parts)))


def _extract_central_timeseries_line(
    p: str,
    *,
    t: float,
    center,
    verbose: bool = False,
) -> str | None:
    import yt

    try:
        ds = yt.load(p)
        point = ds.point([float(center[0]), float(center[1]), float(center[2])])
        rho_req = _field_at_point(ds, point, "rho_req")
        lapse = _field_at_point(ds, point, "lapse")
        activity = _scalar_activity(ds, point)
        phi_re = _field_at_point(ds, point, "phi_re")
        phi_im = _field_at_point(ds, point, "phi_im")
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: central-timeseries failed for {os.path.basename(p)}: {exc}"
            )
        return None

    def _fmt(value: float | None) -> str:
        if value is None or not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"

    if rho_req is None and lapse is None and activity is None:
        return None

    return (
        f"{t:.16e}  {_fmt(rho_req)}  {_fmt(lapse)}  {_fmt(activity)}  "
        f"{_fmt(phi_re)}  {_fmt(phi_im)}"
    )
