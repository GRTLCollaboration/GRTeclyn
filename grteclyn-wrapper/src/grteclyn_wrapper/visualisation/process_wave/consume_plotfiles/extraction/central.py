from __future__ import annotations

import math
import os

import numpy as np


CENTRAL_TIMESERIES_HEADER = (
    "# time  rho_req  lapse  scalar_activity  phi_re  phi_im  "
    "noether_charge  phase_coherence  ham_abs  mom_abs  chi  K  weyl4"
)

_FIELD_CANDIDATES = {
    "rho_req": ("rho_req",),
    "lapse": ("lapse",),
    "phi_re": ("phi", "phi_lump_sum"),
    "phi_im": ("phi_lump0",),
    "pi_re": ("Pi", "Pi_lump_sum"),
    "pi_im": ("Pi_lump0",),
    "scalar_activity": ("scalar_activity",),
    "ham_abs": ("Ham_abs_terms",),
    "mom_abs": ("Mom_abs_terms",),
    # Geometric (spacetime) fields for the spacetime-splash diagnostics:
    #   chi   -> conformal factor (drops as curvature concentrates at center)
    #   K     -> trace of extrinsic curvature (the gravitational "crunch" rate)
    #   weyl4 -> Re(Psi4) Weyl scalar = gravitational-wave content at the center
    "chi": ("chi",),
    "K": ("K",),
    "weyl4": ("Weyl4_Re", "Weyl4", "Weyl4_re"),
    "A11": ("A11",),
    "A12": ("A12",),
    "A22": ("A22",),
}


def _gw_wave_magnitude(a11: float, a12: float, a22: float) -> float:
    """GW strain proxy |h| from conformal traceless A_ij (see visualisation README §8)."""
    h_plus = a11 - a22
    h_cross = 2.0 * a12
    return float(math.sqrt(h_plus * h_plus + h_cross * h_cross))


def _gw_wave_from_aij(ds, point) -> float | None:
    a11 = _field_at_point(ds, point, "A11")
    a12 = _field_at_point(ds, point, "A12")
    a22 = _field_at_point(ds, point, "A22")
    if a11 is None or a12 is None or a22 is None:
        return None
    return _gw_wave_magnitude(a11, a12, a22)


def _resolve_gw_wave_signal(ds, point) -> float | None:
    """Prefer Weyl4_Re; fall back to A_ij GW proxy for splash/RadialRecipe runs."""
    weyl = _field_at_point(ds, point, "weyl4")
    if weyl is not None and math.isfinite(weyl):
        return weyl
    return _gw_wave_from_aij(ds, point)


def _read_scalar(field_obj) -> float | None:
    try:
        value = float(np.asarray(field_obj, dtype=float).reshape(-1)[0])
    except Exception:  # noqa: BLE001
        return None
    if math.isfinite(value):
        return value
    return None


def _resolve_field_name(ds, name: str) -> str | None:
    for candidate in _FIELD_CANDIDATES.get(name, (name,)):
        if ("boxlib", candidate) in ds.field_list:
            return candidate
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


def _scalar_activity_from_values(values: dict[str, float | None]) -> float | None:
    direct = values.get("scalar_activity")
    if direct is not None:
        return direct
    parts: list[float] = []
    for key in ("phi_re", "phi_im", "pi_re", "pi_im"):
        value = values.get(key)
        if value is not None:
            parts.append(value)
    if not parts:
        return None
    return float(math.sqrt(sum(value * value for value in parts)))


def _scalar_activity(ds, point) -> float | None:
    values = {key: _field_at_point(ds, point, key) for key in ("scalar_activity", "phi_re", "phi_im", "pi_re", "pi_im")}
    return _scalar_activity_from_values(values)


def _default_ball_radius(ds) -> float:
    try:
        width = float(np.max(np.asarray(ds.domain_width, dtype=float)))
        dims = np.asarray(ds.domain_dimensions, dtype=float)
        base_n = float(np.max(dims))
        max_level = int(getattr(ds, "max_level", 0) or 0)
        if base_n <= 0.0 or width <= 0.0:
            return 0.5
        dx_finest = width / (base_n * (2**max_level))
        return float(max(2.0 * dx_finest, 1.0e-6))
    except Exception:  # noqa: BLE001
        return 0.5


def _sphere_sample_points(center, radius: float, n_points: int) -> np.ndarray:
    center = np.asarray(center, dtype=float)
    theta = np.linspace(0.0, math.pi, n_points)
    phi = np.linspace(0.0, 2.0 * math.pi, n_points, endpoint=False)
    theta_grid, phi_grid = np.meshgrid(theta, phi, indexing="ij")
    sin_t = np.sin(theta_grid)
    sx = float(radius) * sin_t * np.cos(phi_grid) + center[0]
    sy = float(radius) * sin_t * np.sin(phi_grid) + center[1]
    sz = float(radius) * np.cos(theta_grid) + center[2]
    return np.column_stack((sx.ravel(), sy.ravel(), sz.ravel()))


def _mean_on_sphere(
    ds,
    center,
    radius: float,
    *,
    n_points: int = 16,
) -> dict[str, float | None]:
    left = np.asarray(ds.domain_left_edge, dtype=float)
    right = np.asarray(ds.domain_right_edge, dtype=float)
    pts = _sphere_sample_points(center, radius, n_points)
    in_domain = (
        (pts[:, 0] >= left[0])
        & (pts[:, 0] <= right[0])
        & (pts[:, 1] >= left[1])
        & (pts[:, 1] <= right[1])
        & (pts[:, 2] >= left[2])
        & (pts[:, 2] <= right[2])
    )
    idxs = np.where(in_domain)[0]
    if idxs.size < 4:
        return {}

    sample_pts = pts[idxs]
    field_names = [
        name
        for name in (
            "rho_req", "lapse", "phi_re", "phi_im", "pi_re", "pi_im",
            "scalar_activity", "ham_abs", "mom_abs", "chi", "K", "weyl4",
            "A11", "A12", "A22",
        )
        if _resolve_field_name(ds, name) is not None
    ]
    yt_fields = [("boxlib", _resolve_field_name(ds, name)) for name in field_names]
    vals = ds.find_field_values_at_points(yt_fields, sample_pts)

    out: dict[str, float | None] = {}
    for name, samples in zip(field_names, vals):
        arr = np.asarray(samples, dtype=float).reshape(-1)
        arr = arr[np.isfinite(arr)]
        if arr.size == 0:
            out[name] = None
        else:
            out[name] = float(np.mean(arr))

    if out.get("scalar_activity") is None:
        out["scalar_activity"] = _scalar_activity_from_values(out)

    phi_re_name = _resolve_field_name(ds, "phi_re")
    phi_im_name = _resolve_field_name(ds, "phi_im")
    pi_re_name = _resolve_field_name(ds, "pi_re")
    pi_im_name = _resolve_field_name(ds, "pi_im")

    if phi_re_name and phi_im_name and pi_re_name and pi_im_name:
        phi_re_arr = np.asarray(
            vals[yt_fields.index(("boxlib", phi_re_name))], dtype=float
        ).reshape(-1)
        phi_im_arr = np.asarray(
            vals[yt_fields.index(("boxlib", phi_im_name))], dtype=float
        ).reshape(-1)
        pi_re_arr = np.asarray(
            vals[yt_fields.index(("boxlib", pi_re_name))], dtype=float
        ).reshape(-1)
        pi_im_arr = np.asarray(
            vals[yt_fields.index(("boxlib", pi_im_name))], dtype=float
        ).reshape(-1)
        mask = (
            np.isfinite(phi_re_arr)
            & np.isfinite(phi_im_arr)
            & np.isfinite(pi_re_arr)
            & np.isfinite(pi_im_arr)
        )
        if np.any(mask):
            out["noether_charge"] = float(
                np.mean(phi_re_arr[mask] * pi_im_arr[mask] - phi_im_arr[mask] * pi_re_arr[mask])
            )
        else:
            out["noether_charge"] = 0.0
        if np.any(mask) and np.any(np.abs(phi_re_arr[mask]) + np.abs(phi_im_arr[mask]) > 0.0):
            phases = np.arctan2(phi_im_arr[mask], phi_re_arr[mask])
            out["phase_coherence"] = float(max(0.0, 1.0 - float(np.std(phases)) / math.pi))
        else:
            out["phase_coherence"] = 0.0
    else:
        out["noether_charge"] = 0.0
        out["phase_coherence"] = 0.0

    if out.get("weyl4") is None and all(
        _resolve_field_name(ds, name) is not None for name in ("A11", "A12", "A22")
    ):
        def _samples(logical: str) -> np.ndarray:
            resolved = _resolve_field_name(ds, logical)
            assert resolved is not None
            return np.asarray(
                vals[yt_fields.index(("boxlib", resolved))], dtype=float
            ).reshape(-1)

        a11_vals = _samples("A11")
        a12_vals = _samples("A12")
        a22_vals = _samples("A22")
        mask = np.isfinite(a11_vals) & np.isfinite(a12_vals) & np.isfinite(a22_vals)
        if np.any(mask):
            gw = np.sqrt(
                (a11_vals[mask] - a22_vals[mask]) ** 2 + (2.0 * a12_vals[mask]) ** 2
            )
            out["weyl4"] = float(np.mean(gw))

    return out


def _extract_at_center(
    ds,
    center,
    *,
    central_ball: bool,
    central_ball_radius: float | None,
) -> dict[str, float | None]:
    if central_ball:
        radius = central_ball_radius if central_ball_radius is not None else _default_ball_radius(ds)
        values = _mean_on_sphere(ds, center, radius)
        if values:
            return values

    point = ds.point([float(center[0]), float(center[1]), float(center[2])])
    point_values = {
        "rho_req": _field_at_point(ds, point, "rho_req"),
        "lapse": _field_at_point(ds, point, "lapse"),
        "scalar_activity": _scalar_activity(ds, point),
        "phi_re": _field_at_point(ds, point, "phi_re"),
        "phi_im": _field_at_point(ds, point, "phi_im"),
        "ham_abs": _field_at_point(ds, point, "ham_abs"),
        "mom_abs": _field_at_point(ds, point, "mom_abs"),
        "noether_charge": 0.0,
        "phase_coherence": 0.0,
        "chi": _field_at_point(ds, point, "chi"),
        "K": _field_at_point(ds, point, "K"),
        "weyl4": _resolve_gw_wave_signal(ds, point),
    }

    radius = central_ball_radius if central_ball_radius is not None else _default_ball_radius(ds)
    try:
        sphere_values = _mean_on_sphere(ds, center, radius)
    except Exception:  # noqa: BLE001
        sphere_values = {}
    if not sphere_values:
        return point_values

    for key in ("chi", "K", "weyl4", "lapse", "rho_req"):
        point_val = point_values.get(key)
        if point_val is None or not math.isfinite(point_val):
            fallback = sphere_values.get(key)
            if fallback is not None and math.isfinite(fallback):
                point_values[key] = fallback
    return point_values


def _extract_central_timeseries_line(
    p: str,
    *,
    t: float,
    center,
    central_ball: bool = False,
    central_ball_radius: float | None = None,
    verbose: bool = False,
) -> str | None:
    import yt

    try:
        ds = yt.load(p)
        values = _extract_at_center(
            ds,
            center,
            central_ball=central_ball,
            central_ball_radius=central_ball_radius,
        )
    except Exception as exc:  # noqa: BLE001
        if verbose:
            print(
                f"WARNING: central-timeseries failed for {os.path.basename(p)}: {exc}"
            )
        return None

    rho_req = values.get("rho_req")
    lapse = values.get("lapse")
    activity = values.get("scalar_activity")

    def _fmt(value: float | None) -> str:
        if value is None or not math.isfinite(value):
            return "nan"
        return f"{value:.16e}"

    if rho_req is None and lapse is None and activity is None:
        return None

    return (
        f"{t:.16e}  {_fmt(rho_req)}  {_fmt(lapse)}  {_fmt(activity)}  "
        f"{_fmt(values.get('phi_re'))}  {_fmt(values.get('phi_im'))}  "
        f"{_fmt(values.get('noether_charge', 0.0))}  {_fmt(values.get('phase_coherence', 0.0))}  "
        f"{_fmt(values.get('ham_abs'))}  {_fmt(values.get('mom_abs'))}  "
        f"{_fmt(values.get('chi'))}  {_fmt(values.get('K'))}  {_fmt(values.get('weyl4'))}"
    )
