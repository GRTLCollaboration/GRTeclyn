"""Co-moving stability metrics from shell profiles and shift data."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..io.dat import numeric_rows
from ..probes.ftl.analytic import _axis_profiles
from ..types.diagnostics import STATIONARY_BETA_EPS, ComovingMetrics


def _parse_shell_profiles(path: Path) -> tuple[list[float], dict[str, list[float]]] | None:
    if not path.exists():
        return None
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    if not lines:
        return None
    header = lines[0].lstrip("#").split()
    if not header or header[0] != "time":
        return None
    cols: dict[str, list[float]] = {name: [] for name in header[1:]}
    times: list[float] = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) != len(header):
            continue
        times.append(float(parts[0]))
        for name, value in zip(header[1:], parts[1:]):
            cols[name].append(float(value))
    if not times:
        return None
    return times, cols


def _mean_beta_in_bubble(
    overrides: dict[str, object],
    *,
    ftl_L: float | None = None,
) -> float | None:
    L = ftl_L if ftl_L is not None else float(overrides.get("recipe_basis_radius_max", 8.0))
    _x, _r, _chi, _alpha, beta_x = _axis_profiles(overrides, L=L, n_points=512)
    beta_max = float(max(abs(float(v)) for v in beta_x))
    if beta_max < 1.0e-8:
        return 0.0
    bubble_mask = abs(beta_x) >= 0.1 * beta_max
    if not bubble_mask.any():
        return float(beta_x.mean())
    return float(beta_x[bubble_mask].mean())


def _peak_shift_in_bubble_from_gridinit(
    episode_dir: Path,
    *,
    ftl_L: float | None = None,
) -> float | None:
    gridinit_path = episode_dir / "initial_data.gridinit"
    if not gridinit_path.exists():
        return None
    try:
        from ...grtresna.io import read_gridinit

        grid = read_gridinit(gridinit_path)
    except Exception:
        return None

    try:
        i1 = grid.comp_names.index("shift1")
        i2 = grid.comp_names.index("shift2")
        i3 = grid.comp_names.index("shift3")
    except ValueError:
        return None

    data = grid.data
    beta_mag = np.sqrt(
        data[..., i1] ** 2 + data[..., i2] ** 2 + data[..., i3] ** 2
    )

    if ftl_L is not None and ftl_L > 0.0:
        nz, ny, nx, _ = data.shape
        dx_x, dx_y, dx_z = (float(v) for v in grid.dx_xyz)
        ox, oy, oz = (float(v) for v in grid.origin)
        xs = ox + (np.arange(nx) + 0.5) * dx_x - (ox + 0.5 * nx * dx_x)
        ys = oy + (np.arange(ny) + 0.5) * dx_y - (oy + 0.5 * ny * dx_y)
        zs = oz + (np.arange(nz) + 0.5) * dx_z - (oz + 0.5 * nz * dx_z)
        zz, yy, xx = np.meshgrid(zs, ys, xs, indexing="ij")
        bubble = (xx * xx + yy * yy + zz * zz) <= ftl_L * ftl_L
        if bubble.any():
            beta_mag = beta_mag[bubble]

    return float(beta_mag.max()) if beta_mag.size else 0.0


def read_comoving_metrics(
    episode_dir: Path,
    overrides: dict[str, object] | None,
    *,
    ftl_L: float | None = None,
) -> ComovingMetrics | None:
    if overrides is None:
        return None

    beta_mean = _mean_beta_in_bubble(overrides, ftl_L=ftl_L)
    if beta_mean is None or abs(beta_mean) < 1.0e-8:
        peak_shift = _peak_shift_in_bubble_from_gridinit(episode_dir, ftl_L=ftl_L)
        if peak_shift is not None:
            beta_mean = peak_shift
    if beta_mean is not None and abs(beta_mean) < STATIONARY_BETA_EPS:
        return ComovingMetrics(
            beta_mean=beta_mean,
            delta_comoving=None,
            score=None,
            stationary=True,
        )

    shell_path = episode_dir / "small_data" / "shell_profiles.dat"
    if not shell_path.exists():
        shell_path = episode_dir / "shell_profiles.dat"
    parsed = _parse_shell_profiles(shell_path)

    chi_series: list[float] | None = None
    times: list[float] | None = None
    if parsed is not None:
        times, cols = parsed
        chi_keys = sorted(key for key in cols if key.startswith("chi_mean_"))
        if chi_keys:
            chi_series = cols[chi_keys[0]]

    if chi_series is None or times is None or len(chi_series) < 2:
        collapse_path = episode_dir / "data" / "collapse_diagnostics.dat"
        if not collapse_path.exists():
            collapse_path = episode_dir / "collapse_diagnostics.dat"
        rows = numeric_rows(collapse_path, 3)
        if len(rows) < 2:
            return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)
        times = [row[0] for row in rows]
        chi_series = [row[2] for row in rows]

    if beta_mean is None or len(chi_series) < 2 or len(times) < 2:
        return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)

    times_arr = np.asarray(times, dtype=float)
    chi_arr = np.asarray(chi_series, dtype=float)
    if times_arr[-1] <= times_arr[0]:
        return ComovingMetrics(beta_mean=beta_mean, delta_comoving=None, score=None)

    dchi_dt = np.gradient(chi_arr, times_arr)
    eulerian_rate = float(np.max(np.abs(dchi_dt)))
    shift_scale = max(abs(beta_mean), 1.0e-3)
    delta_comoving = eulerian_rate / shift_scale
    score = 1.0 / (1.0 + delta_comoving)
    return ComovingMetrics(beta_mean=beta_mean, delta_comoving=delta_comoving, score=score)
