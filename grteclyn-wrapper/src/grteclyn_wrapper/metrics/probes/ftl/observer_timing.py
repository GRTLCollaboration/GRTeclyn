"""Freely falling emitter-to-receiver timing on an evolving metric.

This probe is deliberately separate from campaign scoring.  It integrates two
timelike geodesics, emits a null ray at a prescribed emitter proper time, and
shoots the ray until it intersects the receiver worldline.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import trapezoid
from scipy.optimize import least_squares

from .geodesic import _hamiltonian_rhs, _null_relative_drift
from .metric_field import EvolvingMetricField, MetricField
from .metric_stack_cache import evolving_field_from_metric_stack_cache


Vector = NDArray[np.float64]


@dataclass(frozen=True)
class TimelikeWorldline:
    tau: Vector
    x: NDArray[np.float64]
    p_cov: NDArray[np.float64]
    max_mass_shell_drift: float


@dataclass(frozen=True)
class FreeFallTimingReport:
    emission_tau: float
    emission_t: float
    reception_tau: float | None
    reception_t: float | None
    initial_proper_separation: float
    flat_reception_tau: float
    receiver_clock_advance: float | None
    fractional_arrival_advance: float
    miss_distance: float
    reception_tolerance: float
    reached: bool
    max_null_relative_drift: float
    emitter_mass_shell_drift: float
    receiver_mass_shell_drift: float
    emitter_displacement: tuple[float, float, float]
    receiver_displacement: tuple[float, float, float]
    optimizer_success: bool
    optimizer_evaluations: int
    notes: tuple[str, ...] = ()


def _mass_shell(g_inv: NDArray[np.float64], p_cov: Vector) -> float:
    return float(g_inv @ p_cov @ p_cov)


def _project_mass_shell(
    g: NDArray[np.float64],
    g_inv: NDArray[np.float64],
    p_cov: Vector,
    *,
    target: float,
) -> Vector:
    """Adjust p_0 so g^{mu nu} p_mu p_nu equals ``target``."""
    spatial = p_cov[1:]
    a = float(g_inv[0, 0])
    b = 2.0 * float(g_inv[0, 1:] @ spatial)
    c = float(spatial @ g_inv[1:, 1:] @ spatial) - target
    roots = np.roots((a, b, c))
    candidates: list[Vector] = []
    for root in roots:
        if abs(float(np.imag(root))) > 1.0e-10:
            continue
        candidate = p_cov.copy()
        candidate[0] = float(np.real(root))
        if float((g_inv @ candidate)[0]) > 0.0:
            candidates.append(candidate)
    if not candidates:
        raise ValueError("could not construct a future-directed mass-shell covector")
    return min(candidates, key=lambda candidate: abs(candidate[0] - p_cov[0]))


def _eulerian_covector(field: MetricField, x: Vector) -> Vector:
    g, g_inv, _ = field.sample(x)
    alpha = 1.0 / np.sqrt(max(-float(g_inv[0, 0]), 1.0e-30))
    u_contra = -alpha * g_inv[:, 0]
    return g @ u_contra


def _rk4_step(
    field: MetricField,
    x: Vector,
    p_cov: Vector,
    step: float,
) -> tuple[Vector, Vector]:
    def rhs(xp: Vector, pp: Vector) -> tuple[Vector, Vector]:
        _g, g_inv, dg_inv = field.sample(xp)
        return _hamiltonian_rhs(g_inv, dg_inv, pp)

    k1x, k1p = rhs(x, p_cov)
    k2x, k2p = rhs(x + 0.5 * step * k1x, p_cov + 0.5 * step * k1p)
    k3x, k3p = rhs(x + 0.5 * step * k2x, p_cov + 0.5 * step * k2p)
    k4x, k4p = rhs(x + step * k3x, p_cov + step * k3p)
    return (
        x + (step / 6.0) * (k1x + 2.0 * k2x + 2.0 * k3x + k4x),
        p_cov + (step / 6.0) * (k1p + 2.0 * k2p + 2.0 * k3p + k4p),
    )


def integrate_timelike_geodesic(
    field: EvolvingMetricField,
    initial_position: Vector,
    *,
    d_tau: float = 0.025,
) -> TimelikeWorldline:
    """Integrate an initially Eulerian freely falling observer."""
    x = np.asarray(initial_position, dtype=float).copy()
    p_cov = _eulerian_covector(field, x)
    tau = 0.0
    taus = [tau]
    positions = [x.copy()]
    momenta = [p_cov.copy()]
    max_drift = 0.0
    t_end = float(field.times[-1])

    while x[0] < t_end:
        x_new, p_new = _rk4_step(field, x, p_cov, d_tau)
        g, g_inv, _ = field.sample(x_new)
        drift = abs(_mass_shell(g_inv, p_new) + 1.0)
        max_drift = max(max_drift, drift)
        if drift > 1.0e-8:
            p_new = _project_mass_shell(g, g_inv, p_new, target=-1.0)
        if x_new[0] <= x[0]:
            raise ValueError("timelike worldline ceased to be future directed")
        tau += d_tau
        x, p_cov = x_new, p_new
        taus.append(tau)
        positions.append(x.copy())
        momenta.append(p_cov.copy())

    return TimelikeWorldline(
        tau=np.asarray(taus),
        x=np.asarray(positions),
        p_cov=np.asarray(momenta),
        max_mass_shell_drift=max_drift,
    )


def _worldline_at(
    worldline: TimelikeWorldline,
    value: float,
    *,
    parameter: str,
) -> tuple[Vector, Vector, float]:
    grid = worldline.tau if parameter == "tau" else worldline.x[:, 0]
    if value < float(grid[0]) or value > float(grid[-1]):
        raise ValueError(f"{parameter}={value} lies outside observer worldline")
    x = np.array([np.interp(value, grid, worldline.x[:, i]) for i in range(4)])
    p = np.array([np.interp(value, grid, worldline.p_cov[:, i]) for i in range(4)])
    tau = value if parameter == "tau" else float(np.interp(value, grid, worldline.tau))
    return x, p, tau


def _orthonormal_spatial_triad(g: NDArray[np.float64], u: Vector) -> NDArray[np.float64]:
    triad: list[Vector] = []
    for spatial_idx in range(1, 4):
        v = np.zeros(4)
        v[spatial_idx] = 1.0
        v += u * float(u @ g @ v)
        for previous in triad:
            v -= previous * float(previous @ g @ v)
        norm2 = float(v @ g @ v)
        if norm2 <= 1.0e-14:
            raise ValueError("could not build emitter spatial tetrad")
        triad.append(v / np.sqrt(norm2))
    return np.asarray(triad)


def _initial_direction(
    g: NDArray[np.float64],
    u: Vector,
    triad: NDArray[np.float64],
    displacement: Vector,
) -> Vector:
    delta = np.zeros(4)
    delta[1:] = displacement
    components = np.array([float(e @ g @ delta) for e in triad])
    norm = float(np.linalg.norm(components))
    if norm <= 1.0e-14:
        raise ValueError("emitter and receiver are coincident")
    return components / norm


def _direction_basis(direction: Vector) -> tuple[Vector, Vector]:
    trial = np.array([1.0, 0.0, 0.0])
    if abs(float(trial @ direction)) > 0.8:
        trial = np.array([0.0, 1.0, 0.0])
    first = trial - direction * float(trial @ direction)
    first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    return first, second


def _integrate_null_to_affine(
    field: EvolvingMetricField,
    x0: Vector,
    p0: Vector,
    affine_end: float,
    *,
    ds: float,
) -> tuple[Vector, float]:
    x = x0.copy()
    p_cov = p0.copy()
    remaining = float(affine_end)
    max_drift = 0.0
    while remaining > 1.0e-12:
        step = min(ds, remaining)
        x, p_cov = _rk4_step(field, x, p_cov, step)
        g, g_inv, _ = field.sample(x)
        max_drift = max(max_drift, _null_relative_drift(g_inv, p_cov))
        if abs(_mass_shell(g_inv, p_cov)) > 1.0e-8:
            p_cov = _project_mass_shell(g, g_inv, p_cov, target=0.0)
        remaining -= step
    return x, max_drift


def _initial_proper_separation(
    field: EvolvingMetricField,
    start: Vector,
    end: Vector,
    *,
    samples: int = 257,
) -> float:
    delta = end[1:] - start[1:]
    values = []
    for fraction in np.linspace(0.0, 1.0, samples):
        point = start.copy()
        point[1:] += fraction * delta
        gamma = field.sample_g(point)[1:, 1:]
        values.append(np.sqrt(max(float(delta @ gamma @ delta), 0.0)))
    return float(trapezoid(values, dx=1.0 / (samples - 1)))


def compute_freefall_observer_timing(
    field: EvolvingMetricField,
    *,
    emission_tau: float = 4.0,
    observer_step: float = 0.025,
    ray_step: float = 0.05,
    reception_tolerance: float | None = None,
) -> FreeFallTimingReport:
    """Measure arrival at a freely falling receiver's worldline."""
    origin = field.origin
    shape = field.spatial_shape
    spacing = field.spatial_spacing
    center = origin + 0.5 * (np.asarray(shape) - 1.0) * np.asarray(spacing)
    start = center.copy()
    end = center.copy()
    start[0] = origin[0] + 0.05 * (shape[0] - 1) * spacing[0]
    end[0] = origin[0] + 0.95 * (shape[0] - 1) * spacing[0]
    t0 = float(field.times[0])
    emitter_initial = np.array([t0, *start])
    receiver_initial = np.array([t0, *end])

    emitter = integrate_timelike_geodesic(
        field, emitter_initial, d_tau=observer_step
    )
    receiver = integrate_timelike_geodesic(
        field, receiver_initial, d_tau=observer_step
    )
    emission_x, emission_p, _ = _worldline_at(
        emitter, emission_tau, parameter="tau"
    )
    receiver_at_emission, _, _ = _worldline_at(
        receiver, emission_x[0], parameter="t"
    )
    g_emit, g_inv_emit, _ = field.sample(emission_x)
    u_emit = g_inv_emit @ emission_p
    triad = _orthonormal_spatial_triad(g_emit, u_emit)
    direction0 = _initial_direction(
        g_emit,
        u_emit,
        triad,
        receiver_at_emission[1:] - emission_x[1:],
    )
    basis1, basis2 = _direction_basis(direction0)
    coordinate_distance = float(
        np.linalg.norm(receiver_at_emission[1:] - emission_x[1:])
    )
    initial_separation = _initial_proper_separation(
        field, emitter_initial, receiver_initial
    )

    best_null_drift = 0.0

    def trace(parameters: Vector) -> tuple[Vector, float]:
        nonlocal best_null_drift
        direction = (
            direction0 + parameters[0] * basis1 + parameters[1] * basis2
        )
        direction /= np.linalg.norm(direction)
        k_contra = u_emit + direction @ triad
        k_cov = g_emit @ k_contra
        ray_x, drift = _integrate_null_to_affine(
            field,
            emission_x,
            k_cov,
            float(parameters[2]),
            ds=ray_step,
        )
        best_null_drift = max(best_null_drift, drift)
        return ray_x, drift

    scale = np.maximum(np.asarray(spacing), 1.0e-6)

    def residual(parameters: Vector) -> Vector:
        try:
            ray_x, _ = trace(parameters)
            receiver_x, _, _ = _worldline_at(
                receiver, float(ray_x[0]), parameter="t"
            )
            return (ray_x[1:] - receiver_x[1:]) / scale
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            return np.full(3, 1.0e6)

    solution = least_squares(
        residual,
        np.array([0.0, 0.0, max(coordinate_distance, 1.0e-3)]),
        bounds=(
            np.array([-1.5, -1.5, max(0.1, 0.25 * coordinate_distance)]),
            np.array([1.5, 1.5, max(1.0, 2.5 * coordinate_distance)]),
        ),
        xtol=1.0e-9,
        ftol=1.0e-9,
        gtol=1.0e-9,
        max_nfev=80,
    )
    ray_x, _ = trace(solution.x)
    reception_in_window = True
    try:
        receiver_x, _, reception_tau = _worldline_at(
            receiver, float(ray_x[0]), parameter="t"
        )
        separation = ray_x[1:] - receiver_x[1:]
        gamma_reception = field.sample_g(receiver_x)[1:, 1:]
        miss_distance = np.sqrt(
            max(float(separation @ gamma_reception @ separation), 0.0)
        )
    except ValueError:
        reception_in_window = False
        receiver_x = receiver.x[-1]
        reception_tau = float(receiver.tau[-1])
        miss_distance = 1.0e30
    tolerance = (
        float(reception_tolerance)
        if reception_tolerance is not None
        else 0.25 * min(spacing)
    )
    reached = bool(
        solution.success and reception_in_window and miss_distance <= tolerance
    )
    flat_reception_tau = emission_tau + initial_separation
    advance = flat_reception_tau - reception_tau if reached else None
    fraction = (
        max(0.0, advance / flat_reception_tau)
        if advance is not None and flat_reception_tau > 0.0
        else 0.0
    )

    return FreeFallTimingReport(
        emission_tau=emission_tau,
        emission_t=float(emission_x[0]),
        reception_tau=reception_tau if reached else None,
        reception_t=float(ray_x[0]) if reached else None,
        initial_proper_separation=initial_separation,
        flat_reception_tau=flat_reception_tau,
        receiver_clock_advance=advance,
        fractional_arrival_advance=fraction,
        miss_distance=miss_distance,
        reception_tolerance=tolerance,
        reached=reached,
        max_null_relative_drift=best_null_drift,
        emitter_mass_shell_drift=emitter.max_mass_shell_drift,
        receiver_mass_shell_drift=receiver.max_mass_shell_drift,
        emitter_displacement=tuple(
            float(v) for v in emission_x[1:] - emitter_initial[1:]
        ),
        receiver_displacement=tuple(
            float(v) for v in receiver_x[1:] - receiver_initial[1:]
        ),
        optimizer_success=bool(solution.success),
        optimizer_evaluations=int(solution.nfev),
        notes=(
            "initially Eulerian freely falling observers",
            "emission scheduled by emitter proper time",
            "flat reference uses matched initial proper separation",
        ),
    )


def compute_freefall_observer_timing_from_cache(
    cache_dir: Path,
    *,
    emission_tau: float = 4.0,
) -> FreeFallTimingReport | None:
    field = evolving_field_from_metric_stack_cache(cache_dir)
    if field is None:
        return None
    return compute_freefall_observer_timing(field, emission_tau=emission_tau)


def write_freefall_observer_timing_json(
    path: Path, report: FreeFallTimingReport
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(asdict(report), indent=2) + "\n", encoding="utf-8")
