"""Pure 3D geometry for GRTresna shell/ring lump layouts."""

from __future__ import annotations

import math

def _vec_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def _vec_scale(a: tuple[float, float, float], s: float) -> tuple[float, float, float]:
    return (a[0] * s, a[1] * s, a[2] * s)


def _vec_cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _vec_norm(a: tuple[float, float, float]) -> tuple[float, float, float]:
    mag = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    if mag < 1e-12:
        return (0.0, 0.0, 0.0)
    return (a[0] / mag, a[1] / mag, a[2] / mag)


def _shell_lump_velocity(
    axis: tuple[float, float, float],
    world: tuple[float, float, float],
    *,
    v_rad: float,
    v_tor: float,
    v_pol: float,
) -> tuple[float, float, float]:
    """Decompose shell currents into radial / toroidal / poloidal components."""
    tor = _vec_norm(_vec_cross(axis, world))
    pol = _vec_cross(tor, world)
    return _vec_add(
        _vec_add(_vec_scale(world, v_rad), _vec_scale(tor, v_tor)),
        _vec_scale(pol, v_pol),
    )


def _shell_exotic_ids(
    num_lumps: int,
    exotic_fraction: float,
    exotic_phase: float,
    azimuths: list[float],
) -> set[int]:
    exotic_count = int(round(exotic_fraction * num_lumps))
    ranked = sorted(
        range(num_lumps),
        key=lambda k: math.cos(azimuths[k] - exotic_phase),
        reverse=True,
    )
    return set(ranked[:exotic_count])


def _expand_shell_sphere(
    *,
    num_lumps: int,
    amp0: float,
    width0: float,
    radius: float,
    thickness: float,
    axis: tuple[float, float, float],
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    dipole: float,
    quadrupole: float,
    radial_jitter: float,
    v_rad: float,
    v_tor: float,
    v_pol: float,
    omega: float,
    mode: int,
    exotic_ids: set[int],
    profile_ids: set[int],
) -> list[dict]:
    golden = math.pi * (3.0 - math.sqrt(5.0))
    lumps: list[dict] = []
    for k in range(num_lumps):
        uz = 1.0 - 2.0 * (k + 0.5) / num_lumps
        rho = math.sqrt(max(0.0, 1.0 - uz * uz))
        phi_k = golden * k
        ux = rho * math.cos(phi_k)
        uy = rho * math.sin(phi_k)
        world = _vec_norm(
            _vec_add(
                _vec_add(_vec_scale(e1, ux), _vec_scale(e2, uy)),
                _vec_scale(axis, uz),
            )
        )
        mu = uz
        legendre = 1.0 + dipole * mu + quadrupole * (1.5 * mu * mu - 0.5)
        amp = max(0.0, min(0.35, amp0 * max(0.15, legendre)))
        width = max(1.0, min(6.0, width0 * max(0.5, 1.0 + 0.25 * quadrupole * mu)))
        r_k = radius + thickness * mu + radial_jitter * radius * 0.3 * math.sin(golden * k)
        center = _vec_scale(world, r_k)
        lumps.append({
            "amp": amp,
            "width": width,
            "center": center,
            "velocity": _shell_lump_velocity(axis, world, v_rad=v_rad, v_tor=v_tor, v_pol=v_pol),
            "omega": omega,
            "mode": max(0, mode),
            "exotic": 1 if k in exotic_ids else 0,
            "profile": 1 if k in profile_ids else 0,
        })
    return lumps


def _expand_shell_channel(
    *,
    num_lumps: int,
    amp0: float,
    width0: float,
    radius: float,
    thickness: float,
    axis: tuple[float, float, float],
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    dipole: float,
    quadrupole: float,
    v_rad: float,
    v_tor: float,
    v_pol: float,
    omega: float,
    mode: int,
    exotic_ids: set[int],
    profile_ids: set[int],
    azimuths: list[float],
) -> list[dict]:
    """Place lumps along the polar axis -- a directed corridor for null rays."""
    span = max(radius, thickness, 1.0)
    lumps: list[dict] = []
    for k in range(num_lumps):
        t = (k + 0.5) / num_lumps - 0.5
        along = span * t * 2.0
        transverse = thickness * math.sin(azimuths[k])
        offset = _vec_add(
            _vec_scale(axis, along),
            _vec_add(_vec_scale(e1, transverse), _vec_scale(e2, 0.5 * thickness * math.cos(azimuths[k]))),
        )
        center = offset
        world = _vec_norm(center) if _vec_norm(center) != (0.0, 0.0, 0.0) else axis
        mu = along / span if span > 0 else 0.0
        legendre = 1.0 + dipole * mu + quadrupole * (1.5 * mu * mu - 0.5)
        amp = max(0.0, min(0.35, amp0 * max(0.15, legendre)))
        width = max(1.0, min(6.0, width0 * max(0.5, 1.0 + 0.25 * quadrupole * mu)))
        lumps.append({
            "amp": amp,
            "width": width,
            "center": center,
            "velocity": _vec_add(_vec_scale(axis, v_rad), _vec_scale(e1, v_tor)),
            "omega": omega,
            "mode": max(0, mode),
            "exotic": 1 if k in exotic_ids else 0,
            "profile": 1 if k in profile_ids else 0,
        })
    return lumps


def _expand_shell_bipolar(
    *,
    num_lumps: int,
    amp0: float,
    width0: float,
    radius: float,
    thickness: float,
    axis: tuple[float, float, float],
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    dipole: float,
    quadrupole: float,
    radial_jitter: float,
    v_rad: float,
    v_tor: float,
    v_pol: float,
    omega: float,
    mode: int,
    exotic_ids: set[int],
    profile_ids: set[int],
    azimuths: list[float],
) -> list[dict]:
    """Two matter clusters on opposite sides of the origin along the axis."""
    half = max(1, num_lumps // 2)
    lumps: list[dict] = []
    for k in range(num_lumps):
        sign = 1.0 if k < half else -1.0
        mu = 1.0 - 2.0 * abs((k % half) + 0.5) / half
        r_k = radius + thickness * abs(mu) + radial_jitter * radius * 0.2 * math.sin(azimuths[k])
        offset = _vec_add(
            _vec_scale(axis, sign * r_k),
            _vec_add(
                _vec_scale(e1, 0.3 * thickness * math.cos(azimuths[k])),
                _vec_scale(e2, 0.3 * thickness * math.sin(azimuths[k])),
            ),
        )
        world = _vec_norm(_vec_scale(axis, sign)) if sign != 0 else axis
        legendre = 1.0 + dipole * sign * mu + quadrupole * (1.5 * mu * mu - 0.5)
        amp = max(0.0, min(0.35, amp0 * max(0.15, legendre)))
        width = max(1.0, min(6.0, width0 * max(0.5, 1.0 + 0.25 * quadrupole * mu)))
        lumps.append({
            "amp": amp,
            "width": width,
            "center": offset,
            "velocity": _shell_lump_velocity(axis, world, v_rad=v_rad * sign, v_tor=v_tor, v_pol=v_pol),
            "omega": omega,
            "mode": max(0, mode),
            "exotic": 1 if k in exotic_ids else 0,
            "profile": 1 if k in profile_ids else 0,
        })
    return lumps


def _expand_shell_ring(
    *,
    num_lumps: int,
    amp0: float,
    width0: float,
    radius: float,
    thickness: float,
    axis: tuple[float, float, float],
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    dipole: float,
    quadrupole: float,
    v_rad: float,
    v_tor: float,
    v_pol: float,
    omega: float,
    mode: int,
    exotic_ids: set[int],
    profile_ids: set[int],
    azimuths: list[float],
) -> list[dict]:
    """Planar ring in the plane orthogonal to the polar axis."""
    lumps: list[dict] = []
    for k, theta in enumerate(azimuths):
        c = math.cos(theta)
        s = math.sin(theta)
        amp_mod = 1.0 + dipole * c + quadrupole * math.cos(2.0 * theta)
        width_mod = 1.0 + 0.25 * quadrupole * math.sin(2.0 * theta)
        amp = max(0.0, min(0.35, amp0 * max(0.15, amp_mod)))
        width = max(1.0, min(6.0, width0 * max(0.5, width_mod)))
        center = _vec_add(
            _vec_add(_vec_scale(e1, radius * c), _vec_scale(e2, radius * s)),
            _vec_scale(axis, thickness * math.sin(theta)),
        )
        world = _vec_norm(center) if _vec_norm(center) != (0.0, 0.0, 0.0) else e1
        tor = _vec_norm(_vec_cross(axis, world))
        pol = _vec_cross(tor, world)
        velocity = _vec_add(
            _vec_add(_vec_scale(world, v_rad), _vec_scale(tor, v_tor)),
            _vec_scale(pol, v_pol),
        )
        lumps.append({
            "amp": amp,
            "width": width,
            "center": center,
            "velocity": velocity,
            "omega": omega,
            "mode": max(0, mode),
            "exotic": 1 if k in exotic_ids else 0,
            "profile": 1 if k in profile_ids else 0,
        })
    return lumps


# Low-discrepancy R3 (Roberts 2018) additive-recurrence basis: deterministic,
# parameter-free quasi-random offsets in [0,1)^3. The 'cloud' layout uses it to
# scatter lumps irregularly (asymmetrically) yet reproducibly -- the same
# searched params always yield the same cloud, so the optimizer sees a smooth,
# learnable knob rather than per-candidate noise.
_R3_ALPHA = (0.8191725133961644, 0.6710436067037892, 0.5497004779019703)


def _r3_point(k: int) -> tuple[float, float, float]:
    return (
        (0.5 + _R3_ALPHA[0] * (k + 1)) % 1.0,
        (0.5 + _R3_ALPHA[1] * (k + 1)) % 1.0,
        (0.5 + _R3_ALPHA[2] * (k + 1)) % 1.0,
    )


def _expand_shell_cloud(
    *,
    num_lumps: int,
    amp0: float,
    width0: float,
    radius: float,
    thickness: float,
    axis: tuple[float, float, float],
    e1: tuple[float, float, float],
    e2: tuple[float, float, float],
    dipole: float,
    quadrupole: float,
    radial_jitter: float,
    v_rad: float,
    v_tor: float,
    v_pol: float,
    omega: float,
    mode: int,
    exotic_ids: set[int],
    profile_ids: set[int],
) -> list[dict]:
    """Deterministic quasi-random cloud: lumps scattered irregularly in a ball.

    Unlike the symmetric sphere/ring/bipolar layouts, the cloud breaks all
    spatial symmetry while staying bounded and reproducible (fixed R3 sequence),
    letting the optimizer explore fully asymmetric matter configurations.
    """
    span = max(radius + thickness, 1.0)
    lumps: list[dict] = []
    for k in range(num_lumps):
        u1, u2, u3 = _r3_point(k)
        mu = 2.0 * u2 - 1.0  # cos(theta)
        sin_theta = math.sqrt(max(0.0, 1.0 - mu * mu))
        azimuth = 2.0 * math.pi * u3
        world = _vec_norm(
            _vec_add(
                _vec_add(
                    _vec_scale(e1, sin_theta * math.cos(azimuth)),
                    _vec_scale(e2, sin_theta * math.sin(azimuth)),
                ),
                _vec_scale(axis, mu),
            )
        )
        if world == (0.0, 0.0, 0.0):
            world = axis
        r_k = span * (u1 ** (1.0 / 3.0)) * (1.0 + radial_jitter * 0.3 * (u3 - 0.5))
        center = _vec_scale(world, r_k)
        legendre = 1.0 + dipole * mu + quadrupole * (1.5 * mu * mu - 0.5)
        amp = max(0.0, min(0.35, amp0 * max(0.15, legendre)))
        width = max(1.0, min(6.0, width0 * max(0.5, 1.0 + 0.25 * quadrupole * mu)))
        lumps.append({
            "amp": amp,
            "width": width,
            "center": center,
            "velocity": _shell_lump_velocity(axis, world, v_rad=v_rad, v_tor=v_tor, v_pol=v_pol),
            "omega": omega,
            "mode": max(0, mode),
            "exotic": 1 if k in exotic_ids else 0,
            "profile": 1 if k in profile_ids else 0,
        })
    return lumps


def _orthonormal_frame(
    axis_theta: float, axis_phi: float
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    """Return an orthonormal frame ``(a, e1, e2)`` with ``a`` the polar axis.

    ``a`` points along ``(theta, phi)``; ``e1, e2`` span the plane orthogonal
    to it.  Placing canonical Fibonacci points in this frame rotates the whole
    sphere to the searched orientation for free.
    """
    sa = math.sin(axis_theta)
    a = (sa * math.cos(axis_phi), sa * math.sin(axis_phi), math.cos(axis_theta))
    # Reference axis least parallel to a, to keep the cross product well-conditioned.
    ref = (0.0, 0.0, 1.0) if abs(a[2]) < 0.9 else (1.0, 0.0, 0.0)
    e1 = _vec_norm(_vec_cross(ref, a))
    if e1 == (0.0, 0.0, 0.0):
        e1 = _vec_norm(_vec_cross((0.0, 1.0, 0.0), a))
    e2 = _vec_cross(a, e1)
    return a, e1, e2

