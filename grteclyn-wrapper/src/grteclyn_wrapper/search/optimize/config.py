"""Build GRTresnaConfig from optimizer override dicts."""

from __future__ import annotations

import logging
import math
import re
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping

from .geometry import (
    _expand_shell_bipolar,
    _expand_shell_channel,
    _expand_shell_cloud,
    _expand_shell_ring,
    _expand_shell_sphere,
    _orthonormal_frame,
    _shell_exotic_ids,
)
from .spaces import GRTRESNA_DEFAULT_NUM_LUMPS

if TYPE_CHECKING:
    from ...grtresna.solver import GRTresnaConfig

logger = logging.getLogger(__name__)

_LUMP_KEY_RE = re.compile(r"^grtresna_lump(\d+)_(\w+)$")
# Trajectory campaigns name their lumps differently, and every consumer that
# matched only the line above treated them as having no lumps at all.
_TRAJECTORY_LUMP_KEY_RE = re.compile(r"^trajectory_lump(\d+)_(\w+)$")

# The four knobs that give the shell's matter a velocity, with the defaults the
# decoder falls back to when a campaign does not supply them.  The toroidal
# default is 0.35c -- a real, large frame-dragging current that no search
# dimension and no trajectory record ever mentioned.  See DebugPreGPU.md PG-2.
SHELL_CURRENT_DEFAULTS: dict[str, float] = {
    "grtresna_shell_toroidal_velocity": 0.35,
    "grtresna_shell_poloidal_velocity": 0.0,
    "grtresna_shell_radial_velocity": 0.0,
    "grtresna_shell_omega": 0.0,
}


def resolved_shell_currents(overrides: Mapping[str, Any]) -> dict[str, float]:
    """The matter currents the decoder will actually use, defaults included.

    Two rules, and they have to agree or the campaign is describing something it
    is not running:

    * ``grtresna_shell_static >= 1`` means static matter -- every current zero.
    * A campaign that supplies **no** current key is not searching motion, so
      it gets static matter too.  Previously such a campaign still picked up the
      0.35c toroidal default for any candidate whose ``shell_static`` bit
      happened to decode to 0, which is half of them, since that dimension stays
      in the space with its centre at 0 = moving.
    """
    supplied = [key for key in SHELL_CURRENT_DEFAULTS if key in overrides]

    static = False
    if not supplied:
        static = True
    else:
        try:
            static = int(round(float(overrides.get("grtresna_shell_static", 0.0)))) >= 1
        except (TypeError, ValueError):
            static = False

    resolved: dict[str, float] = {}
    for key, default in SHELL_CURRENT_DEFAULTS.items():
        if static:
            resolved[key] = 0.0
            continue
        try:
            resolved[key] = float(overrides.get(key, default))
        except (TypeError, ValueError):
            resolved[key] = default
    return resolved


def unrecorded_shell_currents(
    overrides: Mapping[str, Any],
    searched_keys: frozenset[str] | set[str],
) -> dict[str, float]:
    """Non-zero currents that no search dimension and no override declares.

    A trajectory record lists the searched dimensions only, so a current that
    arrives as a decode default leaves no trace anywhere -- which is how a fixed
    0.35c toroidal flow rode along through a whole campaign unnoticed.
    """
    resolved = resolved_shell_currents(overrides)
    return {
        key: value
        for key, value in resolved.items()
        if value != 0.0 and key not in searched_keys and key not in overrides
    }


def _expand_shell_lumps_from_overrides(
    overrides: Mapping[str, Any],
    get_float: Any,
    *,
    canonical_boson: bool = False,
) -> list[dict]:
    """Expand ``grtresna_shell_*`` search knobs into a lump list."""
    num_lumps = max(3, int(round(get_float("grtresna_shell_lumps", GRTRESNA_DEFAULT_NUM_LUMPS))))
    amp0 = get_float("grtresna_shell_amp", 0.15)
    width0 = get_float("grtresna_shell_width", 3.0)
    radius = get_float("grtresna_shell_radius", 3.5)
    thickness = get_float("grtresna_shell_thickness", 0.5)
    axis_theta = get_float("grtresna_shell_axis_theta", 0.5 * math.pi)
    axis_phi = get_float("grtresna_shell_axis_phi", 0.0)
    # Single source of truth for the static/moving decision -- the space builder
    # strips the velocity dimensions in the static families, and this is what
    # makes "stripped" mean "no motion" instead of "0.35c you cannot see".
    currents = resolved_shell_currents(overrides)
    v_tor = currents["grtresna_shell_toroidal_velocity"]
    v_pol = currents["grtresna_shell_poloidal_velocity"]
    v_rad = currents["grtresna_shell_radial_velocity"]
    omega = currents["grtresna_shell_omega"]
    dipole = get_float("grtresna_shell_dipole_amp", 0.0)
    quadrupole = get_float("grtresna_shell_quadrupole_amp", 0.0)
    if canonical_boson:
        exotic_fraction = 0.0
        exotic_phase = 0.0
    else:
        exotic_fraction = min(1.0, max(0.0, get_float("grtresna_shell_exotic_fraction", 0.4)))
        exotic_phase = get_float("grtresna_shell_exotic_phase", 0.0)
    profile_fraction = min(1.0, max(0.0, get_float("grtresna_shell_profile_fraction", 0.0)))
    profile_phase = get_float("grtresna_shell_profile_phase", 0.0)
    mode = int(round(get_float("grtresna_shell_mode", 1.0)))
    radial_jitter = min(1.0, max(0.0, get_float("grtresna_shell_radial_jitter", 0.0)))

    layout = int(round(get_float("grtresna_matter_layout", 0.0)))
    layout = max(0, min(4, layout))
    axis, e1, e2 = _orthonormal_frame(axis_theta, axis_phi)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    azimuths = [golden * k for k in range(num_lumps)]
    exotic_ids = set() if canonical_boson else _shell_exotic_ids(
        num_lumps, exotic_fraction, exotic_phase, azimuths
    )
    profile_ids = _shell_exotic_ids(num_lumps, profile_fraction, profile_phase, azimuths)

    common = dict(
        num_lumps=num_lumps,
        amp0=amp0,
        width0=width0,
        radius=radius,
        thickness=thickness,
        axis=axis,
        e1=e1,
        e2=e2,
        dipole=dipole,
        quadrupole=quadrupole,
        v_rad=v_rad,
        v_tor=v_tor,
        v_pol=v_pol,
        omega=omega,
        mode=mode,
        exotic_ids=exotic_ids,
        profile_ids=profile_ids,
        azimuths=azimuths,
    )
    if layout == 1:
        lumps = _expand_shell_channel(**common)
    elif layout == 2:
        lumps = _expand_shell_bipolar(**common, radial_jitter=radial_jitter)
    elif layout == 3:
        ring_az = [2.0 * math.pi * k / num_lumps for k in range(num_lumps)]
        ring_common = {k: v for k, v in common.items() if k != "azimuths"}
        lumps = _expand_shell_ring(**ring_common, azimuths=ring_az)
    elif layout == 4:
        cloud_common = {k: v for k, v in common.items() if k != "azimuths"}
        lumps = _expand_shell_cloud(**cloud_common, radial_jitter=radial_jitter)
    else:
        sphere_common = {k: v for k, v in common.items() if k != "azimuths"}
        lumps = _expand_shell_sphere(**sphere_common, radial_jitter=radial_jitter)
    if canonical_boson:
        for lump in lumps:
            lump["exotic"] = 0
    return lumps


def _expand_sh_lumps_from_overrides(
    overrides: Mapping[str, Any],
    get_float: Any,
) -> list[dict]:
    """Expand ``grtresna_sh_*`` search knobs into a lump list.

    Uses spherical-harmonic coefficients to modulate per-lump amplitudes
    on a Fibonacci-lattice sphere, giving genuine per-lump variation
    through a smooth angular expansion.
    """
    from ...grtresna.fields.sh import (
        cartesian_to_spherical,
        eval_sh_modulation,
    )
    from .spaces import SH_DEFAULT_ELL_MAX

    num_lumps = max(3, int(round(get_float("grtresna_sh_lumps", GRTRESNA_DEFAULT_NUM_LUMPS))))
    amp0 = get_float("grtresna_sh_amp", 0.13)
    width0 = get_float("grtresna_sh_width", 2.4)
    radius = get_float("grtresna_sh_radius", 3.5)
    thickness = get_float("grtresna_sh_thickness", 0.5)

    v_tor = get_float("grtresna_sh_toroidal_velocity", 0.4)
    v_pol = get_float("grtresna_sh_poloidal_velocity", 0.0)
    v_rad = get_float("grtresna_sh_radial_velocity", 0.0)
    omega = get_float("grtresna_sh_omega", 0.0)
    if int(round(get_float("grtresna_sh_static", 0.0))) >= 1:
        v_tor = v_pol = v_rad = omega = 0.0

    exotic_fraction = min(1.0, max(0.0, get_float("grtresna_sh_exotic_fraction", 0.4)))
    exotic_phase = get_float("grtresna_sh_exotic_phase", 0.0)

    ell_max = int(round(get_float("grtresna_sh_ell_max", float(SH_DEFAULT_ELL_MAX))))
    n_sh = (ell_max + 1) ** 2

    # Collect SH modulation coefficients (idx 1..N-1; monopole excluded).
    sh_coeffs: list[float] = [0.0] * n_sh
    for idx in range(1, n_sh):
        sh_coeffs[idx] = get_float(f"grtresna_sh_c{idx}", 0.0)

    # Fibonacci lattice on the unit sphere.
    golden = math.pi * (3.0 - math.sqrt(5.0))
    axis = (0.0, 0.0, 1.0)

    # Exotic wedge selection (same idiom as shell ansatz).
    azimuths = [golden * k for k in range(num_lumps)]
    exotic_ids = _shell_exotic_ids(num_lumps, exotic_fraction, exotic_phase, azimuths)

    lumps: list[dict] = []
    for k in range(num_lumps):
        # Fibonacci sphere position.
        uz = 1.0 - 2.0 * (k + 0.5) / num_lumps
        rho = math.sqrt(max(0.0, 1.0 - uz * uz))
        phi_k = golden * k
        ux = rho * math.cos(phi_k)
        uy = rho * math.sin(phi_k)
        mag = math.sqrt(ux * ux + uy * uy + uz * uz)
        if mag < 1e-12:
            world = (0.0, 0.0, 1.0)
        else:
            world = (ux / mag, uy / mag, uz / mag)

        # SH-modulated amplitude.
        _, theta_k, phi_ang = cartesian_to_spherical(*world)
        mod = eval_sh_modulation(sh_coeffs, theta_k, phi_ang, ell_max)
        amp = max(0.0, min(0.35, amp0 * mod))

        r_k = radius + thickness * uz
        center = (world[0] * r_k, world[1] * r_k, world[2] * r_k)

        # Velocity decomposition (toroidal/poloidal/radial in local frame).
        velocity = _shell_lump_velocity(
            axis, world, v_rad=v_rad, v_tor=v_tor, v_pol=v_pol
        )

        lumps.append({
            "amp": amp,
            "width": width0,
            "center": center,
            "velocity": velocity,
            "omega": omega,
            "mode": 0,
            "exotic": 1 if k in exotic_ids else 0,
        })
    return lumps


def _shell_lump_velocity(
    axis: tuple[float, float, float],
    world: tuple[float, float, float],
    *,
    v_rad: float,
    v_tor: float,
    v_pol: float,
) -> tuple[float, float, float]:
    """Decompose shell currents into radial / toroidal / poloidal.

    Thin wrapper so SH expansion can reuse the shell velocity decomposition
    without importing geometry directly (avoids circular deps).
    """
    from .geometry import _shell_lump_velocity as _slv

    return _slv(axis, world, v_rad=v_rad, v_tor=v_tor, v_pol=v_pol)


def _expand_trajectory_lumps_from_overrides(
    overrides: Mapping[str, Any],
    get_float: Any,
) -> list[dict]:
    """Generate GRTresna initial lumps at the t=0 per-lump trajectory positions.

    Each lump has its own orbit defined by ``trajectory_lump{k}_*`` keys.
    This function evaluates the per-lump trajectory at t=0 (pure Python,
    mirrors the C++ TrajectoryEvaluator) and creates lumps at those positions.
    GRTresna solves constraints for this static t=0 snapshot; the C++ side
    then drives the lumps along their trajectories during evolution.
    """
    from .spaces import TRAJECTORY_DEFAULT_NUM_LUMPS

    num_lumps = max(1, int(round(get_float("trajectory_num_lumps", float(TRAJECTORY_DEFAULT_NUM_LUMPS)))))
    well_width = get_float("trajectory_well_width", 1.5)

    lumps: list[dict] = []
    for k in range(num_lumps):
        pfx = f"trajectory_lump{k}_"
        R0 = get_float(f"{pfx}R0", 5.0)
        phase0 = get_float(f"{pfx}phase0", 2.0 * math.pi * k / num_lumps)
        tilt_theta = get_float(f"{pfx}tilt_theta", 0.0)
        tilt_phi = get_float(f"{pfx}tilt_phi", 0.0)
        well_depth = get_float(f"{pfx}well_depth", 0.05)

        # Position at t=0: orbit in plane, then rotate by (tilt_theta, tilt_phi).
        x_orb = R0 * math.cos(phase0)
        y_orb = R0 * math.sin(phase0)

        ct = math.cos(tilt_theta)
        st = math.sin(tilt_theta)
        cp = math.cos(tilt_phi)
        sp = math.sin(tilt_phi)

        cx = cp * ct * x_orb - sp * y_orb
        cy = sp * ct * x_orb + cp * y_orb
        cz = -st * x_orb

        exotic = int(round(get_float(f"{pfx}exotic", 0.0)))

        lumps.append({
            "amp": min(0.15, max(0.0, well_depth)),
            "width": max(1.5, well_width),
            "center": (cx, cy, cz),
            "velocity": (0.0, 0.0, 0.0),
            "omega": 0.0,
            "mode": 0,
            "exotic": exotic,
        })
    return lumps


def _expand_trajectory_boson_lumps_from_overrides(
    overrides: Mapping[str, Any],
    get_float: Any,
) -> list[dict]:
    """Generate GRTresna boson-star lumps at t=0 trajectory positions.

    Like ``_expand_trajectory_lumps_from_overrides`` but each lump is a compact
    boson star soliton.  If ``grtresna_boost_lumps`` is enabled (default True
    when trajectory_mode=1), each lump receives an initial velocity matching its
    trajectory angular velocity at t=0, capped at ``grtresna_boost_v_max``
    (default 0.8c) to stay sub-luminal on the flat initial slice.

    The Lorentz-boosted initial data lets the Q-ball start already moving along
    its orbit so the PD trap only needs to *maintain* the trajectory rather than
    *create* the entire motion from rest.

    The pump amplitude (well_depth) is capped at 0.15, matching the real-scalar
    trajectory cap.  The search space bounds can restrict this further (e.g. the
    ``discovery`` profile defaults to [0.001, 0.02] for corrective-only pumps)
    but the expansion itself does not impose a tighter limit — replay of
    real-scalar elites with bosonic matter needs the original amplitude range.
    """
    from ...grtresna.profiles.boson_star import (
        PROFILE_ODE_BOUND,
        PROFILE_SECH_BOUND,
        PROFILE_SELFGRAV_BOUND,
        bound_width,
    )
    from ...grtresna.profiles.qball_couplings import QBallCouplings
    from .spaces import TRAJECTORY_DEFAULT_NUM_LUMPS

    num_lumps = max(1, int(round(get_float("trajectory_num_lumps", float(TRAJECTORY_DEFAULT_NUM_LUMPS)))))

    # Bound boson lumps: width = 1/sqrt(m^2 - omega^2), the physical decay scale.
    mass = get_float("grtresna_scalar_mass", 0.1)
    lam = get_float("grtresna_scalar_lambda", 0.0)
    mu = get_float("grtresna_scalar_mu", 0.0)
    omega = get_float("grtresna_bs_omega", 0.0)
    width = bound_width(mass, omega)

    couplings = QBallCouplings(mass=mass, lam=lam, mu=mu, omega=omega)
    use_equilibrium = bool(int(round(get_float("grtresna_qball_equilibrium_amplitude", 0.0))))
    use_ode = bool(int(round(get_float("grtresna_qball_ode_profile", 0.0))))
    # Both painters rescale the ODE table by amp/phi_c, so a stationary seed
    # requires amp == the table's own phi_c.  cap_well_depth instead clamps to
    # the thin-wall estimate sqrt(3*lam/4*mu), ~4.5% under phi_c at the Bondi
    # couplings -- an off-shell seed that breathes and sheds.  Opt-in so existing
    # campaigns keep bit-identical initial data.
    use_exact_amp = bool(int(round(get_float("grtresna_qball_exact_amplitude", 0.0))))
    # Self-gravitating boson star: gravity binds the lump, so its SEED is a true
    # equilibrium (profile 4).  This replaces only the dispersing sech seed -- the
    # closed-loop PD trap pump still transports the lump along its trajectory.  The
    # star's central amplitude (which fixes its ADM mass via the ODE eigenvalue) is
    # decoupled from the pump well_depth: it comes from grtresna_bs_phi_c, while
    # well_depth stays a searchable transport-corrective amplitude.
    use_selfgrav = bool(int(round(get_float("grtresna_bs_selfgrav", 0.0))))
    bs_phi_c = get_float("grtresna_bs_phi_c", 0.08)
    if use_selfgrav:
        profile = int(PROFILE_SELFGRAV_BOUND)
    elif use_ode:
        profile = int(PROFILE_ODE_BOUND)
    else:
        profile = int(PROFILE_SECH_BOUND)

    # Velocity boost: match trajectory angular velocity at t=0.
    boost_lumps = bool(int(round(get_float("grtresna_boost_lumps", 1.0))))
    boost_v_max = get_float("grtresna_boost_v_max", 0.8)

    lumps: list[dict] = []
    for k in range(num_lumps):
        pfx = f"trajectory_lump{k}_"
        R0 = get_float(f"{pfx}R0", 5.0)
        phase0 = get_float(f"{pfx}phase0", 2.0 * math.pi * k / num_lumps)
        tilt_theta = get_float(f"{pfx}tilt_theta", 0.0)
        tilt_phi = get_float(f"{pfx}tilt_phi", 0.0)
        omega_rot = get_float(f"{pfx}omega_rot", 0.0)
        well_depth = get_float(f"{pfx}well_depth", 0.005)

        # Tangential speed of the orbit centre the co-moving trap must track.
        # v_t = R0 * |omega_rot| (geometric units, c = 1).  Should already be
        # clamped sub-luminal by _clamp_trajectory_speed; warn loudly if a lump
        # still asks for a relativistic/superluminal trajectory the soliton
        # cannot follow (the dominant dispersal mechanism).
        v_t = abs(omega_rot) * R0
        if v_t > 1.0:
            logger.warning(
                "trajectory lump %d: v_t = R0*|omega_rot| = %.3fc is SUPERLUMINAL "
                "-- the co-moving trap target outruns any soliton; expect dispersal. "
                "Cap omega_rot via trajectory_v_max.",
                k, v_t,
            )
        elif v_t > 0.5:
            logger.info(
                "trajectory lump %d: v_t = %.3fc (relativistic; tight curves will radiate)",
                k, v_t,
            )
        else:
            logger.debug("trajectory lump %d: v_t = %.3fc", k, v_t)

        # Position at t=0: orbit in plane, then rotate by (tilt_theta, tilt_phi).
        x_orb = R0 * math.cos(phase0)
        y_orb = R0 * math.sin(phase0)

        ct = math.cos(tilt_theta)
        st = math.sin(tilt_theta)
        cp = math.cos(tilt_phi)
        sp = math.sin(tilt_phi)

        cx = cp * ct * x_orb - sp * y_orb
        cy = sp * ct * x_orb + cp * y_orb
        cz = -st * x_orb

        # Velocity at t=0: derivative of the spiral orbit.
        # Spiral: r(t) = R0 + v_rad * t, phi(t) = phase0 + omega_rot * t.
        # vx_orb = v_rad * cos(phase0) - R0 * omega_rot * sin(phase0)
        # vy_orb = v_rad * sin(phase0) + R0 * omega_rot * cos(phase0)
        v_rad = get_float(f"{pfx}v_rad", 0.0)
        vx_orb = v_rad * math.cos(phase0) - R0 * omega_rot * math.sin(phase0)
        vy_orb = v_rad * math.sin(phase0) + R0 * omega_rot * math.cos(phase0)

        # Same rotation matrix as position.
        vx = cp * ct * vx_orb - sp * vy_orb
        vy = sp * ct * vx_orb + cp * vy_orb
        vz = -st * vx_orb

        # Cap at v_max to stay sub-luminal on the flat initial slice.
        if boost_lumps and boost_v_max > 0.0:
            v_mag = math.sqrt(vx * vx + vy * vy + vz * vz)
            if v_mag > boost_v_max:
                scale = boost_v_max / v_mag
                vx *= scale
                vy *= scale
                vz *= scale
            velocity = (vx, vy, vz)
        else:
            velocity = (0.0, 0.0, 0.0)

        exotic = int(round(get_float(f"{pfx}exotic", 0.0)))
        if use_selfgrav and exotic and not (lam > 0.0 and mu > 0.0 and omega > 0.0):
            # A gravity-BOUND phantom star cannot exist (its self-gravity is
            # repulsive), so the mini-star path (lam=mu=0, phi_c-parameterized)
            # has nothing to seed.  With sextic couplings + a target frequency
            # the binding is the scalar interaction and the dressed phantom
            # star DOES exist (solved with gravity_sign=-1); let it through.
            logger.warning(
                "trajectory lump %d: exotic flag ignored in self-gravitating "
                "mini-star mode (no gravity-bound phantom equilibrium); using "
                "canonical matter.  Sextic couplings + bs_omega enable it.",
                k,
            )
            exotic = 0

        if use_selfgrav:
            # Initial-data central amplitude sets the star on its mass-radius
            # branch via the gravitational eigenvalue; it is NOT the pump depth.
            amp = bs_phi_c
        elif use_ode and use_exact_amp and lam > 0.0 and mu > 0.0:
            from ...grtresna.profiles.qball_ode import cached_qball_radial_profile

            amp = float(cached_qball_radial_profile(mass, lam, mu, omega).phi_c)
        elif use_equilibrium and lam > 0.0 and mu > 0.0:
            amp = couplings.cap_well_depth(well_depth)
        else:
            amp = min(0.15, max(0.0, well_depth))

        lump: dict[str, Any] = {
            "amp": amp,
            "width": width,
            "center": (cx, cy, cz),
            "velocity": velocity,
            "omega": 0.0,
            "mode": 0,
            "profile": profile,
            "exotic": exotic,
        }
        if (use_ode or use_selfgrav) and lam > 0.0 and mu > 0.0:
            lump.update(
                {
                    "qball_mass": mass,
                    "qball_lam": lam,
                    "qball_mu": mu,
                    "qball_omega": omega,
                }
            )
        lumps.append(lump)
    return lumps


def _is_boson_sector(overrides: Mapping[str, Any], matter_model: str) -> bool:
    from ...grtresna.matter.models import (
        MATTER_SECTOR_BOSON_STAR,
        normalize_matter_sector,
        uses_boson_solve,
    )

    if uses_boson_solve(matter_model):
        return True
    if any(str(k).startswith("grtresna_bs_") for k in overrides):
        return True
    sector = str(overrides.get("grtresna_matter_sector", "")).strip()
    if not sector:
        return False
    try:
        return normalize_matter_sector(sector) == MATTER_SECTOR_BOSON_STAR
    except ValueError:
        return False


def build_grtresna_config(
    overrides: Mapping[str, Any], base: "GRTresnaConfig | None" = None
) -> "GRTresnaConfig":
    """Build a GRTresnaConfig from a candidate's ``grtresna_*`` overrides.

    ``base`` provides the fixed solver settings (grid, MPI ranks, BH content,
    cleanup); the searched ``grtresna_lump{k}_*`` keys are assembled into the
    ``lumps`` basis.  The legacy un-indexed ``grtresna_lump_*`` keys are also
    honoured (assembled as a single lump) for backward compatibility.

    Boson-star (complex scalar) campaigns set ``grtresna_matter_model`` in
    base overrides or search ``grtresna_bs_*`` / scalar-mass knobs directly.
    """
    import dataclasses

    from ...grtresna.matter.models import (
        apply_boson_star_overrides,
        finalize_independent_scalar_config,
    )
    from ...grtresna.solver import GRTresnaConfig, apply_exotic_safe_solver

    cfg = dataclasses.replace(base) if base is not None else GRTresnaConfig()

    def _enable_exotic_safe_solver() -> None:
        apply_exotic_safe_solver(cfg)

    def _gate_signs(out: "GRTresnaConfig") -> "GRTresnaConfig":
        # Fail-closed ID ↔ evolution gravitational-sign consistency (eval-118
        # bug class). Opt out only with GRTRESNA_ALLOW_SIGN_MISMATCH=1.
        from ...grtresna.matter.sign_consistency import (
            enforce_id_evolution_sign_consistency,
        )

        enforce_id_evolution_sign_consistency(out)
        return out

    def _return_scalar() -> "GRTresnaConfig":
        return _gate_signs(
            finalize_independent_scalar_config(
                cfg,
                overrides,
                enable_exotic_safe_solver=_enable_exotic_safe_solver,
            )
        )

    def _get_float(key: str, default: float) -> float:
        try:
            return float(overrides.get(key, default))
        except (TypeError, ValueError):
            return default

    # Matter-model tag from campaign base overrides (set by grtresna_context).
    matter_model = str(overrides.get("grtresna_matter_model", "")).strip()
    is_boson = _is_boson_sector(overrides, matter_model)

    if is_boson:
        allow_exotic = float(overrides.get("grtresna_boson_allow_exotic", 0.0)) >= 1.0
        has_trajectory = any(str(k).startswith("trajectory_") for k in overrides)
        if has_trajectory:
            # Trajectory + boson star: compact solitons on independent orbits.
            # Lumps get bulk velocity from omega_rot for physical momentum.
            cfg.lumps = _expand_trajectory_boson_lumps_from_overrides(
                overrides, _get_float
            )
        elif any(str(k).startswith("grtresna_shell_") for k in overrides):
            cfg.lumps = _expand_shell_lumps_from_overrides(
                overrides, _get_float, canonical_boson=not allow_exotic
            )
        # Read boson physics params (scalar_mass, bs_omega, etc.) and set
        # matter_model / example to complex scalar / BosonStarBH.
        apply_boson_star_overrides(
            cfg,
            overrides,
            enable_exotic_safe_solver=_enable_exotic_safe_solver,
        )
        if cfg.lumps and any(int(lump.get("exotic", 0)) for lump in cfg.lumps):
            _enable_exotic_safe_solver()
        return _gate_signs(cfg)

    if apply_boson_star_overrides(
        cfg,
        overrides,
        enable_exotic_safe_solver=_enable_exotic_safe_solver,
    ):
        return _gate_signs(cfg)

    # Global (ansatz-independent) initial-shift seed: aligns a non-zero beta^i
    # with the matter momentum flux so the warp/channel mechanism is reachable.
    cfg.shift_seed = _get_float("grtresna_shift_seed", cfg.shift_seed)

    # Scalar field mass (matter "weight"), ansatz-independent. Heavier mass
    # binds the lumps (smaller Compton wavelength 1/m) so they persist instead
    # of dispersing/flying away. Propagated to both the GRTresna constraint
    # solve and the GRTeclyn evolution, keeping the two T_ab identical.
    cfg.scalar_mass = _get_float("grtresna_scalar_mass", cfg.scalar_mass)
    cfg.scalar_lambda = _get_float("grtresna_scalar_lambda", cfg.scalar_lambda)
    # Sextic self-interaction (V += mu/6 |Phi|^6). With lambda>0, mu>0 the lump
    # is a stable Q-ball that self-binds, so it stays confined instead of
    # dispersing like a free field. Evolution-only (GRTeclyn); the GRTresna solve
    # stays quartic since gravity is negligible at these amplitudes.
    cfg.scalar_mu = _get_float("grtresna_scalar_mu", getattr(cfg, "scalar_mu", 0.0))

    # Trajectory-guided geometry survey: lumps at t=0 positions from trajectory
    # equations.  The trajectory_* keys flow through to GRTeclyn's params.txt
    # (they are NOT grtresna_ prefixed) and the C++ side drives the pump during
    # evolution.  GRTresna just solves constraints for the static t=0 snapshot.
    if any(str(k).startswith("trajectory_") for k in overrides):
        cfg.lumps = _expand_trajectory_lumps_from_overrides(overrides, _get_float)
        # Enable exotic-safe solver when ANY lump is exotic (per-lump flag).
        # The global coupling check in _return_scalar handles the case where
        # ALL lumps are forced exotic; this handles mixed canonical+exotic.
        if any(int(lump.get("exotic", 0)) for lump in cfg.lumps):
            _enable_exotic_safe_solver()
        return _return_scalar()

    if any(str(k).startswith("grtresna_ring_") for k in overrides):
        num_lumps = max(3, int(round(_get_float("grtresna_ring_lumps", GRTRESNA_DEFAULT_NUM_LUMPS))))
        amp0 = _get_float("grtresna_ring_amp", 0.15)
        width0 = _get_float("grtresna_ring_width", 3.0)
        radius = _get_float("grtresna_ring_radius", 3.5)
        z_scale = _get_float("grtresna_ring_z_scale", 0.4)
        phase = _get_float("grtresna_ring_phase", 0.0)
        v_tan = _get_float("grtresna_ring_tangential_velocity", 0.35)
        v_rad = _get_float("grtresna_ring_radial_velocity", 0.0)
        v_z = _get_float("grtresna_ring_vertical_velocity", 0.0)
        omega = _get_float("grtresna_ring_omega", 0.0)
        exotic_fraction = min(1.0, max(0.0, _get_float("grtresna_ring_exotic_fraction", 0.4)))
        exotic_phase = _get_float("grtresna_ring_exotic_phase", 0.0)
        dipole = _get_float("grtresna_ring_dipole_amp", 0.0)
        quadrupole = _get_float("grtresna_ring_quadrupole_amp", 0.0)
        mode = int(round(_get_float("grtresna_ring_mode", 1.0)))

        angles = [phase + 2.0 * math.pi * k / num_lumps for k in range(num_lumps)]
        exotic_count = int(round(exotic_fraction * num_lumps))
        ranked = sorted(
            range(num_lumps),
            key=lambda k: math.cos(angles[k] - exotic_phase),
            reverse=True,
        )
        exotic_ids = set(ranked[:exotic_count])

        lumps: list[dict] = []
        for k, theta in enumerate(angles):
            c = math.cos(theta)
            s = math.sin(theta)
            amp_mod = 1.0 + dipole * c + quadrupole * math.cos(2.0 * theta)
            width_mod = 1.0 + 0.25 * quadrupole * math.sin(2.0 * theta)
            amp = max(0.0, min(0.35, amp0 * max(0.15, amp_mod)))
            width = max(1.0, min(6.0, width0 * max(0.5, width_mod)))
            center = (radius * c, radius * s, z_scale * math.sin(theta - phase))
            r_hat = (c, s)
            t_hat = (-s, c)
            velocity = (
                v_rad * r_hat[0] + v_tan * t_hat[0],
                v_rad * r_hat[1] + v_tan * t_hat[1],
                v_z * math.cos(theta - phase),
            )
            lumps.append({
                "amp": amp,
                "width": width,
                "center": center,
                "velocity": velocity,
                "omega": omega,
                "mode": max(0, mode),
                "exotic": 1 if k in exotic_ids else 0,
            })
        cfg.lumps = lumps
        if any(lump["exotic"] for lump in lumps):
            _enable_exotic_safe_solver()
        return _return_scalar()

    if any(str(k).startswith("grtresna_sh_") for k in overrides):
        cfg.lumps = _expand_sh_lumps_from_overrides(overrides, _get_float)
        if any(lump["exotic"] for lump in cfg.lumps):
            _enable_exotic_safe_solver()
        return _return_scalar()

    if any(str(k).startswith("grtresna_shell_") for k in overrides):
        cfg.lumps = _expand_shell_lumps_from_overrides(overrides, _get_float)
        if any(lump["exotic"] for lump in cfg.lumps):
            _enable_exotic_safe_solver()
        return _return_scalar()

    # Group indexed lump keys by index.
    by_index: dict[int, dict[str, float]] = {}
    for key, val in overrides.items():
        m = _LUMP_KEY_RE.match(str(key))
        if m:
            by_index.setdefault(int(m.group(1)), {})[m.group(2)] = float(val)

    if by_index:
        lumps: list[dict] = []
        for k in sorted(by_index):
            f = by_index[k]
            lumps.append({
                "amp": f.get("amp", 0.0),
                "width": f.get("width", 5.0),
                "center": (f.get("center_x", 0.0), f.get("center_y", 0.0),
                           f.get("center_z", 0.0)),
                "velocity": (f.get("velocity_x", 0.0), f.get("velocity_y", 0.0),
                             f.get("velocity_z", 0.0)),
                "omega": f.get("omega", 0.0),
                "mode": int(round(f.get("mode", 0.0))),
                "exotic": int(round(f.get("exotic", 0.0))),
                "profile": int(round(f.get("profile", 0.0))),
            })
        cfg.lumps = lumps
        if any(lump["exotic"] for lump in lumps):
            _enable_exotic_safe_solver()
        return _return_scalar()

    # Backward-compatible single (un-indexed) lump.
    def _get(key: str, default: float) -> float:
        return float(overrides.get(key, default))

    if any(str(k).startswith("grtresna_lump_") for k in overrides):
        cfg.lump_amp = _get("grtresna_lump_amp", cfg.lump_amp)
        cfg.lump_width = _get("grtresna_lump_width", cfg.lump_width)
        cfg.lump_velocity = (
            _get("grtresna_lump_velocity_x", cfg.lump_velocity[0]),
            _get("grtresna_lump_velocity_y", cfg.lump_velocity[1]),
            _get("grtresna_lump_velocity_z", cfg.lump_velocity[2]),
        )
        cfg.lump_omega = _get("grtresna_lump_omega", cfg.lump_omega)
        if "grtresna_lump_mode" in overrides:
            cfg.lump_mode = int(round(float(overrides["grtresna_lump_mode"])))
        if "grtresna_lump_exotic" in overrides:
            cfg.lump_exotic = int(round(float(overrides["grtresna_lump_exotic"])))
    if cfg.lump_exotic:
        _enable_exotic_safe_solver()
    return _return_scalar()


def parse_convergence_safe(work_dir: Path) -> dict[str, float] | None:
    """Best-effort read of GRTresna Ham/Mom convergence (never raises)."""
    try:
        from ...grtresna.solver import parse_convergence

        return parse_convergence(work_dir)
    except Exception:
        return None
