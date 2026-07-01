"""Search-space definitions for CMA-ES optimization."""

from __future__ import annotations

import math
from typing import Any

from .dimension import SearchDimension

DEFAULT_SEARCH_SPACE: list[SearchDimension] = [
    SearchDimension("recipe_chi_coeff_0", -0.5, 0.1, -0.1),
    SearchDimension("recipe_chi_coeff_1", -0.3, 0.3, 0.0),
    SearchDimension("recipe_chi_coeff_2", -0.2, 0.2, 0.0),
    SearchDimension("recipe_chi_coeff_3", -0.2, 0.2, 0.0),
    SearchDimension("recipe_basis_width", 0.3, 3.0, 1.0),
    SearchDimension("recipe_alpha_coeff_0", -0.3, 0.3, 0.0),
    SearchDimension("recipe_alpha_coeff_1", -0.3, 0.3, 0.0),
    SearchDimension("recipe_beta_coeff_0", -0.5, 0.5, 0.0),
    SearchDimension("recipe_beta_coeff_1", -0.5, 0.5, 0.0),
]


# --- Non-spherical (angular) search ----------------------------------------
#
# The radial recipe above is spherically symmetric, so the optimizer can only
# ever discover breathing shells / radial channels.  To open up *directional*
# FTL geometries we add axisymmetric Legendre angular modes -- but ONLY on the
# lapse (alpha) and shift (beta) fields.
#
# This is a deliberate, constraint-preserving choice:
#   * alpha and beta are pure GAUGE: they do not appear in the t=0 Hamiltonian
#     or momentum constraints, so making them aspherical adds NO exotic matter
#     and does NOT break the 1D constrained (phi-from-chi) solve.
#   * Angular structure on the shift tilts/reshapes the local light cones (the
#     same mechanism Alcubierre uses), creating a *directional* operational-FTL
#     channel that the search can sculpt.
#   * Angular modes on chi or K, by contrast, WOULD violate the spherical 1D
#     constraint solve and surface as constraint error / forced exotic content,
#     so they are intentionally left out of this stage.
#
# Only the mode amplitudes are searched; the angular order (ell), radial center
# and width are fixed constants (ANGULAR_BASE_OVERRIDES) to keep the dimension
# count low.  ell=1 gives a fore/aft dipole (a genuine "direction"), ell=2 a
# quadrupolar pinch.
# Amplitude AND radial placement (center rc, width rw) of every angular mode
# are searched, so the optimizer can decide not just how strong each lobe is
# but *where along the radius* the directional deformation sits and how sharp
# it is -- a much richer directional-channel family than fixed-placement modes.
ANGULAR_SEARCH_SPACE: list[SearchDimension] = [
    # lapse mode 0 (ell=1 dipole)
    SearchDimension("recipe_lapse_mode_amp_0", -0.25, 0.25, 0.0),
    SearchDimension("recipe_lapse_mode_rc_0", 0.5, 6.0, 2.5),
    SearchDimension("recipe_lapse_mode_rw_0", 0.8, 4.0, 2.0),
    # lapse mode 1 (ell=2 quadrupole)
    SearchDimension("recipe_lapse_mode_amp_1", -0.25, 0.25, 0.0),
    SearchDimension("recipe_lapse_mode_rc_1", 0.5, 6.0, 2.5),
    SearchDimension("recipe_lapse_mode_rw_1", 0.8, 4.0, 2.0),
    # shift mode 0 (ell=1 dipole)
    SearchDimension("recipe_beta_mode_amp_0", -0.6, 0.6, 0.0),
    SearchDimension("recipe_beta_mode_rc_0", 0.5, 6.0, 2.5),
    SearchDimension("recipe_beta_mode_rw_0", 0.8, 4.0, 2.0),
    # shift mode 1 (ell=2 quadrupole)
    SearchDimension("recipe_beta_mode_amp_1", -0.6, 0.6, 0.0),
    SearchDimension("recipe_beta_mode_rc_1", 0.5, 6.0, 2.5),
    SearchDimension("recipe_beta_mode_rw_1", 0.8, 4.0, 2.0),
]

# Fixed (non-searched) parameters that activate the angular modes above: the
# mode counts and the angular order (ell).  Radial center/width are searched
# (see ANGULAR_SEARCH_SPACE) and therefore intentionally NOT pinned here.
ANGULAR_BASE_OVERRIDES: dict[str, Any] = {
    "recipe_num_lapse_angular_modes": 2,
    "recipe_lapse_mode_ell_0": 1,
    "recipe_lapse_mode_ell_1": 2,
    "recipe_num_beta_angular_modes": 2,
    "recipe_beta_mode_ell_0": 1,
    "recipe_beta_mode_ell_1": 2,
}


# --- GRTresna (constraint-satisfying, momentum-carrying matter) search -------
#
# When GRTresna is in the loop the initial data is produced by a full 3D
# elliptic constraint solve rather than the 1D radial recipe.  This unlocks a
# genuinely new FTL axis the recipe cannot represent: a scalar-field cloud that
# carries net linear and/or angular MOMENTUM, with the momentum constraint
# solved (S_i = -Pi d_i phi != 0 -> matter-momentum-driven frame dragging).
#
# These dimensions use the ``grtresna_`` prefix so the objective can route them
# to a GRTresnaConfig (which writes GRTresna's params.txt and runs the solver)
# instead of writing them into GRTeclyn's params.txt.  Everything downstream
# (GPU evolution + scoring) is unchanged: GRTresna just supplies the .gridinit.
#
# The matter source is a basis of K scalar "lumps"; each lump k contributes 11
# searched dimensions (grtresna_lump{k}_*).  A superposition paints an arbitrary
# energy/momentum distribution, which the elliptic solve maps to a rich family
# of constraint-satisfying geometries (chi wells + momentum-driven A_ij).
GRTRESNA_DEFAULT_NUM_LUMPS = 5


def grtresna_search_space(num_lumps: int = GRTRESNA_DEFAULT_NUM_LUMPS) -> list[SearchDimension]:
    """Build the K-lump GRTresna matter search space.

    Initial lump centres are staggered near the origin so the lumps start
    distinct without painting broad matter across the entire visual frame.
    """
    dims: list[SearchDimension] = []
    for k in range(num_lumps):
        # staggered initial x-centre, symmetric about 0
        cx0 = (k - (num_lumps - 1) / 2.0) * 2.5
        dims += [
            # The FTL channel here is shift-driven frame dragging: lump
            # velocity/omega -> conjugate momentum Pi -> momentum density S_i ->
            # A_ij -> evolved shift, and the cones only tilt superluminal once
            # the momentum density is strong.  The right lever for that is
            # VELOCITY/OMEGA, which we widen aggressively -- NOT amplitude.  The
            # elliptic Hamiltonian solve only converges for modest matter
            # amplitude (the prior search lived at amp~0.1-0.15 with clean
            # constraints; amp>~0.4 degrades to a ~99% Hamiltonian residual =
            # unphysical initial data), so amplitude is kept close to that proven
            # range while the momentum knobs are opened up.
            SearchDimension(f"grtresna_lump{k}_amp", 0.0, 0.35, 0.15),
            SearchDimension(f"grtresna_lump{k}_width", 1.5, 4.5, 3.0),
            SearchDimension(f"grtresna_lump{k}_center_x", -6.0, 6.0, cx0),
            SearchDimension(f"grtresna_lump{k}_center_y", -6.0, 6.0, 0.0),
            SearchDimension(f"grtresna_lump{k}_center_z", -3.0, 3.0, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_x", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_y", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_velocity_z", -0.9, 0.9, 0.0),
            SearchDimension(f"grtresna_lump{k}_omega", -0.6, 0.6, 0.0),
            # azimuthal mode m (rounded to int): 0 axisymmetric, >=1 enables L_z
            SearchDimension(f"grtresna_lump{k}_mode", 0.0, 2.0, 1.0),
            # phantom/ghost flag (rounded to int): 0 canonical (+rho), 1 EXOTIC
            # (-rho).  The search can flip individual lumps exotic to source the
            # NEC violation a persistent FTL channel needs, while keeping others
            # canonical.  Starts canonical (0.0) so the optimizer departs from
            # known-good positive matter and discovers exotic mixes from there.
            SearchDimension(f"grtresna_lump{k}_exotic", 0.0, 1.0, 0.0),
        ]
    return dims


def grtresna_ring_search_space() -> list[SearchDimension]:
    """Reduced rotating-ring GRTresna matter search space.

    This is a low-dimensional template for the same lump backend.  CMA-ES
    chooses global ring/momentum/exotic-pattern knobs; build_grtresna_config()
    expands them into K concrete ``cfg.lumps`` entries.  It keeps the physically
    relevant levers (exotic support, angular momentum, counterflow/shear) while
    reducing the optimizer from 5 * 11 = 55 dimensions to 14.
    """
    return [
        SearchDimension("grtresna_ring_amp", 0.04, 0.30, 0.15),
        SearchDimension("grtresna_ring_width", 1.5, 4.5, 3.0),
        SearchDimension("grtresna_ring_radius", 1.5, 6.0, 3.5),
        SearchDimension("grtresna_ring_z_scale", 0.0, 2.5, 0.4),
        SearchDimension("grtresna_ring_phase", 0.0, 2.0 * math.pi, 0.0),
        SearchDimension("grtresna_ring_tangential_velocity", -0.9, 0.9, 0.35),
        SearchDimension("grtresna_ring_radial_velocity", -0.5, 0.5, 0.0),
        SearchDimension("grtresna_ring_vertical_velocity", -0.3, 0.3, 0.0),
        SearchDimension("grtresna_ring_omega", -0.5, 0.5, 0.0),
        SearchDimension("grtresna_ring_exotic_fraction", 0.0, 1.0, 0.4),
        SearchDimension("grtresna_ring_exotic_phase", 0.0, 2.0 * math.pi, 0.0),
        SearchDimension("grtresna_ring_dipole_amp", -0.6, 0.6, 0.0),
        SearchDimension("grtresna_ring_quadrupole_amp", -0.5, 0.5, 0.0),
        SearchDimension("grtresna_ring_mode", 0.0, 2.0, 1.0),
    ]


SHELL_PROFILE_CHOICES = ("middle", "compact", "outer_precursor", "inner_shift")


def grtresna_shell_search_space(profile: str = "compact") -> list[SearchDimension]:
    """Full-sphere ``shell`` GRTresna matter search space (discovery ansatz).

    The ``ring`` ansatz is a planar (equatorial) loop: it cannot place matter
    or currents anywhere on the 2-sphere, only on one circle plus a small
    out-of-plane wobble.  This ``shell`` ansatz lifts that restriction.  It
    distributes the ``K`` lumps over the *whole* unit sphere (a deterministic
    Fibonacci lattice), gives the configuration an arbitrary orientation axis
    ``(theta, phi)``, and decomposes the matter current into the two physically
    distinct families a sphere supports:

    * a *toroidal* current circulating about the axis (net angular momentum /
      gravitomagnetic dipole along the axis -- the warp "motor"), and
    * a *poloidal* current flowing over the poles (a current topology the
      equatorial ring simply does not contain).

    Amplitude/width are modulated by low Legendre orders along the axis
    (``dipole`` = ell=1 hemisphere asymmetry, ``quadrupole`` = ell=2), and the
    exotic sector is a longitudinal wedge selected about the axis.  The ring is
    recovered as the thin-equatorial-band limit, so this is a strict superset
    used for *discovery* (broad coverage) rather than *refinement* (the ring).

    18 searched dimensions for any ``K``: like the ring, the lump count is a
    mesh-resolution knob, not an optimizer dimension.

    ``profile`` selects preset bounds.  ``compact`` caps lump width to reduce
    diffuse-blob basins (optimize_20260604T170329Z eval 128 sat at width~3.9).
    """
    if profile not in SHELL_PROFILE_CHOICES:
        raise ValueError(
            f"unknown shell profile: {profile!r} "
            f"(choices: {', '.join(SHELL_PROFILE_CHOICES)})"
        )

    # Scalar momentum source ~ amp^2 * velocity / width.  Diffuse amp~0.15
    # shells converged but shift~3e-2; amp~0.25 width~1.5 v~3 blew up GRTresna.
    # Bounds tightened to the basin where the elliptic solve actually converges
    # (every prior success had amp<=0.145, thickness>0, |toroidal|<=~1.2,
    # |omega|<=~0.5); the old wide extremes accounted for the bulk of the ~82%
    # pre-GPU rejections.
    amp = (0.08, 0.16, 0.13)
    width = (1.8, 4.0, 2.4)
    radius = (1.5, 6.0, 3.5)
    thickness = (0.1, 2.5, 0.6)
    toroidal_v = (-1.2, 1.2, 0.6)
    exotic_frac = (0.0, 1.0, 0.4)
    # Scalar mass = matter "weight".  At m=0.1 the Compton wavelength 1/m=10 is
    # the box size, so a lump disperses (free-streams) even at zero velocity --
    # the "fly away" the campaign keeps hitting.  Searching m in [0.3, 1.5]
    # (1/m in [3.3, 0.67]) keeps the field bound within the shell width (~2.4)
    # so the search can settle on persistent, near-static matter.  Honoured
    # identically by the GRTresna solve and the GRTeclyn evolution.
    scalar_mass = (0.3, 1.5, 0.6)

    if profile == "compact":
        width = (1.8, 3.0, 2.4)
    elif profile == "outer_precursor":
        width = (2.8, 4.0, 3.5)
        radius = (4.0, 6.0, 5.0)
        toroidal_v = (-2.0, 2.0, -1.2)
        exotic_frac = (0.6, 1.0, 0.8)
    elif profile == "inner_shift":
        width = (1.8, 3.0, 2.4)
        radius = (1.5, 3.0, 2.0)
        toroidal_v = (-1.0, 1.0, 0.4)
        exotic_frac = (0.0, 0.6, 0.4)

    return [
        SearchDimension("grtresna_shell_amp", *amp),
        SearchDimension("grtresna_shell_width", *width),
        SearchDimension("grtresna_shell_radius", *radius),
        SearchDimension("grtresna_shell_thickness", *thickness),
        SearchDimension("grtresna_shell_axis_theta", 0.0, math.pi, 0.5 * math.pi),
        SearchDimension("grtresna_shell_axis_phi", 0.0, 2.0 * math.pi, 0.0),
        SearchDimension("grtresna_scalar_mass", *scalar_mass),
        # Self-interaction on the shared field sum: V = 1/2(m s)^2 - (lambda/4)s^4.
        # lambda=0 recovers the mass-only potential (v12 behaviour).
        SearchDimension("grtresna_scalar_lambda", 0.0, 0.1, 0.0),
        # Matter placement topology: 0=sphere, 1=channel, 2=bipolar, 3=ring,
        # 4=quasi-random cloud (deterministic asymmetric scatter).
        SearchDimension("grtresna_matter_layout", 0.0, 4.0, 0.0),
        # Toroidal current is the warp/frame-drag motor -- kept at full range.
        SearchDimension("grtresna_shell_toroidal_velocity", *toroidal_v),
        # Poloidal/radial currents carry the lumps' net outflow (the "fly away"),
        # not the warp drive, so they are capped tighter than before
        # (poloidal +-1.5 -> +-0.8, radial +-0.8 -> +-0.3) to keep matter bound.
        SearchDimension("grtresna_shell_poloidal_velocity", -0.8, 0.8, 0.0),
        SearchDimension("grtresna_shell_radial_velocity", -0.3, 0.3, 0.0),
        SearchDimension("grtresna_shell_omega", -0.5, 0.5, 0.2),
        SearchDimension("grtresna_shell_dipole_amp", -0.6, 0.6, 0.0),
        SearchDimension("grtresna_shell_quadrupole_amp", -0.5, 0.5, 0.0),
        SearchDimension("grtresna_shell_exotic_fraction", *exotic_frac),
        SearchDimension("grtresna_shell_exotic_phase", 0.0, 2.0 * math.pi, 0.0),
        # Per-lump matter PROFILE: a searchable fraction of the lumps get the
        # smoothed top-hat ("ball") profile (profile=1); the rest stay Gaussian
        # (profile=0). fraction=0 (default) keeps every lump Gaussian (v13).
        # The phase selects WHICH lumps (by azimuth) become top-hat, mirroring
        # the exotic_fraction/exotic_phase idiom.
        SearchDimension("grtresna_shell_profile_fraction", 0.0, 1.0, 0.0),
        SearchDimension("grtresna_shell_profile_phase", 0.0, 2.0 * math.pi, 0.0),
        SearchDimension("grtresna_shell_mode", 0.0, 2.0, 1.0),
        SearchDimension("grtresna_shell_radial_jitter", 0.0, 1.0, 0.0),
        # Initial-shift seed (peak |beta| along the matter momentum flux). This
        # unlocks the warp/channel mechanism that a zero-shift gauge cannot
        # reach; sign selects the drive direction. Default 0 = native gauge.
        SearchDimension("grtresna_shift_seed", -0.6, 0.6, 0.0),
        # STATIC-matter toggle (rounded to int): 0 = momentum-carrying matter
        # (toroidal/poloidal/radial currents + spin as searched above), 1 =
        # fully static matter (all lump velocities and omega forced to zero).
        # Lets the search test purely static lumps -- where the only available
        # FTL channel is the gauge/shift seed and the geometry itself, not
        # frame-dragging -- alongside the moving-matter family. Starts at 0 so
        # the optimizer departs from the established moving-matter behaviour.
        SearchDimension("grtresna_shell_static", 0.0, 1.0, 0.0),
    ]


def grtresna_boson_shell_search_space(
    profile: str = "compact",
    *,
    static: bool = True,
    allow_exotic: bool = False,
    matched_bounds: bool = False,
) -> list[SearchDimension]:
    """Canonical bosonic shell: shell geometry + boson-star convergence bounds.

    Multiple Gaussian sites superpose into one U(1) complex scalar solved by
    BosonStarBH (``num_lumps`` + ``lump{k}_*`` params).  Amplitudes and mass
    are capped tighter than the scalar shell so BosonStarBH stays inside the
    Ham/Mom gate.

    When ``static=True`` (default), velocity dims are excluded and the campaign
    pins ``grtresna_shell_static=1``.  When ``static=False``, conservative
    velocity bounds are included for the moving-matter frame-dragging splash.

    When ``allow_exotic=True`` (FTL wormhole / RL campaigns), ``shell_exotic_*``
    search dims are kept and lump exotic flags are wired into BosonStarBH paint.

    When ``matched_bounds=True`` (v2 fair comparison), mass/amplitude/lambda
    bounds are raised to match the scalar shell space and geometric dims
    (radius, thickness, layout, profile, dipole, quadrupole) are kept.  This
    allows a fair physics comparison — any FTL difference is attributable to
    matter type, not search handicap.
    """
    shell = grtresna_shell_search_space(profile=profile)
    # Dims always replaced with boson-specific versions.
    skip = {
        "grtresna_scalar_mass",
        "grtresna_scalar_lambda",
        "grtresna_shell_amp",
        "grtresna_shell_width",
    }
    if not matched_bounds:
        # v1 behaviour: also filter out geometric dims that the boson solver
        # historically didn't need.  Velocity dims filtered for static.
        skip |= {
            "grtresna_shell_toroidal_velocity",
            "grtresna_shell_poloidal_velocity",
            "grtresna_shell_radial_velocity",
            "grtresna_shell_omega",
        }
    else:
        # v2 matched: keep geometric dims from scalar shell but still replace
        # velocity dims when static (they get pinned anyway).
        if static:
            skip |= {
                "grtresna_shell_toroidal_velocity",
                "grtresna_shell_poloidal_velocity",
                "grtresna_shell_radial_velocity",
                "grtresna_shell_omega",
            }
    if not allow_exotic:
        skip |= {
            "grtresna_shell_exotic_fraction",
            "grtresna_shell_exotic_phase",
        }
    dims = [d for d in shell if d.param_key not in skip]

    if matched_bounds:
        # v2: bounds widened toward scalar shell for fairer comparison.
        # Mass [0.15, 0.5]: BosonStarBH convergence ceiling is ~0.5 at 64^3
        # resolution; v2-attempt-1 showed 84% grtresna_rejected at [0.3, 1.0].
        # The [0.15, 0.5] range keeps the Compton wavelength <= 6.7 (within
        # shell width) while staying inside the solver's convergence basin.
        # Amplitude [0.08, 0.15]: raised from v1's 0.12 cap; chi-tagger
        # blowup at 0.15 will postload_reject naturally.
        # bs_omega [0.1, 0.45]: must satisfy omega < mass for bound states.
        # Capped below mass lower bound (0.15) is fine; the upper 0.45 < 0.5
        # prevents unbound-field divergence that killed v2-attempt-1.
        # Lambda [0, 0.1]: matched to scalar.
        dims.extend([
            SearchDimension("grtresna_shell_amp", 0.08, 0.15, 0.12),
            SearchDimension("grtresna_shell_width", 1.8, 4.0, 2.4),
            SearchDimension("grtresna_scalar_mass", 0.15, 0.5, 0.25),
            SearchDimension("grtresna_scalar_lambda", 0.0, 0.1, 0.0),
            SearchDimension("grtresna_bs_omega", 0.1, 0.45, 0.2),
        ])
    else:
        # v1: conservative bounds for BosonStarBH stability.
        # Amplitude capped at 0.12: pushing to 0.15 deepens the chi well so
        # broadly that the GPU evolution's chi-tagger refines ~100% of the
        # domain (256^3 Level 1), which blows up RK4 storeRKCoarseData and
        # segfaults.  0.12 keeps the matter compact enough that refinement
        # stays on the central region.
        dims.extend([
            SearchDimension("grtresna_shell_amp", 0.04, 0.12, 0.08),
            SearchDimension("grtresna_shell_width", 2.0, 4.0, 3.0),
            SearchDimension("grtresna_scalar_mass", 0.05, 0.35, 0.1),
            SearchDimension("grtresna_scalar_lambda", 0.0, 0.05, 0.0),
            SearchDimension("grtresna_bs_omega", 0.05, 0.35, 0.15),
        ])

    if not allow_exotic:
        dims.append(SearchDimension("grtresna_scalar_sign", 1.0, 1.0, 1.0))
    if not static:
        dims.extend([
            # Spacetime-splash velocities.  RADIAL is the primary driver: a
            # negative v_rad sends the shell inward, sourcing the momentum
            # constraint with converging momentum -> an inward-propagating
            # gravitational wave that focuses into a curvature splash at the
            # center.  Bounds kept modest so the NL Mom solve stays convergent
            # at canonical resolution (fast shells need finer grids).
            SearchDimension("grtresna_shell_radial_velocity", -0.2, 0.2, -0.2),
            SearchDimension("grtresna_shell_toroidal_velocity", -0.2, 0.2, 0.0),
            SearchDimension("grtresna_shell_poloidal_velocity", -0.2, 0.2, 0.0),
            SearchDimension("grtresna_shell_omega", -0.15, 0.15, 0.0),
        ])
    return dims


def grtresna_boson_splash_search_space() -> list[SearchDimension]:
    """Boson-star splash search space with unpinned phase velocity ``omega``.

    Bounds are tighter than the generic boson-star space to keep canonical
    GRTresna solves in the Ham/Mom convergence window for early splash sweeps.
    Profile width defaults small (Gaussian σ ≈ 3.5 on a 64-wide domain) so
    frames show a compact core with room for inward wave structure.
    """
    return [
        SearchDimension("grtresna_scalar_mass", 0.05, 0.35, 0.1),
        SearchDimension("grtresna_scalar_lambda", 0.0, 0.05, 0.0),
        SearchDimension("grtresna_bs_phi_c", 0.03, 0.12, 0.06),
        SearchDimension("grtresna_bs_profile_width", 2.0, 8.0, 3.5),
        SearchDimension("grtresna_bs_omega", 0.05, 0.4, 0.2),
        SearchDimension("grtresna_scalar_sign", 1.0, 1.0, 1.0),
        SearchDimension("grtresna_shift_seed", -0.6, 0.6, 0.0),
    ]


def grtresna_boson_star_search_space() -> list[SearchDimension]:
    """Complex-scalar (mini boson star) search space.

    A single U(1)-charged Gaussian profile at the origin replaces the
    multi-lump real-scalar shell.  Matter persists under evolution because
    the phase velocity Pi2 = -omega phi / alpha conserves global charge,
    unlike dispersing independent real scalars (MapElites action plan step 3).
    """
    return [
        SearchDimension("grtresna_scalar_mass", 0.05, 0.35, 0.1),
        SearchDimension("grtresna_scalar_lambda", 0.0, 0.05, 0.0),
        SearchDimension("grtresna_bs_phi_c", 0.03, 0.15, 0.08),
        SearchDimension("grtresna_bs_profile_width", 4.0, 16.0, 8.0),
        SearchDimension("grtresna_bs_omega", 0.0, 0.0, 0.0),
        SearchDimension("grtresna_scalar_sign", 1.0, 1.0, 1.0),
        SearchDimension("grtresna_shift_seed", -0.6, 0.6, 0.0),
    ]


SH_DEFAULT_ELL_MAX = 4
"""Default maximum spherical-harmonic degree for the ``sh`` ansatz."""

SH_PROFILE_CHOICES = ("discovery", "conservative")


def grtresna_sh_search_space(
    ell_max: int = SH_DEFAULT_ELL_MAX,
    profile: str = "discovery",
) -> list[SearchDimension]:
    """Spherical-harmonic matter ansatz — angular freedom via SH expansion.

    Replaces the shell ansatz's fixed dipole+quadrupole (ell<=2, 2 DOF)
    amplitude modulation with a full real-SH expansion up to ``ell_max``,
    giving ``(ell_max+1)^2 - 1`` angular modulation coefficients.  Each
    coefficient scales a real spherical harmonic that modulates the per-lump
    amplitude around a searched base value, giving genuine per-lump
    variation through a smooth, low-dimensional parameterization.

    The lump positions remain a Fibonacci lattice on the sphere (as in the
    shell ansatz); only their amplitudes (and optionally radii) vary.

    Dimensions:
      - 1 base amplitude + ``(ell_max+1)^2 - 1`` SH modulation = angular DOF
      - Global geometry: radius, width, thickness
      - Physics: scalar_mass, scalar_lambda
      - Kinematics: toroidal/poloidal/radial velocity, omega, shift_seed
      - Exotic sector: exotic_fraction, exotic_phase
      - Static toggle

    For ell_max=4: 24 modulation + ~13 global = ~37 D.
    For ell_max=3: 15 modulation + ~13 global = ~28 D.
    """
    if profile not in SH_PROFILE_CHOICES:
        raise ValueError(
            f"unknown SH profile: {profile!r} "
            f"(choices: {', '.join(SH_PROFILE_CHOICES)})"
        )

    n_sh = (ell_max + 1) ** 2  # total SH modes including monopole
    # Modulation coefficients (skip monopole idx=0, absorbed into sh_amp).
    # With 24 coefficients and |Y_lm| up to ~0.85, the modulation stdev
    # is ~coeff_bound * sqrt(N_sh / 4pi).  We need 2-sigma excursions
    # to keep effective amp below ~0.20 so GRTresna converges (v1 used
    # +-0.8 which pushed amps to 0.35 → 100% Ham rejection).
    coeff_bound = 0.35 if profile == "discovery" else 0.25

    dims: list[SearchDimension] = []

    # Base amplitude (the monopole-equivalent).  Upper raised to 0.22
    # (from 0.16 in v22) together with looser postload gate (ham_l2
    # 0.03→0.05) to let stronger-matter configs reach the GPU where
    # they may actually curve light cones enough for f_geo > 0.
    dims.append(SearchDimension("grtresna_sh_amp", 0.04, 0.22, 0.12))

    # SH modulation coefficients c_1 .. c_{N-1}.
    for idx in range(1, n_sh):
        dims.append(
            SearchDimension(
                f"grtresna_sh_c{idx}", -coeff_bound, coeff_bound, 0.0
            )
        )

    # Shell geometry (shared with shell ansatz).
    dims.append(SearchDimension("grtresna_sh_radius", 1.5, 6.0, 3.5))
    dims.append(SearchDimension("grtresna_sh_width", 1.8, 3.5, 2.4))
    dims.append(SearchDimension("grtresna_sh_thickness", 0.0, 2.0, 0.5))

    # Scalar field physics.
    dims.append(SearchDimension("grtresna_scalar_mass", 0.3, 1.5, 0.6))
    dims.append(SearchDimension("grtresna_scalar_lambda", 0.0, 0.1, 0.0))

    # Kinematics — toroidal is the warp motor.
    if profile == "discovery":
        dims.append(
            SearchDimension("grtresna_sh_toroidal_velocity", -1.2, 1.2, 0.4)
        )
    else:
        dims.append(
            SearchDimension("grtresna_sh_toroidal_velocity", -0.6, 0.6, 0.3)
        )
    dims.append(
        SearchDimension("grtresna_sh_poloidal_velocity", -0.8, 0.8, 0.0)
    )
    dims.append(
        SearchDimension("grtresna_sh_radial_velocity", -0.3, 0.3, 0.0)
    )
    dims.append(SearchDimension("grtresna_sh_omega", -0.5, 0.5, 0.2))

    # Exotic sector.
    dims.append(
        SearchDimension("grtresna_sh_exotic_fraction", 0.0, 1.0, 0.4)
    )
    dims.append(
        SearchDimension("grtresna_sh_exotic_phase", 0.0, 2.0 * math.pi, 0.0)
    )

    # Initial-shift seed (peak |beta| along matter momentum flux).
    dims.append(SearchDimension("grtresna_shift_seed", -0.6, 0.6, 0.0))

    # Static toggle — starts at 0 (dynamics ON) unlike the shell default.
    dims.append(SearchDimension("grtresna_sh_static", 0.0, 1.0, 0.0))

    return dims


TRAJECTORY_DEFAULT_NUM_LUMPS = 5

TRAJECTORY_PROFILE_CHOICES = ("discovery", "rotation_only", "breathing_only")


def grtresna_trajectory_search_space(
    num_lumps: int = TRAJECTORY_DEFAULT_NUM_LUMPS,
    profile: str = "discovery",
) -> list[SearchDimension]:
    """Trajectory-guided FTL geometry survey — per-lump independent orbits.

    Each lump gets its own tilted orbit (circle or spiral) defined by 8 searched
    parameters:
      - R0: orbital radius at t=0
      - omega_rot: angular velocity (sign = direction; counter-rotating lumps
        create shear -> frame dragging, the primary FTL mechanism)
      - phase0: initial azimuthal position
      - tilt_theta: orbital plane tilt from z-axis
      - tilt_phi: orbital plane azimuth (full 3D orientation)
      - well_depth: per-lump pump amplitude
      - v_rad: radial drift speed (v_rad=0 recovers a circle; signed spiral)
      - exotic: per-lump canonical vs phantom matter flag

    Shared parameters control collective phenomena:
      - A_breath, omega_breath: radial pulsation (modulates all radii)
      - z_amp, omega_z: confined axial oscillation (always bounded)
      - well_width: Gaussian spotlight sigma

    The C++ evaluator runs on the CPU each coarse step — no GPU trig.
    GRTresna provides initial data from the t=0 lump positions.

    Dimensionality: 8 * num_lumps + 5 shared.
      - 5 lumps: 45 D
      - 3 lumps: 29 D

    ``profile`` controls which per-lump DOFs are searched:
      - ``discovery``: full per-lump freedom (R0, omega, phase, tilt, amp, v_rad)
      - ``rotation_only``: fix R0/tilt, search only omega_rot + phase0 + amp
      - ``breathing_only``: fix omega_rot/tilt, search R0 + amplitude
    """
    if profile not in TRAJECTORY_PROFILE_CHOICES:
        raise ValueError(
            f"unknown trajectory profile: {profile!r} "
            f"(choices: {', '.join(TRAJECTORY_PROFILE_CHOICES)})"
        )

    dims: list[SearchDimension] = []

    # --- Shared parameters ---
    dims.append(SearchDimension("trajectory_A_breath", 0.0, 2.0, 0.0))
    dims.append(SearchDimension("trajectory_omega_breath", 0.0, 3.0, 0.5))
    dims.append(SearchDimension("trajectory_z_amp", 0.0, 3.0, 0.0))
    dims.append(SearchDimension("trajectory_omega_z", 0.0, 2.0, 0.0))
    dims.append(SearchDimension("trajectory_well_width", 0.8, 3.0, 1.5))

    # --- Per-lump parameters ---
    # Stagger initial phases and radii so lumps start distinct.
    for k in range(num_lumps):
        phase0_default = 2.0 * math.pi * k / num_lumps
        R0_default = 4.0 + (k - (num_lumps - 1) / 2.0) * 0.5

        dims.append(
            SearchDimension(f"trajectory_lump{k}_R0", 1.5, 8.0, R0_default)
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_omega_rot", -1.0, 1.0, 0.0
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_phase0",
                0.0,
                2.0 * math.pi,
                phase0_default,
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_tilt_theta", 0.0, math.pi, 0.0
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_tilt_phi",
                0.0,
                2.0 * math.pi,
                0.0,
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_well_depth", 0.01, 0.15, 0.05
            )
        )
        # Radial drift for spiral orbits (signed: inward < 0, outward > 0).
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_v_rad", -0.3, 0.3, 0.0
            )
        )
        # Per-lump exotic flag: continuous [0,1], rounded to 0/1 in config
        # expansion.  Lets the search discover which lumps should carry
        # negative energy density (exotic) vs positive (canonical).
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_exotic", 0.0, 1.0, 0.0
            )
        )

    # Profile-specific simplifications.
    if profile == "rotation_only":
        # Fix R0, tilt — only rotation and amplitude searched per lump.
        remove = {
            f"trajectory_lump{k}_{p}"
            for k in range(num_lumps)
            for p in ("R0", "tilt_theta", "tilt_phi")
        }
        # Also remove breathing shared dims.
        remove |= {"trajectory_A_breath", "trajectory_omega_breath"}
        dims = [d for d in dims if d.param_key not in remove]
    elif profile == "breathing_only":
        # Fix rotation and tilt — only R0 and amplitude searched per lump.
        remove = {
            f"trajectory_lump{k}_{p}"
            for k in range(num_lumps)
            for p in ("omega_rot", "phase0", "tilt_theta", "tilt_phi")
        }
        dims = [d for d in dims if d.param_key not in remove]

    return dims


# ---------------------------------------------------------------------------
# Trajectory + boson star (compact soliton matter on trajectory orbits)
# ---------------------------------------------------------------------------

TRAJECTORY_BOSON_PROFILE_CHOICES = ("discovery", "rotation_only")


def grtresna_trajectory_boson_search_space(
    num_lumps: int = TRAJECTORY_DEFAULT_NUM_LUMPS,
    profile: str = "discovery",
) -> list[SearchDimension]:
    """Trajectory-guided FTL with **boson star** matter — non-dispersing solitons.

    Each lump is a compact U(1)-charged boson star placed on an independent
    tilted circular orbit.  Unlike the real-scalar trajectory (where the pump
    IS the matter source), here the boson stars are self-bound via charge
    conservation and the pump merely provides gentle orbit corrections.

    Physics:
      - Each boson star has internal phase velocity ``omega_bs`` (shared across
        lumps for now) satisfying ``omega_bs < scalar_mass`` for a bound state.
      - At t=0 each star gets bulk tangential velocity from ``omega_rot × r``
        so GRTresna solves the *full* momentum constraint (non-trivial).
      - Pump amplitude is reduced: ``well_depth ∈ [0.001, 0.02]`` (corrective,
        not generative).

    Search dimensions:
      Per-lump (8 each):
        R0, omega_rot, phase0, tilt_theta, tilt_phi, well_depth, v_rad, exotic
      Shared trajectory (5):
        A_breath, omega_breath, z_amp, omega_z, well_width
      Shared boson physics (3):
        grtresna_scalar_mass, grtresna_scalar_lambda, grtresna_bs_omega

    Total: 8*N + 8   (5 lumps → 48 D).
    """
    if profile not in TRAJECTORY_BOSON_PROFILE_CHOICES:
        raise ValueError(
            f"unknown trajectory-boson profile: {profile!r} "
            f"(choices: {', '.join(TRAJECTORY_BOSON_PROFILE_CHOICES)})"
        )

    dims: list[SearchDimension] = []

    # --- Shared trajectory parameters ---
    dims.append(SearchDimension("trajectory_A_breath", 0.0, 2.0, 0.0))
    dims.append(SearchDimension("trajectory_omega_breath", 0.0, 3.0, 0.5))
    dims.append(SearchDimension("trajectory_z_amp", 0.0, 3.0, 0.0))
    dims.append(SearchDimension("trajectory_omega_z", 0.0, 2.0, 0.0))
    dims.append(SearchDimension("trajectory_well_width", 0.8, 3.0, 1.5))

    # --- Shared boson star physics ---
    # scalar_mass: binding energy.  Compton wavelength 1/m should be comparable
    # to the Gaussian profile width so the star is self-bound.
    dims.append(SearchDimension("grtresna_scalar_mass", 0.05, 0.5, 0.2))
    # scalar_lambda: quartic self-interaction (attractive for lambda>0).
    dims.append(SearchDimension("grtresna_scalar_lambda", 0.0, 0.05, 0.0))
    # bs_omega: shared internal phase velocity (must satisfy omega < mass for
    # bound state; the config builder enforces this).
    dims.append(SearchDimension("grtresna_bs_omega", 0.05, 0.45, 0.15))

    # --- Per-lump trajectory parameters ---
    for k in range(num_lumps):
        phase0_default = 2.0 * math.pi * k / num_lumps
        R0_default = 4.0 + (k - (num_lumps - 1) / 2.0) * 0.5

        dims.append(
            SearchDimension(f"trajectory_lump{k}_R0", 1.5, 8.0, R0_default)
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_omega_rot", -1.0, 1.0, 0.0
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_phase0",
                0.0,
                2.0 * math.pi,
                phase0_default,
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_tilt_theta", 0.0, math.pi, 0.0
            )
        )
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_tilt_phi",
                0.0,
                2.0 * math.pi,
                0.0,
            )
        )
        # Corrective pump — much weaker than generative real-scalar well_depth.
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_well_depth", 0.001, 0.02, 0.005
            )
        )
        # Radial drift speed for spiral orbits.  v_rad > 0 spirals outward,
        # v_rad < 0 spirals inward.  Magnitude is clamped by the trajectory
        # speed cap so the pump target remains sub-luminal.
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_v_rad", -0.3, 0.3, 0.0
            )
        )
        # Per-lump exotic flag (0=canonical, 1=phantom).
        dims.append(
            SearchDimension(
                f"trajectory_lump{k}_exotic", 0.0, 1.0, 0.0
            )
        )

    # Profile simplifications.
    if profile == "rotation_only":
        remove = {
            f"trajectory_lump{k}_{p}"
            for k in range(num_lumps)
            for p in ("R0", "tilt_theta", "tilt_phi")
        }
        remove |= {"trajectory_A_breath", "trajectory_omega_breath"}
        dims = [d for d in dims if d.param_key not in remove]

    return dims


def build_search_space(
    nonspherical: bool = False,
    grtresna: bool = False,
    grtresna_lumps: int = GRTRESNA_DEFAULT_NUM_LUMPS,
    grtresna_ansatz: str = "free",
    grtresna_shell_profile: str = "compact",
    grtresna_matter_sector: str = "scalar",
    grtresna_shell_static: bool = True,
    grtresna_boson_allow_exotic: bool = False,
    grtresna_boson_matched_bounds: bool = False,
) -> list[SearchDimension]:
    """Return the optimizer search space.

    When ``nonspherical`` is True the spherically-symmetric radial space is
    extended with axisymmetric Legendre angular modes on the lapse and shift
    (see ANGULAR_SEARCH_SPACE).  The caller is responsible for also merging
    ANGULAR_BASE_OVERRIDES into the simulation base overrides so the modes are
    actually activated in params.txt.

    When ``grtresna`` is True the GRTresna momentum-carrying-matter dimensions
    REPLACE the radial recipe space, because the initial data is then produced
    entirely by the GRTresna solve.

    ``grtresna_matter_sector`` selects scalar lumps (shell/ring/free ansatz)
    vs boson-star (7-D complex scalar); independent of geometry ansatz.

    ``grtresna_shell_static`` controls whether velocity dims are included in
    the boson shell space (False = include momentum-carrying dims).

    ``grtresna_boson_allow_exotic`` keeps ``shell_exotic_fraction/phase`` in the
    boson-shell space (FTL campaigns); splash pins canonical-only.

    ``grtresna_boson_matched_bounds`` (v2 fair comparison) raises boson mass/amp/
    lambda bounds to match the scalar shell and inherits full geometric dims.
    """
    if grtresna:
        if grtresna_matter_sector == "boson_star":
            if grtresna_ansatz == "trajectory":
                traj_profile = (
                    grtresna_shell_profile
                    if grtresna_shell_profile in TRAJECTORY_BOSON_PROFILE_CHOICES
                    else "discovery"
                )
                return grtresna_trajectory_boson_search_space(
                    num_lumps=grtresna_lumps, profile=traj_profile
                )
            if grtresna_ansatz == "shell":
                return grtresna_boson_shell_search_space(
                    profile=grtresna_shell_profile,
                    static=grtresna_shell_static,
                    allow_exotic=grtresna_boson_allow_exotic,
                    matched_bounds=grtresna_boson_matched_bounds,
                )
            if grtresna_ansatz == "splash":
                return grtresna_boson_splash_search_space()
            return grtresna_boson_star_search_space()
        if grtresna_ansatz == "ring":
            return grtresna_ring_search_space()
        if grtresna_ansatz == "shell":
            return grtresna_shell_search_space(profile=grtresna_shell_profile)
        if grtresna_ansatz == "sh":
            sh_profile = (
                grtresna_shell_profile
                if grtresna_shell_profile in SH_PROFILE_CHOICES
                else "discovery"
            )
            return grtresna_sh_search_space(profile=sh_profile)
        if grtresna_ansatz == "trajectory":
            traj_profile = (
                grtresna_shell_profile
                if grtresna_shell_profile in TRAJECTORY_PROFILE_CHOICES
                else "discovery"
            )
            return grtresna_trajectory_search_space(
                num_lumps=grtresna_lumps, profile=traj_profile
            )
        if grtresna_ansatz == "boson_star":
            return grtresna_boson_star_search_space()
        if grtresna_ansatz == "splash":
            return grtresna_boson_splash_search_space()
        if grtresna_ansatz != "free":
            raise ValueError(f"unknown GRTresna ansatz: {grtresna_ansatz}")
        return grtresna_search_space(grtresna_lumps)
    space = list(DEFAULT_SEARCH_SPACE)
    if nonspherical:
        space += ANGULAR_SEARCH_SPACE
    return space

# Default GRTresna search space (used as a fallback when none is supplied).
GRTRESNA_SEARCH_SPACE: list[SearchDimension] = grtresna_search_space()
