"""Write GRTresna params.txt from solver configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..profiles.boson_star import (
    PROFILE_ODE_BOUND,
    PROFILE_SELFGRAV_BOUND,
    grtresna_lump_profile,
)
from .config import GRTresnaConfig


def _fmt(v: Any) -> str:
    if isinstance(v, (tuple, list)):
        return " ".join(str(x) for x in v)
    return str(v)


def _lump_lines(cfg: GRTresnaConfig) -> list[str]:
    if cfg.lumps:
        lines = [f"num_lumps = {len(cfg.lumps)}"]
        for k, lump in enumerate(cfg.lumps):
            lines.extend([
                f"lump{k}_amp = {lump.get('amp', 0.0)}",
                f"lump{k}_width = {lump.get('width', 5.0)}",
                f"lump{k}_center = {_fmt(tuple(lump.get('center', (0.0, 0.0, 0.0))))}",
                f"lump{k}_velocity = {_fmt(tuple(lump.get('velocity', (0.0, 0.0, 0.0))))}",
                f"lump{k}_omega = {lump.get('omega', 0.0)}",
                f"lump{k}_mode = {int(lump.get('mode', 0))}",
                f"lump{k}_exotic = {int(lump.get('exotic', 0))}",
                f"lump{k}_winding = {int(lump.get('winding', 0))}",
                f"lump{k}_profile = {grtresna_lump_profile(int(lump.get('profile', 0)))}",
            ])
            # Per-lump table (mixed-sector selfgrav pairs); C++ falls back to
            # the global qball_profile_path when this key is absent.
            if lump.get("profile_path"):
                lines.append(f"lump{k}_profile_path = {lump['profile_path']}")
        return lines
    return [
        f"lump_amp = {cfg.lump_amp}",
        f"lump_width = {cfg.lump_width}",
        f"lump_center = {_fmt(cfg.lump_center)}",
        f"lump_velocity = {_fmt(cfg.lump_velocity)}",
        f"lump_omega = {cfg.lump_omega}",
        f"lump_mode = {cfg.lump_mode}",
        f"lump_exotic = {int(cfg.lump_exotic)}",
    ]


def _maybe_write_qball_profile(cfg: GRTresnaConfig, params_path: Path) -> str | None:
    """Emit the tabulated radial profile phi0(r) for ODE-seeded lumps.

    Handles two seed families that both load through the coupling-agnostic C++
    profile==3 tabulated loader:

      * ``PROFILE_ODE_BOUND`` (3): flat-space Q-ball, bound by self-interaction.
      * ``PROFILE_SELFGRAV_BOUND`` (4): self-gravitating boson star, bound by
        gravity.  Its frequency is the gravitational eigenvalue, which we write
        back onto ``cfg.bs_omega`` (and each self-grav lump's ``omega``) so the
        painter's ``Pi2 = -omega*phi0/alpha`` uses the correct value.

    Self-gravitating lumps take priority when both are present (they should not
    be mixed in one campaign).
    """
    from ..profiles.qball_ode import cached_qball_radial_profile

    selfgrav_lumps = [
        lump for lump in cfg.lumps
        if int(lump.get("profile", 0)) == PROFILE_SELFGRAV_BOUND
    ]
    if selfgrav_lumps:
        return _write_selfgrav_profile(cfg, selfgrav_lumps, params_path)

    ode_lumps = [
        lump for lump in cfg.lumps
        if int(lump.get("profile", 0)) == PROFILE_ODE_BOUND
    ]
    if not ode_lumps:
        return None

    src = ode_lumps[0]
    mass = float(src.get("qball_mass", cfg.scalar_mass))
    lam = float(src.get("qball_lam", cfg.scalar_lambda))
    mu = float(src.get("qball_mu", cfg.scalar_mu))
    omega = float(src.get("qball_omega", cfg.bs_omega))

    profile = cached_qball_radial_profile(mass, lam, mu, omega)
    dat_path = params_path.parent / "qball_profile.dat"
    with dat_path.open("w", encoding="utf-8") as fh:
        fh.write(f"# flat-space Q-ball phi0(r): m={mass} lam={lam} mu={mu} omega={omega}\n")
        fh.write(f"# phi_c(0) = {profile.phi_c}\n")
        for r_i, phi_i in zip(profile.r, profile.phi0):
            fh.write(f"{float(r_i):.10g} {float(phi_i):.10g}\n")
    return str(dat_path.resolve())


def _emit_star_table(dat_path: Path, profile, *, mass: float, lam: float, mu: float, sector: str) -> str:
    with dat_path.open("w", encoding="utf-8") as fh:
        fh.write(
            f"# self-gravitating boson star phi0(r) [{sector}]: m={mass} lam={lam} "
            f"mu={mu} omega_eigenvalue={profile.omega} ADM_mass={profile.adm_mass}\n"
        )
        # Three columns: r  phi0(r)  alpha(r).  The lapse alpha(r) lets GRTresna
        # set the stationary momentum Pi_im = -(omega/alpha) phi0 (a boson star is
        # stationary, NOT static: painting Pi with the flat-space alpha=1 gives
        # the wrong kinetic energy and the seed radiates).
        fh.write(f"# phi_c(0) = {profile.phi_c}\n")
        fh.write("# columns: r  phi0  alpha\n")
        for r_i, phi_i, alpha_i in zip(profile.r, profile.phi0, profile.alpha):
            fh.write(f"{float(r_i):.10g} {float(phi_i):.10g} {float(alpha_i):.10g}\n")
    return str(dat_path.resolve())


def _write_selfgrav_profile(
    cfg: GRTresnaConfig, selfgrav_lumps: list[dict], params_path: Path
) -> str:
    """Solve the self-gravitating star ODE(s) and emit tabulated phi0(r).

    The gravitational eigenvalue ``omega`` is written back onto ``cfg.bs_omega``
    and every self-grav lump's ``omega`` so the painter and the constraint solve
    agree on the phase velocity.

    Mixed-sector pairs (canonical + exotic lumps) get ONE TABLE PER SECTOR: the
    phantom star's repulsive self-gravity gives it a slightly puffed profile,
    its own lapse (alpha > 1 inside) and a negative ADM mass, so sharing the
    canonical table would seed it off-equilibrium.  Each lump gets its own
    ``profile_path`` (read per-lump by GRTresna, with the canonical table as the
    global fallback).  Both sectors solve to the SAME requested frequency, so
    the C++ painter's global-omega momentum paint stays correct; per-lump omega
    support in C++ is only needed for mixed-frequency pairs.
    """
    from ..profiles.boson_star_ode import cached_selfgrav_profile

    src = selfgrav_lumps[0]
    mass = float(src.get("qball_mass", cfg.scalar_mass))
    lam = float(src.get("qball_lam", cfg.scalar_lambda))
    mu = float(src.get("qball_mu", cfg.scalar_mu))
    # A requested frequency (qball_omega, sextic couplings) pins the star on
    # the intended eigenvalue branch and makes phi_c a solved OUTPUT;
    # otherwise phi_c parameterizes the star (mini-star path, omega solved).
    omega_target = float(src.get("qball_omega", 0.0))
    if omega_target > 0.0 and lam > 0.0 and mu > 0.0:
        from ..profiles.boson_star_ode import cached_selfgrav_at_omega

        global_path: str | None = None
        omega_out = 0.0
        for sector, sign in (("canonical", 1.0), ("exotic", -1.0)):
            sector_lumps = [
                lump for lump in selfgrav_lumps
                if bool(int(lump.get("exotic", 0))) == (sign < 0)
            ]
            if not sector_lumps:
                continue
            profile = cached_selfgrav_at_omega(mass, lam, mu, omega_target, sign)
            name = "qball_profile.dat" if sign > 0 else "qball_profile_exotic.dat"
            path = _emit_star_table(
                params_path.parent / name, profile,
                mass=mass, lam=lam, mu=mu, sector=sector,
            )
            for lump in sector_lumps:
                # The C++/Python painters rescale the table by amp/phi_c; only
                # amp == the solved phi_c keeps the seed on the eigenstate.
                lump["amp"] = float(profile.phi_c)
                lump["omega"] = float(profile.omega)
                lump["profile_path"] = path
            if global_path is None or sign > 0:
                global_path = path
                omega_out = float(profile.omega)
        cfg.bs_omega = omega_out
        assert global_path is not None
        return global_path

    phi_c = float(src.get("amp", 0.0)) or float(cfg.bs_phi_c)
    profile = cached_selfgrav_profile(mass, lam, mu, phi_c)
    # Propagate the eigenvalue so the painter's Pi2 = -omega*phi0/alpha and the
    # constraint solve share the phase velocity.
    cfg.bs_omega = float(profile.omega)
    for lump in selfgrav_lumps:
        lump["omega"] = float(profile.omega)
    return _emit_star_table(
        params_path.parent / "qball_profile.dat", profile,
        mass=mass, lam=lam, mu=mu, sector="canonical",
    )


def write_grtresna_params(cfg: GRTresnaConfig, path: Path) -> None:
    """Write a GRTresna params.txt."""
    # Emit the tabulated ODE/self-grav profile FIRST: for a self-gravitating
    # star this solves the eigenvalue and writes it back onto cfg.bs_omega, so
    # it must run before the bs_omega / lump lines below are built.
    qball_path = _maybe_write_qball_profile(cfg, path)

    lines = [
        "output_path = Outputs/",
        "output_filename = InitialDataFinal.3d.hdf5",
        "error_filename = Ham_and_Mom_errors",
        f"write_diagnostics = {1 if cfg.write_diagnostics else 0}",
        "",
        f"N = {_fmt(cfg.N)}",
        f"L = {cfg.L}",
        f"max_level = {cfg.max_level}",
        f"refine_threshold = {cfg.refine_threshold}",
        f"block_factor = {cfg.block_factor}",
        f"max_grid_size = {cfg.max_grid_size}",
        f"regrid_radius = {cfg.regrid_radius}",
        f"coefficient_average_type = {cfg.coefficient_average_type}",
        "",
        "is_periodic = 0 0 0",
        f"use_compact_Vi_ansatz = {cfg.use_compact_Vi_ansatz}",
        f"hi_boundary = {_fmt(cfg.hi_boundary)}",
        f"lo_boundary = {_fmt(cfg.lo_boundary)}",
        "",
        f"G_Newton = {cfg.G_Newton}",
        f"phi_0 = {cfg.phi_0}",
        f"dphi = {cfg.dphi}",
        f"dphi_length = {cfg.dphi_length}",
        f"pi_0 = {cfg.pi_0}",
        f"dpi = {cfg.dpi}",
        f"dpi_length = {cfg.dpi_length}",
        f"scalar_mass = {cfg.scalar_mass}",
        f"scalar_lambda = {cfg.scalar_lambda}",
        f"scalar_mu = {getattr(cfg, 'scalar_mu', 0.0)}",
    ]
    if cfg.matter_model in (
        "grtresna_complex_scalar",
        "grtresna_bicomplex_scalar",
    ):
        lines.extend([
            f"scalar_sign = {cfg.scalar_sign}",
            f"bs_phi_c = {cfg.bs_phi_c}",
            f"bs_profile_width = {cfg.bs_profile_width}",
            f"bs_omega = {cfg.bs_omega}",
        ])
    if qball_path is not None:
        lines.append(f"qball_profile_path = {qball_path}")
    lines.extend(_lump_lines(cfg))
    lines.extend([
        f"regularised_part_psi = {cfg.regularised_part_psi}",
        f"sign_of_K = {cfg.sign_of_K}",
        f"maximal_slicing = {1 if cfg.maximal_slicing else 0}",
        f"psi_relaxation = {cfg.psi_relaxation}",
        f"psi_floor = {cfg.psi_floor}",
        f"maximal_jacobian_cap = {cfg.maximal_jacobian_cap}",
        "",
        f"bh1_bare_mass = {cfg.bh1_bare_mass}",
        f"bh1_spin = {_fmt(cfg.bh1_spin)}",
        f"bh1_momentum = {_fmt(cfg.bh1_momentum)}",
        f"bh1_offset = {_fmt(cfg.bh1_offset)}",
        f"bh2_bare_mass = {cfg.bh2_bare_mass}",
        f"bh2_spin = {_fmt(cfg.bh2_spin)}",
        f"bh2_momentum = {_fmt(cfg.bh2_momentum)}",
        f"bh2_offset = {_fmt(cfg.bh2_offset)}",
        "",
        f"max_NL_iterations = {cfg.max_NL_iterations}",
        f"NL_exit_tolerance = {cfg.nl_exit_tolerance}",
        f"NL_stall_tolerance = {cfg.nl_stall_tolerance}",
        "deactivate_zero_mode = 0",
    ])
    path.write_text("\n".join(lines) + "\n")
