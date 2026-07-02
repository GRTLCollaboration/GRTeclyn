"""Write GRTresna params.txt from solver configuration."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..profiles.boson_star import PROFILE_ODE_BOUND, grtresna_lump_profile
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
                f"lump{k}_profile = {grtresna_lump_profile(int(lump.get('profile', 0)))}",
            ])
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
    from ..profiles.qball_ode import cached_qball_radial_profile

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


def write_grtresna_params(cfg: GRTresnaConfig, path: Path) -> None:
    """Write a GRTresna params.txt."""
    lines = [
        f"output_path = Outputs/",
        f"output_filename = InitialDataFinal.3d.hdf5",
        f"error_filename = Ham_and_Mom_errors",
        f"write_diagnostics = {1 if cfg.write_diagnostics else 0}",
        f"",
        f"N = {_fmt(cfg.N)}",
        f"L = {cfg.L}",
        f"max_level = {cfg.max_level}",
        f"refine_threshold = {cfg.refine_threshold}",
        f"block_factor = {cfg.block_factor}",
        f"max_grid_size = {cfg.max_grid_size}",
        f"regrid_radius = {cfg.regrid_radius}",
        f"coefficient_average_type = {cfg.coefficient_average_type}",
        f"",
        f"is_periodic = 0 0 0",
        f"use_compact_Vi_ansatz = {cfg.use_compact_Vi_ansatz}",
        f"hi_boundary = {_fmt(cfg.hi_boundary)}",
        f"lo_boundary = {_fmt(cfg.lo_boundary)}",
        f"",
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
    qball_path = _maybe_write_qball_profile(cfg, path)
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
        f"",
        f"bh1_bare_mass = {cfg.bh1_bare_mass}",
        f"bh1_spin = {_fmt(cfg.bh1_spin)}",
        f"bh1_momentum = {_fmt(cfg.bh1_momentum)}",
        f"bh1_offset = {_fmt(cfg.bh1_offset)}",
        f"bh2_bare_mass = {cfg.bh2_bare_mass}",
        f"bh2_spin = {_fmt(cfg.bh2_spin)}",
        f"bh2_momentum = {_fmt(cfg.bh2_momentum)}",
        f"bh2_offset = {_fmt(cfg.bh2_offset)}",
        f"",
        f"max_NL_iterations = {cfg.max_NL_iterations}",
        f"NL_exit_tolerance = {cfg.nl_exit_tolerance}",
        f"NL_stall_tolerance = {cfg.nl_stall_tolerance}",
        f"deactivate_zero_mode = 0",
    ])
    path.write_text("\n".join(lines) + "\n")
