#!/usr/bin/env python3
"""Isolated rotating Q-torus initial data (stationary m=1 eigenstate).

This is the "torus alone first" validation of the stationary-eigenstate fix
(research/rotatingwormhole/OrbitalPumpPlan.md, "make mass not disperse").
Instead of twisting a spherical Q-ball into an m=1 torus (a non-stationary
ansatz that drains its Noether charge with a t~13-16 half-life), we solve the
GENUINE 2D spinning Q-ball eigenstate

    Phi = f(rho, z) e^{i(m phi_az - omega t)}

via the Newton/amplitude-continuation solver in
``grteclyn_wrapper.grtresna.profiles.qball_torus``, tabulate f(rho, z), and paint
it into GRTresna as a single ``profile == 4`` winding lump on a THROAT-FREE,
CENTERED, full-box background (bh1_bare_mass = 0, all-Sommerfeld boundaries so
GRTresna centres the grid at L/2).  GRTresna re-solves the Hamiltonian + momentum
constraints on top, and the solved Chombo HDF5 is flattened to a ``.gridinit`` at
the evolution's level-0 dx (lesson L2).

The resulting gridinit is a flat-space + torus configuration; evolving it
passively (no throat, no pump) and measuring Q_sphere / confined_frac retention
is the direct test of whether the eigenstate holds together (it should, unlike
the twisted-sphere baseline).

Env knobs (all optional; defaults are the validated compact couplings):
    MASS, LAMBDA, MU6        scalar potential V = 1/2 m^2|Phi|^2 - 1/4 lam|Phi|^4
                             + 1/6 mu6|Phi|^6   (default 0.5 / 170 / 14450)
    TORUS_OMEGA              target internal U(1) frequency (default 0.25;
                             compact/thick-wall.  Deep thin-wall omega<~0.15 makes
                             the torus too large for the box.)
    M_AZ                     azimuthal winding number m (default 1)
    KAPPA                    amplitude scale f -> kappa*f_max (default 1.0 =
                             true eigenstate; <1 = sub-critical for later use)
    EXOTIC                   0 = normal positive-energy torus (default, cleanest
                             eigenstate test); 1 = phantom (wormhole support)
    EVO_L, RES_N             solve/evolution box size and resolution
                             (default L=64, N=128 -> dx=0.5)
    TORUS_NRHO, TORUS_NZ     solver grid (default 160 x 160)

Usage:
    solve_torus.sh [NRANKS]
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

GRTECLYN_ROOT = Path(
    os.environ.get("GRTECLYN_ROOT", "/home/jovyan/nachevsky/test/simulation/GRTeclyn")
).resolve()
GRTRESNA_ROOT = Path(
    os.environ.get("GRTRESNA_ROOT", str(GRTECLYN_ROOT.parent / "GRTresna"))
).resolve()

BOSONSTAR_DIR = GRTRESNA_ROOT / "Examples" / "BosonStarBH"
SOLVER_EXE = BOSONSTAR_DIR / "Main_BosonStarBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"
ID_ROOT = GRTECLYN_ROOT / "runs" / "rotating_torus_id"

MASS = float(os.environ.get("MASS", "0.5"))
LAMBDA = float(os.environ.get("LAMBDA", "170.0"))
MU6 = float(os.environ.get("MU6", "14450.0"))
TORUS_OMEGA = float(os.environ.get("TORUS_OMEGA", "0.25"))
M_AZ = int(os.environ.get("M_AZ", "1"))
KAPPA = float(os.environ.get("KAPPA", "1.0"))
EXOTIC = int(os.environ.get("EXOTIC", "0"))
EVO_L = float(os.environ.get("EVO_L", "64.0"))
RES_N = int(os.environ.get("RES_N", os.environ.get("EVO_N", "128")))
TORUS_NRHO = int(os.environ.get("TORUS_NRHO", "160"))
TORUS_NZ = int(os.environ.get("TORUS_NZ", "160"))
# THROAT_MASS > 0 turns the isolated torus into a torus-SUPPORTED rotating
# wormhole: bh1_bare_mass = THROAT_MASS creates the Bowen-York throat psi ~
# (m/2)/r (the positive bare mass that implodes when the exotic support is cut).
# The exotic Q-torus (EXOTIC=1) hugs the throat and holds it open; ramping the
# evolution's wormhole_support_strength -> 0 triggers collapse (Phase 6/8).
THROAT_MASS = float(os.environ.get("THROAT_MASS", "0.0"))


def _tag() -> str:
    om = f"{TORUS_OMEGA:.3f}".replace(".", "p")
    ka = f"{KAPPA:.2f}".replace(".", "p")
    dx = f"{EVO_L / RES_N:.3g}".replace(".", "p")
    ex = "_exotic" if EXOTIC else ""
    th = f"_throat{THROAT_MASS:g}" if THROAT_MASS > 0 else ""
    return (f"torus_m{M_AZ}_om{om}_kappa{ka}_dx{dx}_L{int(EVO_L)}"
            f"_lam{LAMBDA:g}_mu6{MU6:g}{ex}{th}").replace(".", "p")


def _solve_and_write_torus(run_dir: Path) -> tuple[Path, float, float]:
    """Solve the 2D Q-torus eigenstate, tabulate it, return (path, f_max, omega)."""
    from grteclyn_wrapper.grtresna.profiles.qball_couplings import QBallCouplings
    from grteclyn_wrapper.grtresna.profiles.qball_torus import (
        solve_qball_torus_profile,
        write_torus_profile,
    )

    c = QBallCouplings(mass=MASS, lam=LAMBDA, mu=MU6, omega=TORUS_OMEGA)
    print(f"[torus] solving m={M_AZ} eigenstate: mass={MASS} lam={LAMBDA} "
          f"mu6={MU6} target_omega={TORUS_OMEGA} grid={TORUS_NRHO}x{TORUS_NZ}",
          flush=True)
    prof = solve_qball_torus_profile(
        c, m_az=M_AZ, n_rho=TORUS_NRHO, n_z=TORUS_NZ, verbose=True
    )
    dat = run_dir / "qball_torus.dat"
    write_torus_profile(prof, dat)
    print(f"[torus] solved: omega={prof.omega:.5f} f_max={prof.f_max:.6f} "
          f"Q={prof.noether_charge:.4f} -> {dat}", flush=True)
    return dat.resolve(), prof.f_max, prof.omega


def _write_params(dst: Path, torus_path: Path, f_max: float, omega: float) -> None:
    """Write a throat-free, centered, full-box GRTresna params for a single
    profile==4 winding torus lump.  All-Sommerfeld boundaries => GRTresna centres
    the grid at L/2, so lump_center 0 0 0 sits at the physical box centre."""
    amp = KAPPA * f_max
    throat_desc = (f"torus-supported wormhole (throat bh1_bare_mass={THROAT_MASS:g})"
                   if THROAT_MASS > 0 else "throat-free flat background")
    text = f"""# Rotating Q-torus (stationary m={M_AZ} eigenstate) initial data.
# Genuine 2D spinning Q-ball f(rho, z) e^(i(m phi - omega t)) painted from the
# tabulated eigenstate solve (profile == 4), {throat_desc}.

output_path = Outputs/
output_filename = InitialDataFinal.3d.hdf5
error_filename = Ham_and_Mom_errors
write_diagnostics = 0

N = {RES_N} {RES_N} {RES_N}
L = {EVO_L:g}
max_level = 0
refine_threshold = 0.5
block_factor = 16
max_grid_size = 32
regrid_radius = 10
coefficient_average_type = arithmetic

is_periodic = 0 0 0
use_compact_Vi_ansatz = 1
hi_boundary = 0 0 0
lo_boundary = 0 0 0

G_Newton = 1.0
scalar_mass = {MASS:g}
scalar_lambda = {LAMBDA:g}
scalar_mu = {MU6:g}
bs_phi_c = {amp:.10g}
bs_profile_width = 8.0
bs_omega = {omega:.10g}

# Single winding torus lump (profile 4 = tabulated 2D Q-torus eigenstate).
num_lumps = 1
lump0_amp = {amp:.10g}
lump0_width = 8.0
lump0_center = 0.0 0.0 0.0
lump0_velocity = 0.0 0.0 0.0
lump0_omega = {omega:.10g}
lump0_mode = {M_AZ}
lump0_exotic = {EXOTIC}
lump0_winding = 1
lump0_profile = 4
qball_profile_path = {torus_path}

# Throat: Bowen-York puncture psi ~ (bh1_bare_mass/2)/r (=0 when THROAT_MASS=0,
# i.e. flat regular background for the isolated torus).
regularised_part_psi = 1.0
sign_of_K = 1
maximal_slicing = 1
psi_floor = 0.1
maximal_jacobian_cap = 1.0e6
psi_relaxation = 0.6

bh1_bare_mass = {THROAT_MASS:g}
bh1_spin = 0.0 0.0 0.0
bh1_momentum = 0.0 0.0 0.0
bh1_offset = 0.0 0.0 0.0
bh2_bare_mass = 0.0
bh2_spin = 0.0 0.0 0.0
bh2_momentum = 0.0 0.0 0.0
bh2_offset = 0.0 0.0 0.0

max_NL_iterations = 30
NL_exit_tolerance = 1.0
NL_stall_tolerance = 0.02
deactivate_zero_mode = 0
"""
    dst.write_text(text)


def _parse_convergence(err_file: Path):
    if not err_file.exists():
        return None
    last = None
    for line in err_file.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3:
            try:
                last = (int(parts[0]), float(parts[1]), float(parts[2]))
            except ValueError:
                continue
    return (last[1], last[2]) if last else None


def main() -> int:
    nranks = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    if not SOLVER_EXE.exists():
        print(f"error: solver binary not found: {SOLVER_EXE}\n"
              "build it with scripts/wormhole/build/build_grtresna_bosonstar.sh",
              file=sys.stderr)
        return 2

    run_dir = ID_ROOT / _tag()
    (run_dir / "Outputs").mkdir(parents=True, exist_ok=True)
    (run_dir / "pout").mkdir(exist_ok=True)

    torus_path, f_max, omega = _solve_and_write_torus(run_dir)

    params_path = run_dir / "params.txt"
    _write_params(params_path, torus_path, f_max, omega)

    print(f"\n=== isolated torus ID -> {run_dir} ===", flush=True)
    mpirun = os.environ.get("MPIRUN", "mpirun")
    cmd = [mpirun, "--oversubscribe", "-np", str(nranks),
           str(SOLVER_EXE), str(params_path)]
    subprocess.run(cmd, cwd=str(run_dir), check=True)

    conv = _parse_convergence(run_dir / "Ham_and_Mom_errors.txt")
    hdf5 = run_dir / "Outputs" / "InitialDataFinal.3d.hdf5"
    if not hdf5.exists():
        raise FileNotFoundError(f"solve produced no {hdf5}")

    from grteclyn_wrapper.grtresna.io.conversion import convert_chombo_to_gridinit

    gridinit = run_dir / "initial_data.gridinit"
    convert_chombo_to_gridinit(
        hdf5,
        gridinit,
        N=RES_N,
        L=EVO_L,
        target_center=(EVO_L / 2.0, EVO_L / 2.0, EVO_L / 2.0),
        lumps=None,
        matter_model="grtresna_complex_scalar",
        delete_source=True,
    )

    if not os.environ.get("KEEP_SOLVE_SCRATCH"):
        import shutil
        for junk in (run_dir / "Outputs", run_dir / "pout"):
            shutil.rmtree(junk, ignore_errors=True)

    ham, mom = conv if conv else (float("nan"), float("nan"))
    print("\n================ isolated torus ID summary ================")
    print(f"omega={omega:.5f}  f_max={f_max:.6f}  amp={KAPPA*f_max:.6f}")
    print(f"Ham%={ham:.4g}  Mom%={mom:.4g}")
    print(f"gridinit: {'OK' if gridinit.exists() else 'MISSING'} {gridinit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
