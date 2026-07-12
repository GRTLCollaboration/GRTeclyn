#!/usr/bin/env python3
"""Amplitude-reduction (kappa) initial-data family for the rotating wormhole.

Phase B4/B5 of docs/RotatingWormholePlan.md.  Reuses the *same* runtime machinery
the QD campaign uses under the hood -- the validated GRTresna BosonStarBH solve
binary plus ``convert_chombo_to_gridinit`` -- rather than the search-space CLI
(whose ``scalar``/``boson_star`` sectors do not express the exotic-winding throat).

For each kappa the throat/phantom field amplitude is scaled ``f -> kappa*f`` and
GRTresna re-solves the Hamiltonian + momentum constraints for that configuration,
so every member is an *exact* Einstein solution that is simply out of equilibrium
(kappa=1 is the equilibrium arm; kappa<1 collapses).  The solved Chombo HDF5 is
then flattened to a ``.gridinit`` at the evolution's level-0 dx (lesson L2: never
evolve finer than the ID file's native dx).

Usage:
    solve_kappa_family.sh [KAPPAS] [NRANKS]
      KAPPAS  comma-separated list (default: 1.0,0.9,0.7,0.5)
      NRANKS  MPI ranks for the solve (default: 2)
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

# --- Paths -------------------------------------------------------------------
GRTECLYN_ROOT = Path(
    os.environ.get(
        "GRTECLYN_ROOT",
        "/home/jovyan/nachevsky/test/simulation/GRTeclyn",
    )
).resolve()
GRTRESNA_ROOT = Path(
    os.environ.get(
        "GRTRESNA_ROOT",
        str(GRTECLYN_ROOT.parent / "GRTresna"),
    )
).resolve()

BOSONSTAR_DIR = GRTRESNA_ROOT / "Examples" / "BosonStarBH"
BASE_PARAMS = BOSONSTAR_DIR / "params_rotating_wormhole_test.txt"
SOLVER_EXE = BOSONSTAR_DIR / "Main_BosonStarBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"

# Output ID tree (matches the convention referenced by the evolution params,
# e.g. params_rotating_grtresna_exotic.txt -> runs/rotating_wormhole_id/...).
ID_ROOT = GRTECLYN_ROOT / "runs" / "rotating_wormhole_id"

# --- Evolution target grid (RotatingWormholeCollapse, half-z RadialRecipe) ----
# Matches params_rotating_grtresna_exotic.txt: L=64, N=64 (dx=1.0), center at
# (32,32,0).  The GRTresna solve is also L=64/N=64, so the gridinit export is at
# the native solve dx == evolution level-0 dx (no interpolation; L2 satisfied).
# RES_N sets BOTH the GRTresna solve grid and the gridinit-export grid so the
# solved data is never up-sampled and the exported .gridinit is at the true
# solved dx.  N=64 -> dx=1.0 (coarse); N=128 -> dx=0.5 (high-res, L2-safe for a
# dx=0.5 level-0 evolution).
DEFAULT_L = 64.0
EVO_L = float(os.environ.get("EVO_L", str(DEFAULT_L)))
RES_N = int(os.environ.get("RES_N", os.environ.get("EVO_N", "64")))
EVO_N = RES_N
EVO_CENTER = (EVO_L / 2.0, EVO_L / 2.0, 0.0)
# Field mass mu for the confining potential V = 1/2 mu^2 |Phi|^2.  MASS=0 is the
# massless ghost (dispersive, the default studied so far); MASS>0 binds the
# phantom cloud into a soliton so the throat keeps its support.  This value must
# equal the evolution's `phantom_mass` (wormhole_case.py --mass) for the solved
# ID to stay in equilibrium at t=0.
MASS = float(os.environ.get("MASS", "0.0"))
# Q-ball self-interaction couplings.  V = 1/2 m^2 |Phi|^2 - 1/4 lam |Phi|^4 +
# 1/6 mu6 |Phi|^6.  With LAMBDA>0 and MU6>0 the lump is painted from the SOLVED
# flat-space Q-ball radial eigenstate phi0(r) (lump0_profile = 3) instead of the
# analytic Gaussian, so the field is a genuine bound state whose confinement
# survives the phantom sign flip (Rung 1a).  These MUST equal the evolution's
# phantom_lambda / phantom_mu6 (wormhole_case.py --lambda/--mu6).  Physics
# constraint for a bound profile: omega < mass (MASS).
LAMBDA = float(os.environ.get("LAMBDA", "0.0"))
MU6 = float(os.environ.get("MU6", "0.0"))
# Phase-winding rotation rate of the throat lump (must match lump0_omega in the
# base params and the evolution --omega).  Used as the Q-ball ODE frequency.
LUMP_OMEGA = float(os.environ.get("LUMP_OMEGA", "0.05"))


def _read_base_amp(params_text: str) -> float:
    m = re.search(r"^\s*lump0_amp\s*=\s*([0-9eE.+-]+)", params_text, re.MULTILINE)
    if not m:
        raise RuntimeError("lump0_amp not found in base params")
    return float(m.group(1))


def _write_qball_profile(dst_dir: Path) -> Path:
    """Solve the flat-space Q-ball radial eigenstate phi0(r) and tabulate it.

    Reuses the wrapper's validated ODE solver (same one the QD/boson-star
    campaigns use).  Returns the path GRTresna's ``qball_profile_path`` should
    point at; the winding painter loads it for profile==3 lumps and normalises
    to phi0(0)=1, so ``lump0_amp`` sets the throat's central amplitude.
    """
    from grteclyn_wrapper.grtresna.profiles.qball_ode import (
        cached_qball_radial_profile,
    )

    profile = cached_qball_radial_profile(MASS, LAMBDA, MU6, LUMP_OMEGA)
    dat_path = dst_dir / "qball_profile.dat"
    with dat_path.open("w", encoding="utf-8") as fh:
        fh.write(f"# flat-space Q-ball phi0(r): m={MASS} lam={LAMBDA} "
                 f"mu6={MU6} omega={LUMP_OMEGA}\n")
        fh.write(f"# phi_c(0) = {profile.phi_c}\n")
        for r_i, phi_i in zip(profile.r, profile.phi0):
            fh.write(f"{float(r_i):.10g} {float(phi_i):.10g}\n")
    return dat_path.resolve()


def _write_scaled_params(dst: Path, base_text: str, amp: float) -> None:
    """Copy base params with lump0_amp -> amp and Outputs redirected locally."""
    text = re.sub(
        r"^(\s*lump0_amp\s*=\s*)[0-9eE.+-]+",
        rf"\g<1>{amp:.10g}",
        base_text,
        flags=re.MULTILINE,
    )
    # Normalise the output dir so the driver always finds Outputs/... regardless
    # of the base file's output_path (base uses Outputs_rotwh/).
    text = re.sub(
        r"^(\s*output_path\s*=\s*).*$",
        r"\g<1>Outputs/",
        text,
        flags=re.MULTILINE,
    )
    # Solve at the requested resolution so the exported gridinit carries real
    # (not up-sampled) data at its dx (RES_N/L). Rewrite the "N = nx ny nz" line.
    text = re.sub(
        r"^(\s*N\s*=\s*).*$",
        rf"\g<1>{RES_N} {RES_N} {RES_N}",
        text,
        flags=re.MULTILINE,
    )
    # Rewrite L so the GRTresna solve domain matches the evolution box size.
    # Without this, EVO_L != base-params L would solve at the wrong domain size.
    text = re.sub(
        r"^(\s*L\s*=\s*).*$",
        rf"\g<1>{EVO_L}",
        text,
        flags=re.MULTILINE,
    )
    # Rewrite scalar_mass so the confining potential V = 1/2 mu^2 |Phi|^2 is
    # solved into the ID.  Must match the evolution's phantom_mass (same
    # convention V = 0.5 m^2 phi^2 in both codes) or the lump disperses at t=0.
    text = re.sub(
        r"^(\s*scalar_mass\s*=\s*).*$",
        rf"\g<1>{MASS:.10g}",
        text,
        flags=re.MULTILINE,
    )
    # Q-ball (Rung 1a): a self-interacting bound eigenstate instead of the
    # analytic Gaussian.  Set the self-interaction couplings, switch the throat
    # lump to the tabulated ODE profile (profile 3), and point the shared
    # qball_profile_path at the solved phi0(r) table.  The C++ potential_value
    # matches (V = 1/2 m^2|Phi|^2 - 1/4 lam|Phi|^4 + 1/6 mu|Phi|^6), so the
    # solved ID stays in equilibrium under the matching evolution potential.
    if LAMBDA > 0.0:
        text = re.sub(r"^(\s*scalar_lambda\s*=\s*).*$",
                      rf"\g<1>{LAMBDA:.10g}", text, flags=re.MULTILINE)
        text = re.sub(r"^(\s*scalar_mu\s*=\s*).*$",
                      rf"\g<1>{MU6:.10g}", text, flags=re.MULTILINE)
        text = re.sub(r"^(\s*lump0_profile\s*=\s*).*$",
                      r"\g<1>3", text, flags=re.MULTILINE)
        qball_path = _write_qball_profile(dst.parent)
        # Append (or replace) the shared tabulated-profile path.
        if re.search(r"^\s*qball_profile_path\s*=", text, re.MULTILINE):
            text = re.sub(r"^(\s*qball_profile_path\s*=\s*).*$",
                          rf"\g<1>{qball_path}", text, flags=re.MULTILINE)
        else:
            text += f"\nqball_profile_path = {qball_path}\n"
    dst.write_text(text)


def _parse_convergence(err_file: Path) -> tuple[float, float] | None:
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
    if last is None:
        return None
    return last[1], last[2]


def _L_suffix() -> str:
    """Append _L<val> only when L differs from 64 (backward compat with existing runs)."""
    if abs(EVO_L - DEFAULT_L) < 0.1:
        return ""
    return f"_L{int(EVO_L)}"


def _mass_suffix() -> str:
    """Append _mass<val> only when MASS>0 (backward compat with massless runs)."""
    if MASS <= 0.0:
        return ""
    return f"_mass{MASS:g}".replace(".", "p")


def _qball_suffix() -> str:
    """Append _qball_lam<val>_mu6<val> for a self-interacting (Q-ball) profile.
    Must match wormhole_case.py's _qball_suffix so the evolution finds the ID."""
    if LAMBDA <= 0.0:
        return ""
    return f"_qball_lam{LAMBDA:g}_mu6{MU6:g}".replace(".", "p")


def solve_one(kappa: float, base_amp: float, base_text: str, nranks: int) -> dict:
    amp = kappa * base_amp
    # Tag by resolution (and box size / mass when non-default) so families coexist.
    dx = EVO_L / RES_N
    dx_tag = f"dx{dx:.3g}".replace(".", "p")
    tag = (f"rotwh_omega_p0p05_m1_kappa_{kappa:.2f}_{dx_tag}"
           f"{_L_suffix()}{_mass_suffix()}{_qball_suffix()}").replace(".", "p")
    run_dir = ID_ROOT / tag
    outputs = run_dir / "Outputs"
    outputs.mkdir(parents=True, exist_ok=True)
    (run_dir / "pout").mkdir(exist_ok=True)

    params_path = run_dir / "params.txt"
    _write_scaled_params(params_path, base_text, amp)

    print(f"\n=== kappa={kappa:.2f}  amp={amp:.5g}  -> {run_dir} ===", flush=True)

    # Run the (already-built, c_Ham_abs-fixed) solver.  cwd=run_dir so Outputs/
    # and pout/ land inside the run directory.
    mpirun = os.environ.get("MPIRUN", "mpirun")
    cmd = [mpirun, "--oversubscribe", "-np", str(nranks),
           str(SOLVER_EXE), str(params_path)]
    subprocess.run(cmd, cwd=str(run_dir), check=True)

    conv = _parse_convergence(run_dir / "Ham_and_Mom_errors.txt")
    hdf5 = outputs / "InitialDataFinal.3d.hdf5"
    if not hdf5.exists():
        raise FileNotFoundError(f"solve produced no {hdf5}")

    # Flatten Chombo AMR -> uniform .gridinit at the evolution's level-0 dx.
    from grteclyn_wrapper.grtresna.io.conversion import convert_chombo_to_gridinit

    gridinit = run_dir / "initial_data.gridinit"
    convert_chombo_to_gridinit(
        hdf5,
        gridinit,
        N=EVO_N,
        L=EVO_L,
        target_center=EVO_CENTER,
        lumps=None,                       # use the solved fields directly
        matter_model="grtresna_complex_scalar",
        delete_source=True,               # drop the ~0.8 GB solve HDF5 once flattened
    )

    # Reclaim the rest of the heavy solve scratch (NL diagnostics, pout); keep
    # only the gridinit + params + Ham_and_Mom_errors.txt.  Disk discipline (L6).
    if not os.environ.get("KEEP_SOLVE_SCRATCH"):
        import shutil
        for junk in (outputs, run_dir / "pout"):
            shutil.rmtree(junk, ignore_errors=True)

    return {"kappa": kappa, "amp": amp, "run_dir": run_dir,
            "gridinit": gridinit, "conv": conv}


def main() -> int:
    kappas_arg = sys.argv[1] if len(sys.argv) > 1 else "1.0,0.9,0.7,0.5"
    nranks = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    kappas = [float(x) for x in kappas_arg.split(",") if x.strip()]

    if not SOLVER_EXE.exists():
        print(f"error: solver binary not found: {SOLVER_EXE}", file=sys.stderr)
        print("build it with scripts/wormhole/build/build_grtresna_bosonstar.sh",
              file=sys.stderr)
        return 2

    base_text = BASE_PARAMS.read_text()
    base_amp = _read_base_amp(base_text)
    print(f"base amp = {base_amp}  kappas = {kappas}  nranks = {nranks}")
    print(f"evolution target grid: N={EVO_N} L={EVO_L} center={EVO_CENTER}")
    print(f"field mass (scalar_mass) = {MASS}"
          f"{'  (massless ghost)' if MASS <= 0 else '  (confining potential)'}")

    results = [solve_one(k, base_amp, base_text, nranks) for k in kappas]

    print("\n================ kappa family summary ================")
    print(f"{'kappa':>6} {'amp':>8} {'Ham%':>9} {'Mom%':>9}  gridinit")
    for r in results:
        ham, mom = (r["conv"] if r["conv"] else (float("nan"), float("nan")))
        ok = "OK" if r["gridinit"].exists() else "MISSING"
        print(f"{r['kappa']:>6.2f} {r['amp']:>8.4g} {ham:>9.4g} {mom:>9.4g}  "
              f"{ok} {r['gridinit']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
