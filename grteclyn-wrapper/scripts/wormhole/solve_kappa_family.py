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
EVO_L = float(os.environ.get("EVO_L", "64.0"))
RES_N = int(os.environ.get("RES_N", os.environ.get("EVO_N", "64")))
EVO_N = RES_N
EVO_CENTER = (32.0, 32.0, 0.0)


def _read_base_amp(params_text: str) -> float:
    m = re.search(r"^\s*lump0_amp\s*=\s*([0-9eE.+-]+)", params_text, re.MULTILINE)
    if not m:
        raise RuntimeError("lump0_amp not found in base params")
    return float(m.group(1))


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


def solve_one(kappa: float, base_amp: float, base_text: str, nranks: int) -> dict:
    amp = kappa * base_amp
    # Tag by resolution so dx=1.0 (N=64) and dx=0.5 (N=128) families coexist.
    dx = EVO_L / RES_N
    dx_tag = f"dx{dx:.3g}".replace(".", "p")
    tag = f"rotwh_omega_p0p05_m1_kappa_{kappa:.2f}_{dx_tag}".replace(".", "p")
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
        print("build it with scripts/wormhole/build_grtresna_bosonstar.sh",
              file=sys.stderr)
        return 2

    base_text = BASE_PARAMS.read_text()
    base_amp = _read_base_amp(base_text)
    print(f"base amp = {base_amp}  kappas = {kappas}  nranks = {nranks}")
    print(f"evolution target grid: N={EVO_N} L={EVO_L} center={EVO_CENTER}")

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
