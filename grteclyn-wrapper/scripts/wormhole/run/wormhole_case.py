#!/usr/bin/env python3
"""Single-case rotating-wormhole evolution driver (CLI, no per-case param files).

Mirrors the QD-campaign pattern: instead of committing a params_*.txt per case,
this generates the evolution params from one in-code template + CLI knobs, writes
it into the run directory, launches the GRTeclyn CUDA binary, and (optionally)
renders frames offline then prunes the heavy plotfiles.

The GRTresna-solved initial data is located by the same tag scheme the ID driver
(`solve_kappa_family.py`) writes:

    runs/rotating_wormhole_id/rotwh_omega_p<omega>_m<m>_kappa_<kappa>_dx<dx>/initial_data.gridinit

L2 lesson: the evolution level-0 dx must equal the gridinit's native dx, so
--dx selects BOTH (dx=1.0 -> N=64 solve, dx=0.5 -> N=128 solve).

--box-size L controls the domain size (default 64).  A larger box pushes the
Sommerfeld boundary farther from the throat; used for the resolution / boundary
convergence study (is the t~5 dispersal physical or a boundary artifact?).

The stop-time auto-scales with the outermost extraction radius when not
explicitly set: stop_time = r_outer + 6 (light-crossing to the shell plus a
buffer for signal development).  Extraction uses 2 radii (inner fixed at 12,
outer at L/2 - 8) to keep the consumer lightweight.

--mass mu adds a confining potential V = 0.5 mu^2 |Phi|^2 to the phantom scalar
(default 0 = massless ghost, which disperses). A nonzero mass binds the phantom
cloud into a soliton so the throat keeps its support. The ID must be solved with
the SAME mass: `MASS=<mu> solve_kappa_family.sh` (the potential convention
V=0.5 m^2 phi^2 is shared by both codes).

Usage:
    wormhole_case.sh --kappa 1.0 --dx 0.5 --max-level 3 --stop-time 8 --gpu 0
    wormhole_case.sh --kappa 0.7 --dx 0.5 --gpu 1 --no-frames
    wormhole_case.sh --kappa 1.0 --dx 0.5 --box-size 96 --gpu 0   # larger box
    wormhole_case.sh --kappa 1.0 --dx 0.5 --mass 0.1 --gpu 0      # confined
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

GRTECLYN_ROOT = Path(
    os.environ.get("GRTECLYN_ROOT", "/home/jovyan/nachevsky/test/simulation/GRTeclyn")
).resolve()
EXAMPLE_DIR = GRTECLYN_ROOT / "Examples" / "RotatingWormholeCollapse"
ID_ROOT = GRTECLYN_ROOT / "runs" / "rotating_wormhole_id"
RUN_ROOT = GRTECLYN_ROOT / "runs" / "rotating_wormhole"
BIN = EXAMPLE_DIR / "main3d.gnu.MPI.CUDA.ex"
WRAPPER_SRC = Path(__file__).resolve().parents[3] / "src"

DEFAULT_L = 64.0
FRAME_FIELDS = "chi chi_minus_1 K lapse phi Pi phi2 Pi2 Weyl4_Re Weyl4_Im Weyl4_Mag"


def _p(x: float) -> str:
    """GRTeclyn-run tag formatting: 0.05 -> '0p05', 1.0 -> '1p00'."""
    return f"{x:.2f}".replace(".", "p")


def _dx_tag(dx: float) -> str:
    return f"dx{dx:.3g}".replace(".", "p")


def _L_suffix(L: float) -> str:
    """Append _L<val> only when L differs from the default (backward compat)."""
    if abs(L - DEFAULT_L) < 0.1:
        return ""
    return f"_L{int(L)}"


def _mass_suffix(mass: float) -> str:
    """Append _mass<val> only when mass>0 (backward compat with massless runs)."""
    if mass <= 0.0:
        return ""
    return f"_mass{mass:g}".replace(".", "p")


def _qball_suffix(lam: float, mu6: float) -> str:
    """Append _qball_lam<val>_mu6<val> only for a self-interacting (Q-ball)
    profile (lambda>0), keeping massless/mass-only families separate.  Must
    match solve_kappa_family.py's tag."""
    if lam <= 0.0:
        return ""
    return f"_qball_lam{lam:g}_mu6{mu6:g}".replace(".", "p")


def _constellation_suffix(n_lumps: int, orbit_radius: float,
                          orbit_omega: float) -> str:
    """Append _nlump<N>_R<R0>_worb<omega> for a multi-lump constellation ID.
    Must match solve_kappa_family.py's _constellation_suffix."""
    if n_lumps <= 1:
        return ""
    return (f"_nlump{n_lumps}_R{orbit_radius:g}_worb{orbit_omega:g}"
            .replace(".", "p"))


def id_tag(kappa: float, omega: float, m: int, dx: float,
           L: float = DEFAULT_L, mass: float = 0.0,
           lam: float = 0.0, mu6: float = 0.0,
           id_n_lumps: int = 1, id_orbit_radius: float = 6.0,
           id_orbit_omega: float = 0.1) -> str:
    # Must match solve_kappa_family.py's run_dir tag.
    return (f"rotwh_omega_p{_p(omega)}_m{m}_kappa_{_p(kappa)}"
            f"_{_dx_tag(dx)}{_L_suffix(L)}{_mass_suffix(mass)}"
            f"{_qball_suffix(lam, mu6)}"
            f"{_constellation_suffix(id_n_lumps, id_orbit_radius, id_orbit_omega)}")


def case_tag(kappa: float, omega: float, m: int, dx: float,
             max_level: int, L: float = DEFAULT_L, mass: float = 0.0,
             lam: float = 0.0, mu6: float = 0.0,
             id_n_lumps: int = 1, id_orbit_radius: float = 6.0,
             id_orbit_omega: float = 0.1) -> str:
    return (f"evo_omega_p{_p(omega)}_m{m}_kappa_{_p(kappa)}"
            f"_{_dx_tag(dx)}_ml{max_level}{_L_suffix(L)}{_mass_suffix(mass)}"
            f"{_qball_suffix(lam, mu6)}"
            f"{_constellation_suffix(id_n_lumps, id_orbit_radius, id_orbit_omega)}")


def extraction_radii(L: float, override: "list[float] | None" = None) -> tuple[float, ...]:
    """Weyl4 / Psi4 spherical-extraction shells.

    Default is the near-zone/wave-zone detector pair (r=12, r=24): close enough
    that the burst reaches r=24 within a ~t=45 run, while the box (L>=64) keeps
    the sponge + Sommerfeld boundary far outside so neither detector sees a
    reflection during the run.  Pass ``override`` (from ``--extraction-radii``)
    for a custom shell set, e.g. ``12 18 24`` for a denser radial ladder that
    lets you fit the 1/r fall-off and check propagation speed between shells.

    Radii are validated against the box: any shell at or beyond ``L/2 - 2`` (the
    Sommerfeld boundary buffer, i.e. inside/behind the sponge) is dropped so we
    never extract on a contaminated surface.
    """
    if override:
        radii = sorted({float(r) for r in override})
        r_max_ok = L / 2.0 - 2.0
        kept = [r for r in radii if 0.0 < r < r_max_ok]
        if not kept:
            raise SystemExit(
                f"--extraction-radii {override} has no shell inside the clean "
                f"zone (0, {r_max_ok:g}) for box L={L}.")
        if len(kept) < len(radii):
            dropped = [r for r in radii if r not in kept]
            print(f"  [extraction] dropping radii outside clean zone: {dropped}")
        return tuple(kept)

    r_inner = 12.0
    r_outer = 24.0
    # Guard against pathologically small boxes that cannot contain r=24 + sponge.
    if L / 2.0 - 2.0 <= r_outer:
        r_outer = max(L / 2.0 - 8.0, r_inner + 4.0)
    return (r_inner, r_outer)


def default_stop_time(radii: "tuple[float, ...]") -> float:
    """stop_time = r_outer + 6: light-crossing to the outermost extraction
    shell plus a buffer for signal development."""
    return radii[-1] + 6.0


def regrid_interval_str(max_level: int) -> str:
    """AMReX ParmParse wants one regrid_interval entry per level.

    max_level=0 is unigrid (no regridding). For AMR, coarse levels regrid every
    32 steps, finer levels every 16 (matches SupportedWormholeCollapse params).
    """
    if max_level <= 0:
        return "0"
    base = [32, 32, 16, 16, 16, 16]
    vals = (base + [16] * max_level)[:max_level]
    return " ".join(str(v) for v in vals)


PARAMS_TEMPLATE = """\
# Rotating-wormhole evolution (generated by wormhole_case.py; do not hand-edit).
# GRTresna-solved complex phantom-scalar ID, m={m}, omega={omega}, kappa={kappa}.
# Level-0 dx={dx} (N={N}, L={L}) == gridinit native dx (L2). max_level={max_level}.
verbosity = 1

output_path = "{run_out}"
amr.check_file = "{run_out}/RotatingWormholeChk"
amr.plot_file  = "{run_out}/RotatingWormholePlt"

checkpoint_interval = -1
plot_interval       = {plot_interval}

amr.plot_vars = chi K lapse shift1 shift2 shift3 phi Pi phi2 Pi2
amr.derive_plot_vars = Weyl4

hdf5_subpath = "hdf5"
pout_subpath = "pout"
data_subpath = "data"

print_progress_only_to_rank_0 = 1

recipe_initial_data_file = "{gridinit}"
wormhole_matter_model = "complex_scalar"

wormhole_initial_lapse_type = {lapse_type}
wormhole_throat_radius = 0.5
wormhole_centerA = 0.0 0.0 0.0
wormhole_phi_monopole_amplitude     = 0.0
wormhole_phi_perturbation_amplitude = 0.0
wormhole_phi_perturbation_width     = 0.5
wormhole_support_strength = 1.0
phantom_mass = {mass}
phantom_lambda = {lam}
phantom_mu6 = {mu6}

wormhole_azimuthal_m   = {m}
wormhole_rotation_omega = {omega}
{trajectory_block}
# Support-strength ramp (Phase 6 collapse-on-command trigger).  Ramps the
# exotic stress-energy coupling from base -> base*floor over [t_start,t_end];
# t_start<0 => never (support held constant).
support_ramp_t_start = {support_ramp_t_start}
support_ramp_t_end   = {support_ramp_t_end}
support_ramp_floor   = {support_ramp_floor}

# Numerical sponge zone (outer-shell KO dissipation for clean GW extraction).
sponge_enabled      = {sponge_enabled}
sponge_inner_radius = {sponge_inner}
sponge_outer_radius = {sponge_outer}
sponge_strength     = {sponge_strength}
sponge_ramp_power   = 4
sponge_center       = {center_x} {center_y} {center_z}

L = {L}
N1 = {N}
N2 = {N}
N3 = {N}

center = {center_x} {center_y} {center_z}

max_level = {max_level}
regrid_interval = {regrid_interval}
regrid_threshold = 0.01

max_box_size = 32
min_box_size = 16

isPeriodic = 0 0 0
hi_boundary = 1 1 1
lo_boundary = {lo_boundary}

nonzero_asymptotic_vars = chi h11 h22 h33 lapse
nonzero_asymptotic_values = 1.0 1.0 1.0 1.0 1.0

dt_multiplier = 0.02
stop_time = {stop_time}

max_spatial_derivative_order = 4
nan_check = 1

min_chi = 1.0e-8
min_lapse = 1.0e-10

lapse_advec_coeff = 1.0
lapse_coeff = 2.0
lapse_power = 1.0

shift_advec_coeff = 1.0
shift_Gamma_coeff = 0.75
eta = 1.0

formulation = 0
kappa1 = 3.0
kappa2 = 0.
kappa3 = 1.
covariantZ4 = 0

sigma = 2.0

puncture_tracking.enabled = 0
calculate_constraint_norms = 1

extraction_center = {center_x} {center_y} {center_z}
activate_extraction = 1
# write_extraction=1 dumps the full raw Weyl4 surface (~100 KB) as one file per
# extraction step -- thousands of files over a long dense run.  Default 0: keep
# only the compact appended Weyl4_mode_*.dat spherical-harmonic time series
# (the clean GW signal).  Enable with --write-extraction-surfaces if you need
# the raw spheres for e.g. a full 2D angular reconstruction.
write_extraction = {write_extraction}
extraction_subpath = "data"

num_extraction_radii = {num_extraction_radii}
extraction_radii = {extraction_radii_str}
extraction_levels = {extraction_levels_str}

num_points_phi = 24
num_points_theta = 37
num_modes = 3
modes = 2 0 \\
        2 1 \\
        2 2

amrex.the_arena_init_size = 1048576
amrex.use_gpu_aware_mpi = 0
"""


def build_trajectory_block(args, center) -> str:
    """Render the trajectory-pump params (Rung 2 active support).

    Returns an empty string when --num-lumps <= 0 (pump OFF, no-op regression:
    the generated params contain no trajectory_* keys, so the run is identical
    to a passive evolution).

    Places `num_lumps` orbit sites equally spaced on a circle of radius
    `orbit_radius` in the z=0 plane about the throat, each with tangential
    angular velocity `orbit_omega` (so the pump orbit matches the ID lump
    placement).  The PD trap (k_p>0) drives the field toward the moving target
    soliton; the pump frequency should equal the field's U(1) phase rate omega.
    """
    n = args.num_lumps
    if n <= 0:
        return ""
    cx, cy, cz = center
    lines = [
        "",
        "# --- Trajectory-guided matter pump (Rung 2 active support) ---",
        "trajectory_mode = 1",
        f"trajectory_num_lumps = {n}",
        f"trajectory_well_width = {args.well_width}",
        f"trajectory_pump_frequency = {args.pump_frequency}",
        "trajectory_A_breath = 0.0",
        "trajectory_omega_breath = 0.0",
        "trajectory_z_amp = 0.0",
        "trajectory_omega_z = 0.0",
        f"rl_pump_width = {args.well_width}",
        f"rl_pump_kp = {args.pump_kp}",
        f"rl_pump_kd = {args.pump_kd}",
        f"rl_pump_target_profile = {args.pump_target_profile}",
        f"rl_pump_target_width = {args.pump_target_width}",
        f"rl_pump_target_amp = {args.pump_target_amp}",
        f"rl_l2_ham_governor_center = {args.governor_center}",
        f"rl_l2_ham_governor_width = {args.governor_width}",
        f"rl_pump_max_amplitude = {args.well_depth}",
        f"pump_ramp_t_start = {args.pump_ramp_t_start}",
        f"pump_ramp_t_end = {args.pump_ramp_t_end}",
        f"pump_ramp_floor = {args.pump_ramp_floor}",
    ]
    import math
    for k in range(n):
        phase = 2.0 * math.pi * k / n
        lines += [
            f"trajectory_lump{k}_R0 = {args.orbit_radius}",
            f"trajectory_lump{k}_omega_rot = {args.orbit_omega}",
            f"trajectory_lump{k}_phase0 = {phase:.6f}",
            f"trajectory_lump{k}_tilt_theta = 0.0",
            f"trajectory_lump{k}_tilt_phi = 0.0",
            f"trajectory_lump{k}_well_depth = {args.well_depth}",
            f"trajectory_lump{k}_v_rad = 0.0",
        ]
    return "\n".join(lines) + "\n"


def _consumer_env() -> dict:
    env = dict(os.environ)
    env["PYTHONPATH"] = f"{WRAPPER_SRC}:{env.get('PYTHONPATH', '')}"
    # The consumer is CPU-only; keep it off the evolution GPU.
    env.pop("CUDA_VISIBLE_DEVICES", None)
    return env


# Module-level state set by main() before any consumer call.
_run_center: tuple[float, float, float] = (32.0, 32.0, 0.0)
_run_radii: tuple[float, ...] = (12.0, 24.0)
_run_L: float = DEFAULT_L


def _consumer_cmd(run_out: Path, jobs: int, *, watch: bool,
                  delete: bool, keep_last: int, frames: bool = True,
                  keep_existing_frames: bool = False) -> list[str]:
    cx, cy, cz = _run_center
    cmd = [
        sys.executable, "-m",
        "grteclyn_wrapper.visualisation.process_wave.consume_plotfiles",
        "--data", str(run_out), "--out", str(run_out / "small_data"),
        "--center", str(cx), str(cy), str(cz),
        "--radii", *(str(r) for r in _run_radii), "--n-points", "128",
        "--confinement-timeseries", "--confinement-well-width", "2.5",
        "-j", str(jobs),
    ]
    if frames:
        # Frame zoom: cover the inner ~75% of the box so the boundary is cropped.
        zoom = _run_L * 0.75
        # Frame center offset from the lower-left corner (the consumer's
        # --frames-corner mode places origin at (0,0,0) of the plotfile).
        fc_x = cx - zoom / 2.0
        fc_y = cy - zoom / 2.0
        cmd += [
            "--frames-fields", *FRAME_FIELDS.split(),
            "--frames-axis", "z", "--frames-zoom", str(zoom),
            "--frames-center", str(fc_x), str(fc_y), "0", "--frames-corner",
            # Slice at the run center's z-plane (z=cz). For a full-box centered
            # object (Q-torus, cz=L/2) the matter/throat lives at z=L/2, NOT at
            # the z=0 domain edge that --frames-corner defaults to -- slicing at
            # z=0 renders the matter fields (phi, Pi) as blank white because they
            # are ~0 there. For the half-z octant (cz=0) this is a no-op.
            "--frames-coord", str(cz),
            "--frames-global-zlim", "--frames-out", str(run_out / "frames"),
        ]
    if delete:
        cmd += ["--delete", "--keep-last", str(keep_last)]
    if watch:
        cmd += ["--watch", "--poll-seconds", "15"]
    if keep_existing_frames:
        # Consumer deletes existing frames at startup by default; a post-run
        # drain MUST preserve the frames the live sidecar already rendered.
        cmd += ["--keep-existing-frames"]
    return cmd


def start_frame_consumer(run_out: Path, jobs: int,
                         keep_last: int = 3,
                         frames: bool = True) -> subprocess.Popen:
    """Streaming sidecar (L6 disk discipline): extract Psi4 + confinement and
    delete plotfiles *during* the evolution, keeping only the newest `keep_last`
    for safety.  When ``frames`` is True it also renders PNG SlicePlots; when
    False it is a delete-only sidecar (still enforces L6 for --no-frames test
    runs, so plotfiles never accumulate).
    """
    return subprocess.Popen(
        _consumer_cmd(run_out, jobs, watch=True, delete=True,
                      keep_last=keep_last, frames=frames),
        env=_consumer_env(),
    )


def render_frames(run_out: Path, jobs: int) -> None:
    """One-shot render that keeps plotfiles (used only with --keep-plotfiles)."""
    subprocess.run(
        _consumer_cmd(run_out, jobs, watch=False, delete=False, keep_last=0),
        env=_consumer_env(), check=True,
    )


def drain_frames(run_out: Path, jobs: int, frames: bool = True) -> None:
    """One-shot post-run pass: process any stragglers the sidecar didn't reach
    (its keep_last tail) and delete every remaining plotfile.  When ``frames``
    is False it is a delete-only drain (no PNG rendering)."""
    subprocess.run(
        _consumer_cmd(run_out, jobs, watch=False, delete=True, keep_last=0,
                      frames=frames, keep_existing_frames=True),
        env=_consumer_env(), check=False,
    )


def prune_plotfiles(run_out: Path) -> None:
    for d in list(run_out.glob("RotatingWormholePlt*")) + list(
        run_out.glob("RotatingWormholeChk*")
    ):
        if d.is_dir():
            subprocess.run(["rm", "-rf", str(d)], check=False)


def main() -> int:
    global _run_center, _run_radii, _run_L

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--kappa", type=float, required=True)
    ap.add_argument("--dx", type=float, default=0.5, choices=[0.5, 1.0])
    ap.add_argument("--omega", type=float, default=0.05)
    ap.add_argument("--m", type=int, default=1)
    ap.add_argument("--max-level", type=int, default=3)
    ap.add_argument("--stop-time", type=float, default=None,
                    help="evolution stop time (default: r_outer + 6)")
    ap.add_argument("--box-size", type=float, default=DEFAULT_L,
                    help="domain side length L (default 64); larger box for "
                         "boundary convergence study")
    ap.add_argument("--mass", type=float, default=0.0,
                    help="phantom scalar field mass mu (default 0 = massless "
                         "ghost). >0 adds a confining potential V=0.5 mu^2|Phi|^2 "
                         "so the throat keeps its support; must match the ID "
                         "solved with MASS=<same> solve_kappa_family.sh")
    ap.add_argument("--lambda", type=float, default=0.0, dest="lam",
                    help="phantom scalar attractive-quartic coupling lambda "
                         "(default 0). Together with --mu6>0 this binds the "
                         "field into a Q-ball whose confinement survives the "
                         "phantom sign flip; V=0.5 m^2|Phi|^2 - 0.25 lambda|Phi|^4"
                         " + (1/6) mu6|Phi|^6. Must match the ID solved with "
                         "LAMBDA=<same> MU6=<same> solve_kappa_family.sh")
    ap.add_argument("--mu6", type=float, default=0.0,
                    help="phantom scalar sextic stabiliser mu6 (default 0). "
                         "REQUIRED (>0) for a stable 3D Q-ball with --lambda>0.")
    # --- Trajectory-guided matter pump (Rung 2 active support) ----------
    # --num-lumps 0 (default) => pump OFF; the generated params contain no
    # trajectory_* keys (no-op regression, identical to a passive run).
    ap.add_argument("--num-lumps", type=int, default=0,
                    help="number of orbiting pump sites (0 = pump OFF). >0 "
                         "enables the trajectory-guided matter pump: N lumps "
                         "equally spaced on a circle of radius --orbit-radius "
                         "orbiting at --orbit-omega, PD-trapped toward the "
                         "moving target soliton.")
    ap.add_argument("--orbit-radius", type=float, default=6.0,
                    help="pump orbit radius about the throat (default 6).")
    ap.add_argument("--orbit-omega", type=float, default=0.1,
                    help="pump orbital angular velocity (default 0.1).")
    ap.add_argument("--well-depth", type=float, default=0.05,
                    help="per-lump pump amplitude / max amplitude (default 0.05).")
    ap.add_argument("--well-width", type=float, default=2.5,
                    help="pump spotlight Gaussian sigma (default 2.5).")
    ap.add_argument("--pump-kp", type=float, default=0.0,
                    help="PD trap proportional gain (0 => legacy open-loop "
                         "source). >0 drives the field toward the target soliton.")
    ap.add_argument("--pump-kd", type=float, default=0.0,
                    help="PD trap derivative gain (default 0).")
    ap.add_argument("--pump-frequency", type=float, default=0.0,
                    help="pump U(1) phase rate (should match --omega for "
                         "coherent charge injection; default 0).")
    ap.add_argument("--pump-target-profile", type=int, default=0,
                    choices=[0, 2], help="PD target shape: 0 Gaussian, 2 sech.")
    ap.add_argument("--pump-target-width", type=float, default=0.0,
                    help="PD target physical width (<=0 => use well width).")
    ap.add_argument("--pump-target-amp", type=float, default=0.0,
                    help="PD target central amplitude (<=0 => use well depth). "
                         "Set to the ID lump central amplitude (phi_c).")
    ap.add_argument("--governor-center", type=float, default=0.035,
                    help="L2_Ham governor tanh center (pump cuts off above).")
    ap.add_argument("--governor-width", type=float, default=0.003,
                    help="L2_Ham governor tanh width.")
    ap.add_argument("--pump-ramp-t-start", type=float, default=-1.0,
                    help="collapse-trigger ramp start time (<0 => never ramp).")
    ap.add_argument("--pump-ramp-t-end", type=float, default=0.0,
                    help="collapse-trigger ramp end time.")
    ap.add_argument("--pump-ramp-floor", type=float, default=0.0,
                    help="pump amplitude fraction after the ramp (0 = full cut).")
    # --- Support-strength ramp (Phase 6 collapse-on-command trigger) ----
    ap.add_argument("--support-ramp-t-start", type=float, default=-1.0,
                    help="support-strength ramp start time (<0 => never ramp). "
                         "Ramps wormhole_support_strength (exotic stress-energy "
                         "coupling) down to trigger collapse -- the rotating "
                         "analogue of the static Ellis-Bronnikov S_support cut.")
    ap.add_argument("--support-ramp-t-end", type=float, default=0.0,
                    help="support-strength ramp end time.")
    ap.add_argument("--support-ramp-floor", type=float, default=0.0,
                    help="support-strength fraction after the ramp (0 = full "
                         "cut of the exotic support => throat collapses).")
    # --- Numerical sponge zone (clean GW extraction on a large box) -----
    ap.add_argument("--sponge", action="store_true",
                    help="enable the outer-shell sponge (extra KO dissipation) "
                         "to absorb outgoing waves before the boundary.")
    ap.add_argument("--sponge-inner", type=float, default=None,
                    help="sponge inner radius (default: r_outer_extraction + 4).")
    ap.add_argument("--sponge-outer", type=float, default=None,
                    help="sponge outer radius (default: L/2 - 2).")
    ap.add_argument("--sponge-strength", type=float, default=4.0,
                    help="extra sigma at full ramp (default 4.0).")
    # --- Constellation ID selection (Phase 2/3) -------------------------
    # Which .gridinit to load.  Defaults to the pump orbit geometry so a
    # pumped constellation run (--num-lumps N ...) auto-loads the matching
    # constellation ID.  For a PASSIVE constellation run (pump off) set
    # --id-num-lumps N alone (leave --num-lumps 0).  id-num-lumps<=1 loads a
    # single-throat (winding) ID as before.
    ap.add_argument("--id-num-lumps", type=int, default=None,
                    help="number of lumps in the ID gridinit to load "
                         "(default: --num-lumps). Set alone for a passive "
                         "constellation run with the pump off.")
    ap.add_argument("--id-orbit-radius", type=float, default=None,
                    help="ID constellation orbit radius (default: --orbit-radius).")
    ap.add_argument("--id-orbit-omega", type=float, default=None,
                    help="ID constellation orbit omega (default: --orbit-omega).")
    ap.add_argument("--initial-lapse-type", type=int, default=0,
                    choices=[0, 1, 2, 3],
                    help="lapse seeding for the loaded ID (default 0 = use the "
                         "GRTresna lapse as-is). 1=sqrt(chi), 2=1-3log(chi), "
                         "3=chi: a precollapsed lapse that damps the t~0.5 "
                         "max_K gauge transient which kicks the throat off "
                         "equilibrium.")
    ap.add_argument("--plot-interval", type=int, default=40,
                    help="coarse steps between plotfiles (=> frame/movie "
                         "cadence). In-code Psi4 extraction is INDEPENDENT of "
                         "this (it runs every coarse step), so use a moderate "
                         "value for smooth movies without a plotfile flood.")
    ap.add_argument("--extraction-radii", type=float, nargs="+", default=None,
                    metavar="R",
                    help="Weyl4/Psi4 spherical-extraction shell radii in code "
                         "units (e.g. --extraction-radii 12 18 24). Default: "
                         "12 24. Shells at/behind the sponge (>= L/2-2) are "
                         "dropped. Feeds both the in-code extraction and the "
                         "plotfile-sidecar Psi4.")
    ap.add_argument("--write-extraction-surfaces", action="store_true",
                    help="also dump the raw Weyl4 surface data (one file per "
                         "extraction step). Off by default; the compact "
                         "Weyl4_mode_*.dat time series is always written.")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--gridinit", default=None, help="override ID .gridinit path")
    ap.add_argument("--full-box", action="store_true",
                    help="full-box centered evolution (center z=L/2, all-Sommerfeld "
                         "lo_boundary 1 1 1) instead of the default half-z octant "
                         "(center z=0, z-reflective 1 1 2). Use for a centered, "
                         "z-symmetric object like an isolated Q-torus solved with "
                         "solve_torus.py (target_center=(L/2,L/2,L/2)).")
    ap.add_argument("--frames", dest="frames", action="store_true", default=True)
    ap.add_argument("--no-frames", dest="frames", action="store_false")
    # NFS plotfiles are large (~GB at N=128/ml=3); -j>4 thrashes the network
    # and hangs (see scripts/plot/drain_plotfiles_frames.sh). Keep it modest,
    # especially when several cases render concurrently.
    ap.add_argument("--frame-jobs", type=int, default=4)
    ap.add_argument("--keep-plotfiles", action="store_true",
                    help="do not prune plotfiles after frame extraction")
    ap.add_argument("--run-suffix", default="",
                    help="extra suffix appended to the run-dir tag (does NOT "
                         "affect the ID tag). Use to keep parallel pump-gain "
                         "sweeps (same ID, different k_p) in distinct dirs.")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    L = args.box_size
    N = int(round(L / args.dx))
    center = (L / 2.0, L / 2.0, L / 2.0 if args.full_box else 0.0)
    lo_boundary = "1 1 1" if args.full_box else "1 1 2"
    radii = extraction_radii(L, args.extraction_radii)
    stop_time = args.stop_time if args.stop_time is not None else default_stop_time(radii)

    # Publish to module state so the consumer helpers pick them up.
    _run_center = center
    _run_radii = radii
    _run_L = L

    # Constellation ID selection: default to the pump orbit geometry so a
    # pumped constellation auto-loads its matching ID; overridable for a
    # passive (pump-off) constellation run via --id-num-lumps.
    id_n_lumps = args.id_num_lumps if args.id_num_lumps is not None else args.num_lumps
    id_orbit_radius = (args.id_orbit_radius if args.id_orbit_radius is not None
                       else args.orbit_radius)
    id_orbit_omega = (args.id_orbit_omega if args.id_orbit_omega is not None
                      else args.orbit_omega)

    gridinit = (
        Path(args.gridinit).resolve()
        if args.gridinit
        else ID_ROOT / id_tag(args.kappa, args.omega, args.m, args.dx, L,
                              args.mass, args.lam, args.mu6,
                              id_n_lumps, id_orbit_radius, id_orbit_omega)
        / "initial_data.gridinit"
    )
    if not gridinit.is_file():
        hint = " ".join(
            ([f"EVO_L={int(L)} RES_N={N}"] if L != DEFAULT_L else [f"RES_N={N}"])
            + ([f"MASS={args.mass:g}"] if args.mass > 0 else [])
            + ([f"LAMBDA={args.lam:g} MU6={args.mu6:g}"] if args.lam > 0 else [])
        )
        print(f"error: gridinit not found: {gridinit}\n"
              f"       solve it first: {hint} solve_kappa_family.sh {args.kappa}",
              file=sys.stderr)
        return 2
    if not BIN.exists():
        print(f"error: CUDA binary not found: {BIN}", file=sys.stderr)
        return 2

    tag = case_tag(args.kappa, args.omega, args.m, args.dx, args.max_level, L,
                   args.mass, args.lam, args.mu6,
                   id_n_lumps, id_orbit_radius, id_orbit_omega)
    if args.run_suffix:
        tag = f"{tag}_{args.run_suffix}"
    run_dir = RUN_ROOT / tag
    run_out = run_dir / "output"
    run_out.mkdir(parents=True, exist_ok=True)

    cx, cy, cz = center
    radii_str = " ".join(f"{r:.1f}" for r in radii)
    levels_str = " ".join("0" for _ in radii)

    # Sponge shell defaults: inner just outside the outer extraction radius,
    # outer just inside the boundary, so it absorbs waves after they pass the
    # detectors but before they reflect off the Sommerfeld boundary.
    sponge_inner = (args.sponge_inner if args.sponge_inner is not None
                    else radii[-1] + 4.0)
    sponge_outer = (args.sponge_outer if args.sponge_outer is not None
                    else L / 2.0 - 2.0)

    # GRTeclyn resolves relative paths against its cwd (the example dir); use the
    # absolute run path so the generated params are location-independent.
    params_text = PARAMS_TEMPLATE.format(
        m=args.m, omega=args.omega, kappa=args.kappa, dx=args.dx, N=N, L=L,
        mass=args.mass, lam=args.lam, mu6=args.mu6,
        lapse_type=args.initial_lapse_type,
        center_x=cx, center_y=cy, center_z=cz, lo_boundary=lo_boundary,
        max_level=args.max_level, stop_time=stop_time,
        plot_interval=args.plot_interval, run_out=str(run_out),
        gridinit=str(gridinit),
        regrid_interval=regrid_interval_str(args.max_level),
        num_extraction_radii=len(radii),
        extraction_radii_str=radii_str,
        extraction_levels_str=levels_str,
        write_extraction=(1 if args.write_extraction_surfaces else 0),
        trajectory_block=build_trajectory_block(args, center),
        support_ramp_t_start=args.support_ramp_t_start,
        support_ramp_t_end=args.support_ramp_t_end,
        support_ramp_floor=args.support_ramp_floor,
        sponge_enabled=(1 if args.sponge else 0),
        sponge_inner=sponge_inner, sponge_outer=sponge_outer,
        sponge_strength=args.sponge_strength,
    )
    params_path = run_dir / "params.txt"
    params_path.write_text(params_text)
    print(f"[wormhole_case] {tag}")
    print(f"  gridinit: {gridinit}")
    print(f"  params:   {params_path}")
    print(f"  L={L} N={N} dx={args.dx} max_level={args.max_level} "
          f"mass={args.mass} lambda={args.lam} mu6={args.mu6} "
          f"stop_time={stop_time} radii={radii} gpu={args.gpu}")
    if args.dry_run:
        return 0

    env = dict(os.environ)
    env["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    log = run_dir / "evolve.log"

    # Stream Psi4 + confinement and delete plotfiles *during* the run so the
    # heavy HDF5 backlog never accumulates (L6, disk discipline).  A delete
    # sidecar ALWAYS runs unless --keep-plotfiles, regardless of --frames:
    # --no-frames skips only the PNG rendering, NOT the plotfile deletion (a
    # bare --no-frames must never leave tens of GB of Plt* behind).  Only the
    # newest few are retained live; a post-run drain handles that tail.
    consumer = None
    if not args.keep_plotfiles:
        mode = "frames+Psi4" if args.frames else "delete-only (no frames)"
        print(f"  starting plotfile sidecar: {mode}, stream + delete (L6)")
        consumer = start_frame_consumer(run_out, args.frame_jobs,
                                        frames=args.frames)

    print(f"  running -> {log}")
    run_rc = 0
    try:
        with log.open("w") as fh:
            proc = subprocess.run([str(BIN), str(params_path)],
                                  cwd=str(EXAMPLE_DIR), env=env, stdout=fh,
                                  stderr=subprocess.STDOUT, check=False)
        run_rc = proc.returncode
        if run_rc != 0:
            # A NaN/abort is a legitimate physics outcome (e.g. collapse); do
            # NOT let it skip the plotfile cleanup below (L6: a crashed run must
            # never leave tens of GB of Plt* behind).
            print(f"  WARNING: binary exited {run_rc} (NaN/abort?) -- "
                  f"cleaning up plotfiles anyway", file=sys.stderr)
    finally:
        if consumer is not None:
            consumer.terminate()
            try:
                consumer.wait(timeout=120)
            except subprocess.TimeoutExpired:
                consumer.kill()
        # Plotfile disposal ALWAYS runs, success or crash (L6 disk discipline).
        if args.keep_plotfiles:
            if args.frames and run_rc == 0:
                print("  rendering frames (offline, keeping plotfiles)...")
                render_frames(run_out, args.frame_jobs)
        else:
            print(f"  draining remaining plotfiles ("
                  f"{'frames + ' if args.frames else ''}delete)...")
            drain_frames(run_out, args.frame_jobs, frames=args.frames)
            prune_plotfiles(run_out)
            print("  plotfiles streamed + pruned "
                  f"({'frames + ' if args.frames else ''}data kept)")

    print(f"[wormhole_case] done: {run_out}")
    if run_rc != 0:
        return run_rc
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
