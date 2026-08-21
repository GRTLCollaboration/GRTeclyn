#!/usr/bin/env bash
# Bondi dipole -- gravity-dressed PAIR cell, ultraweak rung.
#
# Two self-gravitating omega=0.55 sextic stars at rest, separation SEP along x
# (lump0 at +R0, lump1 at -R0).  Sector per lump via BONDI_S0/BONDI_S1
# (0=canonical, 1=phantom):
#   S0=0 S1=1 -> pair_pm  (the runaway cell: canonical at +x, phantom at -x;
#                          negative mass chases positive -> both drift +x)
#   S0=0 S1=0 -> pair_pp  (control: mutual attraction, no net drift)
#   S0=1 S1=1 -> pair_mm  (control: mutual repulsion, no net drift)
#
# Star data (bondi_dipole_debug.md §8): canonical phi_c=0.019695 ADM=+0.0640,
# phantom phi_c=0.019589 ADM=-0.0770.  |M-| > |M+| by ~20% -- expect imperfect
# momentum cancellation in pair_pm.
#
# Launch-time verification (t=0 row of small_data/sector_barycenters.dat):
#   canonical lump present: total_canon ~= 15.92, rms_canon ~= 5.05
#   phantom  lump present: total_phantom ~= 20.99, rms_phantom ~= 5.43
# Runaway signal: canon x (col 3) and phantom x (col 8) both drifting the
# same direction, canonical side leading.
#
# Default stop 60, not 100: the massive-field radiation bath cannot exit the
# massless-wave boundaries, and chi drifted 0.99 -> 0.86 by t=40 in the
# singles -- beyond t~60 the bath starts to poison the drift measurement.
#
# Usage:
#   BONDI_S0=0 BONDI_S1=1 BONDI_GPU=3 bash scripts/campaigns/bondi_dipole/run_pair_selfgrav.sh
# Overrides: BONDI_GPU (default 3), BONDI_STOP_TIME (default 60),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8), BONDI_S0/BONDI_S1.
#   BONDI_NL_TOL / BONDI_NL_STALL_TOL: elliptic-solve stopping rule.  Scale
#     BONDI_NL_TOL as dx^4 across a ladder to get a convergence order out of it.
#   BONDI_SPONGE=1: switch on the boundary sponge (inner/outer/strength/ramp
#     via BONDI_SPONGE_INNER etc.) to run past t=60.
#   BONDI_DRYRUN=1: resolve and print the parameters, then exit -- no solve, no GPU.
#   BONDI_GRTRESNA_RANKS: MPI ranks for the elliptic solve (default 8).  The
#     evolution is single-GPU regardless.
#   BONDI_S0_OMEGA: base star frequency for BOTH lumps (default 0.55).  0.55
#     is on the UNSTABLE branch -- the canonical star collapses on its own by
#     t~60 (docs/MatterDebugg.md).  The stable branch starts at ~0.67; the
#     design point is 0.80.  A non-default value appends _wNNN to the run name.
#   BONDI_SCALAR_LAMBDA / BONDI_SCALAR_MU: potential rung (defaults are the
#     published 10240 / 21845333).
#   BONDI_S1_OMEGA: per-lump star frequency for lump1 (equal-|ADM| cells --
#     the phantom star at omega=0.56598 weighs the canonical star's 0.0640;
#     its t=0 fingerprint becomes total ~= 17.05 / rms ~= 5.16).  Appends
#     "_eqm" to the run name so the standard cell is never clobbered.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

# Site paths (GRTRESNA_ENV -> mpirun on PATH).  Needed once the solve runs on
# more than one rank: at a single rank the solver binary is exec'd directly and
# mpirun is never looked up, which is why this was not required before and why
# raising the rank count fails with "Could not locate 'mpirun'" without it.
# shellcheck source=../../lib/env.sh
source "${WRAPPER_DIR}/scripts/lib/env.sh"
# env.sh exports a SCRIPT_DIR of its own and cds to the repo root -- put ours
# back, or every sibling path below resolves against scripts/lib/ instead.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
GPU="${BONDI_GPU:-3}"
STOP_TIME="${BONDI_STOP_TIME:-60}"
SEP="${BONDI_SEP:-8}"
S0="${BONDI_S0:-0}"
S1="${BONDI_S1:-1}"
R0="$(python3 -c "print(${SEP}/2)")"
PI="3.141592653589793"

# Rerun knobs (Debug.md §3).  All default to the published campaign values, so
# an un-set environment reproduces the original cell bit-for-bit.
#   BONDI_NFULL   -- item A, convergence series (192, 256)
#   BONDI_LFULL   -- item C, enlarged domain (128, 256)
#   BONDI_RADII   -- item C, wave-zone extraction shells
#   BONDI_DT_MULT -- item I, CFL experiment (0.02 -> 0.2)
NFULL="${BONDI_NFULL:-128}"
LFULL="${BONDI_LFULL:-64}"
RADII="${BONDI_RADII:-8 16}"
DT_MULT="${BONDI_DT_MULT:-0.02}"
#   BONDI_MAXLEVEL -- item A: set 0 for the convergence series.  Two reasons.
#     (1) Richardson extrapolation needs a uniform grid; an AMR level that
#         switches on partway through makes the convergence order meaningless.
#     (2) Memory.  A refined N=128 cell measured 35 GB of FABs against 6 GB
#         unigrid; at N=256 that factor would OOM an 80 GB H100 outright.
#     The published t<=30 science window was unigrid anyway -- tagging did not
#     fire until t~47.5 (mixed) and never in the singles -- so max_level=0
#     reproduces it over the window that matters.
MAXLEVEL="${BONDI_MAXLEVEL:-1}"

# Elliptic-solve stopping rule.  The published campaign exited every rung of
# the resolution ladder on the same criterion, so all three stopped at the same
# 0.0832% Hamiltonian violation -- identical to four digits at N=128, 192, 256.
# The floor is this tolerance, not the grid, which is why the published pack
# quotes no convergence order.  Scale it as dx^4 across a ladder to measure one:
#   N=128 -> 0.1     N=192 -> 0.019     N=256 -> 0.00625
# The stall guard has to come down with it, or the iteration halts on "no
# further progress" before the tighter exit is ever reached.
NL_TOL="${BONDI_NL_TOL:-0.1}"
NL_STALL_TOL="${BONDI_NL_STALL_TOL:-0.002}"
# Tighter tolerances take more Newton iterations and so more wall clock; the
# published solves used 7 of a permitted 50 and finished well inside 7200 s.
GRTRESNA_TIMEOUT="${BONDI_GRTRESNA_TIMEOUT:-7200}"

# Numerical sponge zone (Source/Grids/SpongeZone.hpp): a radially-ramped band
# of extra Kreiss-Oliger dissipation, off by default.  Massive-scalar radiation
# cannot leave through the massless-wave boundaries, washes back over the stars
# and caps every mixed run at t=60; the sponge is what buys t >> 60.  Two
# things to check on the first cell rather than assume: it was validated
# against *massless* reflections elsewhere, and a band placed at 24/32 starts
# eating the canonical star's own halo once its rms radius crosses r=16
# (t ~ 51), which moves the domain-integrated barycentre.  Run BONDI_SCRUTINY=1
# alongside any sponge cell so the halo the sponge removes can be told apart
# from the halo the diagnostic was mis-counting.
SPONGE="${BONDI_SPONGE:-0}"
SPONGE_INNER="${BONDI_SPONGE_INNER:-24}"
SPONGE_OUTER="${BONDI_SPONGE_OUTER:-32}"
SPONGE_STRENGTH="${BONDI_SPONGE_STRENGTH:-4.0}"
SPONGE_RAMP="${BONDI_SPONGE_RAMP:-4}"

# Print the resolved parameter set and exit without solving or touching a GPU.
# Worth spending on a new cell: the alternative is discovering a mis-set knob
# after the solve has already burned an hour of CPU and the card is committed.
dryrun_args=()
if [[ "${BONDI_DRYRUN:-0}" != "0" ]]; then
  dryrun_args=(--dry-run)
fi

sponge_args=()
if [[ "${SPONGE}" != "0" ]]; then
  sponge_args=(
    --extra-override sponge_enabled=1
    --extra-override sponge_inner_radius="${SPONGE_INNER}"
    --extra-override sponge_outer_radius="${SPONGE_OUTER}"
    --extra-override sponge_strength="${SPONGE_STRENGTH}"
    --extra-override sponge_ramp_power="${SPONGE_RAMP}"
  )
fi

S1_OMEGA="${BONDI_S1_OMEGA:-}"

# Matter configuration -- nothing hard-coded below: frequency and potential
# rung come through the environment, defaulting to the published values so an
# un-set environment reproduces every packed cell bit-for-bit.
S0_OMEGA="${BONDI_S0_OMEGA:-0.55}"
SCALAR_LAMBDA="${BONDI_SCALAR_LAMBDA:-10240}"
SCALAR_MU="${BONDI_SCALAR_MU:-21845333}"

sector_tag() { if [[ "$1" == "1" ]]; then echo m; else echo p; fi; }
suffix="$(sector_tag "${S0}")$(sector_tag "${S1}")"
if [[ -n "${S1_OMEGA}" ]]; then
  suffix="${suffix}_eqm"
fi
if [[ "${S0_OMEGA}" != "0.55" ]]; then
  suffix="${suffix}_w$(printf '%s' "${S0_OMEGA}" | tr -d '.')"
fi
out_name="bondi_sg_pair_${suffix}"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi_dipole_selfgrav_${suffix}}"

mkdir -p "${RUNS_DIR}"

# Stop handle for scripts/campaigns/stop_campaign.sh.
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "${RUNS_DIR}"

# Live per-sector centroid stream + radiation extraction.
export GRTECLYN_SECTOR_BARYCENTERS=1
export GRTECLYN_PSI4=1
# Scrutiny stream (sector_dynamics.dat): halo-free core positions, per-sector
# matter momentum (the Bondi momentum-balance check) and a gauge check.  Costs
# ~1.6 s per plotfile, so opt in with BONDI_SCRUTINY=1.
export GRTECLYN_SECTOR_DYNAMICS="${BONDI_SCRUTINY:-0}"
export GRTECLYN_SECTOR_DYNAMICS_LEVEL="${BONDI_SCRUTINY_LEVEL:-0}"
# Higher-multipole Psi4 (l>=3) into its own stream psi4_mode_higher_l.dat.
# Item C: the odd-l/even-l interference carries the momentum-flux beaming along
# the runaway (x) axis.  Only meaningful once the extraction shells are in the
# wave zone (BONDI_RADII), so it is opt-in and off by default.
export GRTECLYN_PSI4_HIGHER_L="${BONDI_PSI4_HIGHER_L:-0}"
export GRTECLYN_PSI4_ELLS="${BONDI_PSI4_ELLS:-3 4}"

# The elliptic solve is the campaign's CPU bottleneck and it parallelises: at one
# rank a single N=256 solve pegs one core for ~70 min while 47 sit idle, and the
# GPU cannot start until it finishes.
#
# The "mpirun segfaults cluster-wide, do not raise" note that used to live here
# was inherited from the OLD node.  This node's MPI is healthy -- verified
# 2026-08-20 by launching 8 ranks, and the solver binary is an MPI build
# (Main_BosonStarBH3d...MPI.ex).  The EVOLUTION is a separate matter and stays
# single-rank: a RadialRecipe-specific AMR crash still blocks multi-GPU there.
#
# Ranks apply to the solve only.  Keep it modest: four cells solve concurrently,
# so the queue's total core draw is 4 x this number.
GRTRESNA_RANKS="${BONDI_GRTRESNA_RANKS:-8}"

if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
  echo "[bondi] ${out_name} already exists -- delete it or set BONDI_RUNS_DIR"
  exit 1
fi

echo "[bondi] === ${out_name} (ultraweak rung, gravity-dressed pair s0=${S0} s1=${S1}) ==="
PYTHONPATH="${WRAPPER_DIR}/src" "${WRAPPER_DIR}/.venv/bin/python" \
  "${WRAPPER_DIR}/scripts/campaigns/hq/replay_eval.py" \
  "${SOURCE_EVAL}" \
  --name "${out_name}" \
  --runs-dir "${RUNS_DIR}" \
  --gpu "${GPU}" \
  --n-full "${NFULL}" --l-full "${LFULL}" \
  --max-level "${MAXLEVEL}" --regrid-threshold 0.02 \
  --stop-time "${STOP_TIME}" --plot-interval 40 \
  --ftl-L 8.0 \
  --grtresna-ranks "${GRTRESNA_RANKS}" \
  --grtresna-iterations 50 \
  --grtresna-nl-exit-tolerance "${NL_TOL}" \
  --grtresna-nl-stall-tolerance "${NL_STALL_TOL}" \
  --grtresna-max-level 3 \
  --grtresna-domain-l 128 \
  --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
  --consumer-radii ${RADII} \
  --consumer-keep-last 2 \
  --objective-mode weighted \
  --extra-override dt_multiplier="${DT_MULT}" \
  --extra-override grtresna_scalar_lambda="${SCALAR_LAMBDA}" \
  --extra-override grtresna_scalar_mu="${SCALAR_MU}" \
  --extra-override grtresna_bs_omega="${S0_OMEGA}" \
  --extra-override grtresna_bs_selfgrav=1 \
  --extra-override trajectory_well_width=1.2 \
  --extra-override rl_pump_stop_time=0 \
  --extra-override grtresna_boost_lumps=0 \
  --extra-override trajectory_A_breath=0 \
  --extra-override trajectory_z_amp=0 \
  --extra-override trajectory_omega_z=0 \
  --extra-override trajectory_num_lumps=2 \
  --extra-override trajectory_lump0_R0="${R0}" \
  --extra-override trajectory_lump0_phase0=0 \
  --extra-override trajectory_lump0_tilt_theta=0 \
  --extra-override trajectory_lump0_tilt_phi=0 \
  --extra-override trajectory_lump0_v_rad=0 \
  --extra-override trajectory_lump0_omega_rot=0 \
  --extra-override trajectory_lump0_well_depth=0.15 \
  --extra-override trajectory_lump0_exotic="${S0}" \
  --extra-override trajectory_lump1_R0="${R0}" \
  --extra-override trajectory_lump1_phase0="${PI}" \
  --extra-override trajectory_lump1_tilt_theta=0 \
  --extra-override trajectory_lump1_tilt_phi=0 \
  --extra-override trajectory_lump1_v_rad=0 \
  --extra-override trajectory_lump1_omega_rot=0 \
  --extra-override trajectory_lump1_well_depth=0.15 \
  --extra-override trajectory_lump1_exotic="${S1}" \
  ${S1_OMEGA:+--extra-override trajectory_lump1_bs_omega="${S1_OMEGA}"} \
  ${sponge_args[@]+"${sponge_args[@]}"} \
  ${dryrun_args[@]+"${dryrun_args[@]}"}

echo "[bondi] pair cell complete: ${RUNS_DIR}/${out_name}"
