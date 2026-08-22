#!/usr/bin/env bash
# Bondi dipole -- gravity-dressed calibration: ONE canonical lump, ultraweak rung.
#
# Why this rung and this seed (bondi_dipole_debug.md §7, 2026-08-06 follow-up):
# the exact-amplitude FLAT Q-ball at the weak rung (lam=2560) still blew off
# (rms 5.22 -> 6.6 by t=20, gate broken) because at that rung's compactness
# (2E/R ~ 0.14) the gravitationally dressed heavy-branch star DOES NOT EXIST --
# the dressed family maps to feather-light stars (ADM ~ 0.05) or ultra-compact
# unstable ones (ADM ~ 2.3), nothing near the flat ball's E=0.28.  The lump had
# no equilibrium to settle into.
#
# Two rungs up (lam=10240, mu=21845333; lam^2/mu = 4.8 so omega_min = 0.316
# unchanged), the dressed omega=0.55 star exists: phi_c = 0.019695, ADM =
# 0.0640, alpha(0) = 0.977.  This run seeds THAT star via the fixed-frequency
# self-grav ODE solve (grtresna_bs_selfgrav=1; phi_c is a solved output, the
# 3-column table carries the lapse, and the painters' amp/phi_c rescale is
# exactly 1).  Expect breathing ~1%, not the +-10-18% of the flat seeds.
#
# Launch-time verification (t=0 row of small_data/sector_barycenters.dat):
#   total_canon ~= 15.92, rms ~= 5.05  -> dressed star painted
# Ground truth earlier: grtresna/params.txt lump0_amp ~= 0.019695 and
# bs_omega ~= 0.55; qball_profile.dat header says "self-gravitating".
#
# Pass gate unchanged: rms flat (+-10%) to t=40 AND min_chi plateau above 0.3.
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_single_selfgrav.sh
# Overrides: BONDI_GPU (default 1), BONDI_STOP_TIME (default 40),
#            BONDI_RUNS_DIR, BONDI_SEP (default 8), BONDI_EXOTIC=1 (lone
#            phantom star single_m -- reads the PHANTOM columns downstream).
#   BONDI_OMEGA: star frequency (default 0.55 -- the UNSTABLE branch; the
#     stable branch starts at ~0.67, design point 0.80).  A non-default value
#     appends _wNNN to the run name so the published cells are never clobbered.
#   BONDI_SCALAR_LAMBDA / BONDI_SCALAR_MU: potential rung (published defaults).
#   BONDI_NL_TOL / BONDI_NL_STALL_TOL / BONDI_GRTRESNA_TIMEOUT: elliptic-solve
#     stopping rule, same semantics as the pair launcher.
#   BONDI_SPONGE=1 (+ BONDI_SPONGE_INNER/_OUTER/_STRENGTH/_RAMP): boundary
#     sponge -- required for any single pushed past t~60, or the trapped
#     massive-scalar bath poisons the tail of the run.
#   BONDI_GRTRESNA_RANKS: MPI ranks for the elliptic solve (default 8).
#   BONDI_DRYRUN=1: resolve and print the parameters, then exit.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_DIR="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
REPO_ROOT="$(cd -- "${WRAPPER_DIR}/.." && pwd)"

# Site paths (GRTRESNA_ENV -> mpirun on PATH) -- required for a multi-rank
# solve; see the identical block in run_pair_selfgrav.sh.
# shellcheck source=../../lib/env.sh
source "${WRAPPER_DIR}/scripts/lib/env.sh"
# env.sh exports a SCRIPT_DIR of its own and cds to the repo root -- put ours
# back, or every sibling path below resolves against scripts/lib/ instead.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

SOURCE_EVAL="${REPO_ROOT}/results/matter-first-automated-discovery-of-transient-spacetime-shortcuts/search/qball-trajectory-evolving-geodesic-shortcut-search/run/eval_000322"
RUNS_DIR="${BONDI_RUNS_DIR:-${REPO_ROOT}/runs/bondi/dipole_selfgrav_v1}"
GPU="${BONDI_GPU:-1}"
#   BONDI_EVO_RANKS -- MPI ranks for the GPU EVOLUTION (default 1).  When >1,
#     BONDI_GPU must be a comma list of the same length ("0,1").  The
#     RadialRecipe MPI+CUDA path is only known-good at max_level=0; under AMR
#     it segfaults at the first RK4 advance (FillPatchIterator) -- smoke-test
#     on a short stop_time before trusting it for a long run.
EVO_RANKS="${BONDI_EVO_RANKS:-1}"
STOP_TIME="${BONDI_STOP_TIME:-40}"
SEP="${BONDI_SEP:-8}"
EXOTIC="${BONDI_EXOTIC:-0}"

# Matter configuration -- nothing hard-coded below: frequency and potential
# rung come through the environment, defaulting to the published values so an
# un-set environment reproduces the original cell bit-for-bit.
OMEGA="${BONDI_OMEGA:-0.55}"
SCALAR_LAMBDA="${BONDI_SCALAR_LAMBDA:-10240}"
SCALAR_MU="${BONDI_SCALAR_MU:-21845333}"

# Elliptic-solve stopping rule (same semantics as the pair launcher).
NL_TOL="${BONDI_NL_TOL:-0.1}"
NL_STALL_TOL="${BONDI_NL_STALL_TOL:-0.002}"
GRTRESNA_TIMEOUT="${BONDI_GRTRESNA_TIMEOUT:-7200}"
R0="$(python3 -c "print(${SEP}/2)")"

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

# Elliptic-solve GRID -- the "solve end point".  Every value below is an
# environment knob defaulting to the published campaign value, so an un-set
# environment reproduces the original cell bit-for-bit.
#   BONDI_GRTRESNA_DOMAIN_L -- solve box width.  Deliberately WIDER than the
#     evolution box (128 vs 64) so the outer boundary condition sits far from
#     the star.
#   BONDI_GRTRESNA_N        -- solve cells per axis.
#   BONDI_GRTRESNA_MAXLEVEL -- refinement depth inside the solve.
#   BONDI_GRTRESNA_MAXIMAL_SLICING -- 1 => build the slice with K=0.
#
# WHY MAXIMAL SLICING IS A KNOB (measured 2026-08-21).  GRTresna picks the
# initial expansion rate with K = sign*sqrt(24 pi G rho), which is imaginary
# wherever rho < 0.  The wrapper therefore switches to the K=0
# York/Lichnerowicz path if and only if EXOTIC matter is present
# (grtresna/fit/motif.py: maximal_slicing = has_exotic or exotic_needed).  The
# consequence is an uncontrolled asymmetry between the two sectors: the phantom
# star is born on a slice at rest (max|K| ~ 3e-05) while the canonical star is
# born on one already contracting at max|K| ~ 1e-01 -- four orders of magnitude
# apart, for a plumbing reason rather than a physical one.  Only the canonical
# star then collapses.  Set this to 1 to build canonical matter the same way,
# so the two sectors differ only in the sign of the energy.
#
# WHY THIS IS A KNOB AND NOT A CONSTANT (measured 2026-08-21, omega=0.80
# canonical).  The solve cell count used to be hardwired to NFULL while the box
# was hardwired to 128, so the solve cell was ALWAYS exactly 2x the evolution
# cell -- at every resolution, with no way to close the gap.  The Hamiltonian
# constraint takes second derivatives, so the interpolation noise that leaves
# is amplified as 1/dx^2: refining the evolution grid 128 -> 192 made the t=0
# violation 1.85x WORSE (L2 5.63e-4 -> 1.02e-3, Linf 2.37e-3 -> 5.18e-3),
# matching 1.5^2 = 2.25 rather than the 0.20 that 4th-order convergence would
# give.  A resolution study on the old wiring is therefore uninterpretable:
# the finer run starts handicapped, so "the finer grid failed too" does NOT
# acquit resolution.
#
# To match cell sizes exactly set N = NFULL * (DOMAIN_L / LFULL); the centres
# then align and the transfer becomes a straight copy rather than an
# interpolation.  MAXLEVEL=0 additionally removes the solve's own refinement
# seams, which are the other candidate noise source.
GRTRESNA_DOMAIN_L="${BONDI_GRTRESNA_DOMAIN_L:-128}"
GRTRESNA_MAXLEVEL="${BONDI_GRTRESNA_MAXLEVEL:-3}"
GRTRESNA_N="${BONDI_GRTRESNA_N:-${NFULL}}"
GRTRESNA_MAXIMAL_SLICING="${BONDI_GRTRESNA_MAXIMAL_SLICING:-0}"
MAXIMAL_SLICING_FLAG=()
if [[ "${GRTRESNA_MAXIMAL_SLICING}" != "0" ]]; then
  MAXIMAL_SLICING_FLAG=(--grtresna-maximal-slicing)
fi

# BONDI_GRIDINIT -- reuse an already-solved initial_data.gridinit and skip the
# elliptic solve entirely.  The point is an evolution-only A/B: two cells that
# start from the SAME BYTES differ only in the evolution binary, so nothing in
# the comparison can be blamed on solver noise or a re-solve landing on a
# slightly different root.  Every BONDI_GRTRESNA_* knob is ignored when this is
# set -- the slice is whatever was solved before, not what those knobs describe.
GRIDINIT_FLAG=()
if [[ -n "${BONDI_GRIDINIT:-}" ]]; then
  if [[ ! -f "${BONDI_GRIDINIT}" ]]; then
    echo "[bondi] BONDI_GRIDINIT set but not a file: ${BONDI_GRIDINIT}" >&2
    exit 1
  fi
  GRIDINIT_FLAG=(--gridinit "${BONDI_GRIDINIT}")
  echo "[bondi] reusing solved initial data, skipping the GRTresna solve:"
  echo "[bondi]   ${BONDI_GRIDINIT}"
fi

mkdir -p "${RUNS_DIR}"

# Stop handle for scripts/campaigns/stop_campaign.sh.
source "${SCRIPT_DIR}/../lib/launcher_common.sh"
campaign_register_launcher "${RUNS_DIR}"

# Live per-sector centroid stream + radiation extraction.
export GRTECLYN_SECTOR_BARYCENTERS=1
export GRTECLYN_PSI4=1
# Scrutiny stream (sector_dynamics.dat) -- same opt-in as the pair launcher.
export GRTECLYN_SECTOR_DYNAMICS="${BONDI_SCRUTINY:-0}"
export GRTECLYN_SECTOR_DYNAMICS_LEVEL="${BONDI_SCRUTINY_LEVEL:-0}"
# Higher-multipole Psi4 (l>=3) into its own stream psi4_mode_higher_l.dat.
# Item C: the odd-l/even-l interference carries the momentum-flux beaming along
# the runaway (x) axis.  Only meaningful once the extraction shells are in the
# wave zone (BONDI_RADII), so it is opt-in and off by default.
export GRTECLYN_PSI4_HIGHER_L="${BONDI_PSI4_HIGHER_L:-0}"
export GRTECLYN_PSI4_ELLS="${BONDI_PSI4_ELLS:-3 4}"

# Frame colour scales.  The renderer's built-in limits are calibrated for the
# search campaign's brighter matter rung and render this campaign's stars as a
# near-black smudge; frame_contrast.sh maps them onto the star this launcher is
# actually evolving.  Disable with BONDI_FRAME_CONTRAST=0.
source "${SCRIPT_DIR}/frame_contrast.sh"
bondi_frame_contrast_env

# The elliptic solve parallelises and this node's MPI is healthy (verified
# 2026-08-20; the old "mpirun segfaults" note came from the previous node).
# The EVOLUTION is single-GPU regardless -- a RadialRecipe-specific AMR crash
# still blocks multi-GPU there.
GRTRESNA_RANKS="${BONDI_GRTRESNA_RANKS:-8}"

# Numerical sponge zone -- same knob set and same caveats as the pair
# launcher (band default 24/32 is sized for L=64).
SPONGE="${BONDI_SPONGE:-0}"
SPONGE_INNER="${BONDI_SPONGE_INNER:-24}"
SPONGE_OUTER="${BONDI_SPONGE_OUTER:-32}"
SPONGE_STRENGTH="${BONDI_SPONGE_STRENGTH:-4.0}"
SPONGE_RAMP="${BONDI_SPONGE_RAMP:-4}"
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

# Print the resolved parameter set and exit without solving or touching a GPU.
dryrun_args=()
if [[ "${BONDI_DRYRUN:-0}" != "0" ]]; then
  dryrun_args=(--dry-run)
fi

out_name="bondi_sg_single_p"
if [[ "${EXOTIC}" == "1" ]]; then out_name="bondi_sg_single_m"; fi
if [[ "${OMEGA}" != "0.55" ]]; then
  out_name="${out_name}_w$(printf '%s' "${OMEGA}" | tr -d '.')"
fi
if [[ -d "${RUNS_DIR}/${out_name}" ]]; then
  echo "[bondi] ${out_name} already exists -- delete it or set BONDI_RUNS_DIR"
  exit 1
fi

echo "[bondi] === ${out_name} (ultraweak rung, gravity-dressed star) ==="
PYTHONPATH="${WRAPPER_DIR}/src" "${WRAPPER_DIR}/.venv/bin/python" \
  "${WRAPPER_DIR}/scripts/campaigns/hq/replay_eval.py" \
  "${SOURCE_EVAL}" \
  --name "${out_name}" \
  --runs-dir "${RUNS_DIR}" \
  --gpu "${GPU}" \
  --evolution-mpi-ranks "${EVO_RANKS}" \
  --n-full "${NFULL}" --l-full "${LFULL}" \
  --max-level "${MAXLEVEL}" --regrid-threshold 0.02 \
  --stop-time "${STOP_TIME}" --plot-interval 80 \
  --ftl-L 8.0 \
  --grtresna-ranks "${GRTRESNA_RANKS}" \
  --grtresna-iterations 50 \
  --grtresna-nl-exit-tolerance "${NL_TOL}" \
  --grtresna-nl-stall-tolerance "${NL_STALL_TOL}" \
  --grtresna-max-level "${GRTRESNA_MAXLEVEL}" \
  --grtresna-domain-l "${GRTRESNA_DOMAIN_L}" \
  --grtresna-n "${GRTRESNA_N}" \
  "${MAXIMAL_SLICING_FLAG[@]}" \
  "${GRIDINIT_FLAG[@]}" \
  --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
  --consumer-radii ${RADII} \
  --consumer-keep-last 2 \
  --objective-mode weighted \
  --extra-override dt_multiplier="${DT_MULT}" \
  --extra-override grtresna_scalar_lambda="${SCALAR_LAMBDA}" \
  --extra-override grtresna_scalar_mu="${SCALAR_MU}" \
  --extra-override grtresna_bs_omega="${OMEGA}" \
  --extra-override grtresna_bs_selfgrav=1 \
  --extra-override trajectory_well_width=1.2 \
  --extra-override rl_pump_stop_time=0 \
  --extra-override grtresna_boost_lumps=0 \
  --extra-override trajectory_A_breath=0 \
  --extra-override trajectory_z_amp=0 \
  --extra-override trajectory_omega_z=0 \
  --extra-override trajectory_num_lumps=1 \
  --extra-override trajectory_lump0_R0="${R0}" \
  --extra-override trajectory_lump0_phase0=0 \
  --extra-override trajectory_lump0_tilt_theta=0 \
  --extra-override trajectory_lump0_tilt_phi=0 \
  --extra-override trajectory_lump0_v_rad=0 \
  --extra-override trajectory_lump0_omega_rot=0 \
  --extra-override trajectory_lump0_well_depth=0.15 \
  --extra-override trajectory_lump0_exotic="${EXOTIC}" \
  ${sponge_args[@]+"${sponge_args[@]}"} \
  ${dryrun_args[@]+"${dryrun_args[@]}"}

echo "[bondi] dressed-star calibration complete: ${RUNS_DIR}/${out_name}"
