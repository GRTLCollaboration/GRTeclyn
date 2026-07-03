#!/usr/bin/env bash
# GW beam MAP-Elites — directional gravitational-wave optimization with normal matter.
#
# Re-uses the compact canonical Q-ball trajectory model from the HQ replay
# qball_traj_spiral_v2_hq_eval000118 (grtresna_complex_scalar, m=1, lambda=640,
# mu=85333, omega=0.8, ODE profile, equilibrium amplitude) but pins every lump
# to canonical matter (exotic=0) and switches the objective from FTL to
# directional Psi4 gravitational-wave emission.
#
# Physics choices:
#   - L=64, N=128 box (dx=0.5) to keep stage-0 QD cost reasonable.
#   - Single extraction radius R=12, well inside the wave-zone and far from the
#     wall at distance 32 from the center.
#   - stop_time=24: light-travel time from orbit (R~5) to R=12 is ~7, leaving a
#     ~17-unit observing window; reflections from the wall need ~54 units, so
#     t=24 is clean.
#
# Objective: gw_beam (total GW power + Z-axis beaming ratio).
# Descriptor: gw_beam (log total power vs beam ratio).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
QD_RUN="${SCRIPT_DIR}/../qd/run.sh"

# --- QD run identity / scale (override via env) ---
export QD_NAME="${QD_NAME:-gw_beam_v1}"
export QD_TARGET_EVALS="${QD_TARGET_EVALS:-200}"
export GPU_IDS="${GPU_IDS:-0 1 2 3}"
export BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
export QD_ITERATIONS="${QD_ITERATIONS:-20}"

# --- Matter model: compact canonical Q-balls on trajectory orbits ---
export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=trajectory
export GRTRESNA_MATTER_COUPLING=canonical
export GRTRESNA_FULL_Z=1

# --- Objective: directional gravitational-wave emission ---
export OBJECTIVE_MODE=gw_beam
export DESCRIPTOR_MODE=gw_beam

# --- No exotic matter: every lump is canonical (+rho) ---
export SCORE_EXOTIC_PENALTY_WEIGHT=0.0

# --- Compact Q-ball physics (same as HQ run qball_traj_spiral_v2_hq_eval000118) ---
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_bs_omega=0.8 \
trajectory_lump0_exotic=0 \
trajectory_lump1_exotic=0 \
trajectory_lump2_exotic=0 \
trajectory_lump3_exotic=0 \
trajectory_lump4_exotic=0 \
trajectory_lump0_well_depth=0.15 \
trajectory_lump1_well_depth=0.15 \
trajectory_lump2_well_depth=0.15 \
trajectory_lump3_well_depth=0.15 \
trajectory_lump4_well_depth=0.15 \
trajectory_well_width=1.667}"

# --- Evolution: long enough for GWs to reach R=12 and pass through the sphere ---
# Light travel from orbit (R~5) to R=12 is ~7; t=24 leaves ~17 units of observing
# window before any wall reflection returns.
export STOP_TIME=24.0
# Coarse plot cadence: only ~5 plotfiles per eval, so the consumer can keep up with
# the GPU and the QD batch turns faster.
export PLOT_INTERVAL=480

# --- GRTresna speed: reduce solve depth/iterations so the CPU does not become
# the bottleneck and evals do not hang for long solves.
export ITERATIONS=25
export GRTRESNA_MAX_LEVEL=1
export GRTRESNA_TIMEOUT=300
export GRTRESNA_NL_EXIT_TOLERANCE=2.0

# Non-search-space overrides (base overrides).  trajectory_r_min matches the
# C++ TRAJECTORY_R_MIN; amr.derive_plot_vars is required for Weyl4_Re/Im output.
# The Q-ball ODE profile, equilibrium amplitude cap, and scalar_mu are not
# optimizer dimensions, so they go here instead of PIN_DIMS.
export EXTRA_SETS="${EXTRA_SETS:-\
grtresna_scalar_mu=85333 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=1 \
trajectory_r_min=0.1 \
stop_time=${STOP_TIME} \
plot_interval=${PLOT_INTERVAL} \
amr.derive_plot_vars=Weyl4}"

# --- GW extraction / no FTL ---
export GRTECLYN_PSI4=1
export GRTECLYN_PSI4_N_POINTS=64
export CONSUMER_RADII="12"
# Drain-minimal: only extract Psi4 and delete plotfiles; skip FTL/geodesic,
# confinement, shell-field, areal-radius, and incremental scoring. The QD
# pipeline still computes the final gw_beam score from psi4_directional.dat.
export GRTECLYN_CONSUMER_DRAIN_MINIMAL=1
# Do not retain plotfiles: the GW signal is summarized in psi4_directional.dat.
export CONSUMER_KEEP_LAST=0
# Run the plotfile consumer with parallel workers so it can keep up with the GPU.
export CONSUMER_JOBS=4
export GRTECLYN_EVOLVING_GEODESIC=0
export GRTECLYN_GEO_DIRECTIONS=""
export GRTECLYN_GEO_EMIT_INTERVAL=""
export GRTECLYN_GEO_MAX_EMISSIONS=""
export GRTECLYN_FRAMES=0

# --- Pipeline ---
export USE_PIPELINE="${USE_PIPELINE:-1}"
_grtresna_cap=6
if (( BATCH_SIZE < _grtresna_cap )); then _grtresna_cap="${BATCH_SIZE}"; fi
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-${_grtresna_cap}}"
# The post-load constraint gate is tuned for FTL warps (rejects weak precursors).
# For GW beam we only need the gridinit to load and the scalar matter to be
# evolved; the GW scorer will separate strong emitters from weak/crashed ones.
export POSTLOAD_GATE=0
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"

echo "== GW beam QD: ${QD_NAME} =="
echo "   GPUs: ${GPU_IDS} (batch=${BATCH_SIZE})  target_evals=${QD_TARGET_EVALS}"
echo "   Box: L=64 N=128  R_ext=12  t_stop=${STOP_TIME}"
echo "   Objective: ${OBJECTIVE_MODE}  descriptor: ${DESCRIPTOR_MODE}"
echo "   Pipeline: max_grtresna=${MAX_CONCURRENT_GRTRESNA}"

if [[ "${PIPELINE_MONITOR:-1}" == "1" ]]; then
  CAMPAIGNS_LIB="${SCRIPT_DIR}/../lib"
  # shellcheck source=../lib/bootstrap.sh
  source "${CAMPAIGNS_LIB}/bootstrap.sh"
  _campaign_bootstrap "${SCRIPT_DIR}"
  # shellcheck source=../lib/pipeline_monitor.sh
  source "${CAMPAIGNS_LIB}/pipeline_monitor.sh"
  ftl_pipeline_monitor_begin "${QD_NAME}" "${GPU_IDS}"
  bash "${QD_RUN}" 2>&1 | tee "${FTL_PIPELINE_LOG}"
  ftl_pipeline_monitor_end
else
  exec bash "${QD_RUN}"
fi
