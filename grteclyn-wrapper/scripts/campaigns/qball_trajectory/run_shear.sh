#!/usr/bin/env bash
# Extreme spacetime-shear MAP-Elites.
#
# Sibling campaign to qball_traj_lentz_v1.  Instead of scoring FTL shortcuts,
# this run directly rewards strong spacetime curvature / frame-dragging shear
# while penalizing collapse to a horizon.  Goal: find the most extreme,
# non-collapsing, positive-energy configurations that GR permits.
#
# Uses the new objective_mode=spacetime_shear and descriptor_mode=spacetime_shear.
#
# Usage:
#   cd grteclyn-wrapper
#   bash scripts/campaigns/qball_trajectory/run_shear.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
QD_RUN="${SCRIPT_DIR}/../qd/run.sh"

# --- QD run identity / scale (override via env) ---
export QD_NAME="${QD_NAME:-qball_traj_shear_v1}"
export QD_TARGET_EVALS="${QD_TARGET_EVALS:-200}"
export GPU_IDS="${GPU_IDS:-0 1 2 3}"
export BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
export QD_ITERATIONS="${QD_ITERATIONS:-20}"

# --- Matter model: canonical boson star only ---
export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=trajectory
export GRTRESNA_MATTER_COUPLING=canonical
export GRTRESNA_FULL_Z=1

# --- Objective: spacetime shear ---
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-spacetime_shear}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-spacetime_shear}"

# Mild exotic penalty.  The geometry of normal matter at extreme speeds may
# effectively violate the NEC; we still prefer cleaner configs but do not
# crush the search.
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-1.0}"

# Dispersion gate: keep so a dissolving cloud cannot bank shear.
export SCORE_FTL_DISPERSION_GATE="${SCORE_FTL_DISPERSION_GATE:-1.0}"

# --- Compact Q-ball physics ---
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_bs_omega=0.8 \
trajectory_lump0_well_depth=0.15 \
trajectory_lump1_well_depth=0.15 \
trajectory_lump2_well_depth=0.15 \
trajectory_lump3_well_depth=0.15 \
trajectory_lump4_well_depth=0.15 \
trajectory_lump0_exotic=0.0 \
trajectory_lump1_exotic=0.0 \
trajectory_lump2_exotic=0.0 \
trajectory_lump3_exotic=0.0 \
trajectory_lump4_exotic=0.0 \
trajectory_well_width=1.667}"

# --- Evolution ---
export STOP_TIME="${STOP_TIME:-16.0}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

# Uncap trajectory speeds to the extreme frame-dragging regime.
export EXTRA_SETS="${EXTRA_SETS:-\
grtresna_scalar_mu=85333 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=1 \
trajectory_r_min=0.1 \
trajectory_v_max=0.99 \
grtresna_boost_v_max=0.99 \
stop_time=${STOP_TIME}}"

# --- Multi-ray emission sweep (kept for diagnostics, not scored) ---
export GRTECLYN_GEO_EMIT_INTERVAL="${GRTECLYN_GEO_EMIT_INTERVAL:-2}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GRTECLYN_GEO_MAX_EMISSIONS:-7}"

# --- Probe / frames ---
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-search}"
export GRTECLYN_GEO_DIRECTIONS="x y z"

# --- Postload gate ---
# Relaxed so extreme 0.99c Lorentz-boosted initial data can pass.
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-0.1}"
export POSTLOAD_MAX_MOM_L2="${POSTLOAD_MAX_MOM_L2:-0.1}"

# --- Pipeline ---
export USE_PIPELINE="${USE_PIPELINE:-1}"
_grtresna_cap=6
if (( BATCH_SIZE < _grtresna_cap )); then _grtresna_cap="${BATCH_SIZE}"; fi
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-${_grtresna_cap}}"

echo "== Extreme Spacetime-Shear QD: ${QD_NAME} =="
echo "   GPUs: ${GPU_IDS} (batch=${BATCH_SIZE})  target_evals=${QD_TARGET_EVALS}"
echo "   Matter: canonical boson star, all exotic flags pinned to 0"
echo "   Speed: trajectory_v_max=0.99, boost_v_max=0.99"
echo "   Objective: ${OBJECTIVE_MODE} / Descriptor: ${DESCRIPTOR_MODE}"
echo "   Score: SCORE_EXOTIC_PENALTY_WEIGHT=${SCORE_EXOTIC_PENALTY_WEIGHT} + SCORE_FTL_DISPERSION_GATE=${SCORE_FTL_DISPERSION_GATE}"
echo "   Postload gate: POSTLOAD_MAX_HAM_L2=${POSTLOAD_MAX_HAM_L2} POSTLOAD_MAX_MOM_L2=${POSTLOAD_MAX_MOM_L2}"
echo "   Pipeline: max_grtresna=${MAX_CONCURRENT_GRTRESNA} scoring_workers=${SCORING_WORKERS:-<#GPUs>}"

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
