#!/usr/bin/env bash
# Q-ball trajectory CMA-ES — local refine of bicomplex QD elites.
#
# Mirrors qball_trajectory/run.sh physics (bicomplex, igniter pump, score
# weights, pins, EXTRA_SETS) so Stage-1 refinement is score-comparable to QD.
# Warm-starts from the QD trajectory.jsonl (top-K gpu_ok genomes).
#
# Usage (after QD finishes):
#   cd grteclyn-wrapper
#   bash scripts/campaigns/qball_trajectory/cmaes_run.sh
#
# Overrides:
#   RUN_NAME=qball_traj_bicomplex_cmaes_v1
#   WARM_START_TRAJECTORY=../runs/grtresna_qd/qball_traj_bicomplex_v1/trajectory.jsonl
#   WARM_START_TOP_K=1 TARGET_EVALS=150 GPU_IDS="0 1 2 3 4 5 6 7"
#   DRY_RUN=1   # smoke without GPU
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CMAES_RUN="${SCRIPT_DIR}/../cmaes/run.sh"
GRTECLYN_ROOT="$(cd -- "${SCRIPT_DIR}/../../.." && pwd)"
# Repo root is parent of grteclyn-wrapper when SCRIPT_DIR is .../scripts/campaigns/qball_trajectory
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"

# --- CMA-ES run identity ---
export RUN_NAME="${RUN_NAME:-qball_traj_bicomplex_cmaes_v1}"
export RUNS_DIR="${RUNS_DIR:-${REPO_ROOT}/runs/grtresna_cmaes}"
export WARM_START_TRAJECTORY="${WARM_START_TRAJECTORY:-${REPO_ROOT}/runs/grtresna_qd/qball_traj_bicomplex_v1/trajectory.jsonl}"
export WARM_START_TOP_K="${WARM_START_TOP_K:-1}"
export WARM_START_JITTER="${WARM_START_JITTER:-0.05}"
export SIGMA0="${SIGMA0:-0.05}"
export TARGET_EVALS="${TARGET_EVALS:-150}"
export MAX_GENERATIONS="${MAX_GENERATIONS:-50}"
# Keep enough elites that HQ/matrix can still find the champion if freeze is late.
export KEEP_TOP_EVAL_DIRS="${KEEP_TOP_EVAL_DIRS:-10}"
export GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
export POPULATION="${POPULATION:-$(wc -w <<< "${GPU_IDS}")}"

# --- Matter model: must match QD (bicomplex) ---
export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=trajectory
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export GRTRESNA_FULL_Z=1
export GRTRESNA_ALLOW_SIGN_MISMATCH="${GRTRESNA_ALLOW_SIGN_MISMATCH:-0}"
export SCORE_PUMP_ENERGY_WEIGHT="${SCORE_PUMP_ENERGY_WEIGHT:-40}"
export RL_PUMP_STOP_TIME="${RL_PUMP_STOP_TIME:-4}"
export FRAMES_FIELDS="${FRAMES_FIELDS:-scalar_activity phi Pi phi_lump0 Pi_lump0 phi_lump1 Pi_lump1 phi_lump2 Pi_lump2 chi chi_minus_1 local_speed shift1 rho_req}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity phi}"

# --- Objective / score parity with QD ---
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"
export SCORE_EXOTIC_PENALTY_WEIGHT="${SCORE_EXOTIC_PENALTY_WEIGHT:-0.2}"
export SCORE_FTL_DISPERSION_GATE="${SCORE_FTL_DISPERSION_GATE:-1.0}"

# --- Pinned physics (same as QD) ---
export PIN_DIMS="${PIN_DIMS:-\
grtresna_scalar_mass=1.0 \
grtresna_scalar_lambda=640 \
grtresna_bs_omega=0.8 \
trajectory_lump0_well_depth=0.15 \
trajectory_lump1_well_depth=0.15 \
trajectory_lump2_well_depth=0.15 \
trajectory_lump3_well_depth=0.15 \
trajectory_lump4_well_depth=0.15 \
trajectory_well_width=1.667}"

export STOP_TIME="${STOP_TIME:-16.0}"
export PLOT_INTERVAL="${PLOT_INTERVAL:-320}"

export EXTRA_SETS="${EXTRA_SETS:-\
grtresna_matter_model=grtresna_bicomplex_scalar \
grtresna_scalar_mu=85333 \
grtresna_qball_ode_profile=1 \
grtresna_qball_equilibrium_amplitude=1 \
trajectory_retrograde_only=1 \
trajectory_r_min=0.1 \
stop_time=${STOP_TIME}}"

# --- Geodesic emit sweep (parity with QD) ---
export GRTECLYN_GEO_EMIT_INTERVAL="${GRTECLYN_GEO_EMIT_INTERVAL:-2}"
export GRTECLYN_GEO_MAX_EMISSIONS="${GRTECLYN_GEO_MAX_EMISSIONS:-7}"
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
export GRTECLYN_EVOLVING_GEODESIC=1
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-search}"
export GRTECLYN_GEO_DIRECTIONS="x y z"

# --- Pipeline ---
export USE_PIPELINE="${USE_PIPELINE:-1}"
_grtresna_cap=6
_pop="${POPULATION}"
if (( _pop < _grtresna_cap )); then _grtresna_cap="${_pop}"; fi
export MAX_CONCURRENT_GRTRESNA="${MAX_CONCURRENT_GRTRESNA:-${_grtresna_cap}}"
export POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-3e-2}"

if [[ ! -f "${WARM_START_TRAJECTORY}" ]]; then
  echo "ERROR: warm-start trajectory missing: ${WARM_START_TRAJECTORY}" >&2
  echo "Finish QD (qball_traj_bicomplex_v1) first, or set WARM_START_TRAJECTORY." >&2
  exit 1
fi

echo "== Q-ball trajectory CMA-ES: ${RUN_NAME} =="
echo "   Warm-start : ${WARM_START_TRAJECTORY} (top_k=${WARM_START_TOP_K})"
echo "   GPUs       : ${GPU_IDS} (pop=${POPULATION})  target_evals=${TARGET_EVALS}"
echo "   Matter     : grtresna_bicomplex_scalar (EXTRA_SETS forwarded)"
echo "   Pump/score : RL_PUMP_STOP_TIME=${RL_PUMP_STOP_TIME}  exotic_w=${SCORE_EXOTIC_PENALTY_WEIGHT}  pump_w=${SCORE_PUMP_ENERGY_WEIGHT}"
echo "   Output     : ${RUNS_DIR}/${RUN_NAME}"

exec bash "${CMAES_RUN}"
