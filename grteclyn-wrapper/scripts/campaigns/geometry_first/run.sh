#!/usr/bin/env bash
# Geometry-first MAP-Elites campaign — search for matter configurations
# that reproduce a target geometry motif via GRTresna constraint solve.
#
# Uses the same GRTresna settings as the QD campaign (N=128, L=64, RANKS=8,
# shell ansatz) but scores by geometry mismatch instead of FTL metrics.
# No GRTeclyn time evolution — scoring is purely from the solved gridinit.
#
# Usage:
#   MOTIF_JSON=/path/to/motif.json bash scripts/campaigns/geometry_first/run.sh
#
# With custom GPU allocation:
#   MOTIF_JSON=/path/to/motif.json GPU_IDS="0 1 2 3" \
#     QD_TARGET_EVALS=80 bash scripts/campaigns/geometry_first/run.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/bootstrap.sh
source "${SCRIPT_DIR}/../lib/bootstrap.sh"
_campaign_bootstrap "${SCRIPT_DIR}"
_campaign_resolve_python

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/geometry_first}"
MOTIF_JSON="${MOTIF_JSON:?MOTIF_JSON is required — path to motif.json}"
QD_ITERATIONS="${QD_ITERATIONS:-10}"
QD_TARGET_EVALS="${QD_TARGET_EVALS:-}"
QD_RESUME="${QD_RESUME:-0}"
QD_NAME="${QD_NAME:-geometry_first_$(date -u +%Y%m%dT%H%M%SZ)}"
DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-geometry_first}"
BINS="${BINS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"

# Use the same GRTresna settings as the QD campaign
export GRTRESNA_ANSATZ="${GRTRESNA_ANSATZ:-shell}"
export GRTRESNA_MATTER_SECTOR="${GRTRESNA_MATTER_SECTOR:-scalar}"
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export SHELL_PROFILE="${SHELL_PROFILE:-compact}"
export LUMPS="${LUMPS:-5}"
export RANKS="${RANKS:-8}"
export ITERATIONS="${ITERATIONS:-50}"
export GRTRESNA_NL_EXIT_TOLERANCE="${GRTRESNA_NL_EXIT_TOLERANCE:-1.0}"
export GRTRESNA_NL_STALL_TOLERANCE="${GRTRESNA_NL_STALL_TOLERANCE:-0.005}"
export GRTRESNA_MAX_LEVEL="${GRTRESNA_MAX_LEVEL:-3}"
export GRTRESNA_REFINE_THRESHOLD="${GRTRESNA_REFINE_THRESHOLD:-0.5}"
export GRTRESNA_JACOBIAN_CAP="${GRTRESNA_JACOBIAN_CAP:-25.0}"
export GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-900}"
export GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-5.0}"
export GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-5.0}"

# Grid: N=128, L=64 (same as QD campaign — dx=0.5)
export GRTRESNA_EVOLUTION_L_FULL="${GRTRESNA_EVOLUTION_L_FULL:-64.0}"
export GRTRESNA_EVOLUTION_N_FULL="${GRTRESNA_EVOLUTION_N_FULL:-128}"
export GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-128.0}"
export GRTRESNA_DOMAIN_NX="${GRTRESNA_DOMAIN_NX:-64}"
export GRTRESNA_DOMAIN_NY="${GRTRESNA_DOMAIN_NY:-64}"
export GRTRESNA_GRIDINIT_NX="${GRTRESNA_GRIDINIT_NX:-128}"
export GRTRESNA_GRIDINIT_NY="${GRTRESNA_GRIDINIT_NY:-128}"
export GRTRESNA_GRIDINIT_NZ="${GRTRESNA_GRIDINIT_NZ:-128}"
export GRTRESNA_FULL_Z="${GRTRESNA_FULL_Z:-1}"

# No time evolution — geometry_first scores from gridinit only
export STOP_TIME="${STOP_TIME:-0.0}"
export GRTECLYN_FRAMES="${GRTECLYN_FRAMES:-0}"
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-0}"

# FTL probe L for solved-geometry analysis
export FTL_L="${FTL_L:-8.0}"

# Pipeline settings
export GPU_SLOTS_PER_DEVICE="${GPU_SLOTS_PER_DEVICE:-1}"
export MAX_CONCURRENT_GRTESNA="${MAX_CONCURRENT_GRTESNA:-4}"
export USE_PIPELINE="${USE_PIPELINE:-1}"
export CLUSTER_CPU_FRACTION="${CLUSTER_CPU_FRACTION:-0.30}"
export PIPELINE_CPU_SHARE="${PIPELINE_CPU_SHARE:-1.0}"

# Objective and descriptor modes for geometry-first
export OBJECTIVE_MODE="geometry_first"
export DESCRIPTOR_MODE="geometry_first"

mkdir -p "${RUNS_DIR}"

NAME_ARGS=(--name "${QD_NAME}")
RESUME_ARGS=()
[[ "${QD_RESUME}" == "1" ]] && RESUME_ARGS+=(--resume)
TARGET_EVALS_ARGS=()
[[ -n "${QD_TARGET_EVALS}" ]] && TARGET_EVALS_ARGS+=(--target-evals "${QD_TARGET_EVALS}")

DRY_RUN_ARGS=()
[[ "${DRY_RUN:-0}" == "1" ]] && DRY_RUN_ARGS+=(--dry-run)

# Shared search args (grid, GRTresna, pipeline)
ftl_search_common_domain_args
ftl_search_common_grtresna_args
ftl_search_common_global_args
ftl_search_common_pipeline_args

# Skip preflight tests for geometry_first (no FTL evolution)
export SKIP_QD_PREFLIGHT_TESTS=1

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  --runs-dir "${RUNS_DIR}" \
  "${NAME_ARGS[@]}" \
  "${FTL_GLOBAL_ARGS[@]}" \
  "${DRY_RUN_ARGS[@]}" \
  qd \
  --motif-json "${MOTIF_JSON}" \
  --descriptor-mode "${DESCRIPTOR_MODE}" \
  --objective-mode "${OBJECTIVE_MODE}" \
  --iterations "${QD_ITERATIONS}" \
  --keep-top-eval-dirs "${QD_KEEP_TOP_EVAL_DIRS:-3}" \
  --no-ftl-retention \
  "${TARGET_EVALS_ARGS[@]}" \
  "${RESUME_ARGS[@]}" \
  --batch-size "${BATCH_SIZE}" \
  --bins "${BINS}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  "${FTL_PIPELINE_ARGS[@]}" \
  "${FTL_GRTRESNA_ARGS[@]}"
