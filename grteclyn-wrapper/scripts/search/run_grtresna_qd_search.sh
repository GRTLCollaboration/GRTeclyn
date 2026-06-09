#!/usr/bin/env bash
# GRTresna-in-the-loop MAP-Elites search.
#
# This is the quality-diversity counterpart to run_grtresna_search.sh: it keeps
# a diverse archive of shell candidates instead of collapsing immediately onto
# the highest scalar CMA-ES score.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  elif command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  else
    PYTHON_BIN="python"
  fi
fi

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/grtresna_qd}"
LUMPS="${LUMPS:-5}"
SHELL_PROFILE="${SHELL_PROFILE:-compact}"
RANKS="${RANKS:-8}"
ITERATIONS="${ITERATIONS:-30}"
GRTRESNA_MAX_LEVEL="${GRTRESNA_MAX_LEVEL:-3}"
GRTRESNA_REFINE_THRESHOLD="${GRTRESNA_REFINE_THRESHOLD:-0.5}"
GRTRESNA_REGRID_RADIUS="${GRTRESNA_REGRID_RADIUS:-0}"
GRTRESNA_JACOBIAN_CAP="${GRTRESNA_JACOBIAN_CAP:-25.0}"
GRTRESNA_TIMEOUT="${GRTRESNA_TIMEOUT:-900}"
GRTRESNA_EVOLUTION_L_FULL="${GRTRESNA_EVOLUTION_L_FULL:-64.0}"
GRTRESNA_EVOLUTION_N_FULL="${GRTRESNA_EVOLUTION_N_FULL:-64}"
GRTRESNA_DOMAIN_L="${GRTRESNA_DOMAIN_L:-128.0}"
GRTRESNA_DOMAIN_NX="${GRTRESNA_DOMAIN_NX:-64}"
GRTRESNA_DOMAIN_NY="${GRTRESNA_DOMAIN_NY:-64}"
GRTRESNA_DOMAIN_NZ="${GRTRESNA_DOMAIN_NZ:-}"
GRTRESNA_GRIDINIT_NX="${GRTRESNA_GRIDINIT_NX:-64}"
GRTRESNA_GRIDINIT_NY="${GRTRESNA_GRIDINIT_NY:-64}"
GRTRESNA_GRIDINIT_NZ="${GRTRESNA_GRIDINIT_NZ:-64}"
QD_ITERATIONS="${QD_ITERATIONS:-10}"
BINS="${BINS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-10}"

SOLVED_FTL_F_OP_FLOOR="${SOLVED_FTL_F_OP_FLOOR:-1.0e-4}"
SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR="${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR:-0.95}"
SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR="${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR:-1.01}"
SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR="${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR:-0.02}"
SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED="${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED:-8.0}"
SOLVED_FTL_MAX_PHYSICAL_F_OP="${SOLVED_FTL_MAX_PHYSICAL_F_OP:-0.85}"
SOLVED_FTL_REJECTION_SPEED_TARGET="${SOLVED_FTL_REJECTION_SPEED_TARGET:-1.01}"
GRTRESNA_MAX_HAM_PCT="${GRTRESNA_MAX_HAM_PCT:-5.0}"
GRTRESNA_MAX_MOM_PCT="${GRTRESNA_MAX_MOM_PCT:-5.0}"

POSTLOAD_GATE="${POSTLOAD_GATE:-1}"
POSTLOAD_MAX_HAM_L2="${POSTLOAD_MAX_HAM_L2:-1e-2}"
POSTLOAD_MAX_MOM_L2="${POSTLOAD_MAX_MOM_L2:-1e-2}"
export POSTLOAD_GATE

CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
FTL_L="${FTL_L:-8.0}"

export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi Pi chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

mkdir -p "${RUNS_DIR}"

DRY_RUN_ARGS=()
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  DRY_RUN_ARGS=(--dry-run)
fi

DOMAIN_ARGS=(
  --grtresna-evolution-l-full "${GRTRESNA_EVOLUTION_L_FULL}"
  --grtresna-evolution-n-full "${GRTRESNA_EVOLUTION_N_FULL}"
  --grtresna-domain-l "${GRTRESNA_DOMAIN_L}"
  --grtresna-domain-nx "${GRTRESNA_DOMAIN_NX}"
  --grtresna-domain-ny "${GRTRESNA_DOMAIN_NY}"
  --grtresna-gridinit-nx "${GRTRESNA_GRIDINIT_NX}"
  --grtresna-gridinit-ny "${GRTRESNA_GRIDINIT_NY}"
  --grtresna-gridinit-nz "${GRTRESNA_GRIDINIT_NZ}"
)
if [[ -n "${GRTRESNA_DOMAIN_NZ}" ]]; then
  DOMAIN_ARGS+=(--grtresna-domain-nz "${GRTRESNA_DOMAIN_NZ}")
fi

# shellcheck disable=SC2086
exec ${PYTHON_BIN} -m grteclyn_wrapper \
  --runs-dir "${RUNS_DIR}" \
  --example RadialRecipe \
  --set stop_time="${STOP_TIME}" \
  --set plot_interval="${PLOT_INTERVAL}" \
  "${DRY_RUN_ARGS[@]}" \
  --consume-plotfiles \
  --consumer-delete \
  --consumer-radii ${CONSUMER_RADII} \
  --ftl-L "${FTL_L}" \
  qd \
  --descriptor-mode channel \
  --objective-mode ftl_first \
  --iterations "${QD_ITERATIONS}" \
  --batch-size "${BATCH_SIZE}" \
  --bins "${BINS}" \
  --seed "${SEED}" \
  --gpu-ids ${GPU_IDS} \
  --grtresna \
  --grtresna-ansatz shell \
  --grtresna-shell-profile "${SHELL_PROFILE}" \
  --grtresna-lumps "${LUMPS}" \
  --grtresna-full-z \
  "${DOMAIN_ARGS[@]}" \
  --grtresna-ranks "${RANKS}" \
  --grtresna-iterations "${ITERATIONS}" \
  --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
  --grtresna-max-level "${GRTRESNA_MAX_LEVEL}" \
  --grtresna-refine-threshold "${GRTRESNA_REFINE_THRESHOLD}" \
  --grtresna-regrid-radius "${GRTRESNA_REGRID_RADIUS}" \
  --grtresna-jacobian-cap "${GRTRESNA_JACOBIAN_CAP}" \
  --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
  --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
  --solved-ftl-f-op-floor "${SOLVED_FTL_F_OP_FLOOR}" \
  --solved-ftl-near-luminal-speed-floor "${SOLVED_FTL_NEAR_LUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-speed-floor "${SOLVED_FTL_SUPERLUMINAL_SPEED_FLOOR}" \
  --solved-ftl-superluminal-fraction-floor "${SOLVED_FTL_SUPERLUMINAL_FRACTION_FLOOR}" \
  --solved-ftl-max-physical-coord-speed "${SOLVED_FTL_MAX_PHYSICAL_COORD_SPEED}" \
  --solved-ftl-max-physical-f-op "${SOLVED_FTL_MAX_PHYSICAL_F_OP}" \
  --solved-ftl-rejection-speed-target "${SOLVED_FTL_REJECTION_SPEED_TARGET}" \
  --grtresna-postload-gate \
  --postload-max-ham-l2 "${POSTLOAD_MAX_HAM_L2}" \
  --postload-max-mom-l2 "${POSTLOAD_MAX_MOM_L2}"
