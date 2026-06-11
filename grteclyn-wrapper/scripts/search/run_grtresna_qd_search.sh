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
# Raised from 30: oscillating near-misses (Ham~7%/Mom~15%) need more iterations
# to settle below the 5% gate instead of being cut off early.
ITERATIONS="${ITERATIONS:-50}"
# GRTresna early-exit: stop once Ham/Mom are good (%%) or improvement stalls.
# Stall tolerance tightened from 0.02 so near-converging solves are not
# abandoned prematurely on a small residual plateau.
GRTRESNA_NL_EXIT_TOLERANCE="${GRTRESNA_NL_EXIT_TOLERANCE:-1.0}"
GRTRESNA_NL_STALL_TOLERANCE="${GRTRESNA_NL_STALL_TOLERANCE:-0.005}"
# Parallel Chombo→gridinit conversion (0 = auto, min(32, cpu_count)).
GRTRESNA_GRIDINIT_WORKERS="${GRTRESNA_GRIDINIT_WORKERS:-0}"
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
QD_TARGET_EVALS="${QD_TARGET_EVALS:-}"
QD_RESUME="${QD_RESUME:-0}"
QD_NAME="${QD_NAME:-}"
# MAP-Elites behaviour grid: channel (path x mechanism, needs shift>0),
# speed_horizon (cone-tilt x horizon-free, the fast-but-not-trapped niche),
# speed_super (recalibrated cone-tilt x superluminal_fraction, stays
# discriminating in the nontrivial regime), legacy.
DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-speed_super}"
# Optional: warm-start the initial population from prior eval dirs (survivors).
SEED_EVAL_DIRS="${SEED_EVAL_DIRS:-}"
# Keep disk bounded: retain only the top-N scored eval_* directories plus the
# current in-flight batch; trajectory/archive metadata stay intact.
QD_KEEP_TOP_EVAL_DIRS="${QD_KEEP_TOP_EVAL_DIRS:-10}"
BINS="${BINS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-7}"
# Long enough for structural dissipation to show: at t=2 a dissipating
# configuration is indistinguishable from a persistent one (both retain ~90% of
# peak energy density), so the persistence-gated survival metric cannot bite.
# By t=8 a genuine survivor has plateaued (~0.7) while dissipators keep falling
# (<0.55), which is the smallest window that cleanly separates them (~4x the
# t=2 GPU cost).
STOP_TIME="${STOP_TIME:-8.0}"
# Scaled with STOP_TIME (4x) so the number of consumed plotfiles -- and hence
# the FTL-probe/plotting cost per eval -- stays roughly constant.
PLOT_INTERVAL="${PLOT_INTERVAL:-40}"

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

# Keep the GRTresna Chombo HDF5 + workdir per eval (conversion validation/debug).
GRTRESNA_KEEP_SOURCE="${GRTRESNA_KEEP_SOURCE:-0}"
KEEP_SOURCE_ARGS=()
if [[ "${GRTRESNA_KEEP_SOURCE}" == "1" ]]; then
  KEEP_SOURCE_ARGS+=(--grtresna-keep-source)
fi

CONSUMER_RADII="${CONSUMER_RADII:-4 8}"
# Retain the last few plotfiles so evolved/geodesic FTL and effective energy
# conditions can be scored (>=3 required); the rest are still consumed+deleted.
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-3}"
FTL_L="${FTL_L:-8.0}"

export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-lump_activity scalar_activity phi_lump_sum Pi_lump_sum chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

mkdir -p "${RUNS_DIR}"

NAME_ARGS=()
if [[ -n "${QD_NAME}" ]]; then
  NAME_ARGS+=(--name "${QD_NAME}")
fi
RESUME_ARGS=()
if [[ "${QD_RESUME}" == "1" ]]; then
  RESUME_ARGS+=(--resume)
fi
TARGET_EVALS_ARGS=()
if [[ -n "${QD_TARGET_EVALS}" ]]; then
  TARGET_EVALS_ARGS+=(--target-evals "${QD_TARGET_EVALS}")
fi
SEED_ARGS=()
if [[ -n "${SEED_EVAL_DIRS}" ]]; then
  # shellcheck disable=SC2206
  SEED_ARGS+=(--seed-eval-dirs ${SEED_EVAL_DIRS})
fi

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
  "${NAME_ARGS[@]}" \
  --example RadialRecipe \
  --set stop_time="${STOP_TIME}" \
  --set plot_interval="${PLOT_INTERVAL}" \
  "${DRY_RUN_ARGS[@]}" \
  --consume-plotfiles \
  --consumer-delete \
  --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
  --consumer-radii ${CONSUMER_RADII} \
  --ftl-L "${FTL_L}" \
  qd \
  --descriptor-mode "${DESCRIPTOR_MODE}" \
  --objective-mode ftl_first \
  --iterations "${QD_ITERATIONS}" \
  --keep-top-eval-dirs "${QD_KEEP_TOP_EVAL_DIRS}" \
  "${TARGET_EVALS_ARGS[@]}" \
  "${RESUME_ARGS[@]}" \
  "${SEED_ARGS[@]}" \
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
  --grtresna-nl-exit-tolerance "${GRTRESNA_NL_EXIT_TOLERANCE}" \
  --grtresna-nl-stall-tolerance "${GRTRESNA_NL_STALL_TOLERANCE}" \
  --grtresna-gridinit-workers "${GRTRESNA_GRIDINIT_WORKERS}" \
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
  --postload-max-mom-l2 "${POSTLOAD_MAX_MOM_L2}" \
  "${KEEP_SOURCE_ARGS[@]}"
