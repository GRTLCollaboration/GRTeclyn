#!/usr/bin/env bash
# Promote one geometry-first scout through GRTresna.
#
# Usage:
#   SOURCE=../runs/ftl_search_nonspherical_01/eval_000123 \
#     bash scripts/search/run_geometry_projection_eval.sh
#
# Or:
#   bash scripts/search/run_geometry_projection_eval.sh <source-eval> <out-dir>
set -euo pipefail

CALLER_CWD="$(pwd)"
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

SOURCE="${1:-${SOURCE:-}}"
if [[ -z "${SOURCE}" ]]; then
  echo "Set SOURCE=/path/to/eval_XXXXXX or pass it as the first argument." >&2
  exit 2
fi

if [[ "${SOURCE}" = /* ]]; then
  SOURCE_ABS="${SOURCE}"
else
  SOURCE_ABS="$(cd -- "${CALLER_CWD}/$(dirname -- "${SOURCE}")" && pwd)/$(basename -- "${SOURCE}")"
fi
OUT="${2:-${OUT:-${GRTECLYN_ROOT}/runs/geometry_projection/$(basename -- "${SOURCE_ABS}")}}"
if [[ "${OUT}" != /* ]]; then
  OUT="${CALLER_CWD}/${OUT}"
fi

MODE="${MODE:-solve-and-evolve}"
RANKS="${RANKS:-8}"
MAX_LUMPS="${MAX_LUMPS:-3}"
FTL_L="${FTL_L:-8.0}"
GRTRESNA_L="${GRTRESNA_L:-128.0}"
GRIDINIT_N="${GRIDINIT_N:-64}"
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-10}"
CUDA_DEVICE="${CUDA_DEVICE:-0}"

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

POSTLOAD_ARGS=()
if [[ "${SKIP_POSTLOAD_GATE:-0}" == "1" ]]; then
  POSTLOAD_ARGS+=(--skip-postload-gate)
fi

DRY_RUN_ARGS=()
if [[ "${DRY_RUN:-0}" == "1" ]]; then
  DRY_RUN_ARGS+=(--dry-run)
fi

CUDA_ARGS=()
if [[ "${MODE}" != "fit-only" && "${SKIP_POSTLOAD_GATE:-0}" != "1" ]]; then
  CUDA_ARGS+=(--cuda-device "${CUDA_DEVICE}")
fi

exec ${PYTHON_BIN} "${WRAPPER_ROOT}/scripts/search/project_geometry_motif.py" \
  "${SOURCE_ABS}" \
  --out-dir "${OUT}" \
  --mode "${MODE}" \
  --mpi-ranks "${RANKS}" \
  --max-lumps "${MAX_LUMPS}" \
  --ftl-L "${FTL_L}" \
  --grtresna-L "${GRTRESNA_L}" \
  --gridinit-n "${GRIDINIT_N}" \
  --stop-time "${STOP_TIME}" \
  --plot-interval "${PLOT_INTERVAL}" \
  "${CUDA_ARGS[@]}" \
  "${POSTLOAD_ARGS[@]}" \
  "${DRY_RUN_ARGS[@]}"
