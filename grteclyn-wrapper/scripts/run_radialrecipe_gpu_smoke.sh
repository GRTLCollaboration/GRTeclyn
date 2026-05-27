#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
if command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
  PYTHON_BIN="uv run python"
fi

OPENMPI_ROOT="${OPENMPI_ROOT:-${GRTECLYN_ROOT}/../local/openmpi-5.0.8}"
if [[ -d "${OPENMPI_ROOT}/bin" ]]; then
  export PATH="${OPENMPI_ROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
fi

EXAMPLE_DIR="${GRTECLYN_ROOT}/Examples/RadialRecipe"
EXECUTABLE="${EXAMPLE_DIR}/main3d.${COMP:-gnu}.CUDA.ex"
CUDA_ARCH="${CUDA_ARCH:-90}"

# Initial data: exactly one of SEED_NAME (seeds.py) or CANDIDATE_ID (validate_guesser).
SEED_NAME="${SEED_NAME:-}"
CANDIDATE_ID="${CANDIDATE_ID:-}"
NONSPHERICAL_ID="${NONSPHERICAL_ID:-}"
VALIDATION_SEED="${VALIDATION_SEED:-42}"

if [[ -n "${SEED_NAME}" && -n "${CANDIDATE_ID}" ]]; then
  echo "Set only one of SEED_NAME or CANDIDATE_ID." >&2
  exit 2
fi
if [[ -n "${SEED_NAME}" && -n "${NONSPHERICAL_ID}" ]]; then
  echo "Set only one of SEED_NAME or NONSPHERICAL_ID." >&2
  exit 2
fi
if [[ -n "${CANDIDATE_ID}" && -n "${NONSPHERICAL_ID}" ]]; then
  echo "Set only one of CANDIDATE_ID or NONSPHERICAL_ID." >&2
  exit 2
fi

SOURCE_LABEL="${SOURCE_LABEL:-}"
if [[ -z "${SOURCE_LABEL}" ]]; then
  if [[ -n "${CANDIDATE_ID}" ]]; then
    SOURCE_LABEL="${CANDIDATE_ID}"
  elif [[ -n "${NONSPHERICAL_ID}" ]]; then
    SOURCE_LABEL="${NONSPHERICAL_ID}"
  else
    SEED_NAME="${SEED_NAME:-ellis_bronnikov}"
    SOURCE_LABEL="${SEED_NAME}"
  fi
fi

FIELD_FLAG="${FIELD_FLAG:---phantom}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/radialrecipe_gpu_smoke}"
DRY_NAME="${DRY_NAME:-${SOURCE_LABEL}_dry_${RUN_STAMP}}"
GPU_NAME="${GPU_NAME:-${SOURCE_LABEL}_gpu_t${STOP_TIME:-2}_${RUN_STAMP}}"
CUDA_DEVICE="${CUDA_VISIBLE_DEVICES_OVERRIDE:-0}"

N_FULL="${N_FULL:-64}"
MAX_LEVEL="${MAX_LEVEL:-0}"
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-1}"
CHECKPOINT_INTERVAL="${CHECKPOINT_INTERVAL:--1}"
DT_MULTIPLIER="${DT_MULTIPLIER:-0.02}"

WRAPPER_ID_ARGS=()
if [[ -n "${CANDIDATE_ID}" ]]; then
  WRAPPER_ID_ARGS=(--candidate-id "${CANDIDATE_ID}" --validation-seed "${VALIDATION_SEED}")
elif [[ -n "${NONSPHERICAL_ID}" ]]; then
  WRAPPER_ID_ARGS=(--nonspherical-id "${NONSPHERICAL_ID}" --validation-seed "${VALIDATION_SEED}")
else
  WRAPPER_ID_ARGS=(--seed-name "${SEED_NAME}")
fi

echo "== RadialRecipe GPU smoke pipeline =="
echo "GRTeclyn root : ${GRTECLYN_ROOT}"
echo "Example       : ${EXAMPLE_DIR}"
echo "Executable    : ${EXECUTABLE}"
echo "Source        : ${SOURCE_LABEL} ${FIELD_FLAG}"
echo "Runs dir      : ${RUNS_DIR}"
echo "CUDA device   : ${CUDA_DEVICE}"
echo

if [[ "${BUILD:-1}" == "1" ]]; then
  echo "== Step 1/5: build RadialRecipe CUDA executable (single-GPU, USE_MPI=FALSE) =="
  make -C "${EXAMPLE_DIR}" COMP="${COMP:-gnu}" USE_CUDA=TRUE USE_MPI=FALSE \
    CUDA_ARCH="${CUDA_ARCH}" -j"${JOBS:-$(nproc)}"
fi

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable not found or not executable: ${EXECUTABLE}" >&2
  echo "Try: BUILD=1 ${BASH_SOURCE[0]}" >&2
  exit 2
fi

echo "== Step 2/5: Python pre-flight validation/dry-run params generation =="
${PYTHON_BIN} -m grteclyn_wrapper \
  --example RadialRecipe \
  "${WRAPPER_ID_ARGS[@]}" \
  --constrained \
  ${FIELD_FLAG} \
  --preflight \
  --dry-run \
  --runs-dir "${RUNS_DIR}" \
  --name "${DRY_NAME}" \
  --set stop_time="${STOP_TIME}" \
  --set plot_interval="${PLOT_INTERVAL}" \
  --set checkpoint_interval="${CHECKPOINT_INTERVAL}" \
  --set N_full="${N_FULL}" \
  --set max_level="${MAX_LEVEL}" \
  --set dt_multiplier="${DT_MULTIPLIER}" \
  reproduce

echo
echo "Dry-run episode: ${RUNS_DIR}/${DRY_NAME}"
echo "Generated params preview:"
grep -E "^(recipe_|stop_time|N_full|max_level|dt_multiplier|plot_interval|checkpoint_interval)" \
  "${RUNS_DIR}/${DRY_NAME}/params.txt" | head -80
echo

echo "== Step 3/5: C++ check_params + short GPU evolution =="
${PYTHON_BIN} -m grteclyn_wrapper \
  --example RadialRecipe \
  "${WRAPPER_ID_ARGS[@]}" \
  --constrained \
  ${FIELD_FLAG} \
  --preflight \
  --cuda-devices "${CUDA_DEVICE}" \
  --runs-dir "${RUNS_DIR}" \
  --name "${GPU_NAME}" \
  --set stop_time="${STOP_TIME}" \
  --set plot_interval="${PLOT_INTERVAL}" \
  --set checkpoint_interval="${CHECKPOINT_INTERVAL}" \
  --set N_full="${N_FULL}" \
  --set max_level="${MAX_LEVEL}" \
  --set dt_multiplier="${DT_MULTIPLIER}" \
  reproduce

EPISODE_DIR="${RUNS_DIR}/${GPU_NAME}"

echo
echo "== Step 4/5: diagnostics summary =="
${PYTHON_BIN} - <<PY
from pathlib import Path
from grteclyn_wrapper.metrics import dataclass_to_dict, read_episode_metrics
import json

episode = Path(${EPISODE_DIR@Q})
metrics = read_episode_metrics(episode)
print(json.dumps(dataclass_to_dict(metrics), indent=2))

constraint_path = episode / "data" / "constraint_norms.dat"
collapse_path = episode / "data" / "collapse_diagnostics.dat"
print()
print(f"constraint_norms.dat: {constraint_path} exists={constraint_path.exists()}")
print(f"collapse_diagnostics.dat: {collapse_path} exists={collapse_path.exists()}")
PY

echo
echo "== Step 5/5: quick pass/fail hints =="
echo "Inspect:"
echo "  ${EPISODE_DIR}/run.log"
echo "  ${EPISODE_DIR}/data/constraint_norms.dat"
echo "  ${EPISODE_DIR}/score.json"
echo
echo "First target: Hamiltonian L2 should remain controlled and not grow exponentially."
echo "For this smoke test, use ||H||_L2 < 1e-1 as an initial rough acceptance target."
echo
echo "Examples:"
echo "  SEED_NAME=ellis_bronnikov ${BASH_SOURCE[0]}"
echo "  CANDIDATE_ID=bubble_wall_016 ${BASH_SOURCE[0]}"
echo "  CANDIDATE_ID=random_000 STOP_TIME=2 ${BASH_SOURCE[0]}"
echo "  NONSPHERICAL_ID=quadrupole_bubble_001 ${BASH_SOURCE[0]}"
