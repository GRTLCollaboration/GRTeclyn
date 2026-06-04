#!/usr/bin/env bash
# CMA-ES closed-loop optimization for RadialRecipe metric discovery.
#
# Runs the full propose→constrain→preflight→simulate→score→update loop
# using all available GPUs in parallel (one candidate per GPU per generation).
#
# Usage:
#   bash grteclyn-wrapper/scripts/run_optimize_loop.sh
#
# Environment overrides:
#   MAX_GENERATIONS=50    Number of CMA-ES generations
#   POPULATION_SIZE=8     Candidates per generation (default: num GPUs)
#   SIGMA0=0.3            Initial CMA-ES step size
#   STOP_TIME=2.0         Simulation time per candidate
#   L_FULL=64.0           Full-domain width
#   N_FULL=128            Base grid resolution
#   MAX_LEVEL=2           AMR levels
#   GPU_IDS="0 1 2 3 4 5 6 7"   GPU indices to use
#   SEED_NAME=ellis_bronnikov    Warm-start from known solution
#   SEED=42               Random seed for reproducibility
#   BUILD=0               Skip building (if already built)
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
if command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
  PYTHON_BIN="uv run python"
fi

EXAMPLE_DIR="${GRTECLYN_ROOT}/Examples/RadialRecipe"
WRAPPER_DIR="${GRTECLYN_ROOT}/grteclyn-wrapper"
EXECUTABLE="${EXAMPLE_DIR}/main3d.${COMP:-gnu}.CUDA.ex"
CUDA_ARCH="${CUDA_ARCH:-90}"

MAX_GENERATIONS="${MAX_GENERATIONS:-50}"
POPULATION_SIZE="${POPULATION_SIZE:-}"
SIGMA0="${SIGMA0:-0.3}"
STOP_TIME="${STOP_TIME:-2.0}"
L_FULL="${L_FULL:-64.0}"
N_FULL="${N_FULL:-128}"
MAX_LEVEL="${MAX_LEVEL:-2}"
SEED="${SEED:-}"
SEED_NAME="${SEED_NAME:-}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/optimize}"

echo "============================================"
echo " CMA-ES Closed-Loop Optimization"
echo "============================================"
echo "GRTeclyn root   : ${GRTECLYN_ROOT}"
echo "Executable      : ${EXECUTABLE}"
echo "GPU IDs         : ${GPU_IDS}"
echo "Max generations : ${MAX_GENERATIONS}"
echo "Population size : ${POPULATION_SIZE:-auto (= num GPUs)}"
echo "Sigma0          : ${SIGMA0}"
echo "Stop time       : ${STOP_TIME}"
echo "L_full          : ${L_FULL}"
echo "N_full          : ${N_FULL}"
echo "Max AMR level   : ${MAX_LEVEL}"
echo "Seed name       : ${SEED_NAME:-none (start from center)}"
echo "Runs dir        : ${RUNS_DIR}"
echo "============================================"
echo

if [[ "${BUILD:-1}" == "1" ]]; then
  echo "== Building RadialRecipe CUDA executable (USE_MPI=FALSE) =="
  make -C "${EXAMPLE_DIR}" COMP="${COMP:-gnu}" USE_CUDA=TRUE USE_MPI=FALSE \
    CUDA_ARCH="${CUDA_ARCH}" -j"$(nproc)"
  echo
fi

if [[ ! -x "${EXECUTABLE}" ]]; then
  echo "Executable not found: ${EXECUTABLE}" >&2
  echo "Try: BUILD=1 $0" >&2
  exit 2
fi

mkdir -p "${RUNS_DIR}"

GLOBAL_ARGS=()
OPTIMIZER_ARGS=()
if [[ -n "${SEED_NAME}" ]]; then
  GLOBAL_ARGS+=(--seed-name "${SEED_NAME}")
fi
if [[ -n "${SEED}" ]]; then
  OPTIMIZER_ARGS+=(--seed "${SEED}")
fi
if [[ -n "${POPULATION_SIZE}" ]]; then
  OPTIMIZER_ARGS+=(--population-size "${POPULATION_SIZE}")
fi

echo "== Starting CMA-ES loop =="
echo "  Each generation evaluates candidates in parallel across GPUs."
echo "  Trajectory: <runs_dir>/optimize_<timestamp>/trajectory.jsonl"
echo

(
cd "${WRAPPER_DIR}"
${PYTHON_BIN} -m grteclyn_wrapper \
  --example RadialRecipe \
  --constrained --phantom --preflight \
  --consume-plotfiles --consumer-delete \
  --consumer-radii 4 8 \
  --runs-dir "${RUNS_DIR}" \
  --set stop_time="${STOP_TIME}" \
  --set L_full="${L_FULL}" \
  --set N_full="${N_FULL}" \
  --set plot_interval=10 \
  --set checkpoint_interval=-1 \
  --set dt_multiplier=0.02 \
  --set max_level="${MAX_LEVEL}" \
  --score-weight ftl_shortcut=5.0 \
  --score-weight expansion_asymmetry=2.0 \
  --score-weight nonflat_geometry=1.0 \
  --score-weight comoving_stability=2.5 \
  --ftl-L "${FTL_L:-${L_FULL}}" \
  "${GLOBAL_ARGS[@]}" \
  -- optimize \
  --max-generations "${MAX_GENERATIONS}" \
  --sigma0 "${SIGMA0}" \
  "${OPTIMIZER_ARGS[@]}" \
  --gpu-ids ${GPU_IDS}
)

echo
echo "== Optimization complete =="
echo "Results in: ${RUNS_DIR}"
echo
echo "Next steps:"
echo "  1. Inspect trajectory.jsonl for score progression"
echo "  2. Run diagnostics on best episode:"
echo "     bash grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh <best_episode_dir>"
echo "  3. Promote best candidate to higher resolution (t=5, N=128):"
echo "     STOP_TIME=5 N_FULL=128 CANDIDATE_DIR=<best_episode> bash grteclyn-wrapper/scripts/radial/run_radialrecipe_gpu_smoke.sh"
