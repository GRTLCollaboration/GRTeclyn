#!/usr/bin/env bash
# Resolution + long-time promotion ladder for smoke survivors.
# Skips REJECTED_LONG_RUN candidates (e.g. random_000).
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

PYTHON_BIN="${PYTHON_BIN:-python}"
if command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
  PYTHON_BIN="uv run python"
fi

STOP_TIME="${STOP_TIME:-5.0}"
RESOLUTIONS="${RESOLUTIONS:-64 96 128}"
RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/radialrecipe_gpu_promote}"
BUILD="${BUILD:-0}"
CUDA_DEVICE_START="${CUDA_DEVICE_START:-0}"

# Comma-separated override, else default promoted survivors.
if [[ -n "${PROMOTE_CANDIDATES:-}" ]]; then
  IFS=',' read -r -a CANDIDATES <<< "${PROMOTE_CANDIDATES}"
else
  CANDIDATES=(ellis_bronnikov bubble_wall_016)
fi

REJECTED=(random_000)

echo "== RadialRecipe GPU promotion ladder =="
echo "stop_time=${STOP_TIME}  N_full in: ${RESOLUTIONS}"
echo "candidates: ${CANDIDATES[*]}"
echo "runs_dir: ${RUNS_DIR}"
echo

gpu="${CUDA_DEVICE_START}"
for candidate in "${CANDIDATES[@]}"; do
  for rejected in "${REJECTED[@]}"; do
    if [[ "${candidate}" == "${rejected}" ]]; then
      echo "SKIP ${candidate}: in REJECTED_LONG_RUN (${rejected})"
      continue 2
    fi
  done

  for n_full in ${RESOLUTIONS}; do
    echo "--- ${candidate} N_full=${n_full} t=${STOP_TIME} (GPU ${gpu}) ---"
    if [[ "${candidate}" == "ellis_bronnikov" ]]; then
      SEED_NAME="${candidate}" CANDIDATE_ID="" \
        STOP_TIME="${STOP_TIME}" N_FULL="${n_full}" BUILD="${BUILD}" \
        RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
        CUDA_VISIBLE_DEVICES_OVERRIDE="${gpu}" \
        bash "${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
    else
      SEED_NAME="" CANDIDATE_ID="${candidate}" \
        STOP_TIME="${STOP_TIME}" N_FULL="${n_full}" BUILD="${BUILD}" \
        RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
        CUDA_VISIBLE_DEVICES_OVERRIDE="${gpu}" \
        bash "${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"
    fi
    gpu=$((gpu + 1))
    BUILD=0
  done
done

echo
echo "== Promotion summary =="
${PYTHON_BIN} - <<PY
import json
from pathlib import Path
from grteclyn_wrapper.metrics import dataclass_to_dict, read_episode_metrics

runs = Path(${RUNS_DIR@Q})
rows = []
for ep in sorted(runs.glob("*_gpu_t*_${RUN_STAMP@Q}")):
    try:
        m = read_episode_metrics(ep)
    except Exception:
        continue
    rows.append({
        "episode": ep.name,
        "max_H_L2": m.constraints.max_hamiltonian_l2,
        "final_H_L2": m.constraints.final_hamiltonian_l2,
        "min_lapse": m.collapse.min_lapse,
        "min_chi": m.collapse.min_chi,
    })
print(json.dumps(rows, indent=2))
PY
