#!/usr/bin/env bash
# GPU smoke batch for accepted non-spherical RadialRecipe candidates.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"
SMOKE="${SCRIPT_DIR}/run_radialrecipe_gpu_smoke.sh"

RUN_STAMP="${RUN_STAMP:-$(date +%Y%m%d_%H%M%S)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/radialrecipe_nonspherical}"
STOP_TIME="${STOP_TIME:-2.0}"
N_FULL="${N_FULL:-64}"
BUILD="${BUILD:-0}"
CUDA_START="${CUDA_START:-0}"
CONSUME_PLOTFILES="${CONSUME_PLOTFILES:-1}"
CONSUMER_DELETE="${CONSUMER_DELETE:-1}"
CONSUMER_RADII="${CONSUMER_RADII:-4 8}"

# Accepted from nonspherical_ray_validation.csv (skip rejected_negative_chi).
CANDIDATES=(
  dipole_lopsided_000
  quadrupole_bubble_001
  mixed_modes_002
  random_angular_004
  random_angular_005
  random_angular_006
  random_angular_007
)

# Optional: include bad control for comparison (preflight may still fail on GPU).
if [[ "${INCLUDE_BAD_CONTROL:-0}" == "1" ]]; then
  CANDIDATES+=(strong_quadrupole_bad_003)
fi

echo "== Non-spherical GPU batch =="
echo "candidates: ${#CANDIDATES[@]}"
echo "runs_dir:   ${RUNS_DIR}"
echo "stop_time:  ${STOP_TIME}  N_full: ${N_FULL}"
echo "consumer:   CONSUME=${CONSUME_PLOTFILES} DELETE=${CONSUMER_DELETE} RADII=${CONSUMER_RADII}"
echo

mkdir -p "${RUNS_DIR}"
PIDS=()
GPU="${CUDA_START}"

for cid in "${CANDIDATES[@]}"; do
  log="${RUNS_DIR}/${cid}_${RUN_STAMP}.log"
  echo "Launch ${cid} on GPU ${GPU} -> ${log}"
  (
    NONSPHERICAL_ID="${cid}" \
    SEED_NAME="" CANDIDATE_ID="" \
    BUILD=0 STOP_TIME="${STOP_TIME}" N_FULL="${N_FULL}" \
    CONSUME_PLOTFILES="${CONSUME_PLOTFILES}" \
    CONSUMER_DELETE="${CONSUMER_DELETE}" \
    CONSUMER_RADII="${CONSUMER_RADII}" \
    RUN_STAMP="${RUN_STAMP}" RUNS_DIR="${RUNS_DIR}" \
    CUDA_VISIBLE_DEVICES_OVERRIDE="${GPU}" \
    bash "${SMOKE}" > "${log}" 2>&1
    echo "DONE ${cid} (GPU ${GPU})" >> "${log}"
  ) &
  PIDS+=($!)
  GPU=$((GPU + 1))
done

fail=0
for pid in "${PIDS[@]}"; do
  wait "${pid}" || fail=1
done

echo
echo "== Batch summary =="
export PYTHONPATH="${SCRIPT_DIR}/../src:${PYTHONPATH:-}"
python3 - <<PY
import json
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics

runs = Path("${RUNS_DIR}")
stamp = "${RUN_STAMP}"
rows = []
for ep in sorted(runs.glob(f"*_gpu_t*_{stamp}")):
    try:
        m = read_episode_metrics(ep)
        rows.append({
            "episode": ep.name,
            "max_H_L2": m.constraints.max_hamiltonian_l2,
            "final_H_L2": m.constraints.final_hamiltonian_l2,
            "min_lapse": m.collapse.min_lapse,
            "min_chi": m.collapse.min_chi,
        })
    except Exception as exc:
        rows.append({"episode": ep.name, "error": str(exc)})
print(json.dumps(rows, indent=2))
PY

exit "${fail}"
