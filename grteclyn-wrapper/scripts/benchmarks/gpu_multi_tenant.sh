#!/usr/bin/env bash
# Benchmark concurrent GRTeclyn evolutions on a single GPU to size slots_per_gpu.
#
# Usage:
#   GPU_ID=0 MAX_CONCURRENCY=5 ./scripts/benchmarks/gpu_multi_tenant.sh [output.json]
#
# Launches N identical short dry-run episodes on CUDA_VISIBLE_DEVICES=$GPU_ID
# for N in 1..MAX_CONCURRENCY and records wall time per job.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
OUTPUT="${1:-${WRAPPER_ROOT}/benchmark_gpu_slots.json}"

GPU_ID="${GPU_ID:-0}"
MAX_CONCURRENCY="${MAX_CONCURRENCY:-5}"
STOP_TIME="${STOP_TIME:-0.05}"

cd "${WRAPPER_ROOT}"

RESULTS=()
BASELINE=""

for n in $(seq 1 "${MAX_CONCURRENCY}"); do
  echo "[benchmark] launching ${n} concurrent job(s) on GPU ${GPU_ID}"
  start_ts=$(date +%s.%N)
  pids=()
  for i in $(seq 1 "${n}"); do
    (
      job_start=$(date +%s.%N)
      CUDA_VISIBLE_DEVICES="${GPU_ID}" python -m grteclyn_wrapper.cli sweep \
        --count 1 \
        --dry-run \
        --cuda-devices "${GPU_ID}" \
        --set "stop_time=${STOP_TIME}" \
        --set "N_full=64" \
        --set "max_level=2" \
        --name "bench_gpu${GPU_ID}_n${n}_job${i}" \
        >/tmp/bench_gpu_${GPU_ID}_n${n}_job${i}.log 2>&1
      job_end=$(date +%s.%N)
      echo "$(echo "${job_end} - ${job_start}" | bc -l)"
    ) &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "${pid}"
  done
  end_ts=$(date +%s.%N)
  wall=$(echo "${end_ts} - ${start_ts}" | bc -l)
  # Per-job wall time at concurrency N equals the batch wall time: all N jobs
  # start together and finish when the slowest completes.
  per_job="${wall}"
  if [[ "${n}" -eq 1 ]]; then
    BASELINE="${per_job}"
  fi
  slowdown="null"
  if [[ -n "${BASELINE}" && "${n}" -gt 1 ]]; then
    slowdown=$(echo "scale=4; ${per_job} / ${BASELINE}" | bc -l)
  fi
  RESULTS+=("{\"concurrency\":${n},\"wall_seconds\":${wall},\"per_job_seconds\":${per_job},\"slowdown_vs_1\":${slowdown}}")
done

RECOMMENDED=1
if [[ -n "${BASELINE}" ]]; then
  for entry in "${RESULTS[@]}"; do
    n=$(echo "${entry}" | python -c 'import json,sys; print(json.loads(sys.stdin.read())["concurrency"])')
    per=$(echo "${entry}" | python -c 'import json,sys; print(json.loads(sys.stdin.read())["per_job_seconds"])')
    ratio=$(python -c "print(${per} / ${BASELINE})")
    ok=$(python -c "print(1 if float('${ratio}') <= 1.25 else 0)")
    if [[ "${ok}" -eq 1 ]]; then
      RECOMMENDED="${n}"
    fi
  done
fi

python - <<PY
import json
from pathlib import Path

data = {
    "gpu_id": int("${GPU_ID}"),
    "max_concurrency_tested": int("${MAX_CONCURRENCY}"),
    "baseline_per_job_seconds": float("${BASELINE}") if "${BASELINE}" else None,
    "recommended_slots_per_gpu": int("${RECOMMENDED}"),
    "decision_rule": "max n where T_gpu[n] <= T_gpu[1] * 1.25",
    "runs": [${RESULTS[0]}$(printf ',%s' "${RESULTS[@]:1}")],
}
Path("${OUTPUT}").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
print(json.dumps(data, indent=2))
PY

echo "[benchmark] wrote ${OUTPUT}; recommended_slots_per_gpu=${RECOMMENDED}"
