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
BENCH_RUNS="${BENCH_RUNS:-${WRAPPER_ROOT}/../runs/_gpu_benchmark_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${BENCH_RUNS}"

cd "${WRAPPER_ROOT}"

WALL_TIMES=()

for n in $(seq 1 "${MAX_CONCURRENCY}"); do
  echo "[benchmark] launching ${n} concurrent job(s) on GPU ${GPU_ID}"
  start_ts=$(date +%s.%N)
  pids=()
  for i in $(seq 1 "${n}"); do
    (
      CUDA_VISIBLE_DEVICES="${GPU_ID}" PYTHONPATH="${WRAPPER_ROOT}/src" python -m grteclyn_wrapper \
        --runs-dir "${BENCH_RUNS}" \
        --dry-run \
        --cuda-devices "${GPU_ID}" \
        --set "stop_time=${STOP_TIME}" \
        --set "N_full=64" \
        --set "max_level=2" \
        --name "bench_gpu${GPU_ID}_n${n}_job${i}" \
        sweep \
        --count 1 \
        >/tmp/bench_gpu_${GPU_ID}_n${n}_job${i}.log 2>&1
    ) &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "${pid}"
  done
  end_ts=$(date +%s.%N)
  wall=$(python -c "print(${end_ts} - ${start_ts})")
  WALL_TIMES+=("${n}:${wall}")
done

python - <<PY
import json
from pathlib import Path

gpu_id = int("${GPU_ID}")
max_concurrency = int("${MAX_CONCURRENCY}")
raw = """$(printf '%s\n' "${WALL_TIMES[@]}")""".strip().splitlines()
runs = []
baseline = None
for line in raw:
    n_s, wall_s = line.split(":", 1)
    n = int(n_s)
    wall = float(wall_s)
    per_job = wall
    if n == 1:
        baseline = per_job
    slowdown = None if baseline is None or n == 1 else round(per_job / baseline, 4)
    runs.append({
        "concurrency": n,
        "wall_seconds": wall,
        "per_job_seconds": per_job,
        "slowdown_vs_1": slowdown,
    })

recommended = 1
if baseline is not None:
    for row in runs:
        if row["per_job_seconds"] / baseline <= 1.25:
            recommended = row["concurrency"]

data = {
    "gpu_id": gpu_id,
    "max_concurrency_tested": max_concurrency,
    "baseline_per_job_seconds": baseline,
    "recommended_slots_per_gpu": recommended,
    "decision_rule": "max n where T_gpu[n] <= T_gpu[1] * 1.25",
    "runs": runs,
}
Path("${OUTPUT}").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
print(json.dumps(data, indent=2))
print(f"[benchmark] recommended_slots_per_gpu={recommended}", flush=True)
PY

echo "[benchmark] wrote ${OUTPUT}"
