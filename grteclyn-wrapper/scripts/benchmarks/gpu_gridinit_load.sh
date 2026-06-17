#!/usr/bin/env bash
# GPU load test: replay duplicated gridinit from a QD eval (GPU-only, no GRTresna).
#
# Usage:
#   SRC_EVAL=runs/grtresna_qd/.../eval_000005 N_JOBS=2 GPU_ID=0 \
#     ./scripts/benchmarks/gpu_gridinit_load.sh
#
# Sweeps concurrency 1..MAX_CONCURRENCY unless N_JOBS is set (single level).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SRC_EVAL="${SRC_EVAL:?SRC_EVAL required (path to eval_00000N)}"
GPU_ID="${GPU_ID:-0}"
MAX_CONCURRENCY="${MAX_CONCURRENCY:-5}"
N_JOBS="${N_JOBS:-}"
STOP_TIME="${STOP_TIME:-4}"
PLOT_INTERVAL="${PLOT_INTERVAL:-40}"
N_FULL="${N_FULL:-128}"
L_FULL="${L_FULL:-64}"
MAX_LEVEL="${MAX_LEVEL:-2}"
BENCH="${BENCH:-${WRAPPER_ROOT}/../runs/_gpu_gridinit_load_$(date +%Y%m%d_%H%M%S)}"
APPEND_SUMMARY="${APPEND_SUMMARY:-0}"

SRC_EVAL="$(cd "$(dirname "${SRC_EVAL}")" && pwd)/$(basename "${SRC_EVAL}")"
GRIDINIT_SRC="${SRC_EVAL}/initial_data.gridinit"
if [[ ! -f "${GRIDINIT_SRC}" ]]; then
  echo "missing gridinit: ${GRIDINIT_SRC}" >&2
  exit 1
fi

mkdir -p "${BENCH}"
echo "src_eval=${SRC_EVAL}"
echo "bench=${BENCH}"
echo "gridinit=$(du -h "${GRIDINIT_SRC}" | awk '{print $1}')"

launch_job() {
  local n=$1 i=$2
  local jobdir="${BENCH}/n${n}/job_${i}"
  mkdir -p "${jobdir}"
  cp "${GRIDINIT_SRC}" "${jobdir}/initial_data.gridinit"
  CUDA_VISIBLE_DEVICES="${GPU_ID}" "${WRAPPER_ROOT}/.venv/bin/python" \
    "${WRAPPER_ROOT}/scripts/campaigns/hq/replay_eval.py" \
    "${SRC_EVAL}" \
    --gridinit "${jobdir}/initial_data.gridinit" \
    --name "dup_${i}" \
    --runs-dir "${BENCH}/n${n}" \
    --gpu "${GPU_ID}" \
    --n-full "${N_FULL}" --l-full "${L_FULL}" \
    --max-level "${MAX_LEVEL}" \
    --stop-time "${STOP_TIME}" \
    --plot-interval "${PLOT_INTERVAL}" \
    --objective-mode general_ftl \
    --consumer-keep-last 2 \
    > "${jobdir}/replay.log" 2>&1 &
  echo $!
}

run_level() {
  local n=$1
  echo ""
  echo "========== N=${n} concurrent gridinit replays on GPU ${GPU_ID} (stop_time=${STOP_TIME}) =========="
  local start_ts end_ts wall_s
  start_ts=$(date +%s.%N)
  local pids=()
  for i in $(seq 1 "${n}"); do
    pids+=("$(launch_job "${n}" "${i}")")
  done
  local max_m3d=0 max_mem=0
  local samples="${BENCH}/n${n}/mem_samples.csv"
  echo "ts,main3d_n,mem_mib,gpu_util" > "${samples}"
  local t=0
  while :; do
    local alive=0
    for pid in "${pids[@]}"; do
      kill -0 "${pid}" 2>/dev/null && alive=$((alive + 1))
    done
    local m3d mem util
    m3d=$(ps aux | grep 'main3d.gnu.CUDA.ex' | grep "${BENCH}" | grep -v grep | wc -l || true)
    mem=$(nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -i "${GPU_ID}" | tr -d ' ' || echo 0)
    util=$(nvidia-smi --query-gpu=utilization.gpu --format=csv,noheader,nounits -i "${GPU_ID}" | tr -d ' ' || echo 0)
    [[ "${m3d}" -gt "${max_m3d}" ]] && max_m3d=${m3d}
    [[ "${mem}" -gt "${max_mem}" ]] && max_mem=${mem}
    echo "$(date -Iseconds),${m3d},${mem},${util}" >> "${samples}"
    printf '\r  t+%ds main3d=%s mem=%sMiB (~%.1fGB) util=%s%% alive=%s' \
      "${t}" "${m3d}" "${mem}" "$(awk "BEGIN{print ${mem}/1024}")" "${util}" "${alive}"
    [[ "${alive}" -eq 0 ]] && break
    sleep 5
    t=$((t + 5))
  done
  end_ts=$(date +%s.%N)
  wall_s=$(awk "BEGIN{printf \"%.1f\", ${end_ts} - ${start_ts}}")
  echo ""
  echo "  MAX main3d=${max_m3d}  PEAK mem=${max_mem}MiB (~$(awk "BEGIN{printf \"%.1f\", ${max_mem}/1024}")GB)"
  echo "  WALL ${wall_s}s (~$(awk "BEGIN{printf \"%.1f\", ${wall_s}/60}") min) for ${n} concurrent jobs"
  echo "${n},${max_m3d},${max_mem},${wall_s}" >> "${BENCH}/summary.csv"
  "${WRAPPER_ROOT}/.venv/bin/python" - <<PY
import json
from pathlib import Path

bench = Path("${BENCH}")
n = int("${n}")
wall_s = float("${wall_s}")
rows = []
for ep in sorted((bench / f"n{n}").glob("dup_*")):
    meta = ep / "metadata.json"
    if not meta.is_file():
        continue
    data = json.loads(meta.read_text(encoding="utf-8"))
    rows.append({
        "episode": ep.name,
        "simulation_elapsed_seconds": data.get("simulation_elapsed_seconds"),
        "exit_code": data.get("simulation_exit_code"),
    })
sim_times = [r["simulation_elapsed_seconds"] for r in rows if r["simulation_elapsed_seconds"] is not None]
timing = {
    "n": n,
    "stop_time": float("${STOP_TIME}"),
    "wall_seconds": wall_s,
    "wall_minutes": round(wall_s / 60.0, 2),
    "peak_mem_mib": int("${max_mem}"),
    "max_main3d": int("${max_m3d}"),
    "per_job_sim_seconds": rows,
    "sim_seconds_mean": round(sum(sim_times) / len(sim_times), 1) if sim_times else None,
    "sim_seconds_max": round(max(sim_times), 1) if sim_times else None,
    "slowdown_vs_solo_wall": None,
}
solo = bench / "timing_n1.json"
if solo.is_file() and n > 1:
    solo_wall = json.loads(solo.read_text(encoding="utf-8")).get("wall_seconds")
    if solo_wall:
        timing["slowdown_vs_solo_wall"] = round(wall_s / float(solo_wall), 3)
out = bench / f"timing_n{n}.json"
out.write_text(json.dumps(timing, indent=2) + "\n", encoding="utf-8")
print(json.dumps(timing, indent=2))
PY
}

if [[ "${APPEND_SUMMARY}" == "1" && -f "${BENCH}/summary.csv" ]]; then
  : # keep existing summary
else
  echo "n,max_main3d,peak_mem_mib,wall_seconds" > "${BENCH}/summary.csv"
fi
if [[ -n "${N_JOBS}" ]]; then
  run_level "${N_JOBS}"
else
  for n in $(seq 1 "${MAX_CONCURRENCY}"); do
    run_level "${n}"
  done
fi

echo ""
echo "summary: ${BENCH}/summary.csv"
cat "${BENCH}/summary.csv"
