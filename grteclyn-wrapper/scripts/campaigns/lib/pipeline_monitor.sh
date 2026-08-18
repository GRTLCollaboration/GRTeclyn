#!/usr/bin/env bash
# Optional nvidia-smi + in-flight pipeline sampling for QD campaigns.
# Enable with PIPELINE_MONITOR=1 before calling qd/run.sh (see general_ftl/run_all.sh).

ftl_pipeline_monitor_begin() {
  local qd_name="$1"
  local gpu_ids="$2"

  FTL_PIPELINE_QD_NAME="${qd_name}"
  FTL_PIPELINE_GPU_IDS="${gpu_ids}"
  FTL_PIPELINE_LOG_DIR="${GRTECLYN_ROOT}/runs/neuralspacetime/_logs"
  FTL_PIPELINE_RUN_DIR="${GRTECLYN_ROOT}/runs/neuralspacetime/search/map_elites/${qd_name}"
  FTL_PIPELINE_LOG="${FTL_PIPELINE_LOG_DIR}/${qd_name}.log"
  FTL_PIPELINE_GPU_SAMPLES="${FTL_PIPELINE_LOG_DIR}/${qd_name}.gpu_samples.csv"
  FTL_PIPELINE_MONITOR_CSV="${FTL_PIPELINE_LOG_DIR}/${qd_name}.pipeline_monitor.csv"
  FTL_PIPELINE_SUMMARY="${FTL_PIPELINE_LOG_DIR}/${qd_name}.pipeline_summary.txt"
  FTL_PIPELINE_DONE="${FTL_PIPELINE_LOG_DIR}/.${qd_name}.done"

  mkdir -p "${FTL_PIPELINE_LOG_DIR}"

  IFS=',' read -r -a _gpu_arr <<< "$(echo "${gpu_ids}" | tr ' ' ',')"
  FTL_PIPELINE_GPU_QUERY="$(IFS=,; echo "${_gpu_arr[*]}")"
  FTL_PIPELINE_NGPU="${#_gpu_arr[@]}"

  echo "Pipeline monitor: ${FTL_PIPELINE_MONITOR_CSV}"
  echo "Campaign log: ${FTL_PIPELINE_LOG}"

  (
    echo "timestamp,gpu_index,util_pct,mem_used_mib"
    while [[ ! -f "${FTL_PIPELINE_DONE}" ]]; do
      ts="$(date -Iseconds)"
      nvidia-smi \
        --query-gpu=index,utilization.gpu,memory.used \
        --format=csv,noheader,nounits \
        -i "${FTL_PIPELINE_GPU_QUERY}" 2>/dev/null | while IFS=, read -r idx util mem; do
        echo "${ts},${idx// /},${util// /},${mem// /}"
      done
      sleep 5
    done
  ) > "${FTL_PIPELINE_GPU_SAMPLES}" &
  FTL_PIPELINE_SAMPLER_PID=$!

  (
    set +eu
    echo "ts,traj_lines,scored,gpu_ok,rejected,grtresna_n,main3d_n,gpu_busy_n,gpu_mem_mib,driver"
    while [[ ! -f "${FTL_PIPELINE_DONE}" ]]; do
      ts="$(date -Iseconds)"
      traj=0 scored=0 gpu_ok=0 rejected=0 grtresna_n=0 main3d_n=0 driver=0
      gpu_busy=0 gpu_mem=0
      if [[ -f "${FTL_PIPELINE_RUN_DIR}/trajectory.jsonl" ]]; then
        traj=$(wc -l < "${FTL_PIPELINE_RUN_DIR}/trajectory.jsonl")
        gpu_ok=$(grep -c '"status": "gpu_ok"' "${FTL_PIPELINE_RUN_DIR}/trajectory.jsonl" 2>/dev/null || true)
        rejected=$(grep -cE '"status": "(solved_ftl_rejected|grtresna_|postload_|convergence_)' "${FTL_PIPELINE_RUN_DIR}/trajectory.jsonl" 2>/dev/null || true)
      fi
      if compgen -G "${FTL_PIPELINE_RUN_DIR}/eval_*/score.json" > /dev/null; then
        scored=$(find "${FTL_PIPELINE_RUN_DIR}" -path '*/eval_*/score.json' 2>/dev/null | wc -l)
      fi
      grtresna_n=$(ps aux | grep 'ScalarFieldBH' | grep "${FTL_PIPELINE_QD_NAME}" | grep -v grep | awk '{print $NF}' | sed 's|.*/eval_||;s|/grtresna.*||' | sort -u | wc -l)
      main3d_n=$(ps aux | grep 'main3d.gnu.CUDA.ex' | grep "${FTL_PIPELINE_QD_NAME}" | grep -v grep | wc -l)
      driver=$(ps aux | grep "grteclyn_wrapper.*${FTL_PIPELINE_QD_NAME}" | grep -v grep | wc -l)

      while IFS=, read -r _idx _u _m; do
        _u="${_u// /}"
        _m="${_m// /}"
        [[ "${_u:-0}" -ge 10 ]] && gpu_busy=$((gpu_busy + 1))
        [[ "${_m:-0}" -gt "${gpu_mem}" ]] && gpu_mem="${_m}"
      done < <(nvidia-smi --query-gpu=index,utilization.gpu,memory.used --format=csv,noheader,nounits -i "${FTL_PIPELINE_GPU_QUERY}" 2>/dev/null | tr -d ' ')

      echo "${ts},${traj},${scored},${gpu_ok},${rejected},${grtresna_n},${main3d_n},${gpu_busy},${gpu_mem},${driver}"
      sleep 10
    done
  ) > "${FTL_PIPELINE_MONITOR_CSV}" &
  FTL_PIPELINE_MONITOR_PID=$!

  ftl_pipeline_monitor_cleanup() {
    touch "${FTL_PIPELINE_DONE}"
    wait "${FTL_PIPELINE_SAMPLER_PID}" 2>/dev/null || true
    wait "${FTL_PIPELINE_MONITOR_PID}" 2>/dev/null || true
  }
  trap ftl_pipeline_monitor_cleanup EXIT
}

ftl_pipeline_monitor_end() {
  ftl_pipeline_monitor_cleanup 2>/dev/null || true
  trap - EXIT

  if [[ -z "${PYTHON_BIN:-}" ]]; then
    return 0
  fi

  # shellcheck disable=SC2016
  ${PYTHON_BIN} - <<PY
import csv
import json
from collections import defaultdict
from pathlib import Path

gpu_path = Path("${FTL_PIPELINE_GPU_SAMPLES}")
mon_path = Path("${FTL_PIPELINE_MONITOR_CSV}")
run_dir = Path("${FTL_PIPELINE_RUN_DIR}")
summary_path = Path("${FTL_PIPELINE_SUMMARY}")

lines = [
    "Pipeline campaign report: ${FTL_PIPELINE_QD_NAME}",
    f"gpus=${FTL_PIPELINE_GPU_IDS}",
    "",
]

if gpu_path.is_file() and gpu_path.stat().st_size > 20:
    by_gpu = defaultdict(list)
    with gpu_path.open(encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            try:
                by_gpu[int(row["gpu_index"])].append(float(row["util_pct"]))
            except (KeyError, ValueError):
                continue
    lines.append("GPU utilization (5s samples):")
    for gpu in sorted(by_gpu):
        vals = by_gpu[gpu]
        busy = sum(1 for v in vals if v >= 10) / len(vals)
        lines.append(
            f"  GPU {gpu}: n={len(vals)} avg={sum(vals)/len(vals):.1f}% "
            f"max={max(vals):.0f}% busy_frac={busy:.2f}"
        )
    lines.append("")

if mon_path.is_file():
    rows = list(csv.DictReader(mon_path.open(encoding="utf-8")))
    if rows:
        n = len(rows)
        dual_main3d = sum(1 for r in rows if int(r.get("main3d_n") or 0) >= 2)
        dual_grt = sum(1 for r in rows if int(r.get("grtresna_n") or 0) >= 2)
        multi_busy = sum(1 for r in rows if int(r.get("gpu_busy_n") or 0) >= 2)
        last = rows[-1]
        lines += [
            "Pipeline overlap (10s samples):",
            f"  samples={n}",
            f"  >=2 GRTresna: {dual_grt} ({100*dual_grt/n:.1f}%)",
            f"  >=2 main3d: {dual_main3d} ({100*dual_main3d/n:.1f}%)",
            f"  >=2 GPUs busy: {multi_busy} ({100*multi_busy/n:.1f}%)",
            f"  last: traj={last.get('traj_lines')} gpu_ok={last.get('gpu_ok')} rejected={last.get('rejected')}",
            "",
        ]

meta = run_dir / "metadata.json"
if meta.is_file():
    m = json.loads(meta.read_text(encoding="utf-8"))
    lines.append(
        f"metadata: use_pipeline={m.get('use_pipeline')} gpu_slots={m.get('gpu_slots')} "
        f"slots_per_gpu={m.get('slots_per_gpu')} last_eval_counter={m.get('last_eval_counter')}"
    )

text = "\\n".join(lines) + "\\n"
summary_path.write_text(text, encoding="utf-8")
print(text, end="")
PY
}
