#!/usr/bin/env bash
# Stage 2 — HQ promotion: fresh GRTresna solve + full-resolution GPU evolution.
#
# Separate from QD/CMA-ES search resolution. Replays elite genomes from a prior
# campaign directory (QD or CMA-ES) at N=256, L=128, ml=3, t=30 with frames,
# incremental scoring, and 4D geodesic in HQ verify mode.
#
# Usage:
#   SOURCE_RUN=runs/grtresna_qd/ftl_4d/ftl_4d_v1 \
#   CANDIDATES="156 0 142 1" \
#   NAME_PREFIX=ftl_4d \
#   bash scripts/campaigns/hq/run_batch.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"
# shellcheck source=../lib/promote_common.sh
source "${CAMPAIGNS_ROOT}/lib/promote_common.sh"

SOURCE_RUN="${SOURCE_RUN:-${QD_RUN:-${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_discovery/qd_ftl_discovery_20260609T162553Z}}"
NAME_PREFIX="${NAME_PREFIX:-ftl_discovery}"
OBJECTIVE_MODE="${OBJECTIVE_MODE:-general_ftl}"
CANDIDATES="${CANDIDATES:-}"
CANDIDATES_FILE="${CANDIDATES_FILE:-}"
TOP_K="${TOP_K:-0}"
MIN_FTL_GEO_EVOL="${MIN_FTL_GEO_EVOL:-${MIN_OPERATIONAL_FTL:-0.0}}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
DRY_RUN="${DRY_RUN:-0}"
REPLAY="${CAMPAIGNS_ROOT}/hq/replay_eval.py"

mkdir -p "${RUNS_DIR}"

declare -a CANDIDATE_ENTRIES=()

if [[ -n "${CANDIDATES}" ]]; then
  # shellcheck disable=SC2206
  pairs=(${CANDIDATES})
  if ((${#pairs[@]} % 2 != 0)); then
    echo "[hq] CANDIDATES must contain eval/gpu pairs (even count), got ${#pairs[@]} tokens" >&2
    exit 2
  fi
  for ((i = 0; i < ${#pairs[@]}; i += 2)); do
    CANDIDATE_ENTRIES+=("${pairs[i]} ${pairs[i + 1]}")
  done
elif [[ -n "${CANDIDATES_FILE}" ]]; then
  if [[ ! -f "${CANDIDATES_FILE}" ]]; then
    echo "[hq] missing CANDIDATES_FILE: ${CANDIDATES_FILE}" >&2
    exit 2
  fi
  while read -r eval_id gpu_id _rest; do
    [[ -z "${eval_id}" || "${eval_id}" == \#* ]] && continue
    CANDIDATE_ENTRIES+=("${eval_id} ${gpu_id}")
  done < "${CANDIDATES_FILE}"
elif [[ "${TOP_K}" -gt 0 ]]; then
  traj="${SOURCE_RUN}/trajectory.jsonl"
  if [[ ! -f "${traj}" ]]; then
    echo "[hq] missing trajectory for TOP_K pick: ${traj}" >&2
    exit 2
  fi
  mapfile -t picked < <(python3 - <<'PY' "${traj}" "${TOP_K}" "${MIN_FTL_GEO_EVOL}" "${GPU_IDS}"
import json
import sys

traj_path, top_k, min_geo, gpu_ids_raw = sys.argv[1:5]
top_k = int(top_k)
min_geo = float(min_geo)
gpus = [g for g in gpu_ids_raw.split() if g]

rows = []
with open(traj_path, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get("status") != "gpu_ok":
            continue
        if rec.get("grtresna_rejected") or rec.get("solved_ftl_rejected") or rec.get("grtresna_failed"):
            continue
        if rec.get("preflight_rejected"):
            continue
        dd = rec.get("descriptor_details") or {}
        comps = rec.get("components") or {}
        geo = float(dd.get("f_geo_evol") or comps.get("ftl_geo_evolving") or 0.0)
        if geo < min_geo:
            continue
        score = float(rec.get("score") or -1.0)
        rows.append((score, int(rec["eval"])))

rows.sort(key=lambda item: item[0], reverse=True)
rows = rows[:top_k]

if not gpus:
    gpus = ["0"]

for idx, (_score, eval_id) in enumerate(rows):
    gpu = gpus[idx % len(gpus)]
    print(f"{eval_id} {gpu}")
PY
)
  if ((${#picked[@]} == 0)); then
    echo "[hq] TOP_K=${TOP_K} found no candidates (min f_geo_evol=${MIN_FTL_GEO_EVOL})" >&2
    exit 2
  fi
  CANDIDATE_ENTRIES=("${picked[@]}")
else
  echo "[hq] set CANDIDATES, CANDIDATES_FILE, or TOP_K>0" >&2
  exit 2
fi

echo "== Stage 2: HQ promotion batch =="
echo "Source run    : ${SOURCE_RUN}"
echo "Runs dir      : ${RUNS_DIR}"
echo "Name prefix   : ${NAME_PREFIX:-<none>}"
echo "Domain        : L=${L_FULL} N=${N_FULL} max_level=${MAX_LEVEL} t=${STOP_TIME}"
echo "Objective     : ${OBJECTIVE_MODE}"
echo "4D geodesic   : mode=${GRTECLYN_EVOLVING_GEODESIC_MODE} dirs=${GRTECLYN_GEO_DIRECTIONS:-x} (full HQ verify)"
echo "Frames        : ${GRTECLYN_FRAMES:-on}"
echo "Candidates    : ${#CANDIDATE_ENTRIES[@]}"
echo

launched=0
for entry in "${CANDIDATE_ENTRIES[@]}"; do
  read -r EVAL_ID GPU_ID <<< "${entry}"
  eval_num="$((10#${EVAL_ID}))"
  eval_padded="$(printf "%06d" "${eval_num}")"
  name="$(promote_build_name "${eval_padded}")"
  source="${SOURCE_RUN}/eval_${eval_padded}"
  out="${RUNS_DIR}/${name}"
  log="${RUNS_DIR}/${name}.log"

  if [[ ! -d "${source}" ]]; then
    echo "[skip] missing source ${source}"
    continue
  fi
  if [[ -e "${out}" && "${FORCE}" != "1" ]]; then
    echo "[skip] output exists ${out} (set FORCE=1 to replace)"
    continue
  fi

  echo "[launch] ${name} GPU=${GPU_ID} source=eval_${eval_padded} log=${log}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "  dry-run: would nohup replay_eval.py ${source} --name ${name} --gpu ${GPU_ID}"
    launched=$((launched + 1))
    continue
  fi

  gridinit_args=()
  [[ -n "${GRIDINIT:-}" ]] && gridinit_args=(--gridinit "${GRIDINIT}")

  geodesic_args=()
  if [[ "${OBJECTIVE_MODE}" != "critical_collapse" ]]; then
    geodesic_args=(--evolving-geodesic)
  fi

  if [[ "${FOREGROUND:-0}" == "1" ]]; then
    echo "  foreground: replay_eval.py ${source} -> ${out}"
    # shellcheck disable=SC2086
    ${PYTHON_BIN} "${REPLAY}" \
      "${source}" \
      --name "${name}" \
      --runs-dir "${RUNS_DIR}" \
      --gpu "${GPU_ID}" \
      --n-full "${N_FULL}" \
      --l-full "${L_FULL}" \
      --grtresna-domain-l "${GRTRESNA_DOMAIN_L}" \
      --max-level "${MAX_LEVEL}" \
      --regrid-threshold "${REGRID_THRESHOLD}" \
      --stop-time "${STOP_TIME}" \
      --plot-interval "${PLOT_INTERVAL}" \
      --grtresna-ranks "${GRTRESNA_RANKS}" \
      --grtresna-iterations "${GRTRESNA_ITERATIONS}" \
      --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
      --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
      --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
      --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
      --objective-mode "${OBJECTIVE_MODE}" \
      "${geodesic_args[@]}" \
      "${gridinit_args[@]}"
  else
    # shellcheck disable=SC2086
    nohup ${PYTHON_BIN} "${REPLAY}" \
      "${source}" \
      --name "${name}" \
      --runs-dir "${RUNS_DIR}" \
      --gpu "${GPU_ID}" \
      --n-full "${N_FULL}" \
      --l-full "${L_FULL}" \
      --grtresna-domain-l "${GRTRESNA_DOMAIN_L}" \
      --max-level "${MAX_LEVEL}" \
      --regrid-threshold "${REGRID_THRESHOLD}" \
      --stop-time "${STOP_TIME}" \
      --plot-interval "${PLOT_INTERVAL}" \
      --grtresna-ranks "${GRTRESNA_RANKS}" \
      --grtresna-iterations "${GRTRESNA_ITERATIONS}" \
      --grtresna-timeout "${GRTRESNA_TIMEOUT}" \
      --grtresna-max-ham-pct "${GRTRESNA_MAX_HAM_PCT}" \
      --grtresna-max-mom-pct "${GRTRESNA_MAX_MOM_PCT}" \
      --consumer-keep-last "${CONSUMER_KEEP_LAST}" \
      --objective-mode "${OBJECTIVE_MODE}" \
      "${geodesic_args[@]}" \
      "${gridinit_args[@]}" \
      > "${log}" 2>&1 &
    echo "  pid=$!"
  fi
  launched=$((launched + 1))
done

echo
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "Dry-run: would launch ${launched} HQ promotions."
elif [[ "${FOREGROUND:-0}" == "1" ]]; then
  echo "Completed ${launched} HQ promotion(s) in foreground."
else
  echo "Launched ${launched} HQ promotions (background). Monitor:"
  echo "  tail -f ${RUNS_DIR}/*_hq_eval*.log"
fi
