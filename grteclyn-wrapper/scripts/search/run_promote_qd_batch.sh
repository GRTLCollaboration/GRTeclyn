#!/usr/bin/env bash
# Parallel HQ promotion of QD elites via fresh GRTresna solves + framed GPU evolution.
#
# REQUIRED: in-flight frame extraction + plotfile deletion (see README.md
# "ALWAYS extract frames on the fly"). replay_grtresna_eval.py enables
# consume_plotfiles; GRTECLYN_FRAMES_* from promote_common.sh must stay set.
#
# Usage:
#   bash scripts/search/run_promote_qd_batch.sh
#
# Candidate selection (one of):
#   CANDIDATES="024 0 055 1 025 2"          # eval_id gpu_id pairs
#   CANDIDATES_FILE=/path/to/candidates.txt # one "eval gpu" pair per line
#   TOP_K=8 MIN_OPERATIONAL_FTL=0.03        # auto-pick from trajectory.jsonl
#
# Env overrides:
#   QD_RUN=.../qd_ftl_discovery_20260609T162553Z
#   NAME_PREFIX=ftl_discovery               # -> l128n256t30_ftl_discovery_qd_eval000055
#   N_FULL=256 L_FULL=128 STOP_TIME=30
#   GRTRESNA_MAX_HAM_PCT=10 GRTRESNA_MAX_MOM_PCT=10
#   FORCE=1 DRY_RUN=1
set -euo pipefail

SEARCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SEARCH_DIR}/../lib/env.sh"
# shellcheck source=../lib/promote_common.sh
source "${SEARCH_DIR}/../lib/promote_common.sh"

QD_RUN="${QD_RUN:-${GRTECLYN_ROOT}/runs/grtresna_qd/ftl_discovery/qd_ftl_discovery_20260609T162553Z}"
NAME_PREFIX="${NAME_PREFIX:-ftl_discovery}"
CANDIDATES="${CANDIDATES:-}"
CANDIDATES_FILE="${CANDIDATES_FILE:-}"
TOP_K="${TOP_K:-0}"
MIN_OPERATIONAL_FTL="${MIN_OPERATIONAL_FTL:-0.03}"
GPU_IDS="${GPU_IDS:-0 1 2 3 4 5 6 7}"
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "${RUNS_DIR}"

declare -a CANDIDATE_ENTRIES=()

if [[ -n "${CANDIDATES}" ]]; then
  # shellcheck disable=SC2206
  pairs=(${CANDIDATES})
  if ((${#pairs[@]} % 2 != 0)); then
    echo "[promote] CANDIDATES must contain eval/gpu pairs (even count), got ${#pairs[@]} tokens" >&2
    exit 2
  fi
  for ((i = 0; i < ${#pairs[@]}; i += 2)); do
    CANDIDATE_ENTRIES+=("${pairs[i]} ${pairs[i + 1]}")
  done
elif [[ -n "${CANDIDATES_FILE}" ]]; then
  if [[ ! -f "${CANDIDATES_FILE}" ]]; then
    echo "[promote] missing CANDIDATES_FILE: ${CANDIDATES_FILE}" >&2
    exit 2
  fi
  while read -r eval_id gpu_id _rest; do
    [[ -z "${eval_id}" || "${eval_id}" == \#* ]] && continue
    CANDIDATE_ENTRIES+=("${eval_id} ${gpu_id}")
  done < "${CANDIDATES_FILE}"
elif [[ "${TOP_K}" -gt 0 ]]; then
  traj="${QD_RUN}/trajectory.jsonl"
  if [[ ! -f "${traj}" ]]; then
    echo "[promote] missing trajectory for TOP_K pick: ${traj}" >&2
    exit 2
  fi
  mapfile -t picked < <(python3 - <<'PY' "${traj}" "${TOP_K}" "${MIN_OPERATIONAL_FTL}" "${GPU_IDS}"
import json
import sys

traj_path, top_k, min_op, gpu_ids_raw = sys.argv[1:5]
top_k = int(top_k)
min_op = float(min_op)
gpus = [g for g in gpu_ids_raw.split() if g]

rows = []
with open(traj_path, encoding="utf-8") as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get("grtresna_rejected") or rec.get("solved_ftl_rejected") or rec.get("grtresna_failed"):
            continue
        if rec.get("preflight_rejected"):
            continue
        comps = rec.get("components") or {}
        op = float(comps.get("operational_ftl") or 0.0)
        if op < min_op:
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
    echo "[promote] TOP_K=${TOP_K} found no candidates (min operational_ftl=${MIN_OPERATIONAL_FTL})" >&2
    exit 2
  fi
  CANDIDATE_ENTRIES=("${picked[@]}")
else
  echo "[promote] set CANDIDATES, CANDIDATES_FILE, or TOP_K>0" >&2
  exit 2
fi

echo "== QD HQ promotion batch =="
echo "QD source     : ${QD_RUN}"
echo "Runs dir      : ${RUNS_DIR}"
echo "Name prefix   : ${NAME_PREFIX:-<none>}"
echo "Domain        : L=${L_FULL} N=${N_FULL} max_level=${MAX_LEVEL} t=${STOP_TIME}"
echo "GRTresna gate : Ham<=${GRTRESNA_MAX_HAM_PCT}% Mom<=${GRTRESNA_MAX_MOM_PCT}% iterations=${GRTRESNA_ITERATIONS}"
echo "Plotting      : frames on the fly, delete plotfiles, keep_last=${CONSUMER_KEEP_LAST}"
echo "Frames fields : ${GRTECLYN_FRAMES_FIELDS}"
echo "Projections   : ${GRTECLYN_PROJECTION_FIELDS} axes=${GRTECLYN_PROJECTION_AXES}"
echo "Candidates    : ${#CANDIDATE_ENTRIES[@]}"
echo

launched=0
for entry in "${CANDIDATE_ENTRIES[@]}"; do
  read -r EVAL_ID GPU_ID <<< "${entry}"
  eval_num="$((10#${EVAL_ID}))"
  eval_padded="$(printf "%06d" "${eval_num}")"
  name="$(promote_build_name "${eval_padded}")"
  source="${QD_RUN}/eval_${eval_padded}"
  out="${RUNS_DIR}/${name}"
  log="${RUNS_DIR}/${name}.log"

  if [[ ! -d "${source}" ]]; then
    echo "[skip] missing source ${source}"
    continue
  fi
  if [[ -e "${out}" && "${FORCE}" != "1" ]]; then
    echo "[skip] output exists ${out} (set FORCE=1 to replace manually after cleanup)"
    continue
  fi

  echo "[launch] ${name} GPU=${GPU_ID} source=eval_${eval_padded} log=${log}"
  if [[ "${DRY_RUN}" == "1" ]]; then
    echo "  dry-run: would nohup replay_grtresna_eval.py ${source} --name ${name} --gpu ${GPU_ID}"
    launched=$((launched + 1))
    continue
  fi

  gridinit_args=()
  if [[ -n "${GRIDINIT:-}" ]]; then
    gridinit_args=(--gridinit "${GRIDINIT}")
  fi

  # shellcheck disable=SC2086
  nohup ${PYTHON_BIN} "${SEARCH_DIR}/replay_grtresna_eval.py" \
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
    --evolving-geodesic \
    "${gridinit_args[@]}" \
    > "${log}" 2>&1 &
  echo "  pid=$!"
  launched=$((launched + 1))
done

echo
if [[ "${DRY_RUN}" == "1" ]]; then
  echo "Dry-run: would launch ${launched} promotions."
else
  echo "Launched ${launched} promotions (background). Monitor:"
  echo "  tail -f ${RUNS_DIR}/l${L_FULL}n${N_FULL}*qd_eval*.log"
fi
