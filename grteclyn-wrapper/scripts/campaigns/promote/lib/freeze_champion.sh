#!/usr/bin/env bash
# Generic: snapshot a search eval so promote/matrix survives keep_top prune.
#
# Required env:
#   SOURCE_RUN   — QD/CMA campaign dir with trajectory.jsonl + eval_*
#   FREEZE_ROOT  — destination tree (writes CHAMPION.json + eval_XXXXXX/)
# Optional:
#   SOURCE_EVAL_ID — pin eval; default = best gpu_ok score
set -euo pipefail

PROMOTE_LIB="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
CAMPAIGNS_ROOT="$(cd -- "${PROMOTE_LIB}/../.." && pwd)"
SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
# shellcheck source=../../../lib/env.sh
source "${SCRIPTS_ROOT}/lib/env.sh"

SOURCE_RUN="${SOURCE_RUN:-}"
FREEZE_ROOT="${FREEZE_ROOT:-}"
SOURCE_EVAL_ID="${SOURCE_EVAL_ID:-}"

if [[ -z "${SOURCE_RUN}" || -z "${FREEZE_ROOT}" ]]; then
  echo "[freeze] set SOURCE_RUN and FREEZE_ROOT" >&2
  exit 2
fi
if [[ "${SOURCE_RUN}" != /* ]]; then
  SOURCE_RUN="${GRTECLYN_ROOT}/${SOURCE_RUN}"
fi
if [[ "${FREEZE_ROOT}" != /* ]]; then
  FREEZE_ROOT="${GRTECLYN_ROOT}/${FREEZE_ROOT}"
fi
if [[ ! -d "${SOURCE_RUN}" ]]; then
  echo "[freeze] missing SOURCE_RUN: ${SOURCE_RUN}" >&2
  exit 2
fi

eval "$(
  python3 - <<'PY' "${SOURCE_RUN}" "${SOURCE_EVAL_ID}"
import json, shlex, sys
from pathlib import Path

src = Path(sys.argv[1])
forced = sys.argv[2].strip()
traj = src / "trajectory.jsonl"
if forced:
    eval_id = int(forced)
    score = None
    if traj.is_file():
        for line in traj.open(encoding="utf-8"):
            rec = json.loads(line)
            if int(rec.get("eval", -1)) == eval_id:
                score = rec.get("score")
                break
else:
    if not traj.is_file():
        raise SystemExit(f"need trajectory or SOURCE_EVAL_ID: {traj}")
    best = None
    for line in traj.open(encoding="utf-8"):
        line = line.strip()
        if not line:
            continue
        rec = json.loads(line)
        if rec.get("status") != "gpu_ok":
            continue
        score = float(rec.get("score") or -1e300)
        if best is None or score > best[0]:
            best = (score, int(rec["eval"]))
    if best is None:
        raise SystemExit("no gpu_ok in trajectory")
    score, eval_id = best

eval_dir = src / f"eval_{eval_id:06d}"
print(f"EVAL_ID={eval_id}")
print(f"SCORE={shlex.quote('' if score is None else str(score))}")
print(f"EVAL_DIR={shlex.quote(str(eval_dir))}")
PY
)"

if [[ ! -d "${EVAL_DIR}" ]]; then
  echo "[freeze] champion eval dir missing (pruned?): ${EVAL_DIR}" >&2
  exit 4
fi
if [[ ! -f "${EVAL_DIR}/metadata.json" ]]; then
  echo "[freeze] missing metadata.json in ${EVAL_DIR}" >&2
  exit 4
fi

DEST="${FREEZE_ROOT}/eval_$(printf '%06d' "${EVAL_ID}")"
mkdir -p "${DEST}"

copy_if() {
  local f="$1"
  if [[ -e "${EVAL_DIR}/${f}" ]]; then
    cp -a "${EVAL_DIR}/${f}" "${DEST}/"
    echo "  + ${f}"
  fi
}

echo "[freeze] ${EVAL_DIR} -> ${DEST}  (score=${SCORE:-n/a})"
copy_if metadata.json
copy_if initial_data.matter.json
copy_if params.txt
copy_if score.json
copy_if postload_gate

python3 - <<PY
import json, pathlib
rec = {
  "source_run": "${SOURCE_RUN}",
  "eval_id": int("${EVAL_ID}"),
  "score": None if "${SCORE}" == "" else float("${SCORE}"),
  "frozen_eval": "${DEST}",
  "matter_model": "${GRTRESNA_MATTER_MODEL:-}",
  "rl_pump_stop_time": float("${RL_PUMP_STOP_TIME:-4}"),
}
path = pathlib.Path("${FREEZE_ROOT}/CHAMPION.json")
path.write_text(json.dumps(rec, indent=2) + "\n")
print(f"[freeze] wrote {path}")
PY
