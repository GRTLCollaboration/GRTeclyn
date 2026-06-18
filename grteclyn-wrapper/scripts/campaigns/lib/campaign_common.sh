#!/usr/bin/env bash
# Shared helpers for scripts/campaigns/run_full_campaign.sh
set -euo pipefail

ftl_campaign_require_name() {
  if [[ -z "${CAMPAIGN_NAME:-}" ]]; then
    echo "[campaign] set CAMPAIGN_NAME (e.g. general_ftl_wormhole_v22)" >&2
    exit 2
  fi
}

ftl_campaign_init_layout() {
  ftl_campaign_require_name
  CAMPAIGN_ROOT="${CAMPAIGN_ROOT:-${GRTECLYN_ROOT}/runs/campaigns/${CAMPAIGN_NAME}}"
  CAMPAIGN_LOG_DIR="${CAMPAIGN_ROOT}/logs"
  QD_DIR="${CAMPAIGN_ROOT}/qd"
  CMAES_DIR="${CAMPAIGN_ROOT}/cmaes"
  HQ_RUNS_DIR="${CAMPAIGN_ROOT}/promote"
  CAMPAIGN_STATE="${CAMPAIGN_ROOT}/campaign_state.json"
  mkdir -p "${CAMPAIGN_ROOT}" "${CAMPAIGN_LOG_DIR}" "${HQ_RUNS_DIR}"
}

ftl_campaign_stage_done() {
  local stage="$1"
  [[ -f "${CAMPAIGN_STATE}" ]] || return 1
  local flag
  flag="$(${PYTHON_BIN} - "${CAMPAIGN_STATE}" "${stage}" <<'PY'
import json
import sys
from pathlib import Path

state_path, stage = sys.argv[1:3]
state = json.loads(Path(state_path).read_text(encoding="utf-8"))
print("1" if state.get(stage, {}).get("status") == "complete" else "0")
PY
)"
  [[ "${flag}" == "1" ]]
}

ftl_campaign_should_run_stage() {
  local stage="$1"
  if [[ "${STAGE:-all}" != "all" && "${STAGE}" != "${stage}" ]]; then
    return 1
  fi
  if [[ "${RESUME:-0}" == "1" ]] && ftl_campaign_stage_done "${stage}"; then
    echo "[campaign] skip ${stage} (already complete; RESUME=1)" >&2
    return 1
  fi
  return 0
}

ftl_campaign_mark_stage() {
  local stage="$1"
  local extra_json="${2-"{}"}"
  CAMPAIGN_STAGE="${stage}" \
  CAMPAIGN_EXTRA="${extra_json}" \
  CAMPAIGN_STATE_PATH="${CAMPAIGN_STATE}" \
  ${PYTHON_BIN} - <<'PY'
import json
import os
from datetime import datetime, timezone
from pathlib import Path

stage = os.environ["CAMPAIGN_STAGE"]
extra_raw = os.environ.get("CAMPAIGN_EXTRA", "{}")
state_path = Path(os.environ["CAMPAIGN_STATE_PATH"])
extra = json.loads(extra_raw) if extra_raw else {}
state = json.loads(state_path.read_text(encoding="utf-8")) if state_path.is_file() else {}
state[stage] = {
    "status": "complete",
    "at": datetime.now(timezone.utc).isoformat(),
    **extra,
}
state_path.parent.mkdir(parents=True, exist_ok=True)
state_path.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")
PY
}

ftl_campaign_pick_top_eval() {
  local trajectory="$1"
  local extra_args=()
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    extra_args=(--include-dry-run)
  fi
  ${PYTHON_BIN} "${CAMPAIGNS_ROOT}/lib/pick_top_eval.py" "${trajectory}" "${extra_args[@]}"
}

ftl_campaign_pick_top_pairs() {
  local trajectory="$1"
  local gpu_id="$2"
  local extra_args=()
  if [[ "${DRY_RUN:-0}" == "1" ]]; then
    extra_args=(--include-dry-run)
  fi
  ${PYTHON_BIN} "${CAMPAIGNS_ROOT}/lib/pick_top_eval.py" \
    "${trajectory}" --format pairs --gpu "${gpu_id}" "${extra_args[@]}"
}

ftl_campaign_write_manifest() {
  ${PYTHON_BIN} - <<'PY' "${CAMPAIGN_ROOT}"
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

root = Path(sys.argv[1])
manifest = {
    "campaign_name": root.name,
    "root": str(root),
    "layout": {
        "qd": str(root / "qd"),
        "cmaes": str(root / "cmaes"),
        "promote": str(root / "promote"),
        "logs": str(root / "logs"),
    },
    "updated_at": datetime.now(timezone.utc).isoformat(),
}
if (root / "campaign_state.json").is_file():
    manifest["state"] = json.loads(
        (root / "campaign_state.json").read_text(encoding="utf-8")
    )
(root / "campaign.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
PY
}

ftl_campaign_first_gpu() {
  # shellcheck disable=SC2206
  local ids=(${GPU_IDS:-0})
  echo "${ids[0]}"
}
