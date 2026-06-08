#!/usr/bin/env bash
# GPU-only diagnostic replay from an existing projected .gridinit.
#
# Usage:
#   GRIDINIT=../runs/geometry_projection/eval_000123/initial_data.gridinit \
#     bash scripts/search/run_geometry_projection_replay_gridinit.sh
set -euo pipefail

CALLER_CWD="$(pwd)"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../lib/env.sh"

PYTHON_BIN="${PYTHON_BIN:-}"
if [[ -z "${PYTHON_BIN}" ]]; then
  if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
    PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
  elif command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
    PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
  else
    PYTHON_BIN="python"
  fi
fi

GRIDINIT="${1:-${GRIDINIT:-}}"
if [[ -z "${GRIDINIT}" ]]; then
  echo "Set GRIDINIT=/path/to/initial_data.gridinit or pass it as the first argument." >&2
  exit 2
fi

if [[ "${GRIDINIT}" = /* ]]; then
  GRIDINIT_ABS="${GRIDINIT}"
else
  GRIDINIT_ABS="$(cd -- "${CALLER_CWD}/$(dirname -- "${GRIDINIT}")" && pwd)/$(basename -- "${GRIDINIT}")"
fi
OUT="${OUT:-$(dirname -- "${GRIDINIT_ABS}")}"
if [[ "${OUT}" != /* ]]; then
  OUT="${CALLER_CWD}/${OUT}"
fi
NAME="${NAME:-projection_evolve}"
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-10}"
CUDA_DEVICE="${CUDA_DEVICE:-0}"
CONSUMER_KEEP_LAST="${CONSUMER_KEEP_LAST:-1}"
FTL_L="${FTL_L:-8.0}"

export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"
export CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-${CUDA_DEVICE}}"
export GRIDINIT_ABS OUT NAME STOP_TIME PLOT_INTERVAL CUDA_DEVICE CONSUMER_KEEP_LAST FTL_L

${PYTHON_BIN} - <<'PY'
import json
import os
from pathlib import Path

from grteclyn_wrapper.core.config import (
    DEFAULT_RADIAL_RECIPE_TEMPLATE,
    resolve_example,
    resolve_executable,
)
from grteclyn_wrapper.core.episode import write_json
from grteclyn_wrapper.core.evaluation import evaluate_overrides

gridinit = Path(os.environ["GRIDINIT_ABS"]).resolve()
out = Path(os.environ["OUT"]).resolve()
name = os.environ["NAME"]
stop_time = float(os.environ["STOP_TIME"])
plot_interval = int(os.environ["PLOT_INTERVAL"])
cuda_device = os.environ.get("CUDA_DEVICE") or None
consumer_keep_last = int(os.environ["CONSUMER_KEEP_LAST"])
ftl_l = float(os.environ["FTL_L"])

example = resolve_example("RadialRecipe")
executable = resolve_executable(example=example)
evaluation = evaluate_overrides(
    {
        "recipe_initial_data_file": str(gridinit),
        "stop_time": stop_time,
        "plot_interval": plot_interval,
    },
    out_dir=out,
    name=name,
    example=example,
    template=DEFAULT_RADIAL_RECIPE_TEMPLATE,
    executable=executable,
    constrained=False,
    phantom=False,
    use_preflight=False,
    cuda_devices=cuda_device,
    dry_run=False,
    target_stop_time=stop_time,
    ftl_L=ftl_l,
    grtresna=False,
    consume_plotfiles=True,
    consumer_radii=(4.0, 8.0),
    consumer_keep_last=consumer_keep_last,
)
score_path = out / f"{name}_score.json"
write_json(
    score_path,
    {
        "score": evaluation.score,
        "components": evaluation.components,
        "episode_path": evaluation.episode_path,
        "exit_code": evaluation.exit_code,
    },
)
print(json.dumps(json.loads(score_path.read_text(encoding="utf-8")), indent=2))
raise SystemExit(int(evaluation.exit_code or 0))
PY
