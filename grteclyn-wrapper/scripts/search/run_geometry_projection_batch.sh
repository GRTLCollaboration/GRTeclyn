#!/usr/bin/env bash
# Geometry-first projection campaign, MAP-Elites style.
#
# Mirrors the QD search philosophy (run_grtresna_qd_search.sh): generate MANY
# geometry-first scout candidates, then gate them through GRTresna projection.
# Only preservation survivors are evolved on GPU -- failures are recorded but
# never consume long replays.
#
# Pipeline:
#   1. SCOUT   geometry-first MAP-Elites campaign (RadialRecipe, no GRTresna)
#   2. RANK    score every scout episode, pick top-K
#   3. PROJECT each top-K scout through GRTresna (solve + motif preservation)
#   4. EVOLVE  short GPU replay only for preservation survivors
#
# Stages can be limited with STAGES="scout rank project evolve" (default all).
set -euo pipefail

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

export GRTRESNA_ROOT="${GRTRESNA_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)/GRTresna}"

NAME="${NAME:-geometry_projection_$(date -u +%Y%m%dT%H%M%SZ)}"
RUNS_DIR="${RUNS_DIR:-${GRTECLYN_ROOT}/runs/geometry_projection}"
PROJ_ROOT="${PROJ_ROOT:-${RUNS_DIR}/${NAME}}"
SCOUT_DIR="${SCOUT_DIR:-${PROJ_ROOT}/scout}"

STAGES="${STAGES:-scout rank project evolve}"

# Scout (geometry-first MAP-Elites) knobs.
SCOUT_MODE="${SCOUT_MODE:-qd}"          # qd or optimize
QD_ITERATIONS="${QD_ITERATIONS:-8}"
BINS="${BINS:-8}"
MAX_GENERATIONS="${MAX_GENERATIONS:-8}"
GPU_IDS="${GPU_IDS:-0 1 2 3}"
BATCH_SIZE="${BATCH_SIZE:-$(wc -w <<< "${GPU_IDS}")}"
SEED="${SEED:-11}"
STOP_TIME="${STOP_TIME:-2.0}"
PLOT_INTERVAL="${PLOT_INTERVAL:-10}"
FTL_L="${FTL_L:-8.0}"
NONSPHERICAL="${NONSPHERICAL:-1}"

# Ranking / projection knobs.
TOP_K="${TOP_K:-3}"
RANKS="${RANKS:-8}"
MAX_LUMPS="${MAX_LUMPS:-3}"
GRTRESNA_L="${GRTRESNA_L:-128.0}"
GRIDINIT_N="${GRIDINIT_N:-64}"
CUDA_DEVICE="${CUDA_DEVICE:-0}"

export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-phi Pi scalar_activity chi chi_minus_1 local_speed shift1 rho_req}"
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"

mkdir -p "${PROJ_ROOT}"

has_stage() { [[ " ${STAGES} " == *" $1 "* ]]; }

# ---------------------------------------------------------------------------
# Stage 1: scout
# ---------------------------------------------------------------------------
if has_stage scout; then
  echo "=== [scout] geometry-first ${SCOUT_MODE} campaign -> ${SCOUT_DIR} ==="
  NONSPHERICAL_ARGS=()
  if [[ "${NONSPHERICAL}" == "1" ]]; then
    NONSPHERICAL_ARGS=(--nonspherical)
  fi
  if [[ "${SCOUT_MODE}" == "qd" ]]; then
    # shellcheck disable=SC2086
    ${PYTHON_BIN} -m grteclyn_wrapper \
      --example RadialRecipe \
      --consume-plotfiles --consumer-delete --consumer-radii 4.0 8.0 \
      --no-phantom \
      --ftl-L "${FTL_L}" \
      --runs-dir "$(dirname -- "${SCOUT_DIR}")" \
      --name "$(basename -- "${SCOUT_DIR}")" \
      --set stop_time="${STOP_TIME}" \
      --set plot_interval="${PLOT_INTERVAL}" \
      qd \
        --descriptor-mode channel \
        --objective-mode ftl_first \
        --iterations "${QD_ITERATIONS}" \
        --batch-size "${BATCH_SIZE}" \
        --bins "${BINS}" \
        --seed "${SEED}" \
        --gpu-ids ${GPU_IDS} \
        "${NONSPHERICAL_ARGS[@]}"
  else
    # shellcheck disable=SC2086
    ${PYTHON_BIN} -m grteclyn_wrapper \
      --example RadialRecipe \
      --consume-plotfiles --consumer-delete --consumer-radii 4.0 8.0 \
      --no-phantom \
      --ftl-L "${FTL_L}" \
      --runs-dir "$(dirname -- "${SCOUT_DIR}")" \
      --name "$(basename -- "${SCOUT_DIR}")" \
      --set stop_time="${STOP_TIME}" \
      --set plot_interval="${PLOT_INTERVAL}" \
      optimize \
        --max-generations "${MAX_GENERATIONS}" \
        --gpu-ids ${GPU_IDS} \
        --seed "${SEED}" \
        "${NONSPHERICAL_ARGS[@]}"
  fi
fi

# ---------------------------------------------------------------------------
# Stage 2: rank
# ---------------------------------------------------------------------------
RANK_FILE="${PROJ_ROOT}/scout_ranking.txt"
if has_stage rank; then
  echo "=== [rank] scoring scouts in ${SCOUT_DIR} (top ${TOP_K}) ==="
  SCOUT_DIR="${SCOUT_DIR}" TOP_K="${TOP_K}" FTL_L="${FTL_L}" \
  STOP_TIME="${STOP_TIME}" RANK_FILE="${RANK_FILE}" \
  ${PYTHON_BIN} - <<'PY'
import os
from pathlib import Path
from grteclyn_wrapper.metrics import read_episode_metrics, score_episode

scout_dir = Path(os.environ["SCOUT_DIR"])
top_k = int(os.environ["TOP_K"])
ftl_l = float(os.environ["FTL_L"])
stop_time = float(os.environ["STOP_TIME"])
rank_file = Path(os.environ["RANK_FILE"])

rows = []
for ep in sorted(scout_dir.glob("eval_*")):
    if not (ep / "metadata.json").exists():
        continue
    try:
        m = read_episode_metrics(ep, ftl_L=ftl_l)
        s = score_episode(m, target_stop_time=stop_time)
        c = s.components
        rows.append((
            s.total,
            ep.name,
            c.get("ftl_shortcut", 0.0) or c.get("ftl_precursor", 0.0),
            c.get("operational_ftl", 0.0),
            c.get("channel_progress", 0.0),
        ))
    except Exception as exc:  # noqa: BLE001
        print(f"skip {ep.name}: {exc}")

rows.sort(reverse=True)
rank_file.parent.mkdir(parents=True, exist_ok=True)
with rank_file.open("w", encoding="utf-8") as fh:
    for total, name, shortcut, op, ch in rows[:top_k]:
        fh.write(f"{name}\n")
        print(f"{name}  score={total:.3f}  shortcut={shortcut:.3f} "
              f"op_ftl={op:.3f}  channel={ch:.3f}")

if not rows:
    print("No scout episodes with metadata found.")
PY
fi

# ---------------------------------------------------------------------------
# Stage 3: project (GRTresna solve + motif preservation), survivors only evolve
# ---------------------------------------------------------------------------
if has_stage project; then
  if [[ ! -f "${RANK_FILE}" ]]; then
    echo "[project] missing ${RANK_FILE}; run the rank stage first." >&2
    exit 2
  fi
  EVOLVE_FLAG=0
  if has_stage evolve; then
    EVOLVE_FLAG=1
  fi
  SURVIVORS="${PROJ_ROOT}/preservation_survivors.txt"
  : > "${SURVIVORS}"

  while IFS= read -r ev; do
    [[ -z "${ev}" ]] && continue
    SRC="${SCOUT_DIR}/${ev}"
    OUT="${PROJ_ROOT}/${ev}"
    echo "=== [project] ${ev} -> ${OUT} ==="

    # Solve + preservation only (CPU/MPI), skip post-load gate here.
    set +e
    MODE=solve-only SKIP_POSTLOAD_GATE=1 \
      SOURCE="${SRC}" OUT="${OUT}" \
      RANKS="${RANKS}" MAX_LUMPS="${MAX_LUMPS}" FTL_L="${FTL_L}" \
      GRTRESNA_L="${GRTRESNA_L}" GRIDINIT_N="${GRIDINIT_N}" \
      bash "${SCRIPT_DIR}/run_geometry_projection_eval.sh"
    rc=$?
    set -e

    PRES="${OUT}/projection_report_preservation.json"
    passed="false"
    if [[ -f "${PRES}" ]]; then
      passed=$(${PYTHON_BIN} -c "import json,sys; print(str(json.load(open('${PRES}'))['passed']).lower())" 2>/dev/null || echo "false")
    fi
    echo "[project] ${ev}: exit=${rc} preservation=${passed}"

    if [[ "${passed}" == "true" ]]; then
      echo "${ev}" >> "${SURVIVORS}"
      if [[ "${EVOLVE_FLAG}" == "1" ]]; then
        echo "=== [evolve] survivor ${ev} short GPU replay ==="
        GRIDINIT="${OUT}/initial_data.gridinit" \
          OUT="${OUT}" NAME="projection_evolve" \
          STOP_TIME="${STOP_TIME}" PLOT_INTERVAL="${PLOT_INTERVAL}" \
          CUDA_DEVICE="${CUDA_DEVICE}" CONSUMER_KEEP_LAST=1 FTL_L="${FTL_L}" \
          bash "${SCRIPT_DIR}/run_geometry_projection_replay_gridinit.sh"
      fi
    else
      echo "[project] ${ev} rejected by preservation; skipping GPU replay."
    fi
  done < "${RANK_FILE}"

  echo "=== [done] survivors written to ${SURVIVORS} ==="
  cat "${SURVIVORS}" || true
fi
