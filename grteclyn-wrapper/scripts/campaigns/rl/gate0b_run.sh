#!/usr/bin/env bash
# Gate 0B — ZMQ bridge smoke: rl_enabled=1 + dummy agent completes without hang.
set -euo pipefail

RL_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${RL_SCRIPT_DIR}/../../lib/env.sh"

GRTECLYN="${GRTECLYN_ROOT}"
EXE="${GRTeclyn_EXE:-${GRTECLYN}/Examples/RadialRecipe/main3d.gnu.CUDA.ex}"
ZMQ_PREFIX="${GRTECLYN}/local/zeromq"
ELITE_EVAL="${ELITE_EVAL:-${GRTECLYN}/runs/grtresna_qd/spacetime_splash_v14_moving/eval_000010}"
RUN_ROOT="${RUN_ROOT:-${GRTECLYN}/runs/rl_gate0b/spacetime_splash_v14_eval010}"
STOP_TIME="${GATE0B_STOP_TIME:-1.0}"
RL_INTERVAL="${GATE0B_RL_INTERVAL:-32}"
RL_PORT="${RL_PORT:-5556}"
GPU_ID="${CUDA_VISIBLE_DEVICES:-0}"
MIN_EXCHANGES="${MIN_EXCHANGES:-2}"

if [[ ! -x "${EXE}" ]]; then
  echo "Missing RL binary: ${EXE}" >&2
  exit 2
fi

pkill -f "dummy_agent.py --port ${RL_PORT}" 2>/dev/null || true
sleep 1
mkdir -p "${RUN_ROOT}/data"

python3 - "${ELITE_EVAL}/params.txt" "${RUN_ROOT}" "${STOP_TIME}" "${RL_INTERVAL}" "${RL_PORT}" "${ELITE_EVAL}" <<'PY'
import sys
from pathlib import Path

src, out_dir, stop_time, interval, port, elite = sys.argv[1:7]
lines = Path(src).read_text(encoding="utf-8").splitlines()
out = Path(out_dir)
elite_dir = Path(elite)

overrides = {
    "output_path": str(out.resolve()),
    "amr.check_file": str((out / "RadialRecipeChk").resolve()),
    "amr.plot_file": str((out / "RadialRecipePlt").resolve()),
    "checkpoint_interval": "-1",
    "plot_interval": "-1",
    "stop_time": stop_time,
    "recipe_initial_data_file": str((elite_dir / "initial_data.gridinit").resolve()),
    "rl_enabled": "1",
    "rl_coarse_step_interval": interval,
    "rl_zmq_port": port,
    "rl_num_lumps": "1",
    "rl_lump_seed_x": "0.0",
    "rl_lump_seed_y": "0.0",
    "rl_lump_seed_z": "0.0",
    "rl_pump_width": "8.0",
    "rl_pump_max_amplitude": "0.05",
}
keys = set(overrides)
rendered = []
seen = set()
for line in lines:
    if "=" in line and not line.lstrip().startswith("#"):
        key = line.split("=", 1)[0].strip()
        if key in overrides:
            rendered.append(f"{key} = {overrides[key]}")
            seen.add(key)
            continue
    rendered.append(line)
for key in sorted(keys - seen):
    rendered.append(f"{key} = {overrides[key]}")
(out / "params.txt").write_text("\n".join(rendered) + "\n", encoding="utf-8")
PY

PARAMS="${RUN_ROOT}/params.txt"
DUMMY_STEPS=100

export LD_LIBRARY_PATH="${ZMQ_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES="${GPU_ID}"
export PYTHONUNBUFFERED=1

echo "== Gate 0B: RL bridge smoke (stop_time=${STOP_TIME}, interval=${RL_INTERVAL}) =="
rm -f "${RUN_ROOT}/data/"*.dat

(
  cd "${WRAPPER_ROOT}"
  uv run python "${RL_SCRIPT_DIR}/dummy_agent.py" --port "${RL_PORT}" --steps "${DUMMY_STEPS}"
) > "${RUN_ROOT}/dummy_agent.log" 2>&1 &
DUMMY_PID=$!
sleep 2
if ! kill -0 "${DUMMY_PID}" 2>/dev/null; then
  echo "FAIL dummy_agent failed to start" >&2
  cat "${RUN_ROOT}/dummy_agent.log" >&2
  exit 2
fi

set +e
"${EXE}" "${PARAMS}" > "${RUN_ROOT}/run.log" 2>&1
SIM_RC=$?
set -e
kill "${DUMMY_PID}" 2>/dev/null || true
wait "${DUMMY_PID}" 2>/dev/null || true

if [[ "${SIM_RC}" -ne 0 ]]; then
  echo "FAIL simulation exit code ${SIM_RC}" >&2
  tail -30 "${RUN_ROOT}/run.log" >&2
  exit "${SIM_RC}"
fi

if ! grep -q "AMReX.*finalized" "${RUN_ROOT}/run.log"; then
  echo "FAIL simulation did not finalize cleanly" >&2
  exit 2
fi

EXCHANGES="$(grep -c '^step=' "${RUN_ROOT}/dummy_agent.log" || true)"
if [[ "${EXCHANGES}" -lt "${MIN_EXCHANGES}" ]]; then
  echo "FAIL expected >= ${MIN_EXCHANGES} ZMQ exchanges, got ${EXCHANGES}" >&2
  cat "${RUN_ROOT}/dummy_agent.log" >&2
  exit 2
fi

OBS_DIM="$(grep '^step=' "${RUN_ROOT}/dummy_agent.log" | head -1 | sed -n 's/.*obs_dim=\([0-9]*\).*/\1/p')"
if [[ "${OBS_DIM}" != "14" ]]; then
  echo "FAIL expected obs_dim=14 for 6+8*1 lump, got ${OBS_DIM:-missing}" >&2
  exit 2
fi

echo "PASS Gate 0B: sim finalized, ${EXCHANGES} ZMQ exchanges, obs_dim=${OBS_DIM}"
