#!/usr/bin/env bash
# Gate 0A — boson chassis IVP neutrality: rl_enabled=0 vs rl_enabled=1 + neutral agent.
set -euo pipefail

RL_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${RL_SCRIPT_DIR}/../../lib/env.sh"

GRTECLYN="${GRTECLYN_ROOT}"
EXE="${GRTeclyn_EXE:-${GRTECLYN}/Examples/RadialRecipe/main3d.gnu.CUDA.ex}"
ZMQ_PREFIX="${GRTECLYN}/local/zeromq"
ELITE_EVAL="${ELITE_EVAL:-${GRTECLYN}/runs/grtresna_qd/spacetime_splash_v14_moving/eval_000010}"
RUN_ROOT="${RUN_ROOT:-${GRTECLYN}/runs/rl_gate0a/spacetime_splash_v14_eval010}"
STOP_TIME="${STOP_TIME:-2.0}"
RL_PORT="${RL_PORT:-5555}"
GPU_ID="${CUDA_VISIBLE_DEVICES:-0}"

if [[ ! -x "${EXE}" ]]; then
  echo "Missing RL binary: ${EXE}" >&2
  echo "Build: cd Examples/RadialRecipe && make USE_RL=TRUE USE_CUDA=TRUE USE_MPI=FALSE CUDA_ARCH=90 NVCC_CCBIN=/usr/bin/g++" >&2
  exit 2
fi

mkdir -p "${RUN_ROOT}/baseline/data" "${RUN_ROOT}/rl_neutral/data"

write_params() {
  local out_dir="$1"
  local rl_enabled="$2"
  local src="${ELITE_EVAL}/params.txt"
  python3 - "${src}" "${out_dir}" "${rl_enabled}" "${STOP_TIME}" "${ELITE_EVAL}" <<'PY'
import sys
from pathlib import Path

src, out_dir, rl_enabled, stop_time, elite = sys.argv[1:6]
lines = Path(src).read_text(encoding="utf-8").splitlines()
out = Path(out_dir)
out.mkdir(parents=True, exist_ok=True)
gridinit = Path(elite) / "initial_data.gridinit"

overrides = {
    "output_path": str(out.resolve()),
    "amr.check_file": str((out / "RadialRecipeChk").resolve()),
    "amr.plot_file": str((out / "RadialRecipePlt").resolve()),
    "checkpoint_interval": "-1",
    "plot_interval": "-1",
    "stop_time": stop_time,
    "recipe_initial_data_file": str(gridinit.resolve()),
    "rl_enabled": "1" if rl_enabled == "1" else "0",
    "rl_coarse_step_interval": "64",
    "rl_zmq_port": "5555",
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
print(out / "params.txt")
PY
}

BASE_PARAMS="$(write_params "${RUN_ROOT}/baseline" 0)"
RL_PARAMS="$(write_params "${RUN_ROOT}/rl_neutral" 1)"

export LD_LIBRARY_PATH="${ZMQ_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES="${GPU_ID}"

echo "== Gate 0A: baseline (rl_enabled=0) =="
rm -f "${RUN_ROOT}/baseline/data/"*.dat
"${EXE}" "${BASE_PARAMS}" > "${RUN_ROOT}/baseline/run.log" 2>&1

echo "== Gate 0A: RL neutral (rl_enabled=1 + dummy agent) =="
rm -f "${RUN_ROOT}/rl_neutral/data/"*.dat
DUMMY_STEPS=500
(
  cd "${WRAPPER_ROOT}"
  uv run python "${RL_SCRIPT_DIR}/dummy_agent.py" --port "${RL_PORT}" --steps "${DUMMY_STEPS}"
) > "${RUN_ROOT}/rl_neutral/dummy_agent.log" 2>&1 &
DUMMY_PID=$!
sleep 2
if ! kill -0 "${DUMMY_PID}" 2>/dev/null; then
  echo "dummy_agent failed to start:" >&2
  cat "${RUN_ROOT}/rl_neutral/dummy_agent.log" >&2
  exit 2
fi
set +e
"${EXE}" "${RL_PARAMS}" > "${RUN_ROOT}/rl_neutral/run.log" 2>&1
SIM_RC=$?
set -e
kill "${DUMMY_PID}" 2>/dev/null || true
wait "${DUMMY_PID}" 2>/dev/null || true
if [[ "${SIM_RC}" -ne 0 ]]; then
  echo "RL simulation failed (rc=${SIM_RC}); tail run.log:" >&2
  tail -40 "${RUN_ROOT}/rl_neutral/run.log" >&2
  exit "${SIM_RC}"
fi

echo "== Gate 0A comparison =="
python3 "${RL_SCRIPT_DIR}/gate0a_compare.py" \
  "${RUN_ROOT}/baseline" \
  "${RUN_ROOT}/rl_neutral" \
  --rtol 1e-10
