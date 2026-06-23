#!/usr/bin/env bash
# Gate 1 — Tax Man T1–T4 + short live SpacetimeFtlEnv episode (plot consumer on).
set -euo pipefail

RL_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${RL_SCRIPT_DIR}/../../lib/env.sh"

GRTECLYN="${GRTECLYN_ROOT}"
EXE="${GRTeclyn_EXE:-${GRTECLYN}/Examples/RadialRecipe/main3d.gnu.CUDA.ex}"
ZMQ_PREFIX="${GRTECLYN}/local/zeromq"
ELITE_EVAL="${ELITE_EVAL:-${GRTECLYN}/runs/grtresna_qd/spacetime_splash_v14_moving/eval_000010}"
RUN_ROOT="${RUN_ROOT:-${GRTECLYN}/runs/rl_gate1/spacetime_splash_v14_eval010}"
GATE1_STOP_TIME="${GATE1_STOP_TIME:-1.0}"
GATE1_PLOT_INTERVAL="${GATE1_PLOT_INTERVAL:-10}"
RL_PORT="${GATE1_ZMQ_PORT:-5557}"
GPU_ID="${CUDA_VISIBLE_DEVICES:-0}"

if [[ ! -x "${EXE}" ]]; then
  echo "Missing RL binary: ${EXE}" >&2
  exit 2
fi

pkill -f "dummy_agent.py --port ${RL_PORT}" 2>/dev/null || true
mkdir -p "${RUN_ROOT}/data" "${RUN_ROOT}/small_data"

python3 - "${ELITE_EVAL}/params.txt" "${RUN_ROOT}" "${GATE1_STOP_TIME}" "${GATE1_PLOT_INTERVAL}" "${RL_PORT}" "${ELITE_EVAL}" <<'PY'
import sys
from pathlib import Path

src, out_dir, stop_time, plot_interval, port, elite = sys.argv[1:7]
lines = Path(src).read_text(encoding="utf-8").splitlines()
out = Path(out_dir)
elite_dir = Path(elite)

overrides = {
    "output_path": str(out.resolve()),
    "amr.check_file": str((out / "RadialRecipeChk").resolve()),
    "amr.plot_file": str((out / "RadialRecipePlt").resolve()),
    "checkpoint_interval": "-1",
    "plot_interval": plot_interval,
    "stop_time": stop_time,
    "recipe_initial_data_file": str((elite_dir / "initial_data.gridinit").resolve()),
    "rl_enabled": "1",
    "rl_coarse_step_interval": "4",
    "rl_zmq_port": port,
    "rl_num_lumps": "1",
    "rl_lump_seed_x": "0.0",
    "rl_lump_seed_y": "0.0",
    "rl_lump_seed_z": "0.0",
    "rl_pump_width": "8.0",
    "rl_pump_max_amplitude": "0.05",
    "rl_zmq_timeout_ms": "300000",
    "N_full": "32",
    "L_full": "16.0",
    "max_level": "0",
    "center": "8 8 8",
    "regrid_interval": "0",
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

export LD_LIBRARY_PATH="${ZMQ_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export CUDA_VISIBLE_DEVICES="${GPU_ID}"

echo "== Gate 1: static T1–T4 + live env smoke =="
cd "${WRAPPER_ROOT}"
uv run python "${RL_SCRIPT_DIR}/gate1_validate.py" \
  --executable "${EXE}" \
  --params "${RUN_ROOT}/params.txt" \
  --episode-path "${RUN_ROOT}" \
  --stop-time "${GATE1_STOP_TIME}" \
  --zmq-port "${RL_PORT}" \
  --gpu-id "${GPU_ID}" \
  --zmq-lib "${ZMQ_PREFIX}/lib" \
  --max-steps 8
