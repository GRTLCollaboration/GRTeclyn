#!/usr/bin/env bash
# Stage 1.5 — RL training launcher (opt-in).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

: "${RL_EPISODE_PATH:?set RL_EPISODE_PATH to eval run directory}"
: "${GRTeclyn_EXE:?set GRTeclyn_EXE to RadialRecipe binary}"
: "${RL_PARAMS:?set RL_PARAMS to params.txt with rl_enabled=1}"

uv run python scripts/campaigns/rl/train_motor.py \
  --executable "${GRTeclyn_EXE}" \
  --params "${RL_PARAMS}" \
  --episode-path "${RL_EPISODE_PATH}"
