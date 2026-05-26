#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

python -m grteclyn_wrapper \
  --name "${ATLAS_NAME:-atlas_gpu0}" \
  --cuda-devices "${CUDA_VISIBLE_DEVICES_OVERRIDE:-0}" \
  --set stop_time="${ATLAS_STOP_TIME:-0.04}" \
  --set N_full="${ATLAS_N_FULL:-32}" \
  --set max_level="${ATLAS_MAX_LEVEL:-0}" \
  atlas --count "${ATLAS_COUNT:-20}" --seed "${ATLAS_SEED:-1}"
