#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

python -m grteclyn_wrapper \
  --name flow_medium_gpu0 \
  --cuda-devices "${CUDA_VISIBLE_DEVICES_OVERRIDE:-0}" \
  --consume-plotfiles \
  --consumer-radii 8 16 \
  --consumer-delete \
  --set N_full=64 \
  --set max_level=2 \
  --set stop_time=2 \
  reproduce
