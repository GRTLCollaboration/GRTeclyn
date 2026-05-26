#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

python -m grteclyn_wrapper \
  --name full_single_gpu0 \
  --cuda-devices "${CUDA_VISIBLE_DEVICES_OVERRIDE:-0}" \
  --consume-plotfiles \
  --consumer-radii 8 16 \
  --consumer-delete \
  reproduce
