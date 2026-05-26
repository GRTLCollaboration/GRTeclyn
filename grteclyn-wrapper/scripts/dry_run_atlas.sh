#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/env.sh"

python -m grteclyn_wrapper \
  --dry-run \
  --name "${ATLAS_NAME:-atlas_dry_run}" \
  atlas --count "${ATLAS_COUNT:-3}" --seed "${ATLAS_SEED:-1}"
