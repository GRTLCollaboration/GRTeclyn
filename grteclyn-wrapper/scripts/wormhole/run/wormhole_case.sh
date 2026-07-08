#!/usr/bin/env bash
# CLI launcher for a single rotating-wormhole evolution case (no per-case param
# files -- params are generated from a template by wormhole_case.py).
#
#   wormhole_case.sh --kappa 1.0 --dx 0.5 --max-level 3 --stop-time 8 --gpu 0
#   wormhole_case.sh --kappa 0.5 --dx 0.5 --gpu 2 --no-frames
#
# The .gridinit for the case must already exist (produce it with
# solve_kappa_family.sh at the matching resolution: RES_N=128 -> dx=0.5).
set -euo pipefail

WH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${WH_DIR}/../../lib/env.sh"

PY="${GRTECLYN_PYTHON:-${WRAPPER_ROOT}/.venv/bin/python}"
if [[ ! -x "${PY}" ]]; then
  echo "error: wrapper venv python not found at ${PY}; set GRTECLYN_PYTHON" >&2
  exit 2
fi

exec "${PY}" "${WH_DIR}/wormhole_case.py" "$@"
