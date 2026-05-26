#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

if [[ -z "${GRTECLYN_ROOT:-}" ]]; then
  if [[ -f "${WRAPPER_ROOT}/../Examples/SupportedWormholeCollapse/params_2gpu.txt" ]]; then
    export GRTECLYN_ROOT="$(cd -- "${WRAPPER_ROOT}/.." && pwd)"
  else
    echo "Set GRTECLYN_ROOT=/path/to/GRTeclyn before running this script." >&2
    exit 2
  fi
fi

export PYTHONPATH="${WRAPPER_ROOT}/src:${PYTHONPATH:-}"
cd "${GRTECLYN_ROOT}"
