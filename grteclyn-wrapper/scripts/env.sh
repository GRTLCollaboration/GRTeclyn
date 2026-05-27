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

# Optional local OpenMPI (multi-GPU builds only; smoke uses USE_MPI=FALSE).
OPENMPI_ROOT="${OPENMPI_ROOT:-${GRTECLYN_ROOT}/../local/openmpi-5.0.8}"
if [[ -d "${OPENMPI_ROOT}/bin" ]]; then
  export PATH="${OPENMPI_ROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
fi

cd "${GRTECLYN_ROOT}"
