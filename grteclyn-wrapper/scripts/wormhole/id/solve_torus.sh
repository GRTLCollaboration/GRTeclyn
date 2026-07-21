#!/usr/bin/env bash
# Isolated rotating Q-torus (stationary m=1 eigenstate) initial data.
# Reuses the GRTresna BosonStarBH solver + convert_chombo_to_gridinit pipeline,
# painting the genuine 2D spinning Q-ball f(rho,z) as a profile==4 winding lump
# on a throat-free centered full box.  See solve_torus.py for env knobs.
#
# Usage:  solve_torus.sh [NRANKS]     (default NRANKS=2)
set -euo pipefail

WH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../lib/env.sh
source "${WH_DIR}/../../lib/env.sh"

if [[ -z "${GRTRESNA_ENV:-}" ]] || [[ ! -d "${GRTRESNA_ENV}/bin" ]]; then
  echo "Set GRTRESNA_ENV (see .env.example) before running solve_torus.sh." >&2
  exit 2
fi
export MPIRUN="${MPIRUN:-${GRTRESNA_ENV}/bin/mpirun}"
export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="${WRAPPER_ROOT}/src:${PYTHONPATH:-}"

PY="${GRTECLYN_PYTHON:-${WRAPPER_ROOT}/.venv/bin/python}"
if [[ ! -x "${PY}" ]]; then
  echo "error: wrapper venv python not found at ${PY}; set GRTECLYN_PYTHON" >&2
  exit 2
fi

exec "${PY}" "${WH_DIR}/solve_torus.py" "${1:-2}"
