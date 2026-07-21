#!/usr/bin/env bash
# Amplitude-reduction (kappa) initial-data family for the rotating wormhole
# (Phase B4/B5).  Reuses the validated GRTresna BosonStarBH solver + the shared
# convert_chombo_to_gridinit pipeline (same machinery the QD campaign uses).
#
# Usage:
#   solve_kappa_family.sh [KAPPAS] [NRANKS]
#     KAPPAS  comma-separated (default: 1.0,0.9,0.7,0.5)
#     NRANKS  MPI ranks for each solve (default: 2)
set -euo pipefail

WH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../../lib/env.sh
source "${WH_DIR}/../../lib/env.sh"

# The GRTresna env (set via .env / GRTRESNA_ENV) supplies mpirun + runtime libs;
# the wrapper's own .venv supplies the conversion deps (scipy/h5py/numpy).
if [[ -z "${GRTRESNA_ENV:-}" ]] || [[ ! -d "${GRTRESNA_ENV}/bin" ]]; then
  echo "Set GRTRESNA_ENV (see .env.example) before running solve_kappa_family.sh." >&2
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

exec "${PY}" "${WH_DIR}/solve_kappa_family.py" "${1:-1.0,0.9,0.7,0.5}" "${2:-2}"
