#!/usr/bin/env bash
# Rebuild Chombo MPI libs + GRTresna ScalarFieldBH MPI executable (fixes SIGILL).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=env.sh
source "${SCRIPT_DIR}/env.sh"

GRTRESNA_ENV="${GRTRESNA_ENV:-${HOME}/.mlspace/envs/grtresna}"
CHOMBO_HOME="${CHOMBO_HOME:-${GRTECLYN_ROOT}/../Chombo/lib}"
GRTRESNA_ROOT="${GRTRESNA_ROOT:-${GRTECLYN_ROOT}/../GRTresna}"
JOBS="${JOBS:-8}"

if [[ ! -d "${GRTRESNA_ENV}/bin" ]]; then
  echo "GRTresna conda env not found: ${GRTRESNA_ENV}" >&2
  exit 1
fi

export PATH="${GRTRESNA_ENV}/bin:${PATH}"
export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="${GRTRESNA_ENV}"

echo "=== Rebuilding Chombo MPI libraries (DIM=3) ==="
cd "${CHOMBO_HOME}"
for lib in BaseTools BoxTools AMRTools AMRElliptic; do
  echo "--- ${lib} ---"
  make --no-print-directory -C "src/${lib}" realclean MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE || true
done
rm -f lib/lib*mpicxx*.a
make lib DIM=3 MPI=TRUE OPT=HIGH DEBUG=FALSE -j"${JOBS}"

echo "=== Rebuilding GRTresna ScalarFieldBH MPI executable ==="
cd "${GRTRESNA_ROOT}/Examples/ScalarFieldBH"
make realclean MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE 2>/dev/null || true
rm -f Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex
rm -rf o/3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI
PATH="${GRTRESNA_ENV}/bin:${PATH}" CONDA_PREFIX="${GRTRESNA_ENV}" \
  make all -j"${JOBS}" CHOMBO_HOME="${CHOMBO_HOME}" MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE

EXE="${GRTRESNA_ROOT}/Examples/ScalarFieldBH/Main_ScalarFieldBH3d.Linux.64.mpicxx.gfortran.OPTHIGH.MPI.ex"
echo "=== Smoke test: mpirun -np 2 (GRTresna is quiet on stdout; watch Ham_and_Mom_errors.txt) ==="
WORK="${GRTRESNA_ROOT}/Examples/ScalarFieldBH/_mpi_smoke_$$"
rm -rf "${WORK}"
mkdir -p "${WORK}/Outputs" "${WORK}/pout"
cp params_canonical_amr_test.txt "${WORK}/params.txt"
stdbuf -oL -eL mpirun -np 2 --oversubscribe "${EXE}" "${WORK}/params.txt" 2>&1 &
MPI_PID=$!
for _ in $(seq 1 60); do
  if ! kill -0 "${MPI_PID}" 2>/dev/null; then break; fi
  if [[ -f "${WORK}/Ham_and_Mom_errors.txt" ]]; then
    tail -1 "${WORK}/Ham_and_Mom_errors.txt"
  else
    echo "  ... waiting for Ham_and_Mom_errors.txt"
  fi
  sleep 5
done
wait "${MPI_PID}" 2>/dev/null || true
if [[ -f "${WORK}/Ham_and_Mom_errors.txt" ]]; then
  echo "Final Ham/Mom:"
  tail -3 "${WORK}/Ham_and_Mom_errors.txt"
else
  echo "No Ham_and_Mom_errors.txt (check pout/ or MPI abort above)" >&2
fi
rm -rf "${WORK}"
echo "Done."
