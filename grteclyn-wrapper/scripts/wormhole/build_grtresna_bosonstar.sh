#!/usr/bin/env bash
# Rebuild the GRTresna BosonStarBH initial-data solver (CTTKHybrid<ComplexScalarField>),
# used to generate rotating-wormhole constraint-clean initial data.
set -euo pipefail

GRTRESNA_ROOT="/home/jovyan/nachevsky/test/simulation/GRTresna"
export CHOMBO_HOME="${CHOMBO_HOME:-/home/jovyan/nachevsky/test/simulation/Chombo/lib}"
GRTRESNA_ENV="${GRTRESNA_ENV:-$HOME/.mlspace/envs/grtresna}"
if [[ -d "${GRTRESNA_ENV}/bin" ]]; then
  export PATH="${GRTRESNA_ENV}/bin:${PATH}"
  export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
fi

echo "[build] CHOMBO_HOME=${CHOMBO_HOME}"
echo "[build] mpicxx=$(command -v mpicxx || echo MISSING) gfortran=$(command -v gfortran || echo MISSING)"

make -C "${GRTRESNA_ROOT}/Examples/BosonStarBH" all \
  -j"$(nproc)" MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE
