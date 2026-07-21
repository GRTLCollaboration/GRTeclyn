#!/usr/bin/env bash
# Rebuild the GRTresna BosonStarBH initial-data solver (CTTKHybrid<ComplexScalarField>),
# used to generate rotating-wormhole constraint-clean initial data.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/../../lib/env.sh"

GRTRESNA_ROOT="${GRTRESNA_ROOT:?GRTRESNA_ROOT unset}"
CHOMBO_HOME="${CHOMBO_HOME:?CHOMBO_HOME unset}"

echo "[build] CHOMBO_HOME=${CHOMBO_HOME}"
echo "[build] GRTRESNA_ROOT=${GRTRESNA_ROOT}"
echo "[build] mpicxx=$(command -v mpicxx || echo MISSING) gfortran=$(command -v gfortran || echo MISSING)"

make -C "${GRTRESNA_ROOT}/Examples/BosonStarBH" all \
  -j"$(nproc)" MPI=TRUE DIM=3 OPT=HIGH DEBUG=FALSE
