#!/usr/bin/env bash
# CMA-ES refinement for boson-star QD elites.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-ftl_first}"
export GRTRESNA_FULL_Z=1

exec bash "${SCRIPT_DIR}/../cmaes/run.sh"
