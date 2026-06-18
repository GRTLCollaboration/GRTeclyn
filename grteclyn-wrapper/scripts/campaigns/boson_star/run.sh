#!/usr/bin/env bash
# Boson-star MAP-Elites survey (complex scalar / U(1) matter, 7-D search space).
#
# Selects matter via GRTRESNA_MATTER_SECTOR (geometry ansatz unchanged from default shell).
#
# Example:
#   cd grteclyn-wrapper
#   QD_NAME=boson_star_v1 QD_TARGET_EVALS=80 \
#     GPU_IDS="0 1 2 3" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/boson_star/run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-ftl_first}"
export GRTRESNA_FULL_Z=1

exec bash "${SCRIPT_DIR}/../qd/run.sh"
