#!/usr/bin/env bash
# HQ promotion for boson-star QD/CMA-ES elites.

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_MATTER_COUPLING="${GRTRESNA_MATTER_COUPLING:-canonical}"
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-ftl_first}"

exec bash "${SCRIPT_DIR}/../hq/run_batch.sh"
