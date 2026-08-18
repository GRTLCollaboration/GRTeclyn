#!/usr/bin/env bash
# Launch one matrix cell for fgeo_max_cmaes_v1 (RM / RC / RF / …).
#
#   bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh --list
#   DRY_RUN=1 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RM
#   GPU_ID=2 bash scripts/campaigns/promote/fgeo_max_cmaes_v1/run.sh FMAX-RI
set -euo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=campaign.env.sh
source "${CAMPAIGN_DIR}/campaign.env.sh"
campaign_resolve_source

echo "== ${CAMPAIGN_NAME} =="
echo "  manifest   : ${MANIFEST}"
echo "  source_run : ${SOURCE_RUN}"
echo "  source_eval: ${SOURCE_EVAL_ID:-auto}"
echo "  pump_stop  : ${RL_PUMP_STOP_TIME}"
echo

export MANIFEST SOURCE_RUN SOURCE_EVAL_ID VALIDATION_LAUNCH_LOG_DIR
export FREEZE_HINT="${CAMPAIGN_DIR}/freeze.sh"

exec bash "${PROMOTE_LIB}/run_matrix.sh" "$@"
