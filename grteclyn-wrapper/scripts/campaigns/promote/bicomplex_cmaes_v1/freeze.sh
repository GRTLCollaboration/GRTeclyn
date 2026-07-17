#!/usr/bin/env bash
# Freeze CMA-ES champion for this campaign.
#
#   bash scripts/campaigns/promote/bicomplex_cmaes_v1/freeze.sh
#   SOURCE_EVAL_ID=63 bash scripts/campaigns/promote/bicomplex_cmaes_v1/freeze.sh
set -euo pipefail

CAMPAIGN_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=campaign.env.sh
source "${CAMPAIGN_DIR}/campaign.env.sh"

export SOURCE_RUN="${SOURCE_RUN:-${LIVE_RUN}}"
export FREEZE_ROOT

bash "${PROMOTE_LIB}/freeze_champion.sh"

echo
echo "Next:"
echo "  bash ${CAMPAIGN_DIR}/run.sh BCMA-RM"
