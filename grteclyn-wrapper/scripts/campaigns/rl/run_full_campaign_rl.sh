#!/usr/bin/env bash
# Full campaign orchestrator: QD -> CMA-ES -> RL -> HQ (manual stage chaining).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "${ROOT}"

echo "Stage 0: MAP-Elites QD (existing launcher)"
echo "  bash scripts/campaigns/general_ftl/run_all.sh"

echo "Stage 1: CMA-ES refine (existing launcher)"
echo "  bash scripts/campaigns/cmaes/run.sh"

echo "Stage 1.5: RL training"
echo "  RL_EPISODE_PATH=... GRTeclyn_EXE=... RL_PARAMS=... bash scripts/campaigns/rl/run.sh"

echo "Stage 2: HQ promotion (existing launcher)"
echo "  bash scripts/campaigns/hq/run_batch.sh"
