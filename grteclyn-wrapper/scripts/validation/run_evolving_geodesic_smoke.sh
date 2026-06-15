#!/usr/bin/env bash
# Fast smoke for the 4D evolving null-geodesic probe (no GPU, no plotfiles).
set -euo pipefail

SEARCH_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SEARCH_DIR}/../lib/env.sh"

cd "${WRAPPER_ROOT}"

echo "== evolving geodesic unit tests =="
uv run pytest -q \
  tests/metrics/ftl/test_evolving_geodesic.py \
  tests/metrics/ftl/test_metric_field.py \
  tests/metrics/aggregation/test_evolving_geodesic_collector.py

echo
echo "== analytic Alcubierre validation =="
uv run python scripts/validation/evolving_geodesic_validation.py --v-s 2.0

echo
echo "Smoke OK."
