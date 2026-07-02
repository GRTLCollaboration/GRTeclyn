#!/usr/bin/env bash
# Fast post-run plotfile drain: parallel Psi4 extraction + HDF5 delete (no PNG frames).
#
# Why the old drain was ~3 min/plotfile:
#   - --evolving-geodesic builds a full metric_stack slice per plotfile (very slow on 256³ AMR)
#   - a bad --ftl-l 128 (L_full) made that stack enormous
#   - 19 slice fields + 6 projections when frames were enabled
#
# This script defaults to minimal mode: Psi4 + delete only, CONSUMER_JOBS=4.
#
# Usage:
#   ./grteclyn-wrapper/scripts/plot/drain_plotfiles_fast.sh RUN_DIR
#   CONSUMER_JOBS=8 ./grteclyn-wrapper/scripts/plot/drain_plotfiles_fast.sh RUN_DIR
#   GRTECLYN_CONSUMER_DRAIN_MINIMAL=0 CONSUMER_JOBS=2 ...  # also FTL/geodesic (slow)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"

RUN_DIR="${1:?RUN_DIR required}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"

export GRTECLYN_FRAMES=0
export GRTECLYN_CONSUMER_DRAIN=1
export GRTECLYN_PSI4="${GRTECLYN_PSI4:-1}"
export GRTECLYN_CONSUMER_DRAIN_MINIMAL="${GRTECLYN_CONSUMER_DRAIN_MINIMAL:-1}"
export CONSUMER_JOBS="${CONSUMER_JOBS:-2}"
export GRTECLYN_PSI4_N_POINTS="${GRTECLYN_PSI4_N_POINTS:-32}"
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-0}"

PYTHON=(python3)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

echo "Fast drain: ${RUN_DIR}"
echo "  minimal=${GRTECLYN_CONSUMER_DRAIN_MINIMAL}  jobs=${CONSUMER_JOBS}  psi4=${GRTECLYN_PSI4}"
echo "Plotfiles on disk: $(find "${RUN_DIR}" -maxdepth 1 -type d -name 'RadialRecipePlt*' 2>/dev/null | wc -l)"

"${PYTHON[@]}" - <<PY
from pathlib import Path
from grteclyn_wrapper.core.episode import Episode
from grteclyn_wrapper.core.plot_consumer import consumer_drain_minimal, consumer_radii_from_env
from grteclyn_wrapper.core.runner import drain_plotfile_backlog

episode = Episode(path=Path("${RUN_DIR}"))
minimal = consumer_drain_minimal()
params = episode.params_path.read_text(encoding="utf-8")
stop_time = 30.0
for line in params.splitlines():
    if line.strip().startswith("stop_time"):
        try:
            stop_time = float(line.split("=", 1)[1].strip().split()[0])
        except ValueError:
            pass
        break

# Match HQ replay_eval default (NOT L_full — that was the stuck-drain bug).
ftl_l = 8.0

result = drain_plotfile_backlog(
    episode,
    example_name="RadialRecipe",
    radii=consumer_radii_from_env(),
    delete=True,
    keep_last=1,
    frames=False,
    ftl_timeseries=not minimal,
    ftl_L=ftl_l,
    incremental_score=not minimal,
    objective_mode="general_ftl",
    target_stop_time=stop_time,
    evolving_geodesic=not minimal,
)
raise SystemExit(result.returncode)
PY

echo "Drain complete. Remaining plotfiles: $(find "${RUN_DIR}" -maxdepth 1 -type d -name 'RadialRecipePlt*' 2>/dev/null | wc -l)"
