#!/usr/bin/env bash
# Post-run frames drain: render PNG slice movies + Psi4 from remaining plotfiles.
#
# Skips the slow paths that made the old consumer ~3 min/plotfile:
#   - no --evolving-geodesic / metric_stack cache
#   - no FTL timeseries / incremental score / shell / areal
#   - no projection panels unless GRTECLYN_PROJECTION_FIELDS is set
#
# Usage:
#   ./grteclyn-wrapper/scripts/plot/drain_plotfiles_frames.sh RUN_DIR
#   CONSUMER_JOBS=2 FRAMES_FIELDS="Weyl4_Re Weyl4_Im Weyl4_Mag chi" ... RUN_DIR
#   GRTECLYN_PSI4=0 ...   # frames only (Psi4 already extracted)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
# shellcheck source=../campaigns/lib/promote_common.sh
source "${SCRIPT_DIR}/../campaigns/lib/promote_common.sh"

RUN_DIR="${1:?RUN_DIR required}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"

export GRTECLYN_FRAMES=1
# Fixed colorbars across time (no per-frame percentile bounce in movies).
export GRTECLYN_FRAMES_STABLE_MOVIE="${GRTECLYN_FRAMES_STABLE_MOVIE:-1}"
export GRTECLYN_CONSUMER_DRAIN=0
export GRTECLYN_CONSUMER_DRAIN_MINIMAL=1
# Psi4 already extracted during fast drain — skip on backlog (saves ~20s/plotfile).
export GRTECLYN_PSI4="${GRTECLYN_PSI4:-0}"
export GRTECLYN_PSI4_N_POINTS="${GRTECLYN_PSI4_N_POINTS:-128}"
# NFS plotfiles (~2.7 GB each): j>2 thrashes network and hangs; j=2 is the sweet spot.
export CONSUMER_JOBS="${CONSUMER_JOBS:-2}"
export GRTECLYN_EVOLVING_GEODESIC=0
# Projections dominate cost (2 fields × 3 axes); off unless caller sets PROJECTION_FIELDS.
export GRTECLYN_PROJECTION_FIELDS="${PROJECTION_FIELDS:-}"
export GRTECLYN_PROJECTION_AXES="${PROJECTION_AXES:-}"
# Full HQ field set (same as live consumer / promote_common boson_star). Override via FRAMES_FIELDS.
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-scalar_activity phi Pi phi_lump0 Pi_lump0 chi chi_minus_1 local_speed shift1 rho_req Weyl4_Re Weyl4_Im Weyl4_Mag}"

PYTHON=(python3)
if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
  PYTHON=("${WRAPPER_ROOT}/.venv/bin/python")
elif command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

echo "Frames drain: ${RUN_DIR}"
echo "  fields=${GRTECLYN_FRAMES_FIELDS}"
echo "  jobs=${CONSUMER_JOBS}  psi4=${GRTECLYN_PSI4}  projections=${GRTECLYN_PROJECTION_FIELDS:-off}"
echo "Plotfiles on disk: $(find "${RUN_DIR}" -maxdepth 1 -type d -name 'RadialRecipePlt*' 2>/dev/null | wc -l)"

"${PYTHON[@]}" - <<PY
from pathlib import Path
from grteclyn_wrapper.core.episode import Episode
from grteclyn_wrapper.core.plot_consumer import consumer_radii_from_env
from grteclyn_wrapper.core.runner import drain_plotfile_backlog

episode = Episode(path=Path("${RUN_DIR}"))
params = episode.params_path.read_text(encoding="utf-8")
stop_time = 30.0
for line in params.splitlines():
    if line.strip().startswith("stop_time"):
        try:
            stop_time = float(line.split("=", 1)[1].strip().split()[0])
        except ValueError:
            pass
        break

result = drain_plotfile_backlog(
    episode,
    example_name="RadialRecipe",
    radii=consumer_radii_from_env(),
    delete=True,
    keep_last=1,
    frames=True,
    ftl_timeseries=False,
    ftl_L=8.0,
    incremental_score=False,
    objective_mode="general_ftl",
    target_stop_time=stop_time,
    evolving_geodesic=False,
)
raise SystemExit(result.returncode)
PY

echo "Frames drain complete. Remaining plotfiles: $(find "${RUN_DIR}" -maxdepth 1 -type d -name 'RadialRecipePlt*' 2>/dev/null | wc -l)"
