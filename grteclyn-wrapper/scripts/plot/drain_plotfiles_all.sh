#!/usr/bin/env bash
# Full post-run (or live) plotfile drain: Psi4 + geodesic + projections + frames.
# This is the ONE-SHOT pass that replaces the old sequence of:
#   drain_plotfiles_fast.sh (psi4)  →  drain_plotfiles_frames.sh (frames)
# Reading each 2.7 GB plotfile twice was killing NFS throughput; this does it
# once and extracts everything.
#
# Use after the GPU stops, OR live as a watcher via --watch flag.
#
# Usage:
#   ./grteclyn-wrapper/scripts/plot/drain_plotfiles_all.sh RUN_DIR
#   CONSUMER_JOBS=1 ./grteclyn-wrapper/scripts/plot/drain_plotfiles_all.sh RUN_DIR
#   GRTECLYN_PSI4=0 ./grteclyn-wrapper/scripts/plot/drain_plotfiles_all.sh RUN_DIR   # skip psi4
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
# shellcheck source=../campaigns/lib/promote_common.sh
source "${SCRIPT_DIR}/../campaigns/lib/promote_common.sh"

RUN_DIR="${1:?RUN_DIR required}"
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"

# Full HQ profile
export GRTECLYN_FRAMES=1
export GRTECLYN_PSI4="${GRTECLYN_PSI4:-1}"
export GRTECLYN_PSI4_N_POINTS="${GRTECLYN_PSI4_N_POINTS:-128}"
# NFS plotfiles (~2.7 GB each): j>1 thrashes network. j=1 is the safe default for
# full-profile drains on NFS. Bump to 2 on local NVMe.
export CONSUMER_JOBS="${CONSUMER_JOBS:-1}"
export GRTECLYN_EVOLVING_GEODESIC="${GRTECLYN_EVOLVING_GEODESIC:-1}"
export GRTECLYN_EVOLVING_GEODESIC_MODE="${GRTECLYN_EVOLVING_GEODESIC_MODE:-hq}"
export GRTECLYN_GEO_DIRECTIONS="${GRTECLYN_GEO_DIRECTIONS:-x y z}"
# FTL metrics + confinement + incremental score
export GRTECLYN_FTL_TIMESERIES="${GRTECLYN_FTL_TIMESERIES:-1}"
export GRTECLYN_CONFINEMENT_TIMESERIES="${GRTECLYN_CONFINEMENT_TIMESERIES:-1}"
export GRTECLYN_INCREMENTAL_SCORE="${GRTECLYN_INCREMENTAL_SCORE:-1}"
# Projections (2 fields × 3 axes). Override PROJECTION_FIELDS="" to disable.
export GRTECLYN_PROJECTION_FIELDS="${GRTECLYN_PROJECTION_FIELDS:-scalar_activity phi}"
export GRTECLYN_PROJECTION_AXES="${GRTECLYN_PROJECTION_AXES:-x y z}"
export GRTECLYN_PROJECTION_METHOD="${GRTECLYN_PROJECTION_METHOD:-mip}"
# Shell profiles + areal radius
export GRTECLYN_SHELL_FIELDS="${GRTECLYN_SHELL_FIELDS:-chi lapse K}"
export GRTECLYN_AREAL_RADIUS="${GRTECLYN_AREAL_RADIUS:-1}"
# Full HQ frame field set (promote_common boson_star branch). Override FRAMES_FIELDS to customize.
export GRTECLYN_FRAMES_FIELDS="${FRAMES_FIELDS:-scalar_activity phi Pi phi_lump0 Pi_lump0 chi chi_minus_1 local_speed shift1 rho_req Weyl4_Re Weyl4_Im Weyl4_Mag}"
export GRTECLYN_FRAMES_ZOOM="${FRAMES_ZOOM:-none}"

PYTHON=(python3)
if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
  PYTHON=("${WRAPPER_ROOT}/.venv/bin/python")
elif command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

echo "Full drain: ${RUN_DIR}"
echo "  fields=${GRTECLYN_FRAMES_FIELDS}"
echo "  jobs=${CONSUMER_JOBS}  psi4=${GRTECLYN_PSI4}(n=${GRTECLYN_PSI4_N_POINTS})  geodesic=${GRTECLYN_EVOLVING_GEODESIC}"
echo "  projections=${GRTECLYN_PROJECTION_FIELDS} axes=${GRTECLYN_PROJECTION_AXES}"
echo "  shell=${GRTECLYN_SHELL_FIELDS}  areal=${GRTECLYN_AREAL_RADIUS}"
echo "  ftl=${GRTECLYN_FTL_TIMESERIES}  confinement=${GRTECLYN_CONFINEMENT_TIMESERIES}  incremental_score=${GRTECLYN_INCREMENTAL_SCORE}"
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

# HQ replay_eval default: ftl_L=8 (NOT L_full=128 — that was the stuck-drain bug).
ftl_l = 8.0

result = drain_plotfile_backlog(
    episode,
    example_name="RadialRecipe",
    radii=consumer_radii_from_env(),
    delete=True,
    keep_last=1,
    frames=True,
    ftl_timeseries=bool(${GRTECLYN_FTL_TIMESERIES}),
    ftl_L=ftl_l,
    incremental_score=bool(${GRTECLYN_INCREMENTAL_SCORE}),
    objective_mode="general_ftl",
    target_stop_time=stop_time,
    evolving_geodesic=bool(${GRTECLYN_EVOLVING_GEODESIC}),
)
raise SystemExit(result.returncode)
PY

echo "Full drain complete. Remaining plotfiles: $(find "${RUN_DIR}" -maxdepth 1 -type d -name 'RadialRecipePlt*' 2>/dev/null | wc -l)"
