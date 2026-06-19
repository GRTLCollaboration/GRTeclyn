#!/usr/bin/env bash
# Gravitational splash MAP-Elites (boson star, critical_collapse objective).
#
# Example:
#   cd grteclyn-wrapper
#   QD_NAME=boson_splash_v1 QD_TARGET_EVALS=12 \
#     GPU_IDS="0 1" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/splash/run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=boson_star
export GRTRESNA_ANSATZ=splash
# Splash studies use positive-energy boson matter only (no phantom / −ρ).
export GRTRESNA_MATTER_COUPLING=canonical
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-critical_collapse}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-wave_focusing}"
export SPLASH_MODE="${SPLASH_MODE:-discovery}"
export GRTRESNA_FULL_Z=1
# Render slice/projection frames (search_common defaults GRTECLYN_FRAMES=0 for speed).
export GRTECLYN_FRAMES=1
export GRTECLYN_FRAMES_ZOOM="${GRTECLYN_FRAMES_ZOOM:-none}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity phi}"
export PROJECTION_AXES="${PROJECTION_AXES:-x y z}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
