#!/usr/bin/env bash
# Gravitational splash MAP-Elites (canonical scalar shell, critical_collapse).
# Uses the validated GRTresna scalar-lump shell ansatz (compact profile) rather
# than boson-star matter. Scoring gates focus on absolute peak rho at origin.
#
# Example:
#   cd grteclyn-wrapper
#   QD_NAME=scalar_splash_gpu_v1 QD_TARGET_EVALS=12 \
#     GPU_IDS="0 1" GPU_SLOTS_PER_DEVICE=1 \
#     bash scripts/campaigns/splash/run.sh

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

export GRTRESNA_MATTER_SECTOR=scalar
export GRTRESNA_ANSATZ=shell
export GRTRESNA_SHELL_PROFILE="${GRTRESNA_SHELL_PROFILE:-compact}"
export GRTRESNA_MATTER_COUPLING=canonical
export OBJECTIVE_MODE="${OBJECTIVE_MODE:-critical_collapse}"
export DESCRIPTOR_MODE="${DESCRIPTOR_MODE:-wave_focusing}"
export SPLASH_MODE="${SPLASH_MODE:-discovery}"
export GRTRESNA_FULL_Z=1
export GRTECLYN_FRAMES=1
export GRTECLYN_FRAMES_ZOOM="${GRTECLYN_FRAMES_ZOOM:-none}"
export PROJECTION_FIELDS="${PROJECTION_FIELDS:-scalar_activity}"
export PROJECTION_AXES="${PROJECTION_AXES:-x y z}"

exec bash "${SCRIPT_DIR}/../qd/run.sh"
