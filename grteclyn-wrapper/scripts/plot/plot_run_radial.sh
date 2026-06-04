#!/usr/bin/env bash
# Stream RadialRecipe plotfiles during a run: extract small-data, optional PNG frames, delete HDF5 plot dirs.
#
# Usage (sidecar while simulation runs):
#   ./grteclyn-wrapper/scripts/plot/plot_run_radial.sh /path/to/episode_dir
#
# Or let grteclyn-wrapper start the same consumer automatically via --consume-plotfiles.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"

RUN_DIR=""
EXTRA_ARGS=()
JOBS="${JOBS:-4}"
RADII=(4 8)
KEEP_LAST="${KEEP_LAST:-1}"
DELETE="${DELETE:-1}"

while [[ $# -gt 0 ]]; do
  case $1 in
    --not-remove)
      EXTRA_ARGS+=(--keep-existing-frames)
      shift
      ;;
    --no-delete)
      DELETE=0
      shift
      ;;
    -j|--jobs)
      JOBS="$2"
      shift 2
      ;;
    --radii)
      shift
      RADII=()
      while [[ $# -gt 0 && "$1" != -* ]]; do
        RADII+=("$1")
        shift
      done
      ;;
    --keep-last)
      KEEP_LAST="$2"
      shift 2
      ;;
    *)
      if [[ -z "${RUN_DIR}" ]]; then
        RUN_DIR="$1"
      fi
      shift
      ;;
  esac
done

if [[ -z "${RUN_DIR}" ]]; then
  echo "Usage: $0 EPISODE_DIR [--radii 4 8 16] [--jobs N] [--keep-last N]" >&2
  exit 2
fi

RUN_DIR="$(cd "${RUN_DIR}" && pwd)"
SMALL_DATA_DIR="${RUN_DIR}/small_data"
FRAMES_DIR="${RUN_DIR}/frames"

if [[ ! -d "${RUN_DIR}" ]]; then
  echo "Episode directory does not exist: ${RUN_DIR}" >&2
  exit 1
fi

mkdir -p "${SMALL_DATA_DIR}" "${FRAMES_DIR}"

echo "=========================================="
echo "RadialRecipe plotfile consumer"
echo "Episode     : ${RUN_DIR}"
echo "Small data  : ${SMALL_DATA_DIR}"
echo "Frames      : ${FRAMES_DIR}"
echo "Radii       : ${RADII[*]}"
echo "Delete HDF5 : ${DELETE} (keep-last=${KEEP_LAST})"
echo "=========================================="

DELETE_ARGS=()
if [[ "${DELETE}" == "1" ]]; then
  DELETE_ARGS=(--delete --keep-last "${KEEP_LAST}")
fi

PYTHON=(python)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

"${PYTHON[@]}" -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data "${RUN_DIR}" \
  --out "${SMALL_DATA_DIR}" \
  --radii "${RADII[@]}" \
  --n-points 64 \
  --center 0 0 0 \
  --no-psi4 \
  --shell-fields chi lapse K \
  --areal-radius \
  --frames-fields chi lapse K \
  --frames-axis z \
  --frames-out "${FRAMES_DIR}" \
  --watch \
  "${DELETE_ARGS[@]}" \
  -j "${JOBS}" \
  --verbose \
  "${EXTRA_ARGS[@]}"
