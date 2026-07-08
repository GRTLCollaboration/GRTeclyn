#!/usr/bin/env bash
#
# Launch a RotatingWormholeCollapse evolution WITH the plotfile-consumer sidecar
# (frame extraction + small-data + plotfile deletion) running concurrently, per
# the disk-discipline rule (lesson L6 / README "ALWAYS extract frames on the
# fly"). Heavy HDF5 plotfiles are streamed to PNG frames + small_data and then
# deleted (--keep-last 2); they never accumulate.
#
# Usage:
#   run_rotating_wormhole.sh [PARAMS] [NGPU]
#
#   PARAMS  params filename inside Examples/RotatingWormholeCollapse
#           (default: params_rotating_complex_smoke.txt)
#   NGPU    number of GPUs / MPI ranks (default: 1)
#
# Env toggles:
#   CONSUME_PLOTFILES=0   disable the consumer sidecar (NOT recommended)
#   FRESH=1               rm -rf the run output dir before launching
#
# Examples:
#   scripts/wormhole/run_rotating_wormhole.sh params_rotating_complex_equilibrium.txt 1
#   scripts/wormhole/run_rotating_wormhole.sh params_rotating_complex_collapse.txt 4
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"

PARAMS="${1:-params_rotating_complex_smoke.txt}"
NGPU="${2:-1}"
CONSUME_PLOTFILES="${CONSUME_PLOTFILES:-1}"
FRESH="${FRESH:-0}"

EXAMPLE_DIR="${GRTECLYN_ROOT}/Examples/RotatingWormholeCollapse"
BIN="${EXAMPLE_DIR}/main3d.gnu.MPI.CUDA.ex"
PARAMS_PATH="${EXAMPLE_DIR}/${PARAMS}"

if [[ ! -x "${BIN}" ]]; then
  echo "error: binary not found: ${BIN}" >&2
  echo "build: make -C '${EXAMPLE_DIR}' -j8 USE_CUDA=TRUE USE_MPI=TRUE COMP=gnu CUDA_ARCH=90" >&2
  exit 2
fi
if [[ ! -f "${PARAMS_PATH}" ]]; then
  echo "error: params file not found: ${PARAMS_PATH}" >&2
  exit 2
fi

# Resolve the run output directory from output_path in the params file.
OUTPUT_REL="$(grep -E '^\s*output_path' "${PARAMS_PATH}" | head -1 | sed -E 's/.*=\s*"?([^"]+)"?.*/\1/')"
if [[ -z "${OUTPUT_REL}" ]]; then
  echo "error: could not parse output_path from ${PARAMS_PATH}" >&2
  exit 2
fi
OUTPUT_DIR="$(cd "${EXAMPLE_DIR}" && mkdir -p "${OUTPUT_REL}" && cd "${OUTPUT_REL}" && pwd)"

if [[ "${FRESH}" == "1" ]]; then
  echo "[run_rotating_wormhole] FRESH=1 -> clearing ${OUTPUT_DIR}"
  rm -rf "${OUTPUT_DIR}"
  mkdir -p "${OUTPUT_DIR}"
fi

echo "[run_rotating_wormhole] params=${PARAMS} ngpu=${NGPU} consumer=${CONSUME_PLOTFILES}"
echo "[run_rotating_wormhole] output=${OUTPUT_DIR}"

CONSUMER_PID=""
cleanup() {
  if [[ -n "${CONSUMER_PID}" ]] && kill -0 "${CONSUMER_PID}" 2>/dev/null; then
    echo "[run_rotating_wormhole] draining consumer (pid ${CONSUMER_PID})"
    # Wait until the sidecar has processed/deleted the trailing plotfiles
    # (yt loading is slow); cap the wait so we never hang forever.
    for _ in $(seq 1 60); do
      local_remaining=$(find "${OUTPUT_DIR}" -maxdepth 1 -type d -name 'RotatingWormholePlt*' 2>/dev/null | wc -l)
      if [[ "${local_remaining}" -le 2 ]]; then break; fi
      sleep 5
    done
    sleep 5
    kill "${CONSUMER_PID}" 2>/dev/null || true
    wait "${CONSUMER_PID}" 2>/dev/null || true
  fi
  # The consumer only manages plotfiles; prune heavy checkpoints too, keeping the
  # most recent KEEP_CHK (default 1) for restart. Set KEEP_CHK=0 to keep all.
  local keep="${KEEP_CHK:-1}"
  mapfile -t _chks < <(find "${OUTPUT_DIR}" -maxdepth 1 -type d -name 'RotatingWormholeChk*' | sort)
  local nchk=${#_chks[@]}
  if [[ "${keep}" -gt 0 && "${nchk}" -gt "${keep}" ]]; then
    local ndel=$(( nchk - keep ))
    echo "[run_rotating_wormhole] pruning ${ndel} old checkpoint(s), keeping last ${keep}"
    for ((i=0; i<ndel; i++)); do rm -rf "${_chks[$i]}"; done
  fi
}
trap cleanup EXIT INT TERM

# Start the consumer sidecar BEFORE the evolution so plotfiles are extracted and
# deleted as they appear (never accumulate). The sidecar's auto-center detection
# needs a plotfile to exist, which it does not at startup, so we set the frame /
# extraction centers explicitly from the params `center` (per Debug.md workflow).
if [[ "${CONSUME_PLOTFILES}" == "1" ]]; then
  CENTER="$(grep -E '^\s*center\s*=' "${PARAMS_PATH}" | head -1 | sed -E 's/.*=\s*//')"
  CENTER="${CENTER:-32 32 0}"
  ZOOM="${GRTECLYN_FRAMES_ZOOM:-48}"
  read -r CX CY CZ <<<"${CENTER}"
  CORNER="$(awk -v x="${CX}" -v y="${CY}" -v z="${CZ}" -v zm="${ZOOM}" \
    'BEGIN{printf "%g %g %g", x-0.5*zm, y-0.5*zm, z}')"
  echo "[run_rotating_wormhole] starting plotfile consumer sidecar (extraction center=${CENTER}, frame corner=${CORNER})"
  GRTECLYN_FRAMES_ZOOM="${ZOOM}" \
  GRTECLYN_FRAMES_CENTER="${CORNER}" \
  GRTECLYN_EXTRACTION_CENTER="${CENTER}" \
    bash "${SCRIPT_DIR}/../plot/plot_run.sh" "${OUTPUT_DIR}" \
    > "${OUTPUT_DIR}/consumer.log" 2>&1 &
  CONSUMER_PID=$!
  sleep 3
fi

cd "${EXAMPLE_DIR}"
mpirun -n "${NGPU}" bash -c \
  'export CUDA_VISIBLE_DEVICES=${OMPI_COMM_WORLD_LOCAL_RANK:-0}; exec ./main3d.gnu.MPI.CUDA.ex '"${PARAMS}"''
