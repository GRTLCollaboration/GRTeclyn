#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
SIM_ROOT="$(cd "${GRTECLYN_ROOT}/.." && pwd)"
VIS_ROOT="${WRAPPER_ROOT}/src/grteclyn_wrapper/visualisation"

DEFAULT_DATA_DIR=""
for candidate in "${SIM_ROOT}/data_2gpu" "${SIM_ROOT}/data_supported" "${SIM_ROOT}/data"; do
  if [[ -d "${candidate}" ]]; then
    DEFAULT_DATA_DIR="${candidate}"
    break
  fi
done
if [[ -z "${DEFAULT_DATA_DIR}" ]]; then
  DEFAULT_DATA_DIR="${SIM_ROOT}/data_2gpu"
fi

REMOVE_STALE=true
DATA_DIR=""
EXTRA_ARGS=""
JOBS=64
N_POINTS="${GRTECLYN_PLOT_N_POINTS:-128}"
FRAME_FIELDS="${GRTECLYN_PLOT_FRAME_FIELDS:-chi chi_minus_1 K lapse phi Pi scalar_activity Weyl4_Re Weyl4_Im Weyl4_Mag}"
FRAME_ZOOM="${GRTECLYN_FRAMES_ZOOM:-48}"
FRAME_CENTER="${GRTECLYN_FRAMES_CENTER:-}"
EXTRACTION_CENTER="${GRTECLYN_EXTRACTION_CENTER:-}"
FRAME_AUTO_ZLIM="${GRTECLYN_FRAMES_AUTO_ZLIM:-0}"
FRAME_GLOBAL_ZLIM="${GRTECLYN_FRAMES_GLOBAL_ZLIM:-1}"

while [[ $# -gt 0 ]]; do
  case $1 in
    --not-remove)
      REMOVE_STALE=false
      EXTRA_ARGS="${EXTRA_ARGS} --keep-existing-frames"
      shift
      ;;
    -j|--jobs)
      if [[ -n "${2:-}" ]] && [[ "$2" != -* ]]; then
        JOBS="$2"
        shift 2
      else
        echo "Error: Argument for $1 is missing" >&2
        exit 1
      fi
      ;;
    *)
      if [[ -z "$DATA_DIR" ]]; then
        DATA_DIR="$1"
      fi
      shift
      ;;
  esac
done

if [[ -z "$DATA_DIR" ]]; then
  DATA_DIR="${DEFAULT_DATA_DIR}"
fi
if [[ -d "${DATA_DIR}" ]]; then
  DATA_DIR="$(cd "${DATA_DIR}" && pwd)"
fi
RUN_DATA_DIR="${DATA_DIR}/data"
SMALL_DATA_DIR="${DATA_DIR}/small_data"
FRAMES_DIR="${DATA_DIR}/frames"

if [[ ! -d "${DATA_DIR}" ]]; then
  echo "Run directory does not exist: ${DATA_DIR}" >&2
  exit 1
fi

if [ "$REMOVE_STALE" = true ]; then
  echo "=========================================="
  echo "Removing stale plotfiles from: ${DATA_DIR}"
  shopt -s nullglob
  STALE_PLOTS=(
    "${DATA_DIR}"/WormholePlt*
    "${DATA_DIR}"/SupportedWormholePlt*
    "${DATA_DIR}"/RotatingWormholePlt*
    "${DATA_DIR}"/plt*
  )
  if (( ${#STALE_PLOTS[@]} > 0 )); then
    rm -rf "${STALE_PLOTS[@]}"
    echo "Removed ${#STALE_PLOTS[@]} stale plotfile directories"
  else
    echo "No stale plotfiles found"
  fi
  shopt -u nullglob

  if [[ -d "${SMALL_DATA_DIR}" ]]; then
    rm -f "${SMALL_DATA_DIR}/psi4_mode_l2m0.dat" \
          "${SMALL_DATA_DIR}/areal_radius.dat" \
          "${SMALL_DATA_DIR}/consume_state.json"
    echo "Reset extracted small-data state in: ${SMALL_DATA_DIR}"
  fi

  if [[ -d "${RUN_DATA_DIR}" ]]; then
    rm -f "${RUN_DATA_DIR}/constraint_norms.dat" \
          "${RUN_DATA_DIR}/collapse_diagnostics.dat"
    echo "Reset run diagnostics in: ${RUN_DATA_DIR}"
  fi

  rm -f "${VIS_ROOT}/constraines/constraints_plot.png" \
        "${VIS_ROOT}/constraines/constraints_plot.eps" \
        "${VIS_ROOT}/constraines/constraints_plot.pdf" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.png" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.eps" \
        "${VIS_ROOT}/diagnostic/collapse_diagnostics_plot.pdf" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.png" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.eps" \
        "${VIS_ROOT}/process_wave/psi4_extracted_R8_R12_R16.pdf" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.png" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.eps" \
        "${VIS_ROOT}/process_wave/psi4_extracted_simulation.pdf" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.png" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.eps" \
        "${VIS_ROOT}/process_wave/psi4_strain_analysis.pdf" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.png" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.eps" \
        "${VIS_ROOT}/process_wave/psi4_propagation_speed.pdf"
  if [[ -d "${FRAMES_DIR}" ]]; then
    rm -rf "${FRAMES_DIR}"
    echo "Reset frames in: ${FRAMES_DIR}"
  fi
  echo "Removed shared generated plot images"
else
  echo "=========================================="
  echo "Keeping existing plotfiles and data in: ${DATA_DIR}"
fi

mkdir -p "${FRAMES_DIR}"

PYTHON=(python)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

if [[ -z "${FRAME_CENTER}" || -z "${EXTRACTION_CENTER}" ]]; then
  read -r AUTO_FRAME_CORNER AUTO_EXTRACTION_CENTER < <(
    "${PYTHON[@]}" - "${DATA_DIR}" "${FRAME_ZOOM}" <<'PY'
import glob
import os
import sys

data_dir = sys.argv[1]
zoom = float(sys.argv[2])
fallback_corner = "8 8 0"
fallback_extract = "32 32 0"
try:
    import yt
except ImportError:
    print(fallback_corner, fallback_extract)
    raise SystemExit(0)

plotfiles = sorted(
    f
    for pattern in ("RotatingWormholePlt*", "WormholePlt*", "*Plt*")
    for f in glob.glob(os.path.join(data_dir, pattern))
    if os.path.isdir(f) and os.path.isfile(os.path.join(f, "Header"))
)
if not plotfiles:
    print(fallback_corner, fallback_extract)
    raise SystemExit(0)

ds = yt.load(plotfiles[-1])
lx = float(ds.domain_right_edge[0] - ds.domain_left_edge[0])
ly = float(ds.domain_right_edge[1] - ds.domain_left_edge[1])
throat_x = 0.5 * lx
throat_y = 0.5 * ly
corner_x = throat_x - 0.5 * zoom
corner_y = throat_y - 0.5 * zoom
print(f"{corner_x:g} {corner_y:g} 0", f"{throat_x:g} {throat_y:g} 0")
PY
  )
  if [[ -z "${FRAME_CENTER}" ]]; then
    FRAME_CENTER="${AUTO_FRAME_CORNER}"
  fi
  if [[ -z "${EXTRACTION_CENTER}" ]]; then
    EXTRACTION_CENTER="${AUTO_EXTRACTION_CENTER}"
  fi
fi

echo "=========================================="
echo "Watching plotfiles in: ${DATA_DIR}"
echo "Writing small-data to: ${SMALL_DATA_DIR}"
echo "Writing frames to:     ${FRAMES_DIR}"
echo "Frame fields:          ${FRAME_FIELDS}"
echo "Frame sample points:   ${N_POINTS}"
echo "Frame zoom:            ${FRAME_ZOOM}"
echo "Frame corner origin:   ${FRAME_CENTER}"
echo "Extraction center:     ${EXTRACTION_CENTER}"
echo "Frame auto colorbar:   ${FRAME_AUTO_ZLIM}"
echo "Frame global colorbar: ${FRAME_GLOBAL_ZLIM}"
echo "=========================================="

AUTO_ZLIM_ARGS=()
if [[ "${FRAME_AUTO_ZLIM}" == "1" ]]; then
  AUTO_ZLIM_ARGS=(--frames-auto-zlim)
fi
GLOBAL_ZLIM_ARGS=(--frames-global-zlim)
if [[ "${FRAME_GLOBAL_ZLIM}" == "0" ]]; then
  GLOBAL_ZLIM_ARGS=(--no-frames-global-zlim)
fi

"${PYTHON[@]}" -m grteclyn_wrapper.visualisation.process_wave.consume_plotfiles \
  --data "${DATA_DIR}" \
  --out "${DATA_DIR}/small_data" \
  --center ${EXTRACTION_CENTER} \
  --radii 12 16 20 24 \
  --n-points "${N_POINTS}" \
  --areal-radius \
  --embedding --embedding-rmax 5.0 \
  --frames-fields ${FRAME_FIELDS} \
  --frames-axis z \
  --frames-zoom "${FRAME_ZOOM}" \
  --frames-center ${FRAME_CENTER} \
  --frames-corner \
  "${AUTO_ZLIM_ARGS[@]}" \
  "${GLOBAL_ZLIM_ARGS[@]}" \
  --frames-out "${FRAMES_DIR}" \
  --watch --delete --keep-last 2 \
  --verbose -j "${JOBS}" ${EXTRA_ARGS}
