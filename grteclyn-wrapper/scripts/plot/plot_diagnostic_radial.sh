#!/usr/bin/env bash
# Post-run RadialRecipe diagnostics (no gravitational-wave extraction).
#
# Usage: ./grteclyn-wrapper/scripts/plot/plot_diagnostic_radial.sh [EPISODE_DIR]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/env.sh
source "${SCRIPT_DIR}/../lib/env.sh"
VIS_DIR="${WRAPPER_ROOT}/src/grteclyn_wrapper/visualisation"
PLOTS_DIR="${VIS_DIR}/plots/radial"

choose_latest_episode() {
  local root="${GRTECLYN_ROOT}/runs"
  local best="" best_mtime=-1
  for dir in "${root}"/radialrecipe_*/* "${root}"/radialrecipe_*; do
    [[ -d "${dir}/data" ]] || continue
    [[ -f "${dir}/data/constraint_norms.dat" ]] || continue
    local mtime
    mtime=$(stat -c %Y "${dir}/data/constraint_norms.dat")
    if (( mtime > best_mtime )); then
      best_mtime="${mtime}"
      best="${dir}"
    fi
  done
  if [[ -n "${best}" ]]; then
    printf '%s\n' "${best}"
  fi
}

RUN_DIR="${1:-$(choose_latest_episode || true)}"
if [[ -z "${RUN_DIR}" ]]; then
  echo "Usage: $0 EPISODE_DIR" >&2
  exit 2
fi
RUN_DIR="$(cd "${RUN_DIR}" && pwd)"

CONSTRAINT_FILE="${RUN_DIR}/data/constraint_norms.dat"
COLLAPSE_FILE="${RUN_DIR}/data/collapse_diagnostics.dat"
SHELL_FILE="${RUN_DIR}/small_data/shell_profiles.dat"
AREAL_FILE="${RUN_DIR}/small_data/areal_radius.dat"

if [[ ! -f "${CONSTRAINT_FILE}" ]]; then
  echo "Missing ${CONSTRAINT_FILE}" >&2
  exit 1
fi
if [[ ! -f "${COLLAPSE_FILE}" ]]; then
  echo "Missing ${COLLAPSE_FILE}" >&2
  exit 1
fi

rm -rf "${PLOTS_DIR}"
mkdir -p "${PLOTS_DIR}"

echo "Using episode: ${RUN_DIR}"
echo "Writing plots to: ${PLOTS_DIR}"

PYTHON=(python3)
if command -v uv >/dev/null 2>&1 && [[ -f "${WRAPPER_ROOT}/pyproject.toml" ]]; then
  PYTHON=(uv run --directory "${WRAPPER_ROOT}" python)
fi

echo "[1/3] Constraint norms..."
"${PYTHON[@]}" -m grteclyn_wrapper.visualisation.constraines \
  "${CONSTRAINT_FILE}" \
  -o "${PLOTS_DIR}/constraints_plot.eps"

echo "[2/3] Collapse diagnostics..."
"${PYTHON[@]}" "${VIS_DIR}/diagnostic/diagnostic.py" \
  "${COLLAPSE_FILE}" \
  --data "${RUN_DIR}" \
  --out "${PLOTS_DIR}" \
  --name collapse_diagnostics_plot.eps

if [[ -f "${SHELL_FILE}" ]]; then
  echo "[3/3] Shell profiles (chi/lapse/K)..."
  "${PYTHON[@]}" -m grteclyn_wrapper.visualisation.radial.plot_shell_profiles \
    "${SHELL_FILE}" \
    --out "${PLOTS_DIR}" \
    --name shell_profiles
else
  echo "[3/3] Skipping shell profiles (missing ${SHELL_FILE})"
fi

if [[ -f "${AREAL_FILE}" ]]; then
  echo "Areal radius series: ${AREAL_FILE}"
fi

echo "Done. Plots in ${PLOTS_DIR}"
ls -1 "${PLOTS_DIR}"
