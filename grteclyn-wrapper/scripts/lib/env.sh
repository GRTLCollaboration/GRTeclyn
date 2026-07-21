#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
export WRAPPER_ROOT="$(cd -- "${SCRIPTS_ROOT}/.." && pwd)"
export SCRIPTS_ROOT

# Gitignored machine overlay. Prefer .env; fall back to legacy site.local.env.
# See .env.example. Already-exported vars are not overwritten (set -a + source
# would overwrite — we only source when the key is unset).
_load_dotenv_file() {
  local file="$1"
  [[ -f "${file}" ]] || return 0
  local line key value
  while IFS= read -r line || [[ -n "${line}" ]]; do
    line="${line%%#*}"
    line="${line#"${line%%[![:space:]]*}"}"
    line="${line%"${line##*[![:space:]]}"}"
    [[ -z "${line}" ]] && continue
    [[ "${line}" == export\ * ]] && line="${line#export }"
    [[ "${line}" == *=* ]] || continue
    key="${line%%=*}"
    value="${line#*=}"
    key="${key%"${key##*[![:space:]]}"}"
    key="${key#"${key%%[![:space:]]*}"}"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    if [[ "${value}" == \"*\" && "${value}" == *\" ]]; then
      value="${value:1:${#value}-2}"
    elif [[ "${value}" == \'*\' && "${value}" == *\' ]]; then
      value="${value:1:${#value}-2}"
    fi
    # Skip if already set in the environment.
    if [[ -n "${key}" ]] && [[ -z "${!key+x}" ]]; then
      # shellcheck disable=SC2086,SC2163
      export "${key}=${value}"
    fi
  done < "${file}"
}

if [[ -f "${WRAPPER_ROOT}/.env" ]]; then
  _load_dotenv_file "${WRAPPER_ROOT}/.env"
elif [[ -f "${WRAPPER_ROOT}/site.local.env" ]]; then
  _load_dotenv_file "${WRAPPER_ROOT}/site.local.env"
fi

# Expand ${VAR} / $VAR in the path knobs after load (bash-native).
_expand_knob() {
  local name="$1"
  local raw="${!name-}"
  [[ -z "${raw}" ]] && return 0
  # shellcheck disable=SC2086
  export "${name}=$(eval "printf '%s' \"${raw}\"")"
}

for _knob in SIM_ROOT GRTECLYN_ROOT GRTRESNA_ROOT CHOMBO_HOME OPENMPI_ROOT GRTRESNA_ENV; do
  _expand_knob "${_knob}"
done
unset _knob

if [[ -z "${GRTECLYN_ROOT:-}" ]]; then
  if [[ -f "${WRAPPER_ROOT}/../Examples/SupportedWormholeCollapse/params_2gpu.txt" ]]; then
    export GRTECLYN_ROOT="$(cd -- "${WRAPPER_ROOT}/.." && pwd)"
  else
    echo "Set GRTECLYN_ROOT=/path/to/GRTeclyn (or create ${WRAPPER_ROOT}/.env from .env.example)." >&2
    exit 2
  fi
fi

export SIM_ROOT="${SIM_ROOT:-$(cd -- "${GRTECLYN_ROOT}/.." && pwd)}"
export GRTRESNA_ROOT="${GRTRESNA_ROOT:-${SIM_ROOT}/GRTresna}"
export CHOMBO_HOME="${CHOMBO_HOME:-${SIM_ROOT}/Chombo/lib}"
export PYTHONPATH="${WRAPPER_ROOT}/src:${PYTHONPATH:-}"

export OPENMPI_ROOT="${OPENMPI_ROOT:-${SIM_ROOT}/local/openmpi}"
if [[ -d "${OPENMPI_ROOT}/bin" ]]; then
  export PATH="${OPENMPI_ROOT}/bin:${PATH}"
  export LD_LIBRARY_PATH="${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH:-}"
fi

# GRTRESNA_ENV must come from the environment or .env — no host-specific defaults.
if [[ -n "${GRTRESNA_ENV:-}" ]] && [[ -d "${GRTRESNA_ENV}/bin" ]]; then
  export PATH="${GRTRESNA_ENV}/bin:${PATH}"
  export LD_LIBRARY_PATH="${GRTRESNA_ENV}/lib:${LD_LIBRARY_PATH:-}"
fi

cd "${GRTECLYN_ROOT}"
