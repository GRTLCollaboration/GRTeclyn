#!/usr/bin/env bash
# Resolve campaign paths and source shared environment. Do not run directly.
set -euo pipefail

_campaign_bootstrap() {
  local stage_dir="$1"
  CAMPAIGNS_ROOT="$(cd -- "${stage_dir}/.." && pwd)"
  SCRIPTS_ROOT="$(cd -- "${CAMPAIGNS_ROOT}/.." && pwd)"
  # shellcheck source=../../lib/env.sh
  source "${SCRIPTS_ROOT}/lib/env.sh"
  # shellcheck source=search_common.sh
  source "${CAMPAIGNS_ROOT}/lib/search_common.sh"
}

_campaign_resolve_python() {
  if [[ -z "${PYTHON_BIN:-}" ]]; then
    if [[ -x "${WRAPPER_ROOT}/.venv/bin/python" ]]; then
      PYTHON_BIN="${WRAPPER_ROOT}/.venv/bin/python"
    elif command -v uv >/dev/null 2>&1 && [[ "${USE_UV:-1}" == "1" ]]; then
      PYTHON_BIN="uv run --project ${WRAPPER_ROOT} python"
    else
      PYTHON_BIN="python"
    fi
  fi
  export PYTHON_BIN
}
