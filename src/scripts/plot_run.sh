#!/usr/bin/env bash
# Deprecated shim — use grteclyn-wrapper/scripts/plot_run.sh
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/../../grteclyn-wrapper/scripts/plot_run.sh" "$@"
