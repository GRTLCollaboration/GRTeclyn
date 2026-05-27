#!/usr/bin/env bash
# Sidecar: run plot_run_radial.sh on a wrapper episode directory.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${SCRIPT_DIR}/../../src/scripts/plot_run_radial.sh" "$@"
