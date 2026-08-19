#!/usr/bin/env bash
# Bondi dipole rerun campaign, wave 1 -- four cells, one per GPU.
#
# Wave 1 deliberately does NOT launch the expensive N=256 convergence cells.
# Item I (the CFL experiment) decides whether every later cell costs 1x or 10x,
# and it answers in ~30 min, so it runs first -- alongside item B, which needs
# no code change and would otherwise leave two H100s idle.
#
#   GPU 0  item I     single_p, dt_multiplier 0.2   -- the CFL test
#   GPU 1  control    single_p, dt_multiplier 0.02  -- CFL baseline AND the
#                                                      first reproduction of a
#                                                      published cell on this
#                                                      machine + rebuilt binary
#   GPU 2  item B     pair_pm     + scrutiny stream (momentum balance)
#   GPU 3  item B     pair_pm_eqm + scrutiny stream (momentum balance)
#
# Every cell writes under runs/bondi_rerun/<cell>/ so the published July tree
# in runs/bondi/ is never touched.
#
# Memory: N=128 measured at ~7 GB arena + 0.5 GB managed (July pair_pm run.log),
# one cell per 80 GB H100 -- no OOM risk at this wave's resolution.
#
# Stop everything with:
#   for d in runs/bondi_rerun/*/; do bash grteclyn-wrapper/scripts/campaigns/stop_campaign.sh "$d"; done
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../../.." && pwd)"
cd "${REPO_ROOT}"

CAMP="grteclyn-wrapper/scripts/campaigns/bondi_dipole"
OUT="${REPO_ROOT}/runs/bondi_rerun"
mkdir -p "${OUT}"

# NOTE: /usr/bin/env, spelled out, is NOT optional on this machine.
# Two shared-workspace directories sit ahead of /usr/bin on PATH and each
# contains an executable named `env` which is really a PATH-setup snippet meant
# to be sourced (they live under other users' home and shared-workspace bin
# directories -- run `type -a env` to see which ones are shadowing it for you).
# Invoked as `env VAR=x cmd` they prepend their own directory to PATH and exit 0
# WITHOUT EVER RUNNING cmd.  A bare `env` here therefore launches nothing,
# silently, with a success exit code and an empty log.
launch() {
  local tag="$1"; shift
  echo "[wave1] launching ${tag}"
  setsid nohup /usr/bin/env "$@" BONDI_RUNS_DIR="${OUT}/${tag}" \
    bash "${CAMP}/${SCRIPT}" \
    > "${OUT}/${tag}.launch.log" 2>&1 < /dev/null &
  disown
}

# --- item I: CFL experiment (dt 0.02 -> 0.2) --------------------------------
SCRIPT=run_single_selfgrav.sh
launch cfl_dt0p2   BONDI_GPU=0 BONDI_DT_MULT=0.2  BONDI_STOP_TIME=40

# --- control: same cell at the published dt, on this machine ---------------
SCRIPT=run_single_selfgrav.sh
launch cfl_base    BONDI_GPU=1 BONDI_DT_MULT=0.02 BONDI_STOP_TIME=40

# --- item B: momentum balance, scrutiny stream on ---------------------------
SCRIPT=run_pair_selfgrav.sh
launch momB_pm     BONDI_GPU=2 BONDI_S0=0 BONDI_S1=1 BONDI_SCRUTINY=1 BONDI_STOP_TIME=60
launch momB_pm_eqm BONDI_GPU=3 BONDI_S0=0 BONDI_S1=1 BONDI_S1_OMEGA=0.56598 \
                   BONDI_SCRUTINY=1 BONDI_STOP_TIME=60

sleep 5
echo
echo "[wave1] launched. Verify with:  pgrep -af bondi_sg_"
