#!/usr/bin/env bash
# Bondi dipole -- STABLE-BRANCH qualification campaign (docs/MatterDebugg.md).
#
# The omega=0.55 matter model is rejected: the canonical star sits past the
# turning point of its family and collapses on its own by t~60 in every cell
# ever run, pair or single.  Before any new pair is booked, the replacement
# star has to prove it can sit still.  Four cells, one per card:
#
#   single_p_w080 / single_m_w080   lone canonical / lone phantom, omega=0.80
#       The design point: the r90 minimum of BOTH sectors (5.11 / 5.18), mass
#       ratio 1.043, and formally on the stable side of the turning point
#       (~0.67) for the canonical star.  Expected t=0 fingerprint from the
#       family scan (results/bondi-dipole-runaway/stars/star_radius.csv):
#       canonical phi_c=0.019126 ADM=+0.011209; phantom phi_c=0.019199
#       ADM=-0.011692.
#   single_p_w075 / single_m_w075   the same pair of tests at omega=0.75
#       The backup frequency: 28% more mass (stronger Bondi signal, t~95
#       instead of t~107 for the same drift) at the price of sitting closer
#       to the turning point.  If 0.80 passes, a pair campaign can still
#       prefer 0.75 for signal -- but only if 0.75 passes here too.
#
# t=120 because the follow-on pair run needs t~110 at omega=0.80 to
# accumulate the drift the old model showed by t=45.
#
# PROTOTYPE GRID, on purpose: N=128 at L=64 (dx=0.5) -- the published rung.
# This is a yes/no stability screen, not a measurement: 128^3 cells run ~8x
# faster than the N=256 design grid and four cells finish in a morning, and
# because dx=0.5 is the resolution every published cell ran at, a result
# here needs no resolution caveat -- the omega=0.55 collapse was seen at
# exactly this grid.  Sponge at the default 24/32 band (sized for L=64);
# these stars' tails end at r ~ 14 from the origin, well inside it, so the
# t=120 tail is protected from the trapped massive-scalar bath without the
# sponge ever touching the star.
# Core tracker on (BONDI_SCRUTINY=1): the domain-integrated barycentre
# rms is halo-contaminated on long runs -- the old phantom's rms grew
# 5.43 -> 17.6 by t=100 while its core never moved 1%.
#
# PASS gate, per cell, over the FULL t=0..120:
#   - core peak amplitude (sector_dynamics.dat col 4 canonical / col 8
#     phantom, 0-indexed) flat to +-2%
#   - min lapse (data/collapse_diagnostics.dat col 1) steady near 1 -- no
#     downward march.  The omega=0.55 canonical went 0.976 -> 0.867 by t=40.
#   - core position (cols 1-3 / 5-7) parked at x=+4 to within a cell (0.5)
# FAIL of the canonical cell at BOTH frequencies = the whole rung
# (lambda=10240, mu=21845333) is dead for pairs, go up the frequency ladder
# or change rung before spending any more GPU time.
#
# Usage:
#   bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_stable_campaign.sh
#   BONDI_STABLE_DRYRUN=1 bash ...        # write the jobs, start nothing
#   BONDI_STABLE_CELLS="single_p_w080" .. # a subset
#   BONDI_STABLE_GPUS="0 1 2" bash ...    # fewer cards
#
# Stop with scripts/campaigns/stop_campaign.sh, or by hand: touch the queue's
# STOP sentinel first so no further cell is dispatched, then kill the runner
# (runner.pid) and sweep the workers.  Killing a running cell on its own only
# frees its slot for the next one.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${REPO_ROOT}"

# New top-level tree, deliberately NOT under runs/bondi/rerun: that tree is
# the closed omega=0.55 archive.  Everything stable-branch lives here.
ROOT="${BONDI_STABLE_ROOT:-${REPO_ROOT}/runs/bondi/correct}"
QDIR="${ROOT}/_queue"
GPUS="${BONDI_STABLE_GPUS:-0 1 2 3}"
DRYRUN="${BONDI_STABLE_DRYRUN:-0}"
QUEUE_RUNNER="${SCRIPT_DIR}/../lib/gpu_queue.sh"
LAUNCHER="${SCRIPT_DIR}/run_single_selfgrav.sh"

ALL_CELLS="single_p_w080 single_m_w080 single_p_w075 single_m_w075"
CELLS="${BONDI_STABLE_CELLS:-${ALL_CELLS}}"

# Shared by every cell: the prototype grid (see header), t=120, sponge on
# (default 24/32 band, sized for L=64), core tracker on.  Solve tolerance
# stays at the published default -- a tightened exit buys nothing on a
# yes/no screen.
COMMON="BONDI_STOP_TIME=120 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 BONDI_SCRUTINY=1 BONDI_SPONGE=1"

cell_env() {
  case "$1" in
    single_p_w080) echo "BONDI_EXOTIC=0 BONDI_OMEGA=0.80" ;;
    single_m_w080) echo "BONDI_EXOTIC=1 BONDI_OMEGA=0.80" ;;
    single_p_w075) echo "BONDI_EXOTIC=0 BONDI_OMEGA=0.75" ;;
    single_m_w075) echo "BONDI_EXOTIC=1 BONDI_OMEGA=0.75" ;;
    *) echo "[stable] unknown cell '$1'" >&2; return 1 ;;
  esac
}

mkdir -p "${QDIR}/pending" "${QDIR}/logs"
rm -f "${QDIR}/STOP"

idx=0
for cell in ${CELLS}; do
  idx=$(( idx + 10 ))
  tag="$(printf '%03d_%s' "${idx}" "${cell}")"
  runs_dir="${ROOT}/${cell}"
  job="${QDIR}/pending/${tag}.job"

  if [[ -d "${runs_dir}" ]]; then
    echo "[stable] ${cell}: run dir already exists -- not enqueued"
    continue
  fi
  for state in running done failed; do
    if compgen -G "${QDIR}/${state}/${tag}.job*" > /dev/null; then
      echo "[stable] ${cell}: already ${state} in the queue -- not enqueued"
      continue 2
    fi
  done

  # The job runs in the FOREGROUND: the queue runner is the thing that is
  # detached, once, and it owns the slot for this cell's whole solve+evolve.
  {
    printf 'cd "%s"\n' "${REPO_ROOT}"
    printf 'BONDI_GPU="${QUEUE_GPU}" %s %s \\\n' "${COMMON}" "$(cell_env "${cell}")"
    printf '  BONDI_RUNS_DIR="%s" \\\n' "${runs_dir}"
    printf '  bash "%s"\n' "${LAUNCHER}"
    printf 'rc=$?\n'
    # Drop the constraint solve once the evolution has consumed it
    # (~0.5 GB per N=128 cell); regenerable from the cell's own params.txt.
    # Only on success, so a cell that died mid-evolution keeps its data.
    printf 'if [ "$rc" -eq 0 ]; then find "%s" -name initial_data.gridinit -delete; fi\n' "${runs_dir}"
    printf 'exit "$rc"\n'
  } > "${job}"
  echo "[stable] enqueued ${tag} -> ${runs_dir}"
done

if [[ "${DRYRUN}" != "0" ]]; then
  echo "[stable] dry run -- queue written to ${QDIR}/pending, runner not started"
  exit 0
fi

if [[ -f "${QDIR}/runner.pid" ]] && kill -0 "$(cat "${QDIR}/runner.pid")" 2>/dev/null; then
  echo "[stable] runner already up (pid $(cat "${QDIR}/runner.pid")); it will pick the new jobs up"
  exit 0
fi

# `/usr/bin/env` spelled out on purpose: on this node other users' bin dirs
# precede /usr/bin on PATH and each holds an `env` that is really a PATH-setup
# snippet meant to be sourced.  Run as `env VAR=x cmd` it prepends its own
# directory to PATH and exits 0 WITHOUT executing cmd -- so a bare `env` here
# launches nothing while looking like a success.
setsid nohup /usr/bin/env QUEUE_IDLE_EXIT=1 \
  bash "${QUEUE_RUNNER}" "${QDIR}" ${GPUS} \
  > "${ROOT}/runner.log" 2>&1 < /dev/null &
disown

echo "[stable] runner starting on GPUs: ${GPUS}"
echo "[stable] verify with: pgrep -af gpu_queue.sh   (a silent no-op looks like success)"
echo "[stable] queue events: ${QDIR}/queue.log ; per-cell logs: ${QDIR}/logs/"
