#!/usr/bin/env bash
# Bondi dipole -- IS THE CANONICAL COLLAPSE A BIRTH ARTEFACT? (docs/MatterDebugg.md)
#
# THE FINDING THIS CAMPAIGN TESTS
# Every canonical cell ever run sinks: min lapse peaks near t=20 and falls
# monotonically to a horizon (t~103 at omega=0.75, t~108 at omega=0.80).  Every
# phantom cell sits perfectly still.  Three explanations have already been
# killed by measurement:
#   - dirty initial data: rebuilt 100-195x cleaner at birth (matched solve grid)
#     and the decline was unchanged to within half a percent
#   - resolution: N=192 tracks N=128 to ~0.1%
#   - boundary/pump: net boundary flux ~1e-6, and rl_pump_stop_time=0 means the
#     matter pump is off for the entire evolution (gains zeroed, num_sites=0)
#
# What is left is an asymmetry in HOW THE TWO SECTORS ARE BUILT.  GRTresna sets
# the initial expansion rate from the CTTK ansatz
#     K = sign_of_K * sqrt(24 pi G rho),
# which is imaginary wherever rho < 0.  The wrapper therefore switches to the
# K=0 York/Lichnerowicz path if and only if exotic matter is present
# (grtresna/fit/motif.py: maximal_slicing = has_exotic or exotic_needed).
# Measured consequence at t=0.01, same rung, same grid, same box:
#     canonical  max|K| = 1.0e-01      (maximal_slicing = 0)
#     phantom    max|K| = 3.0e-05      (maximal_slicing = 1)
# Four orders of magnitude, for a plumbing reason rather than a physical one.
# The canonical star is handed a slice that is already collapsing at ~10% per
# time unit -- about 12% of its own field frequency -- and it is the only one
# that dies.  It cannot push back: the core amplitude is pinned to the
# potential's preferred height (0.01875) at every frequency, so the star has no
# way to stiffen in response.
#
# THE TEST
# Rebuild the canonical star with maximal_slicing=1 -- K=0, matter energy
# carried by the elliptic psi-solve, exactly the phantom's construction.  The
# two sectors then differ ONLY in the sign of the energy.
#
#   flat  -> the "collapse" was an artefact of the initial-data method and this
#            rung is usable for the pair campaign after all
#   sinks -> the birth kick is not the cause; the rung is genuinely unstable to
#            compression and the mass dial (or a new rung) is the next move
#
# Cells 010/020 are direct A/B partners of the already-completed K!=0 runs
# (single_p_w080_n128_matched, single_p_w075_n128_matched) -- same grid, same
# solve box, same tolerance, one flag different.  Cells 030/040 walk up the
# frequency ladder so that a partial cure is still interpretable: higher omega
# is lighter and fluffier (ADM 0.0112 -> 0.0080, r90 5.11 -> 6.85 from 0.80 to
# 0.90), so if K=0 helps but does not fully cure, the ladder shows where the
# stable branch begins.  All four are canonical: the phantom side is already
# proven stable to t=120 and needs no further cells.
#
# Grid: N=128 at L=64 (dx=0.5), the published rung, with the MATCHED solve box
# (L=128/N=256 -> identical 0.5 cell, straight copy instead of interpolation)
# and MAXLEVEL=0 on both sides.  t=120 so a horizon would have formed by ~t=105
# if the star still behaves as before.
#
# PASS gate, per cell, over the full t=0..120:
#   - min lapse steady near 1, no downward march (the K!=0 runs go 0.979 -> 0.21)
#   - core peak amplitude flat to +-2%
#   - core parked at x=+4 to within one cell (0.5)
#
# Usage:
#   bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_slicing_test.sh
#   BONDI_SLICING_DRYRUN=1 bash ...          # write the jobs, start nothing
#   BONDI_SLICING_CELLS="p_w080_k0" bash ... # a subset
#   BONDI_SLICING_GPUS="0 1" bash ...        # fewer cards
#
# Stop with scripts/campaigns/stop_campaign.sh, or by hand: touch the queue's
# STOP sentinel first so no further cell is dispatched, then kill the runner
# (runner.pid) and sweep the workers.  Killing a running cell on its own only
# frees its slot for the next one.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${REPO_ROOT}"

ROOT="${BONDI_SLICING_ROOT:-${REPO_ROOT}/runs/bondi/correct}"
# Own queue dir: the sibling _queue holds the completed K!=0 screen's
# done/failed records and a STOP sentinel.
QDIR="${ROOT}/_queue_k0"
GPUS="${BONDI_SLICING_GPUS:-0 1 2 3}"
DRYRUN="${BONDI_SLICING_DRYRUN:-0}"
QUEUE_RUNNER="${SCRIPT_DIR}/../lib/gpu_queue.sh"
LAUNCHER="${SCRIPT_DIR}/run_single_selfgrav.sh"

ALL_CELLS="p_w080_k0 p_w075_k0 p_w085_k0 p_w090_k0"
CELLS="${BONDI_SLICING_CELLS:-${ALL_CELLS}}"

# Shared by every cell.  The matched solve grid (DOMAIN_L=128 with N=256 gives
# the evolution's own 0.5 cell) and the tightened 0.1% exit are carried over
# from the matched-grid screen so the ONLY difference from those cells is the
# slicing flag.
COMMON="BONDI_STOP_TIME=120 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_EXOTIC=0 \
BONDI_GRTRESNA_DOMAIN_L=128 BONDI_GRTRESNA_N=256 BONDI_GRTRESNA_MAXLEVEL=0 \
BONDI_GRTRESNA_MAXIMAL_SLICING=1 \
BONDI_GRTRESNA_RANKS=32 BONDI_NL_TOL=0.001 BONDI_GRTRESNA_TIMEOUT=21600"

cell_env() {
  case "$1" in
    p_w080_k0) echo "BONDI_OMEGA=0.80" ;;
    p_w075_k0) echo "BONDI_OMEGA=0.75" ;;
    p_w085_k0) echo "BONDI_OMEGA=0.85" ;;
    p_w090_k0) echo "BONDI_OMEGA=0.90" ;;
    *) echo "[slicing] unknown cell '$1'" >&2; return 1 ;;
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
    echo "[slicing] ${cell}: run dir already exists -- not enqueued"
    continue
  fi
  for state in running done failed; do
    if compgen -G "${QDIR}/${state}/${tag}.job*" > /dev/null; then
      echo "[slicing] ${cell}: already ${state} in the queue -- not enqueued"
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
    # Drop the constraint solve once the evolution has consumed it (~0.5 GB per
    # cell); regenerable from the cell's own params.txt.  Only on success, so a
    # cell that died mid-evolution keeps its data.
    printf 'if [ "$rc" -eq 0 ]; then find "%s" -name initial_data.gridinit -delete; fi\n' "${runs_dir}"
    printf 'exit "$rc"\n'
  } > "${job}"
  echo "[slicing] enqueued ${tag} -> ${runs_dir}"
done

if [[ "${DRYRUN}" != "0" ]]; then
  echo "[slicing] dry run -- queue written to ${QDIR}/pending, runner not started"
  exit 0
fi

if [[ -f "${QDIR}/runner.pid" ]] && kill -0 "$(cat "${QDIR}/runner.pid")" 2>/dev/null; then
  echo "[slicing] runner already up (pid $(cat "${QDIR}/runner.pid")); it will pick the new jobs up"
  exit 0
fi

# `/usr/bin/env` spelled out on purpose: on this node other users' bin dirs
# precede /usr/bin on PATH and each holds an `env` that is really a PATH-setup
# snippet meant to be sourced.  Run as `env VAR=x cmd` it prepends its own
# directory to PATH and exits 0 WITHOUT executing cmd -- so a bare `env` here
# launches nothing while looking like a success.
setsid nohup /usr/bin/env QUEUE_IDLE_EXIT=1 \
  bash "${QUEUE_RUNNER}" "${QDIR}" ${GPUS} \
  > "${ROOT}/runner_k0.log" 2>&1 < /dev/null &
disown

echo "[slicing] runner starting on GPUs: ${GPUS}"
echo "[slicing] verify with: pgrep -af gpu_queue.sh   (a silent no-op looks like success)"
echo "[slicing] queue events: ${QDIR}/queue.log ; per-cell logs: ${QDIR}/logs/"
