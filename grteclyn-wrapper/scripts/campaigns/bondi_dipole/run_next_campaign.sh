#!/usr/bin/env bash
# Bondi dipole -- launcher for the follow-up campaign in
# research/bondi_dipole/docs/MatterDebugg.md.  Six cells that answer the open items the
# published pack cannot, enqueued onto the shared work-queue runner
# (../lib/gpu_queue.sh) so whichever card frees first takes the next cell.
#
#   items 1+2  tolB_n128 / tolB_n192 / tolB_n256
#       The d0=12 reference ladder re-solved with the elliptic exit tolerance
#       scaled as dx^4, and with the halo-free peak tracker switched on.  Every
#       rung of the published ladder exited on the SAME tolerance and so floored
#       at the same 0.0832% Hamiltonian violation, identical to four digits --
#       which is why the pack quotes no convergence order.  Scaling the
#       tolerance is what turns that caveat into a number.  The same three cells
#       carry the tracker, so the halo bias is measured at the separation the
#       article actually quotes instead of carried over from d0=8.
#   item 3a    spongeA_boxC
#       boxC_pm_L128_n256 repeated with the boundary sponge on -- the validation
#       against its existing sponge-off twin.  Whether a sponge tuned against
#       massless reflections also clears this campaign's massive scalar bath is
#       an empirical question, and this cell is the experiment.
#   item 3b    spongeB_sep16_long
#       d0=16 out to t=120: the only configuration whose gap stays open long
#       enough to test the textbook constant-gap runaway, and it cannot get
#       past t=60 without the sponge.
#   item 5     lever_omega0615
#       The phantom at omega=0.615, widening the force-scaling lever arm from
#       |M-|/M+ = 0.83 to 0.62.  Measured against convA_pm_n256, already packed.
#
# One cell per card: two 192^3 cells sharing an H100 each ran 2.2x slower, so
# the pair finished later than running them back to back.  One queue slot per
# GPU enforces that.
#
# NOT prestaged, deliberately.  replay_eval.py offers --solve-only/--gridinit to
# run the CPU solve off the card, and for these cells it would be WRONG: the
# --gridinit path re-reads the matter parameters from the frozen champion eval
# rather than from the solve that produced the gridinit, and that champion is a
# five-lump all-phantom configuration at scalar_lambda=640 / mu=85333 against
# this campaign's two-lump canonical+phantom at 10240 / 21845333.  The evolution
# would silently run different matter from its own initial data.  Each cell
# therefore does its own solve on its own card, exactly as every published cell
# did -- the card idles 14-70 min per cell, which is the price of comparability.
#
# Usage:
#   bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_next_campaign.sh
#   BONDI_NEXT_DRYRUN=1 bash ...      # write the jobs, print them, start nothing
#   BONDI_NEXT_CELLS="tolB_n128" ...  # a subset
#   BONDI_NEXT_GPUS="0 1 2" ...       # fewer cards
#
# Stop with scripts/campaigns/stop_campaign.sh, or by hand: touch the queue's
# STOP sentinel first so no further cell is dispatched, then kill the runner
# (runner.pid) and sweep the workers.  Killing a running cell on its own only
# frees its slot for the next one.
# --- SUPERSEDED 2026-08-22 -------------------------------------------------
# This script belongs to a finished campaign and reads/writes runs/bondi/rerun/next/
# which no longer exists.  Running it now would either fail halfway or
# recreate a dead tree and mix superseded cells into the live campaign.
# It is also an auto-chaining launcher, and the campaign no longer works that
# way: cells are launched by hand, one at a time, because an orchestrator
# outlives the session that made it and will relaunch cells behind your back.
# What replaced it: runs/bondi/staging/ -- plan and launch commands in
# research/bondi_dipole/docs/GPU_RUN_PAPER.md.
# Kept for the reasoning in the header above.  Remove this block only if you
# deliberately want the old tree back.
printf '%s: superseded campaign; %s no longer exists.\n' "$(basename "$0")" "runs/bondi/rerun/next/" >&2
printf '  see research/bondi_dipole/docs/GPU_RUN_PAPER.md for the live campaign\n' >&2
exit 2
# ---------------------------------------------------------------------------

set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${REPO_ROOT}"

# Lives under the existing bondi_rerun tree rather than a new top-level dir.
# pack_campaign.sh skips "next" by name -- these are follow-up cells, not part
# of the published pack, and must not be swept into it.
ROOT="${BONDI_NEXT_ROOT:-${REPO_ROOT}/runs/bondi/rerun/next}"
QDIR="${ROOT}/_queue"
GPUS="${BONDI_NEXT_GPUS:-0 1 2 3}"
DRYRUN="${BONDI_NEXT_DRYRUN:-0}"
QUEUE_RUNNER="${SCRIPT_DIR}/../lib/gpu_queue.sh"
LAUNCHER="${SCRIPT_DIR}/run_pair_selfgrav.sh"

# Ordered longest-evolution-first, so the short rungs backfill the cards that
# free up rather than straggling at the end.  Job files dispatch in name order.
ALL_CELLS="tolB_n256 spongeB_sep16_long lever_omega0615 spongeA_boxC tolB_n192 tolB_n128"
CELLS="${BONDI_NEXT_CELLS:-${ALL_CELLS}}"

# Per-cell environment.  Anything not named here keeps the launcher's published
# default, so each cell differs from the packed campaign only in these knobs.
cell_env() {
  case "$1" in
    # --- items 1+2: dx^4 tolerance ladder, tracker on, d0=12 ----------------
    # The stall guard is held at the published 1/50 of the exit tolerance; left
    # at 0.002 it would halt the Newton iteration on "no further progress"
    # before a tightened exit is ever reached.  The solve timeout goes up
    # because a tighter exit costs iterations -- the published solves used 7 of
    # a permitted 50, so there is headroom, but not free headroom.
    tolB_n128) echo "BONDI_SEP=12 BONDI_NFULL=128 BONDI_MAXLEVEL=0 BONDI_SCRUTINY=1 BONDI_NL_TOL=0.1 BONDI_NL_STALL_TOL=0.002" ;;
    tolB_n192) echo "BONDI_SEP=12 BONDI_NFULL=192 BONDI_MAXLEVEL=0 BONDI_SCRUTINY=1 BONDI_NL_TOL=0.019 BONDI_NL_STALL_TOL=0.00038 BONDI_GRTRESNA_TIMEOUT=21600" ;;
    tolB_n256) echo "BONDI_SEP=12 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_SCRUTINY=1 BONDI_NL_TOL=0.00625 BONDI_NL_STALL_TOL=0.000125 BONDI_GRTRESNA_TIMEOUT=21600" ;;

    # --- item 3a: sponge validation against the existing boxC twin ----------
    # Identical to boxC_pm_L128_n256 (L=128, N=256, d0=8, stop 90) except for
    # the sponge.  The band sits at 40/60 rather than the 24/32 default: that
    # default is sized for L=64, and at L=128 it would sit on top of the
    # extraction shells and start eating the canonical star's own halo once its
    # rms radius crosses r=16 (t ~ 51).  The tracker rides along so a halo the
    # sponge removes can be told apart from a halo the domain-integrated
    # barycentre was mis-counting.
    spongeA_boxC) echo "BONDI_SEP=8 BONDI_LFULL=128 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_STOP_TIME=90 BONDI_RADII='16 24 32 40' BONDI_SCRUTINY=1 BONDI_PSI4_HIGHER_L=1 BONDI_SPONGE=1 BONDI_SPONGE_INNER=40 BONDI_SPONGE_OUTER=60" ;;

    # --- item 3b: the constant-gap test ------------------------------------
    spongeB_sep16_long) echo "BONDI_SEP=16 BONDI_LFULL=128 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_STOP_TIME=120 BONDI_RADII='16 24 32 40' BONDI_SCRUTINY=1 BONDI_PSI4_HIGHER_L=1 BONDI_SPONGE=1 BONDI_SPONGE_INNER=40 BONDI_SPONGE_OUTER=60" ;;

    # --- item 5: wider mass lever arm --------------------------------------
    # omega=0.615 interpolates the phantom family to ADM = -0.0397 against the
    # canonical +0.0640, i.e. |M-|/M+ = 0.62, where the published pair sits at
    # 0.83 and the equal-mass control at 1.00.
    lever_omega0615) echo "BONDI_SEP=8 BONDI_NFULL=256 BONDI_MAXLEVEL=0 BONDI_S1_OMEGA=0.615" ;;

    *) echo "[next] unknown cell '$1'" >&2; return 1 ;;
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
    echo "[next] ${cell}: run dir already exists -- not enqueued"
    continue
  fi
  for state in running done failed; do
    if compgen -G "${QDIR}/${state}/${tag}.job*" > /dev/null; then
      echo "[next] ${cell}: already ${state} in the queue -- not enqueued"
      continue 2
    fi
  done

  # The job runs in the FOREGROUND: the queue runner is the thing that is
  # detached, once, and it is what owns the slot for this cell's whole
  # solve+evolve.  QUEUE_GPU is exported by the runner as the slot label.
  {
    printf 'cd "%s"\n' "${REPO_ROOT}"
    printf 'BONDI_S0=0 BONDI_S1=1 BONDI_GPU="${QUEUE_GPU}" %s \\\n' "$(cell_env "${cell}")"
    printf '  BONDI_RUNS_DIR="%s" \\\n' "${runs_dir}"
    printf '  bash "%s"\n' "${LAUNCHER}"
    printf 'rc=$?\n'
    # Drop the constraint solve once the evolution has consumed it: 4.1 GB per
    # N=256 cell, 18.8 GB across this queue, and regenerable from the cell's own
    # params.txt.  Nothing downstream reads it -- each cell solves for itself,
    # and the analysis works off data/ and small_data/.  Only on success, so a
    # cell that died mid-evolution keeps its initial data for inspection.
    printf 'if [ "$rc" -eq 0 ]; then find "%s" -name initial_data.gridinit -delete; fi\n' "${runs_dir}"
    printf 'exit "$rc"\n'
  } > "${job}"
  echo "[next] enqueued ${tag} -> ${runs_dir}"
done

if [[ "${DRYRUN}" != "0" ]]; then
  echo "[next] dry run -- queue written to ${QDIR}/pending, runner not started"
  exit 0
fi

if [[ -f "${QDIR}/runner.pid" ]] && kill -0 "$(cat "${QDIR}/runner.pid")" 2>/dev/null; then
  echo "[next] runner already up (pid $(cat "${QDIR}/runner.pid")); it will pick the new jobs up"
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

echo "[next] runner starting on GPUs: ${GPUS}"
echo "[next] verify with: pgrep -af gpu_queue.sh   (a silent no-op looks like success)"
echo "[next] queue events: ${QDIR}/queue.log ; per-cell logs: ${QDIR}/logs/"
