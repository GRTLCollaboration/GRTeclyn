#!/usr/bin/env bash
# Bondi dipole -- DOES FIXING THE STRESS-TENSOR TRACE STOP THE CANONICAL SINK?
#
# THE BUG THIS CAMPAIGN TESTS (found 2026-08-21, Source/Matter/*.impl.hpp)
# Every matter class built S_ij with the potential already in it,
#     S_ij = d_i phi d_j phi - gamma_ij (Vt/2 + V),
# and then computed the trace as
#     trS = chi * trace(S_ij, h_UU) - 3 V.
# The trace of that S_ij ALREADY contains the -3V, so the potential was
# subtracted twice.  Checked against the covariant T_ab = d_a phi d_b phi -
# g_ab (d phi^2 / 2 + V): rho, j_i and S_ij all come out exact, and trS comes
# out short by exactly -3V.  In the pure vacuum-energy limit (phi constant,
# Pi = 0, so T_ab = -V g_ab) the code returns trS = -6V where the answer is
# -3V -- a clean factor of two.
#
# WHY IT WENT UNSEEN.  trS is consumed in exactly one place,
# CCZ4RHSWithMatter::add_emtensor_rhs -> rhs[c_K]; the A_ij equation re-derives
# its own trace from S_ij, and the Hamiltonian and momentum constraints read
# only rho and j_i.  So the error is invisible in the constraint monitor at
# t=0 and corrupts nothing but the evolution of K.  And every test in Tests/
# instantiates ScalarField<DefaultPotential>, whose V is identically zero --
# the whole suite multiplies the broken term by nought.  A star held together
# by its potential is precisely the case the suite never covers.
#
# WHY THE PHANTOM SECTOR NEVER SHOWED IT.  The sign of the spurious term
# follows the sign of the energy, so canonical and phantom are pushed opposite
# ways -- but that is not why only one sector died.  At omega = 0.80 the star's
# compactness is M/R = 2.3e-03, roughly 1/200 of anything that can form a
# horizon, and the phantom twin is the same size with M = -0.0117.  Negative
# mass gravitates repulsively: there is no collapsed end state for it to fall
# into, so any spurious source disperses and the error saturates.  The phantom
# is the direction with no attractor.  Its stability was never evidence that
# the coupling was right.
#
# THE TEST
# Re-evolve the four K=0 cells from the SAME solved initial data -- the exact
# initial_data.gridinit bytes the killed runs used -- against the rebuilt
# binary.  The elliptic solve is skipped, so the ONLY difference from the runs
# these are compared against is the trS fix.
#
#   flat            -> the double-counted potential was the sink, and this rung
#                      is usable for the pair campaign
#   sinks unchanged -> trS was a real bug but not this one; the remaining
#                      suspect is elsewhere in the geometry sector
#
# WHAT THE OLD RUNS DID, for scoring (min lapse, canonical, K=0 slice):
#   t     20      40      60      80      100     110
#   w075  0.970   0.932   0.876   0.747   ~0.30   horizon
#   w080  0.974   0.943   0.902   0.827   0.578   0.268
#   w085  --      --      --      0.864   --      --
#   w090  --      --      --      0.881   --      --
# Core peak amplitude was FLAT in all of them (0.02419 -> 0.02478, +2.4%) and
# the core never left x=+4, while min chi fell to 0.088.  The matter was fine
# throughout; the geometry was manufacturing several hundred times the star's
# own mass.
#
# PASS gate, per cell, over the full t=0..120:
#   - min lapse steady near 1, no downward march
#   - core peak amplitude flat to +-2%
#   - core parked at x=+4 to within one cell (0.5)
#   - Linf Ham no longer climbing past the physical source scale (16 pi rho
#     ~ 4e-2); the old runs reached 0.24, six times their own matter.
#
# Usage:
#   bash grteclyn-wrapper/scripts/campaigns/bondi_dipole/run_trsfix_ab.sh
#   BONDI_TRSFIX_DRYRUN=1 bash ...          # write the jobs, start nothing
#   BONDI_TRSFIX_CELLS="p_w080_k0" bash ... # a subset
#   BONDI_TRSFIX_GPUS="0 1" bash ...        # fewer cards
#
# Stop with scripts/campaigns/stop_campaign.sh, or by hand: touch the queue's
# STOP sentinel first so no further cell is dispatched, then kill the runner
# (runner.pid) and sweep the workers.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${REPO_ROOT}"

ROOT="${BONDI_TRSFIX_ROOT:-${REPO_ROOT}/runs/bondi_trsfix}"
# Where the already-solved slices live: the cells this campaign is the A/B of.
SRC_ROOT="${BONDI_TRSFIX_SRC:-${REPO_ROOT}/runs/bondi_correct}"
QDIR="${ROOT}/_queue"
GPUS="${BONDI_TRSFIX_GPUS:-0 1 2 3}"
DRYRUN="${BONDI_TRSFIX_DRYRUN:-0}"
QUEUE_RUNNER="${SCRIPT_DIR}/../lib/gpu_queue.sh"
LAUNCHER="${SCRIPT_DIR}/run_single_selfgrav.sh"

ALL_CELLS="p_w080_k0 p_w075_k0 p_w085_k0 p_w090_k0"
CELLS="${BONDI_TRSFIX_CELLS:-${ALL_CELLS}}"

# Evolution knobs only.  Every BONDI_GRTRESNA_* value is deliberately absent:
# the solve is skipped, so quoting solve settings here would only invite the
# reader to believe they still applied.
COMMON="BONDI_STOP_TIME=120 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_EXOTIC=0"

cell_omega() {
  case "$1" in
    p_w080_k0) echo "0.80" ;;
    p_w075_k0) echo "0.75" ;;
    p_w085_k0) echo "0.85" ;;
    p_w090_k0) echo "0.90" ;;
    *) echo "[trsfix] unknown cell '$1'" >&2; return 1 ;;
  esac
}

# The solved slice for a cell, inside the run that produced it.  The launcher
# names a run bondi_sg_single_<tag>, and the K=0 cells were tagged with the
# cell name minus its _k0 suffix.
cell_gridinit() {
  echo "${SRC_ROOT}/$1/bondi_sg_single_${1%_k0}/initial_data.gridinit"
}

mkdir -p "${QDIR}/pending" "${QDIR}/logs"
rm -f "${QDIR}/STOP"

idx=0
for cell in ${CELLS}; do
  idx=$(( idx + 10 ))
  tag="$(printf '%03d_%s' "${idx}" "${cell}")"
  runs_dir="${ROOT}/${cell}"
  job="${QDIR}/pending/${tag}.job"
  omega="$(cell_omega "${cell}")"
  gridinit="$(cell_gridinit "${cell}")"

  if [[ ! -f "${gridinit}" ]]; then
    echo "[trsfix] ${cell}: no solved slice at ${gridinit} -- not enqueued" >&2
    continue
  fi
  if [[ -d "${runs_dir}" ]]; then
    echo "[trsfix] ${cell}: run dir already exists -- not enqueued"
    continue
  fi
  for state in running done failed; do
    if compgen -G "${QDIR}/${state}/${tag}.job*" > /dev/null; then
      echo "[trsfix] ${cell}: already ${state} in the queue -- not enqueued"
      continue 2
    fi
  done

  {
    printf 'cd "%s"\n' "${REPO_ROOT}"
    printf 'BONDI_GPU="${QUEUE_GPU}" %s BONDI_OMEGA=%s \\\n' "${COMMON}" "${omega}"
    printf '  BONDI_GRIDINIT="%s" \\\n' "${gridinit}"
    printf '  BONDI_RUNS_DIR="%s" \\\n' "${runs_dir}"
    printf '  bash "%s"\n' "${LAUNCHER}"
    # No gridinit sweep here: the slice belongs to the cell it was solved in,
    # which is the A/B partner.  Deleting it would destroy the comparison.
  } > "${job}"
  echo "[trsfix] enqueued ${tag} -> ${runs_dir}  (reusing the slice solved in ${cell})"
done

if [[ "${DRYRUN}" != "0" ]]; then
  echo "[trsfix] dry run -- queue written to ${QDIR}/pending, runner not started"
  exit 0
fi

if [[ -f "${QDIR}/runner.pid" ]] && kill -0 "$(cat "${QDIR}/runner.pid")" 2>/dev/null; then
  echo "[trsfix] runner already up (pid $(cat "${QDIR}/runner.pid")); it will pick the new jobs up"
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

echo "[trsfix] runner starting on GPUs: ${GPUS}"
echo "[trsfix] verify with: pgrep -af gpu_queue.sh   (a silent no-op looks like success)"
echo "[trsfix] queue events: ${QDIR}/queue.log ; per-cell logs: ${QDIR}/logs/"
