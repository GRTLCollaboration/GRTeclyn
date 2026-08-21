#!/usr/bin/env bash
# Bondi dipole runaway -- the paper's main result, on lumps that are actually
# stable.
#
# WHY NOW.  Every earlier pair cell was evolved with the double-counted
# potential in trS (Source/Matter/*.impl.hpp, fixed 2026-08-21) and on
# omega=0.55 stars that sit on the unstable branch.  Both are gone: the trS
# term now reads
#     out.trS = chi * compute_trace(out.S, h_UU)
# and this campaign is built on the omega that survived the fix cleanest.
#
# WHICH STAR.  The trs-fix A/B ran omega = 0.75 / 0.80 / 0.85 / 0.90 as
# isolated canonical stars, unigrid N=128, L=64, to t=120.  omega=0.75 won on
# every axis at once:
#   - flattest lapse of the four (drift +1e-05 over t=76, vs -9e-04 at 0.90)
#   - smallest amplitude creep (+0.6%, vs +11.9% at 0.90)
#   - most compact (r90 = 5.24), so a close pair still has distinct cores
#   - heaviest (M = +0.014350), so the runaway acceleration is the largest
#     available on the stable branch
# It is the most stable AND the loudest signal, which does not usually happen.
#
# GEOMETRY.  Separation 10 along x: canonical lump at x=+5, phantom at x=-5.
# r90 = 5.24 each, so the 90%-mass surfaces just kiss -- they cross by 0.48,
# under 5% of a radius.  That is deliberate: the runaway acceleration goes as
# 1/d^2, so pulling in from 12 to 10 buys 44% more signal, and 10.5 is where
# the cores would actually begin to share matter.  Closer than this and the
# "two clean stable stars" premise gets hard to defend; the haloes (r99 = 8.99)
# overlap either way, as they always do for a soliton.
#
# THE MASS MATCH.  A textbook Bondi runaway needs |M-| = |M+| exactly, so the
# pair carries zero total mass and there is no ordinary two-body force mixed
# into the drift.  The phantom twin at the SAME omega is 5.4% heavier
# (-0.015131 vs +0.014350), which would leave a net -0.00078 of mass sitting
# in the box.  Interpolating the phantom branch between omega 0.750 (-0.015131)
# and 0.775 (-0.013226) puts |M-| = 0.014350 at omega = 0.7603.  Hence two
# runaway cells, not one:
#   pm_eqm  -- phantom retuned to 0.7603, masses equal and opposite.  PRIMARY.
#   pm      -- both lumps at 0.75, 5.4% mass excess on the phantom.  Shows the
#              runaway is not an artefact of the retuning.
#
# THE CONTROLS.  A drift measurement on a code that was manufacturing mass
# three days ago is worth nothing without cells that must NOT drift:
#   pp -- two canonical stars.  Must attract; centre of mass must stay put.
#   mm -- two phantom stars.  Must repel; centre of mass must stay put.
# Any centre-of-mass motion in pp or mm is numerical, and its size is the error
# bar on the pm result.
#
# pp NEEDS FORCED MAXIMAL SLICING, and the first attempt at this campaign is why.
# GRTresna switches the K=0 York solve on automatically ONLY when an exotic lump
# is present, because the CTTK ansatz K = sign*sqrt(24 pi G rho) is imaginary for
# rho<0.  So pm / pm_eqm / mm are all born at K=0 while pp -- the only cell with
# no phantom in it -- silently gets K = 1e-01, i.e. a slice that is already
# collapsing before a single step is taken.  Measured on the first attempt: pp
# ran max|K| = 3.1e-02 with the lapse swinging 0.99 -> 0.73 -> 0.88 inside t=6,
# against max|K| = 1.2e-04 and a flat lapse in pm.  Its centre of mass had
# wandered -0.025 by t=8, eight times the real runaway signal at that time.  The
# control was noisier than the measurement, for a reason that has nothing to do
# with the physics under test.  BONDI_GRTRESNA_MAXIMAL_SLICING=1 removes it.  This is the check the phantom-vs-canonical comparison
# could never provide -- a lone negative mass has no state to fall into, so it
# looks stable no matter what is broken.
#
# WHAT TO EXPECT.  Newtonian estimate for the equal-mass pair:
#     a = G M / d^2 = 0.014350 / 100 = 1.435e-04
#     displacement = a t^2 / 2  ->  0.72 at t=100, 2.87 at t=200
# Both lumps drift the same way, +x (phantom chases canonical), separation held.
# The isolated-star runs put the core-position noise floor at 2e-04 grid units
# over t=73, so a 0.5 displacement is already four orders of magnitude clear of
# it.  Stopping at t=200 rather than t=120 buys a 4x larger signal for 1.7x the
# wall clock, and the star only reaches x~8 by then -- still 16 clear of the
# sponge.
#
# SOLVE PRECISION.  The launcher's default exit gate (0.1% Hamiltonian) is not
# good enough here, and the first attempt at this campaign proved it.  The two
# sectors converge by different paths: a pure-canonical cell alternates Ham and
# Mom and falls straight through the gate to 0.002%, while any cell carrying a
# phantom lump converges geometrically at ~0.4 per iteration and stops the
# moment it crosses the gate -- at 0.082%, still improving.  That left pp with
# initial data 40x cleaner than pm/pm_eqm/mm, i.e. the CONTROL better resolved
# than the measurement, which makes the comparison worthless.  The gate below
# is set where the canonical path lands anyway, so all four cells enter the
# evolution at the same quality.  It costs ~4 extra Newton iterations (7 -> 11
# of the 50 permitted), a couple of CPU minutes.  The stall guard has to come
# down with the gate or the iteration halts on "no further progress" first.
#
# READING THE RESULT.  small_data/sector_dynamics.dat, columns 2-4 (canonical
# core xyz) and the phantom block; these are halo-free core positions, immune to
# the sponge eating the outer halo, which the domain-integrated barycentres in
# sector_barycenters.dat are not.  PASS, per cell, over t = 0..200:
#   pm / pm_eqm -- both cores drift +x together, separation held to ~1 cell,
#                  displacement growing as t^2 not t, cores intact (peak
#                  amplitude flat to +-2%), min lapse steady near its birth
#                  value with no downward march
#   pp          -- cores approach each other, centre of mass fixed to <0.01
#   mm          -- cores separate, centre of mass fixed to <0.01
#
# Usage:
#   bash scripts/campaigns/bondi_dipole/run_runaway_stable.sh
# Overrides:
#   BONDI_RUNAWAY_ROOT   run root (default runs/bondi/runaway)
#   BONDI_RUNAWAY_SEP    separation along x (default 10)
#   BONDI_RUNAWAY_CELLS  subset of "pm_eqm pm pp mm"
#   BONDI_RUNAWAY_GPUS   default "0 1 2 3"
#   BONDI_RUNAWAY_DRYRUN 1 = write the queue, do not start the runner
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/../../../.." && pwd)"
cd "${REPO_ROOT}"

ROOT="${BONDI_RUNAWAY_ROOT:-${REPO_ROOT}/runs/bondi/runaway}"
QDIR="${ROOT}/_queue"
GPUS="${BONDI_RUNAWAY_GPUS:-0 1 2 3}"
DRYRUN="${BONDI_RUNAWAY_DRYRUN:-0}"
QUEUE_RUNNER="${SCRIPT_DIR}/../lib/gpu_queue.sh"
LAUNCHER="${SCRIPT_DIR}/run_pair_selfgrav.sh"

ALL_CELLS="pm_eqm pm pp mm"
CELLS="${BONDI_RUNAWAY_CELLS:-${ALL_CELLS}}"

# The stable rung, and the separation that keeps the cores just clear.
OMEGA="${BONDI_RUNAWAY_OMEGA:-0.75}"
OMEGA_PHANTOM_EQM="${BONDI_RUNAWAY_OMEGA_EQM:-0.7603}"
SEP="${BONDI_RUNAWAY_SEP:-10}"

# Elliptic-solve exit gate, in percent -- see SOLVE PRECISION above.  Both
# numbers move together: the stall guard scales with the gate.
NL_TOL="${BONDI_RUNAWAY_NL_TOL:-0.002}"
NL_STALL_TOL="${BONDI_RUNAWAY_NL_STALL_TOL:-0.00004}"

# Grid and evolution.  Unigrid on purpose: an AMR level switching on partway
# through moves the core-position diagnostic, which IS the measurement here.
# Plot cadence 80 matches the singles -- at 40 the consumer falls behind a
# 200-time-unit run and the scratch dir grows without bound.
COMMON="BONDI_STOP_TIME=200 BONDI_NFULL=128 BONDI_LFULL=64 BONDI_MAXLEVEL=0 \
BONDI_PLOT_INTERVAL=80 BONDI_SCRUTINY=1 BONDI_SPONGE=1 BONDI_SEP=${SEP} \
BONDI_GRTRESNA_RANKS=8 BONDI_NL_TOL=${NL_TOL} BONDI_NL_STALL_TOL=${NL_STALL_TOL}"

# sectors: 0 = canonical, 1 = phantom.  s0 sits at +x, s1 at -x.
cell_s0() { case "$1" in pm_eqm|pm|pp) echo 0 ;; mm) echo 1 ;; esac; }
cell_s1() { case "$1" in pm_eqm|pm|mm) echo 1 ;; pp) echo 0 ;; esac; }
# Only the equal-mass cell retunes lump1; everywhere else both lumps share omega.
cell_s1_omega() { case "$1" in pm_eqm) echo "${OMEGA_PHANTOM_EQM}" ;; *) echo "" ;; esac; }
# pp is the one cell with no exotic lump, so it is the one cell that does not get
# maximal slicing for free -- see THE CONTROLS above.
cell_extra() { case "$1" in pp) echo "BONDI_GRTRESNA_MAXIMAL_SLICING=1" ;; *) echo "" ;; esac; }

for cell in ${CELLS}; do
  case " ${ALL_CELLS} " in
    *" ${cell} "*) ;;
    *) echo "[runaway] unknown cell '${cell}' -- expected one of: ${ALL_CELLS}" >&2; exit 1 ;;
  esac
done

mkdir -p "${QDIR}/pending" "${QDIR}/logs"
rm -f "${QDIR}/STOP"

idx=0
for cell in ${CELLS}; do
  idx=$(( idx + 10 ))
  tag="$(printf '%03d_%s' "${idx}" "${cell}")"
  runs_dir="${ROOT}/${cell}"
  job="${QDIR}/pending/${tag}.job"
  s0="$(cell_s0 "${cell}")"
  s1="$(cell_s1 "${cell}")"
  s1_omega="$(cell_s1_omega "${cell}")"
  extra="$(cell_extra "${cell}")"

  if [[ -d "${runs_dir}" ]]; then
    echo "[runaway] ${cell}: run dir already exists -- not enqueued"
    continue
  fi
  for state in running done failed; do
    if compgen -G "${QDIR}/${state}/${tag}.job*" > /dev/null; then
      echo "[runaway] ${cell}: already ${state} in the queue -- not enqueued"
      continue 2
    fi
  done

  {
    printf 'cd "%s"\n' "${REPO_ROOT}"
    printf 'BONDI_GPU="${QUEUE_GPU}" %s %s\\\n' "${COMMON}" "${extra:+${extra} }"
    printf '  BONDI_S0=%s BONDI_S1=%s BONDI_S0_OMEGA=%s \\\n' "${s0}" "${s1}" "${OMEGA}"
    if [[ -n "${s1_omega}" ]]; then
      printf '  BONDI_S1_OMEGA=%s \\\n' "${s1_omega}"
    fi
    printf '  BONDI_RUNS_DIR="%s" \\\n' "${runs_dir}"
    printf '  bash "%s"\n' "${LAUNCHER}"
  } > "${job}"
  echo "[runaway] enqueued ${tag}  s0=${s0} s1=${s1} omega=${OMEGA}${s1_omega:+ / ${s1_omega}} sep=${SEP}${extra:+ ${extra}} -> ${runs_dir}"
done

if [[ "${DRYRUN}" != "0" ]]; then
  echo "[runaway] dry run -- queue written to ${QDIR}/pending, runner not started"
  exit 0
fi

if [[ -f "${QDIR}/runner.pid" ]] && kill -0 "$(cat "${QDIR}/runner.pid")" 2>/dev/null; then
  echo "[runaway] runner already up (pid $(cat "${QDIR}/runner.pid")); it will pick the new jobs up"
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

echo "[runaway] runner starting on GPUs: ${GPUS}"
echo "[runaway] verify with: pgrep -af gpu_queue.sh   (a silent no-op looks like success)"
echo "[runaway] queue events: ${QDIR}/queue.log ; per-cell logs: ${QDIR}/logs/"
