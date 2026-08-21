#!/usr/bin/env bash
# Pack light, GitHub-friendly extracts of the omega = 0.75 runaway campaign
# (runs/bondi/runaway) into results/bondi-dipole-runaway/campaign/.
#
# WHY A THIRD PACKER
# pack_results.sh and pack_campaign.sh both read runs/bondi/rerun, which is the
# tree built by the binary that double-counted the potential in trS.  Those
# extracts now live under campaign/_corrupted_runs/ and the two old scripts are
# kept only as the provenance record for them.  This script owns the campaign
# run after that fix -- a separate source tree, a separate destination, and no
# shared state with either of them, so refreshing one can never rewrite the
# other.
#
# WHAT IS COPIED
# Per-cell time series (small_data/*.dat), the four evolution diagnostic
# streams downsampled to dt = 0.5, the solve and evolution parameter files, the
# elliptic residual history, and the launch configuration that produced the
# cell.  Deliberately NOT copied: PNG frames, movies, plotfiles and the
# gridinit -- they are large, they carry no number an analysis needs, and they
# stay in the gitignored runs/ tree.
#
# Safe to re-run at any time, including while cells are still evolving: each
# cell folder is rebuilt from scratch, so a partial extract is replaced
# cleanly, and cells with no time series yet are skipped with a note.
#
# Machine identity is never hard-coded: absolute paths are scrubbed at runtime
# by grteclyn_wrapper.packaging.scrub_paths using tokens derived from the
# environment.
#
# Usage: bash research/bondi_dipole/pack_runaway.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
RUNS="${ROOT}/runs/bondi/runaway"
DEST="${ROOT}/results/bondi-dipole-runaway"

export ROOT SIM_ROOT
scrub() { python3 "${ROOT}/grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py" "$@"; }

if [[ ! -d "${RUNS}" ]]; then
  echo "[pack-runaway] no run tree at ${RUNS#"${ROOT}"/} -- nothing to do"
  exit 0
fi

mkdir -p "${DEST}/campaign" "${DEST}/stars"

# Cells are every top-level run directory.  Names beginning with '_' are
# bookkeeping (the job queue, the abandoned tree), never cells.
cells=()
for celldir in "${RUNS}"/*/; do
  cell="$(basename "${celldir}")"
  [[ "${cell}" == _* ]] && continue
  cells+=("${cell}|${celldir%/}")
done

for entry in "${cells[@]}"; do
  cell="${entry%%|*}"
  celldir="${entry#*|}"
  run="$(find "${celldir}" -maxdepth 1 -type d -name 'bondi_sg_*' | head -n 1)"
  if [[ -z "${run}" ]] || ! compgen -G "${run}/small_data/*.dat" > /dev/null; then
    echo "[pack-runaway] ${cell}: no time series yet -- skipping"
    continue
  fi
  out="${DEST}/campaign/${cell}"
  rm -rf "${out}"
  mkdir -p "${out}"

  # Time series: every stream is self-documenting (# header line).
  cp "${run}"/small_data/*.dat "${out}/"

  # Provenance: what was solved, how well, and what was evolved.
  [[ -f "${run}/params.txt" ]] && cp "${run}/params.txt" "${out}/evolution_params.txt"
  [[ -f "${run}/grtresna/params.txt" ]] && cp "${run}/grtresna/params.txt" "${out}/grtresna_params.txt"
  [[ -f "${run}/grtresna/Ham_and_Mom_errors.txt" ]] && cp "${run}/grtresna/Ham_and_Mom_errors.txt" "${out}/"
  for f in metadata.json initial_data.matter.json; do
    [[ -f "${celldir}/${f}" ]] && cp "${celldir}/${f}" "${out}/"
    [[ -f "${run}/${f}" ]] && cp "${run}/${f}" "${out}/"
  done

  # The launch configuration.  Cells launched by hand carry their own
  # launch.sh; the first four came off the job queue, where the equivalent
  # record is the .job file.  Either way the packed copy is what lets a reader
  # rebuild the cell without reading a params file line by line.
  if [[ -f "${celldir}/launch.sh" ]]; then
    cp "${celldir}/launch.sh" "${out}/launch_config.sh"
  else
    job="$(ls -t "${RUNS}"/_queue/{done,failed,running}/*_"${cell}".job 2>/dev/null | head -n 1 || true)"
    [[ -n "${job}" ]] && cp "${job}" "${out}/launch_config.sh"
  fi

  # Evolution-side diagnostics arrive at every-step cadence (~6000 rows for a
  # t = 200 cell); downsample to dt = 0.5 exactly as the two older packers do.
  # None of these four streams writes its own header line, so supply the column
  # names the evolution code emits (RadialRecipeLevel.cpp write_header_line).
  # Header-only addition: the columns themselves are untouched, and every
  # reader already skips '#' lines.
  for stream in constraint_norms energy_conditions curvature_invariants collapse_diagnostics; do
    [[ -f "${run}/data/${stream}.dat" ]] || continue
    case "${stream}" in
      constraint_norms)
        header="time L2_Ham L2_Mom min_rho_req max_rho_req integral_neg_rho L2_Ham_rel L2_Mom_rel pump_force_L2 governor pump_fi_L2 L2_Ham_amr L2_Mom_amr L2_Ham_amr_rel L2_Mom_amr_rel Linf_Ham_amr L2_Ham_amr_ref" ;;
      energy_conditions)
        header="time matter_min_NEC matter_min_WEC matter_min_SEC matter_min_DEC matter_integral_NEC_violation" ;;
      curvature_invariants)
        header="time max_abs_ricci_scalar max_ricci_tensor_sq max_Kij_sq L2_ricci_scalar" ;;
      collapse_diagnostics)
        header="time min_lapse min_chi max_abs_K min_lapse_x min_lapse_y min_lapse_z max_ah_r min_theta_plus r_at_min_theta_plus min_phi max_phi min_Pi max_Pi pump_work" ;;
      *) header="" ;;
    esac
    HEADER="${header}" python3 - "${run}/data/${stream}.dat" "${out}/${stream}.dat" <<'PY'
import os, sys

src, dst = sys.argv[1], sys.argv[2]
step, last = 0.5, None
with open(src, encoding="utf-8") as fh, open(dst, "w", encoding="utf-8") as out:
    out.write("# downsampled to dt=0.5 from the every-step stream\n")
    wrote_header = False
    for line in fh:
        if line.startswith("#") or not line.strip():
            out.write(line)
            wrote_header = True
            continue
        if not wrote_header:
            if os.environ.get("HEADER"):
                out.write("# " + os.environ["HEADER"] + "\n")
            wrote_header = True
        t = float(line.split()[0])
        if last is None or t - last >= step - 1e-9:
            out.write(line)
            last = t
PY
  done

  found_txt=$(find "${out}" -maxdepth 1 \( -name '*.txt' -o -name '*.json' -o -name '*.sh' \))
  [[ -n "${found_txt}" ]] && scrub ${found_txt}

  echo "[pack-runaway] campaign/${cell}: $(find "${out}" -type f | wc -l) files, $(du -sh "${out}" | cut -f1)"
done

# ---------------------------------------------------------------------------
# Dressed-star profile tables for this rung -- one copy each, not per cell.
# The masses these tables integrate to are what every Newtonian comparison in
# the campaign is measured against, so they are part of the result, not a
# by-product.  omega = 0.7603 is the phantom frequency chosen to match the
# canonical star's |ADM| mass; the mixed pair at 0.75/0.75 is the unmatched one.
# ---------------------------------------------------------------------------
canon_src="${RUNS}/pm/bondi_sg_pair_pm_w075/grtresna/qball_profile.dat"
phant_src="${RUNS}/pm/bondi_sg_pair_pm_w075/grtresna/qball_profile_exotic.dat"
eqm_src="${RUNS}/pm_eqm/bondi_sg_pair_pm_eqm_w075/grtresna/qball_profile_exotic.dat"
[[ -f "${canon_src}" ]] && cp "${canon_src}" "${DEST}/stars/canonical_omega0.750.dat"
[[ -f "${phant_src}" ]] && cp "${phant_src}" "${DEST}/stars/phantom_omega0.750.dat"
[[ -f "${eqm_src}"   ]] && cp "${eqm_src}"   "${DEST}/stars/phantom_omega0.7603.dat"
scrub "${DEST}"/stars/*.dat
echo "[pack-runaway] stars: profile tables for the omega = 0.75 rung"

echo "[pack-runaway] total campaign size: $(du -sh "${DEST}/campaign" | cut -f1)"
echo "[pack-runaway] done"
