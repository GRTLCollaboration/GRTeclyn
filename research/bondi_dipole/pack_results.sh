#!/usr/bin/env bash
# Pack light, GitHub-friendly extracts of the Bondi dipole runaway campaign
# into results/bondi-dipole-runaway/.
#
# Copies ONLY small artefacts: the per-cell time series (small_data/*.dat),
# the solve + evolution parameter files, the dressed-star profile tables, a
# and the code patches.  Frames en masse (~120 MB per cell), plotfiles and
# gridinit stay in the gitignored /runs tree.
#
# Machine identity is never hard-coded: absolute paths are scrubbed at runtime
# by grteclyn_wrapper.packaging.scrub_paths using tokens derived from the
# environment.
#
# Usage: bash research/bondi_dipole/pack_results.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
RUNS="${ROOT}/runs/bondi_rerun/published"
DEST="${ROOT}/results/bondi-dipole-runaway"
GRTRESNA_ROOT="${GRTRESNA_ROOT:-${SIM_ROOT}/GRTresna}"

export ROOT SIM_ROOT
scrub() { python3 "${ROOT}/grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py" "$@"; }

# cell | run subdir | evolution stop time
CELLS=(
  "single_p|single_p/bondi_sg_single_p|40"
  "single_m|single_m/bondi_sg_single_m|40"
  "pair_pm|pair_pm/bondi_sg_pair_pm|60"
  # The mirror control: sectors swapped in space, so the drift must reverse.
  "pair_mp_mirror|pair_mp_mirror/bondi_sg_pair_mp|60"
)

echo "[pack] dest: ${DEST}"
mkdir -p "${DEST}"/{campaign,stars,patches,analysis}

# ---------------------------------------------------------------------------
# 1. Per-cell time series + provenance (params, residuals, metadata)
# ---------------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  IFS='|' read -r cell sub _stop <<<"${spec}"
  src="${RUNS}/${sub}"
  out="${DEST}/campaign/${cell}"
  if [[ ! -d "${src}" ]]; then
    echo "[pack] WARNING: missing run dir for ${cell} (${src}) -- skipping"
    continue
  fi
  mkdir -p "${out}"
  # Time series: every stream is self-documenting (# header line).
  cp "${src}"/small_data/*.dat "${out}/"
  # Provenance: what was solved, how well, and what was evolved.
  cp "${src}/params.txt" "${out}/evolution_params.txt"
  cp "${src}/grtresna/params.txt" "${out}/grtresna_params.txt"
  cp "${src}/grtresna/Ham_and_Mom_errors.txt" "${out}/"
  cp "${src}/metadata.json" "${src}/initial_data.matter.json" "${out}/"
  # Evolution-side diagnostics live in data/ at every-step cadence (6001 rows,
  # ~2 MB each).  Downsample to dt = 0.5 -- plenty for a drift/violation curve.
  for stream in constraint_norms energy_conditions curvature_invariants collapse_diagnostics; do
    [[ -f "${src}/data/${stream}.dat" ]] || continue
    # None of these four streams writes its own header line, so supply the
    # column names the evolution code emits (RadialRecipeLevel.cpp
    # write_header_line).  Header-only addition: the columns are untouched and
    # every reader already skips '#' lines.  Kept in step with pack_campaign.sh.
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
    HEADER="${header}" python3 - "${src}/data/${stream}.dat" "${out}/${stream}.dat" <<'PY'
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
            # These streams ship without a header line; supply the column
            # names the evolution code writes.
            if os.environ.get("HEADER"):
                out.write("# " + os.environ["HEADER"] + "\n")
            wrote_header = True
        t = float(line.split()[0])
        if last is None or t - last >= step - 1e-9:
            out.write(line)
            last = t
PY
  done
  scrub "${out}"/*.txt "${out}"/*.json
  echo "[pack] campaign/${cell}: $(ls "${out}" | wc -l) files"
done

# ---------------------------------------------------------------------------
# 2. Dressed-star profile tables (r, phi0, alpha) -- one copy each, not per cell
# ---------------------------------------------------------------------------
cp "${RUNS}/pair_pm/bondi_sg_pair_pm/grtresna/qball_profile.dat" \
   "${DEST}/stars/canonical_omega0.550.dat"
cp "${RUNS}/pair_pm/bondi_sg_pair_pm/grtresna/qball_profile_exotic.dat" \
   "${DEST}/stars/phantom_omega0.550.dat"
cp "${RUNS}/pair_pm_eqm/bondi_sg_pair_pm_eqm/grtresna/qball_profile_exotic.dat" \
   "${DEST}/stars/phantom_omega0.566.dat"
scrub "${DEST}"/stars/*.dat
echo "[pack] stars: 3 profile tables"


# ---------------------------------------------------------------------------
# 3b. Movies -- the views that carry the result (60-180 kB each, not the full
#     19-field set).  Stitch first if missing:
#       bash grteclyn-wrapper/scripts/plot/make_movies.sh runs/bondi_rerun/published/<cell>/<run> \
#            --only scalar_activity_proj_z chi_minus_1_z rho_req_z
# ---------------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  IFS='|' read -r cell sub _stop <<<"${spec}"
  mdir="${RUNS}/${sub}/movies"
  [[ -d "${mdir}" ]] || continue
  out="${DEST}/movies/${cell}"
  mkdir -p "${out}"
  # Matter motion (projection = the clearest view of the drift) + geometry sign.
  for view in scalar_activity_proj_z chi_minus_1_z; do
    [[ -f "${mdir}/movie_${view}.mp4" ]] && cp "${mdir}/movie_${view}.mp4" "${out}/${view}.mp4"
  done
  # Signed energy density: only meaningful where both sectors are present.
  if [[ "${cell}" == pair_pm* && -f "${mdir}/movie_rho_req_z.mp4" ]]; then
    cp "${mdir}/movie_rho_req_z.mp4" "${out}/rho_req_z.mp4"
  fi
done
echo "[pack] movies: $(find "${DEST}/movies" -name '*.mp4' | wc -l) clips ($(du -sh "${DEST}/movies" | cut -f1))"

# ---------------------------------------------------------------------------
# 4. Code patches -- the matter-model modifications this campaign required
# ---------------------------------------------------------------------------
git -C "${GRTRESNA_ROOT}" diff -- \
  Source/Matter/BosonStarParams.hpp Source/Matter/ComplexScalarField.cpp \
  > "${DEST}/patches/grtresna-matter.diff" || true
git -C "${ROOT}" show a83ad2d7 -- \
  grteclyn-wrapper/src grteclyn-wrapper/tests grteclyn-wrapper/scripts \
  > "${DEST}/patches/wrapper-per-lump-omega.diff" || true
git -C "${ROOT}" show e4349d35 -- grteclyn-wrapper \
  > "${DEST}/patches/wrapper-phantom-star-solver.diff" || true
scrub "${DEST}"/patches/*.diff
echo "[pack] patches: $(ls "${DEST}/patches" | wc -l) diffs"

# ---------------------------------------------------------------------------
# 5. Full campaign journal (the debugging record the article's method cites)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 6. Derived tables (trajectories, per-cell summary) from the packed data
# ---------------------------------------------------------------------------
python3 "${DEST}/analysis/make_tables.py" "${DEST}"
python3 "${DEST}/analysis/newtonian_reference.py" "${DEST}" >/dev/null
python3 "${DEST}/analysis/separation_scaling.py" "${DEST}" >/dev/null
echo "[pack] analysis: summary.{csv,md}, trajectories.csv, newtonian_reference.csv"
# stars/star_family.csv is NOT regenerated here: the family scan needs the
# wrapper venv (numpy/scipy).  Refresh it with
#   PYTHONPATH=grteclyn-wrapper/src grteclyn-wrapper/.venv/bin/python \
#     results/bondi-dipole-runaway/analysis/star_family_scan.py

echo "[pack] total size: $(du -sh "${DEST}" | cut -f1)"
echo "[pack] done"
