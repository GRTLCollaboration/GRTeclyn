#!/usr/bin/env bash
# Pack light, GitHub-friendly extracts of the Bondi dipole runaway campaign
# into results/bondi-dipole-runaway/.
#
# Copies ONLY small artefacts: the per-cell time series (small_data/*.dat),
# the solve + evolution parameter files, the dressed-star profile tables, a
# curated set of PNG frames, and the code patches.  Frames en masse (~120 MB
# per cell), plotfiles, gridinit and movies stay in the gitignored /runs tree.
#
# Machine identity is never hard-coded: absolute paths are scrubbed at runtime
# by scrub_paths.py using tokens derived from the environment.
#
# Usage: bash research/bondi_dipole/pack_results.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
RUNS="${ROOT}/runs/bondi"
DEST="${ROOT}/results/bondi-dipole-runaway"
GRTRESNA_ROOT="${GRTRESNA_ROOT:-${SIM_ROOT}/GRTresna}"

export ROOT SIM_ROOT
scrub() { python3 "${SCRIPT_DIR}/scrub_paths.py" "$@"; }

# cell | run subdir | evolution stop time
CELLS=(
  "single_p|single_p/bondi_sg_single_p|40"
  "single_m|single_m/bondi_sg_single_m|40"
  "pair_pm|pair_pm/bondi_sg_pair_pm|60"
  "pair_pp|pair_pp/bondi_sg_pair_pp|60"
  "pair_mm|pair_mm/bondi_sg_pair_mm|60"
  "pair_pm_eqm|pair_pm_eqm/bondi_sg_pair_pm_eqm|60"
)

echo "[pack] dest: ${DEST}"
mkdir -p "${DEST}"/{data,stars,figures,patches,debug_log,analysis}

# ---------------------------------------------------------------------------
# 1. Per-cell time series + provenance (params, residuals, metadata)
# ---------------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  IFS='|' read -r cell sub _stop <<<"${spec}"
  src="${RUNS}/${sub}"
  out="${DEST}/data/${cell}"
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
  scrub "${out}"/*.txt "${out}"/*.json
  echo "[pack] data/${cell}: $(ls "${out}" | wc -l) files"
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
# 3. Curated frames -- t = 0, 1/3, 2/3, end for the two article-critical views
# ---------------------------------------------------------------------------
for spec in "${CELLS[@]}"; do
  IFS='|' read -r cell sub stop <<<"${spec}"
  src="${RUNS}/${sub}/frames"
  [[ -d "${src}" ]] || continue
  for field in scalar_activity_z chi_minus_1_z; do
    fdir="${src}/${field}/frames"
    [[ -d "${fdir}" ]] || continue
    CELL="${cell}" FIELD="${field}" STOP="${stop}" DEST="${DEST}" python3 - "${fdir}" <<'PY'
import os, pathlib, re, shutil, sys

frames = sorted(
    pathlib.Path(sys.argv[1]).glob("*.png"),
    key=lambda p: int(re.search(r"(\d+)\s*$", p.stem).group(1)),
)
if not frames:
    raise SystemExit(0)
cell, field = os.environ["CELL"], os.environ["FIELD"]
stop, dest = float(os.environ["STOP"]), pathlib.Path(os.environ["DEST"])
out = dest / "figures" / cell
out.mkdir(parents=True, exist_ok=True)
n = len(frames)
# Singles are calibration runs (start/end suffice); pairs carry the trajectory.
fractions = (0.0, 1.0) if cell.startswith("single") else (0.0, 1 / 3, 2 / 3, 1.0)
for frac in fractions:
    idx = min(n - 1, round(frac * (n - 1)))
    t = stop * idx / (n - 1)
    shutil.copy2(frames[idx], out / f"{field}_t{t:04.0f}.png")
PY
  done
done
echo "[pack] figures: $(find "${DEST}/figures" -name '*.png' | wc -l) frames"

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
cp "${ROOT}/research/bondi_dipole_debug.md" "${DEST}/debug_log/"
scrub "${DEST}/debug_log/bondi_dipole_debug.md"

# ---------------------------------------------------------------------------
# 6. Derived tables (trajectories, per-cell summary) from the packed data
# ---------------------------------------------------------------------------
python3 "${DEST}/analysis/make_tables.py" "${DEST}"

echo "[pack] total size: $(du -sh "${DEST}" | cut -f1)"
echo "[pack] done"
