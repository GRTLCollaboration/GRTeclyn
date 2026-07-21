#!/usr/bin/env bash
# Pack light, GitHub-friendly extracts for the article
# "Matter-first automated discovery of transient spacetime shortcuts"
# into results/<article-slug>/. Skips frames, plotfiles, gridinit, metric stacks.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
SRC_PAPER="${ROOT}/research/neuralspacetime"
ARTICLE_SLUG="matter-first-automated-discovery-of-transient-spacetime-shortcuts"
DEST="${ROOT}/results/${ARTICLE_SLUG}"
RUNS="${ROOT}/runs"

# Rewrite absolute machine paths to repo-relative / env placeholders.
scrub_machine_paths() {
  local f="$1"
  [[ -f "${f}" ]] || return 0
  case "${f}" in
    *.png|*.jpg|*.jpeg|*.pdf|*.npz|*.bin|*.ex|*.so) return 0 ;;
  esac
  ROOT="${ROOT}" SIM_ROOT="${SIM_ROOT}" python3 - "${f}" <<'PY'
import os, pathlib, re, sys
path = pathlib.Path(sys.argv[1])
try:
    text = path.read_text(encoding="utf-8")
except (UnicodeDecodeError, OSError):
    raise SystemExit(0)
root = os.environ["ROOT"]
sim = os.environ["SIM_ROOT"]
for old, new in (
    (root + "/", ""),
    (root, "$GRTECLYN_ROOT"),
    (sim + "/", ""),
    (sim, "$SIM_ROOT"),
):
    text = text.replace(old, new)
text = re.sub(r"/home/jovyan/[^\s\"']*", r"$HOME/<redacted>", text)
text = text.replace("nachevsky", "<user>")
text = text.replace(".mlspace", "<env>")
path.write_text(text, encoding="utf-8")
PY
}

copy_file() {
  local src="$1"
  local dest="$2"
  if [[ ! -e "${src}" ]]; then
    echo "[skip missing] ${src}" >&2
    return 0
  fi
  mkdir -p "$(dirname -- "${dest}")"
  cp -a "${src}" "${dest}"
  scrub_machine_paths "${dest}"
  echo "[ok] ${dest#"${ROOT}/"}"
}

copy_hq_light() {
  local name="$1"
  local src="${RUNS}/grtresna_promote/${name}"
  local dest="${DEST}/runs/grtresna_promote/${name}"
  if [[ ! -d "${src}" ]]; then
    echo "[skip missing run] ${src}" >&2
    return 0
  fi
  mkdir -p "${dest}/small_data" "${dest}/data"
  for f in score.json params.txt metadata.json; do
    [[ -e "${src}/${f}" ]] && copy_file "${src}/${f}" "${dest}/${f}"
  done
  for f in evolving_geodesic.json confinement.dat; do
    [[ -e "${src}/small_data/${f}" ]] && copy_file "${src}/small_data/${f}" "${dest}/small_data/${f}"
  done
  for f in constraint_norms.dat collapse_diagnostics.dat; do
    [[ -e "${src}/data/${f}" ]] && copy_file "${src}/data/${f}" "${dest}/data/${f}"
  done
}

copy_eval_light() {
  local src_eval="$1"
  local dest_eval="$2"
  if [[ ! -d "${src_eval}" ]]; then
    echo "[skip missing eval] ${src_eval}" >&2
    return 0
  fi
  mkdir -p "${dest_eval}"
  for f in score.json metadata.json params.txt initial_data.matter.json; do
    [[ -e "${src_eval}/${f}" ]] && copy_file "${src_eval}/${f}" "${dest_eval}/${f}"
  done
}

echo "Packing publishable results -> ${DEST}"
mkdir -p "${DEST}"
# Refresh regenerable trees; keep README.md if already present.
for d in article figures manifests validation runs; do
  rm -rf "${DEST}/${d}"
done
rm -f "${DEST}/PLOT_DATA_SOURCES.txt"

# ---- article extracts + provenance + used figures ----
mkdir -p "${DEST}/article/data" "${DEST}/figures" "${DEST}/manifests"
copy_file "${SRC_PAPER}/article/PLOT_DATA_SOURCES.txt" "${DEST}/article/PLOT_DATA_SOURCES.txt"
for f in \
  campaign_progress.txt \
  cmaes_progress.txt \
  fgeo_sweeps.txt \
  constraints_rm.txt \
  constraints_resolution.txt \
  constraint_order.txt \
  pump_work_budget.txt \
  casimir_bounds.txt \
  null_constraint_ray.txt \
  null_constraint_ray_rm.txt \
  endpoint_gauge_rm.txt
do
  copy_file "${SRC_PAPER}/article/data/${f}" "${DEST}/article/data/${f}"
done

for f in \
  fig_eval146_lump_trajectories.png \
  fig_rf_rho_t0.png fig_rf_rho_t4.png fig_rf_rho_t16.png fig_rf_rho_t30.png \
  fig_rf_shift_t0.png fig_rf_shift_t4.png fig_rf_shift_t16.png fig_rf_shift_t30.png
do
  copy_file "${SRC_PAPER}/figures/${f}" "${DEST}/figures/${f}"
done

# ---- validation certificates (already light JSON) ----
if [[ -d "${SRC_PAPER}/validation" ]]; then
  mkdir -p "${DEST}/validation"
  cp -a "${SRC_PAPER}/validation/." "${DEST}/validation/"
  while IFS= read -r -d '' f; do
    scrub_machine_paths "${f}"
  done < <(find "${DEST}/validation" -type f -print0)
  echo "[ok] results/${ARTICLE_SLUG}/validation/ (from research/.../validation)"
fi

copy_file \
  "${ROOT}/grteclyn-wrapper/scripts/campaigns/promote/bicomplex_cmaes_v1/manifest_pumpfree.json" \
  "${DEST}/manifests/manifest_pumpfree.json"

# ---- ME / CMA / canonical campaigns ----
for f in trajectory.jsonl archive.json ftl_champions.json metadata.json; do
  copy_file \
    "${RUNS}/grtresna_qd/qball_traj_bicomplex_v1/${f}" \
    "${DEST}/runs/grtresna_qd/qball_traj_bicomplex_v1/${f}"
done
copy_eval_light \
  "${RUNS}/grtresna_qd/qball_traj_bicomplex_v1/eval_000087" \
  "${DEST}/runs/grtresna_qd/qball_traj_bicomplex_v1/eval_000087"

copy_file \
  "${RUNS}/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/trajectory.jsonl" \
  "${DEST}/runs/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/trajectory.jsonl"
copy_eval_light \
  "${RUNS}/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/eval_000146" \
  "${DEST}/runs/grtresna_cmaes/qball_traj_bicomplex_cmaes_v1/eval_000146"

for f in trajectory.jsonl archive.json ftl_champions.json metadata.json; do
  copy_file \
    "${RUNS}/grtresna_qd/qball_traj_canonical_v1/${f}" \
    "${DEST}/runs/grtresna_qd/qball_traj_canonical_v1/${f}"
done
copy_eval_light \
  "${RUNS}/grtresna_qd/qball_traj_canonical_v1/eval_000130" \
  "${DEST}/runs/grtresna_qd/qball_traj_canonical_v1/eval_000130"

# ---- frozen champion source + HQ ladder ----
copy_file \
  "${RUNS}/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/CHAMPION.json" \
  "${DEST}/runs/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/CHAMPION.json"
copy_eval_light \
  "${RUNS}/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/eval_000146" \
  "${DEST}/runs/grtresna_promote/sources/qball_traj_bicomplex_cmaes_v1/eval_000146"

for name in \
  bcma_rc_L128_N192_t30_hq_eval000146 \
  bcma_rm_L128_N256_t30_hq_eval000146 \
  bcma_rf_L128_N384_t30_hq_eval000146 \
  bcma_ds_L96_N192_t30_hq_eval000146 \
  bcma_dl_L160_N320_t30_hq_eval000146 \
  bcma_pfrm_L128_N256_t30_hq_eval000146
do
  copy_hq_light "${name}"
done

# ---- SOURCE map (audit pointer; full provenance stays in article/) ----
copy_file "${SRC_PAPER}/article/PLOT_DATA_SOURCES.txt" "${DEST}/PLOT_DATA_SOURCES.txt"

# Audit README (source of truth next to the packer).
copy_file "${SRC_PAPER}/publishable_results_README.md" "${DEST}/README.md"

# Final sweep + gate: publishable pack must not mention this machine.
while IFS= read -r -d '' f; do
  scrub_machine_paths "${f}"
done < <(find "${DEST}" -type f -print0)

if grep -RIlE '/home/jovyan|nachevsky|mlspace' "${DEST}" >/dev/null 2>&1; then
  echo "ERROR: machine-identifying paths remain in ${DEST}:" >&2
  grep -RIlE '/home/jovyan|nachevsky|mlspace' "${DEST}" >&2 || true
  exit 1
fi

du -sh "${DEST}"
echo "Done. See ${DEST}/README.md for the audit layout."
