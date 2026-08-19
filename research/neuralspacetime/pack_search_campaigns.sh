#!/usr/bin/env bash
# Pack light extracts of the trajectory-search campaigns behind the article
# "Matter-first automated discovery of transient spacetime shortcuts" into
# results/<article-slug>/search/<pack>/.
#
# Companion to pack_publishable_results.sh, which owns the article's figures,
# tables, validation tree and promotion (HQ) runs — this script never touches
# those; it owns the search/ subtree and nothing else.
#
# Safe to re-run at any time, including while a campaign is live:
#   - campaign-level files, analysis tables and progress logs are rebuilt from
#     scratch on every run, so a partial extract from a still-running campaign
#     is replaced cleanly;
#   - campaigns are listed in the CAMPAIGNS table below and skipped with a note
#     (never an error) when their run directory is absent;
#   - per-eval directories are ADDITIVE: an eval already packed here is never
#     deleted just because the live tree no longer has it. The live tree keeps
#     only the top few evals (KEEP_TOP_EVAL_DIRS) and gets trimmed by hand to
#     reclaim disk, so for a pruned eval this pack is the last surviving copy.
#     Set PACK_REFRESH_EVALS=1 to re-pull evals that are already packed.
#
# Heavy files are never copied: initial_data.gridinit (~530 MB/eval),
# small_data/metric_stack (~14 MB/eval), plotfiles and frames all stay in the
# gitignored /runs tree on the machine that produced them.
#
# Hand-written prose (README.md, CAMPAIGN_RESULTS.md) inside a pack is never
# overwritten. Analysis tables are only generated when missing, so a
# hand-curated table with campaign-specific columns survives a re-pack; set
# PACK_REFRESH_ANALYSIS=1 to regenerate them from trajectory.jsonl.
#
# Machine identity is scrubbed at runtime by the shared scrubber in the wrapper
# package (grteclyn_wrapper.packaging.scrub_paths), per the no-identity-in-git
# rule, and the pack is then gated: any surviving host, user or home-path token
# fails the run.
#
# Usage: bash research/neuralspacetime/pack_search_campaigns.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
SEARCH_RUNS="${ROOT}/runs/neuralspacetime/search"
ARTICLE_SLUG="matter-first-automated-discovery-of-transient-spacetime-shortcuts"
DEST="${ROOT}/results/${ARTICLE_SLUG}/search"

: "${PACK_REFRESH_EVALS:=0}"
: "${PACK_REFRESH_ANALYSIS:=0}"

# Optional site knobs (may be unset); used only to learn path tokens to strip.
: "${GRTRESNA_ENV:=}"
: "${OPENMPI_ROOT:=}"
: "${CHOMBO_HOME:=}"
if [[ -f "${ROOT}/grteclyn-wrapper/.env" ]]; then
  set -a
  # shellcheck disable=SC1090
  source <(grep -E '^(GRTRESNA_ENV|OPENMPI_ROOT|CHOMBO_HOME)=' \
    "${ROOT}/grteclyn-wrapper/.env" | sed 's/\r$//' || true)
  set +a
fi

export ROOT SIM_ROOT GRTRESNA_ENV OPENMPI_ROOT CHOMBO_HOME
scrub() { python3 "${ROOT}/grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py" "$@"; }

# pack slug | run dir under runs/neuralspacetime/search | launch-log glob
# Order = lineage order, oldest first; the pack README narrates it.
CAMPAIGNS=(
  "qball-trajectory-pump-seed-search|map_elites/qball_traj_pump_v2|qball_traj_pump_v2*.log"
  "qball-trajectory-evolving-geodesic-shortcut-search|map_elites/qball_traj_fgeo_v1|qball_traj_fgeo_v1*.log"
  "qball-trajectory-geodesic-depth-search|map_elites/qball_traj_fgeo_depth_v1|qball_traj_fgeo_depth_v1*.log"
  "qball-trajectory-cmaes-refinement|cma_es/qball_traj_fgeo_depth_cmaes_v1|qball_traj_fgeo_depth_cmaes_v1*.log"
  "qball-trajectory-fgeo-max-refinement|cma_es/qball_traj_fgeo_max_cmaes_v1|qball_traj_fgeo_max_cmaes_v1*.log"
)

# Campaign-level state. cmaes_state.pkl is deliberately excluded: it is a
# pickle of the live optimiser, only meaningful to the exact library version
# that wrote it, and trajectory.jsonl already carries every genome it holds.
CAMPAIGN_FILES=(
  metadata.json trajectory.jsonl result.json
  ftl_champions.json ftl_retention.jsonl
  archive.json pre_gpu_archive.json validation.json
  warm_start_gen1_seed.jsonl
)

# Per-eval light payload. small_data/metric_stack and initial_data.gridinit are
# the two heavy items and are never listed here.
EVAL_FILES=(metadata.json params.txt score.json initial_data.matter.json run.log)

copy_eval_light() {
  local src="$1" dst="$2"
  mkdir -p "${dst}"
  local f
  for f in "${EVAL_FILES[@]}"; do
    [[ -f "${src}/${f}" ]] && cp -a "${src}/${f}" "${dst}/${f}"
  done
  # Diagnostics: every stream here is a few kB of ASCII.
  if [[ -d "${src}/small_data" ]]; then
    mkdir -p "${dst}/small_data"
    find "${src}/small_data" -maxdepth 1 -type f \
      \( -name '*.dat' -o -name '*.json' \) -exec cp -a {} "${dst}/small_data/" \;
  fi
  if [[ -d "${src}/data" ]]; then
    mkdir -p "${dst}/data"
    find "${src}/data" -maxdepth 1 -type f -name '*.dat' -exec cp -a {} "${dst}/data/" \;
  fi
  # Solver provenance: what was solved and how well it converged.
  if [[ -d "${src}/grtresna" ]]; then
    mkdir -p "${dst}/grtresna"
    for f in params.txt Ham_and_Mom_errors.txt; do
      [[ -f "${src}/grtresna/${f}" ]] && cp -a "${src}/grtresna/${f}" "${dst}/grtresna/${f}"
    done
  fi
  # Post-load constraint gate: the evidence the eval was admitted honestly.
  if [[ -d "${src}/postload_gate" ]]; then
    rm -rf "${dst}/postload_gate"
    mkdir -p "${dst}/postload_gate"
    (cd "${src}/postload_gate" && find . -type f \
      \( -name '*.json' -o -name '*.txt' -o -name '*.dat' -o -name '*.log' \) -print0) \
      | while IFS= read -r -d '' rel; do
          mkdir -p "${dst}/postload_gate/$(dirname -- "${rel}")"
          cp -a "${src}/postload_gate/${rel}" "${dst}/postload_gate/${rel}"
        done
  fi
  find "${dst}" -type f -exec chmod 0644 {} +
}

echo "Packing search campaigns -> results/${ARTICLE_SLUG}/search/"
mkdir -p "${DEST}"

for entry in "${CAMPAIGNS[@]}"; do
  IFS='|' read -r slug relpath logglob <<< "${entry}"
  src="${SEARCH_RUNS}/${relpath}"
  out="${DEST}/${slug}"

  if [[ ! -d "${src}" ]]; then
    echo "[pack-search] ${slug}: no run directory yet -- skipping"
    continue
  fi
  if [[ ! -f "${src}/trajectory.jsonl" ]]; then
    echo "[pack-search] ${slug}: no trajectory yet -- skipping"
    continue
  fi

  mkdir -p "${out}/run" "${out}/analysis" "${out}/logs"

  # Campaign-level files are cheap and always refreshed.
  for f in "${CAMPAIGN_FILES[@]}"; do
    [[ -f "${src}/${f}" ]] && cp -a "${src}/${f}" "${out}/run/${f}"
  done
  # A warm-start seed trajectory is written next to the campaign, not inside.
  for seed in "${src}"_gen1_seed.jsonl "${src}"_warm_start.jsonl; do
    [[ -f "${seed}" ]] && cp -a "${seed}" "${out}/run/warm_start_gen1_seed.jsonl"
  done

  n_new=0 n_kept=0
  while IFS= read -r evdir; do
    name="$(basename -- "${evdir}")"
    if [[ -d "${out}/run/${name}" && "${PACK_REFRESH_EVALS}" != "1" ]]; then
      n_kept=$((n_kept + 1))
      continue
    fi
    copy_eval_light "${evdir}" "${out}/run/${name}"
    n_new=$((n_new + 1))
  done < <(find "${src}" -maxdepth 1 -type d -name 'eval_*' | sort)

  # Evals that live only here now (pruned upstream) are reported, never removed.
  n_orphan=0
  while IFS= read -r packed; do
    [[ -d "${src}/$(basename -- "${packed}")" ]] || n_orphan=$((n_orphan + 1))
  done < <(find "${out}/run" -maxdepth 1 -type d -name 'eval_*' | sort)

  # Progress lines lifted out of the multi-hundred-MB launch log. A pack that
  # already ships a curated progress log keeps it (PACK_REFRESH_LOGS=1 to redo).
  shopt -s nullglob
  existing_logs=("${out}"/logs/*.log)
  if [[ ${#existing_logs[@]} -eq 0 || "${PACK_REFRESH_LOGS:-0}" == "1" ]]; then
    log_out="${out}/logs/progress.log"
    : > "${log_out}"
    for lg in "${SEARCH_RUNS}/$(dirname -- "${relpath}")"/${logglob}; do
      grep -E '^\[(optimize|qd|search|cmaes)\]|generation .* best' "${lg}" \
        >> "${log_out}" 2>/dev/null || true
    done
    if [[ ! -s "${log_out}" ]]; then
      rm -f "${log_out}"
      rmdir "${out}/logs" 2>/dev/null || true
    fi
  fi
  shopt -u nullglob

  # Analysis tables: generated only when absent (see header).
  if [[ "${PACK_REFRESH_ANALYSIS}" == "1" || ! -f "${out}/analysis/per_eval.csv" ]]; then
    OUT_DIR="${out}/analysis" python3 "${SCRIPT_DIR}/search_campaign_tables.py" \
      "${src}/trajectory.jsonl" "${src}/metadata.json"
  fi

  # Scrub every text artefact this pack owns.
  while IFS= read -r -d '' f; do
    scrub "${f}"
  done < <(find "${out}" -type f \
    \( -name '*.json' -o -name '*.jsonl' -o -name '*.txt' -o -name '*.log' -o -name '*.csv' \) -print0)

  echo "[pack-search] ${slug}: $(find "${out}" -type f | wc -l) files, $(du -sh "${out}" | cut -f1)" \
       "(evals: ${n_new} packed, ${n_kept} already present, ${n_orphan} pruned upstream)"
done

# Identity gate over the whole search subtree.
DEST="${DEST}" python3 - <<'PY'
import os, re, socket, sys
from pathlib import Path

dest = Path(os.environ["DEST"])
skip = {".png", ".jpg", ".jpeg", ".pdf", ".npz", ".bin", ".ex", ".so", ".pyc", ".pkl"}
GENERIC = {
    "", "home", "usr", "tmp", "var", "opt", "etc", "lib", "bin", "local",
    "path", "to", "users", "user", "envs", "env", "conda", "mamba",
    "test", "simulation", "src", "lib64", "share", "include", "runs",
}
KEEP = {
    "grteclyn", "grtresna", "chombo", "grteclyn-wrapper", "amrex", "openmpi",
    "radialrecipe", "bondi", "neuralspacetime",
}
tokens: set[str] = set()

def add(value: str) -> None:
    value = (value or "").strip().strip("/.")
    if len(value) < 3 or value.lower() in GENERIC or value.lower() in KEEP:
        return
    if re.split(r"[-_]", value, maxsplit=1)[0].lower() in KEEP:
        return
    tokens.add(value)

add(os.environ.get("USER", ""))
add(os.environ.get("LOGNAME", ""))
add(Path(os.environ.get("HOME", "")).name)
host = socket.gethostname().split(".")[0]
add(host)
for part in re.split(r"[-_.]", host):
    if len(part) >= 4:
        add(part)
for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
    for part in Path(os.environ.get(key) or "").parts:
        add(part)

patterns = [re.compile(rf"(?<![A-Za-z0-9_]){re.escape(t)}(?![A-Za-z0-9_])", re.I)
            for t in sorted(tokens, key=len, reverse=True)]
patterns.append(re.compile(re.escape("/" + "home" + "/"), re.I))

hits, n_files = [], 0
for f in dest.rglob("*"):
    if not f.is_file() or f.suffix.lower() in skip:
        continue
    n_files += 1
    try:
        text = f.read_bytes().decode("utf-8", errors="replace")
    except OSError:
        continue
    for pat in patterns:
        m = pat.search(text)
        if m:
            hits.append(f"{f.relative_to(dest)}:{text.count(chr(10), 0, m.start()) + 1}")
            break

if hits:
    print("ERROR: machine-identifying tokens remain in the search pack:", file=sys.stderr)
    for h in hits[:50]:
        print(f"  {h}", file=sys.stderr)
    if len(hits) > 50:
        print(f"  ... and {len(hits) - 50} more", file=sys.stderr)
    sys.exit(1)
print(f"machine-token gate: clean ({n_files} text files; {len(tokens)} runtime tokens)")
PY

echo "[pack-search] total: $(du -sh "${DEST}" | cut -f1)"
echo "[pack-search] done -- see ${DEST#"${ROOT}/"}/README.md"
