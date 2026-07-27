#!/usr/bin/env bash
# Pack light, GitHub-friendly extracts for the article
# "Matter-first automated discovery of transient spacetime shortcuts"
# into results/<article-slug>/. Skips frames, plotfiles, gridinit, metric stacks.
#
# Machine identity is NEVER hard-coded here. Scrub tokens are derived at runtime
# from $ROOT / $HOME / $USER / hostname / optional env path knobs.
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
SRC_PAPER="${ROOT}/research/neuralspacetime"
ARTICLE_SLUG="matter-first-automated-discovery-of-transient-spacetime-shortcuts"
DEST="${ROOT}/results/${ARTICLE_SLUG}"
RUNS="${ROOT}/runs"

# Optional site knobs (may be unset); used only to learn path tokens to strip.
: "${GRTRESNA_ENV:=}"
: "${OPENMPI_ROOT:=}"
: "${CHOMBO_HOME:=}"

# Rewrite absolute machine paths to repo-relative / env placeholders.
# Tokens come from the live environment — no host/user literals in this file.
scrub_machine_paths() {
  local f="$1"
  [[ -f "${f}" ]] || return 0
  case "${f}" in
    *.png|*.jpg|*.jpeg|*.pdf|*.npz|*.bin|*.ex|*.so) return 0 ;;
  esac
  ROOT="${ROOT}" SIM_ROOT="${SIM_ROOT}" \
  GRTRESNA_ENV="${GRTRESNA_ENV}" OPENMPI_ROOT="${OPENMPI_ROOT}" CHOMBO_HOME="${CHOMBO_HOME}" \
  python3 - "${f}" <<'PY'
import os, pathlib, re, socket, sys

path = pathlib.Path(sys.argv[1])
try:
    text = path.read_text(encoding="utf-8")
except (UnicodeDecodeError, OSError):
    raise SystemExit(0)

root = os.environ["ROOT"]
sim = os.environ["SIM_ROOT"]

# Longest absolute prefixes first.
abs_prefixes = []
for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
    val = os.environ.get(key) or ""
    if val and os.path.isabs(val):
        abs_prefixes.append(os.path.abspath(val))
abs_prefixes = sorted(set(abs_prefixes), key=len, reverse=True)

for prefix in abs_prefixes:
    if prefix == root:
        text = text.replace(prefix + "/", "")
        text = text.replace(prefix, "$GRTECLYN_ROOT")
    elif prefix == sim:
        text = text.replace(prefix + "/", "")
        text = text.replace(prefix, "$SIM_ROOT")
    else:
        text = text.replace(prefix + "/", "$SITE/")
        text = text.replace(prefix, "$SITE")

# Any remaining absolute home-style paths (pattern built at runtime).
_home = "/" + "home" + "/"
text = re.sub(rf"(?i){_home}[^/\s\"']+/[^\s\"']*", r"$HOME/<redacted>", text)
text = re.sub(rf"(?i){_home}[^/\s\"']+", r"$HOME/<redacted>", text)

# Identity tokens derived from this machine (never hard-coded host/user names).
GENERIC = {
    "", "home", "usr", "tmp", "var", "opt", "etc", "lib", "bin", "local",
    "path", "to", "users", "user", "envs", "env", "conda", "mamba",
    "test", "simulation", "src", "lib64", "share", "include",
}
# Product / campaign path segments that must survive scrubbing.
KEEP = {
    "grteclyn", "grtresna", "chombo", "grteclyn-wrapper", "neuralspacetime",
    "amrex", "openmpi", "radialrecipe",
}
tokens: set[str] = set()

def add_token(s: str) -> None:
    s = (s or "").strip().strip("/.")
    if len(s) < 3 or s.lower() in GENERIC or s.lower() in KEEP:
        return
    # Drop versioned product dirs like openmpi-5.0.8
    base = re.split(r"[-_]", s, maxsplit=1)[0].lower()
    if base in KEEP:
        return
    tokens.add(s)

add_token(os.environ.get("USER", ""))
add_token(os.environ.get("LOGNAME", ""))
add_token(pathlib.Path(os.environ.get("HOME", "")).name)
add_token(socket.gethostname())
add_token(socket.gethostname().split(".")[0])
for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
    val = os.environ.get(key) or ""
    if not val:
        continue
    for part in pathlib.Path(val).parts:
        add_token(part)

# Hostname fragments that appear inside MPI process names (host-xxx-yyy).
host = socket.gethostname().split(".")[0]
for part in re.split(r"[-_.]", host):
    if len(part) >= 4:
        add_token(part)

# Word-boundary replace so campaign names like grtresna_promote stay intact.
for token in sorted(tokens, key=len, reverse=True):
    text = re.sub(
        rf"(?<![A-Za-z0-9_]){re.escape(token)}(?![A-Za-z0-9_])",
        "<redacted>",
        text,
        flags=re.IGNORECASE,
    )

path.write_text(text, encoding="utf-8")
PY
}

assert_pack_has_no_machine_identity() {
  ROOT="${ROOT}" SIM_ROOT="${SIM_ROOT}" DEST="${DEST}" \
  GRTRESNA_ENV="${GRTRESNA_ENV}" OPENMPI_ROOT="${OPENMPI_ROOT}" CHOMBO_HOME="${CHOMBO_HOME}" \
  python3 <<'PY'
from pathlib import Path
import os, re, socket, sys

dest = Path(os.environ["DEST"])
skip_suffix = {".png", ".jpg", ".jpeg", ".pdf", ".npz", ".bin", ".ex", ".so", ".pyc"}

GENERIC = {
    "", "home", "usr", "tmp", "var", "opt", "etc", "lib", "bin", "local",
    "path", "to", "users", "user", "envs", "env", "conda", "mamba",
    "test", "simulation", "src", "lib64", "share", "include",
}
KEEP = {
    "grteclyn", "grtresna", "chombo", "grteclyn-wrapper", "neuralspacetime",
    "amrex", "openmpi", "radialrecipe",
}
tokens: set[str] = set()

def add_token(s: str) -> None:
    s = (s or "").strip().strip("/.")
    if len(s) < 3 or s.lower() in GENERIC or s.lower() in KEEP:
        return
    base = re.split(r"[-_]", s, maxsplit=1)[0].lower()
    if base in KEEP:
        return
    tokens.add(s)

add_token(os.environ.get("USER", ""))
add_token(os.environ.get("LOGNAME", ""))
add_token(Path(os.environ.get("HOME", "")).name)
host = socket.gethostname().split(".")[0]
add_token(host)
for part in re.split(r"[-_.]", host):
    if len(part) >= 4:
        add_token(part)
for key in ("ROOT", "SIM_ROOT", "GRTRESNA_ENV", "OPENMPI_ROOT", "CHOMBO_HOME", "HOME"):
    val = os.environ.get(key) or ""
    if not val:
        continue
    for part in Path(val).parts:
        add_token(part)

patterns = [
    re.compile(rf"(?<![A-Za-z0-9_]){re.escape(t)}(?![A-Za-z0-9_])", re.I)
    for t in sorted(tokens, key=len, reverse=True)
]
patterns.append(re.compile(re.escape("/" + "home" + "/"), re.I))

hits = []
n_files = 0
for f in dest.rglob("*"):
    if not f.is_file() or f.suffix.lower() in skip_suffix:
        continue
    n_files += 1
    try:
        raw = f.read_bytes()
    except OSError:
        continue
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("utf-8", errors="replace")
    for pat in patterns:
        m = pat.search(text)
        if m:
            line = text.count("\n", 0, m.start()) + 1
            hits.append(f"{f.relative_to(dest)}:{line}")
            break

if hits:
    print("ERROR: machine-identifying tokens remain in results pack:", file=sys.stderr)
    for h in hits[:50]:
        print(f"  {h}", file=sys.stderr)
    if len(hits) > 50:
        print(f"  ... and {len(hits) - 50} more", file=sys.stderr)
    sys.exit(1)
print(f"machine-token gate: clean ({n_files} text files; {len(tokens)} runtime tokens)")
PY
}

# The runs/ tree is append-only by default. A run directory under /runs may have
# been re-run since the article figures were generated (RC/RM were), so blindly
# re-copying would silently replace the data the published plots were made from.
# Set PACK_REFRESH_RUNS=1 to deliberately re-pull existing run extracts.
: "${PACK_REFRESH_RUNS:=0}"

copy_file() {
  local src="$1"
  local dest="$2"
  if [[ ! -e "${src}" ]]; then
    echo "[skip missing] ${src}" >&2
    return 0
  fi
  if [[ -e "${dest}" && "${dest}" == "${DEST}/runs/"* && "${PACK_REFRESH_RUNS}" != "1" ]]; then
    echo "[keep existing] ${dest#"${ROOT}/"}"
    return 0
  fi
  mkdir -p "$(dirname -- "${dest}")"
  cp -a "${src}" "${dest}"
  chmod 0644 "${dest}"   # data only; avoids spurious 100644->100755 churn in git
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
  for f in evolving_geodesic.json confinement.dat \
           freefall_observer_timing.json freefall_observer_convergence.json; do
    [[ -e "${src}/small_data/${f}" ]] && copy_file "${src}/small_data/${f}" "${dest}/small_data/${f}"
  done
  for f in constraint_norms.dat collapse_diagnostics.dat; do
    [[ -e "${src}/data/${f}" ]] && copy_file "${src}/data/${f}" "${dest}/data/${f}"
  done
}

# Stationary geometry-first atlas (Sec. VI.D, Table VII).
# trajectory.jsonl carries per-eval f_geo, descriptors and E_-, which reproduces
# every quoted atlas statistic; elites/*.json carry the winning genomes.
# The 74 MB-per-eval *.gridinit grids and the eval-level dumps are skipped.
copy_atlas_light() {
  local name="$1"
  local src="${RUNS}/geometry_atlas/${name}"
  local dest="${DEST}/runs/geometry_atlas/${name}"
  if [[ ! -d "${src}" ]]; then
    echo "[skip missing atlas campaign] ${src}" >&2
    return 0
  fi
  mkdir -p "${dest}"
  for f in archive.json metadata.json state.json summary.json trajectory.jsonl; do
    [[ -e "${src}/${f}" ]] && copy_file "${src}/${f}" "${dest}/${f}"
  done
  if [[ -d "${src}/elites" ]]; then
    mkdir -p "${dest}/elites"
    while IFS= read -r -d '' f; do
      copy_file "${f}" "${dest}/elites/$(basename -- "${f}")"
    done < <(find "${src}/elites" -maxdepth 1 -name '*.json' -print0)
  fi
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

# Load private site knobs if present (gitignored) so env-path tokens are known.
if [[ -f "${ROOT}/grteclyn-wrapper/.env" ]]; then
  # shellcheck disable=SC1091
  set -a
  # Only export the path knobs we care about; do not print them.
  # shellcheck disable=SC1090
  source <(grep -E '^(GRTRESNA_ENV|OPENMPI_ROOT|CHOMBO_HOME|SIM_ROOT|GRTECLYN_ROOT)=' \
    "${ROOT}/grteclyn-wrapper/.env" | sed 's/\r$//' || true)
  set +a
fi

echo "Packing publishable results -> ${DEST}"
mkdir -p "${DEST}"
# Only wipe trees that are always fully regenerable from the repo itself.
# runs/ is deliberately NOT wiped: source runs under /runs get pruned to reclaim
# disk, and for those the committed extract here is the last surviving copy
# (e.g. the RF/DS/DL/pump-free ladders behind Table V). Stale run extracts must
# be removed by hand, never by a re-pack.
for d in article figures manifests validation; do
  rm -rf "${DEST}/${d}"
done
rm -f "${DEST}/PLOT_DATA_SOURCES.txt"

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
  mechanism_ablation.txt \
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

if [[ -d "${SRC_PAPER}/validation" ]]; then
  mkdir -p "${DEST}/validation"
  cp -a "${SRC_PAPER}/validation/." "${DEST}/validation/"
  while IFS= read -r -d '' f; do
    chmod 0644 "${f}"
    scrub_machine_paths "${f}"
  done < <(find "${DEST}/validation" -type f -print0)
  echo "[ok] results/${ARTICLE_SLUG}/validation/"
fi

copy_file \
  "${ROOT}/grteclyn-wrapper/scripts/campaigns/promote/bicomplex_cmaes_v1/manifest_pumpfree.json" \
  "${DEST}/manifests/manifest_pumpfree.json"

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
  bcma_pfrm_L128_N256_t30_hq_eval000146 \
  bcma_rc_freefall_corrected_L128_N192_t30_hq_eval000146 \
  bcma_rm_freefall_corrected_L128_N256_t30_hq_eval000146 \
  bcma_rf_freefall_corrected_L128_N384_t30_hq_eval000146
do
  copy_hq_light "${name}"
done

for name in \
  geometry_atlas_topologies_n64_L32_e260_20260723 \
  geometry_atlas_maxfgeo_n64_L32_e500_20260723 \
  geometry_atlas_lowexotic_n64_L32_e500_20260723 \
  geometry_atlas_breadth_n64_L32_e220_20260723
do
  copy_atlas_light "${name}"
done

copy_file "${SRC_PAPER}/article/PLOT_DATA_SOURCES.txt" "${DEST}/PLOT_DATA_SOURCES.txt"
copy_file "${SRC_PAPER}/publishable_results_README.md" "${DEST}/README.md"

while IFS= read -r -d '' f; do
  scrub_machine_paths "${f}"
done < <(find "${DEST}" -type f -print0)

assert_pack_has_no_machine_identity

du -sh "${DEST}"
echo "Done. See ${DEST}/README.md for the audit layout."
