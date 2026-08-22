#!/usr/bin/env bash
# Pack light extracts of the convergence / big-box campaign into
# results/bondi-dipole-runaway/campaign/.
#
# Companion to pack_results.sh, which owns the dx = 0.5 production cells.
# Both now write into campaign/, under disjoint cell names: this script owns
# convA_*/boxC_*, pack_results.sh owns the unprefixed pair_*/single_* cells.
# Neither touches stars/, figures/ or movies/.
#
# Safe to re-run at any time, including while runs are live:
#   - each campaign cell's output folder is rebuilt from scratch on every run,
#     so partial extracts from a still-running cell are replaced cleanly;
#   - cells are discovered dynamically (every top-level dir in
#     runs/bondi/rerun except published/, experiments/, logs/), so the queued
#     equal-mass series is picked up automatically when it appears;
#   - the follow-up campaign in runs/bondi/rerun/next/ is packed too, under a
#     "next_" prefix, so it stays visibly separate from the published series;
#   - cells with no data yet are skipped with a note, never an error.
#
# Machine identity is scrubbed at runtime by the shared scrubber in the
# wrapper package (grteclyn_wrapper.packaging.scrub_paths) -- see the
# no-identity-in-git rule); PNG frames carry no paths.
#
# Usage: bash research/bondi_dipole/pack_campaign.sh
# --- SUPERSEDED 2026-08-22 -------------------------------------------------
# This script belongs to a finished campaign and reads/writes runs/bondi/rerun/
# which no longer exists.  Running it now would either fail halfway or
# recreate a dead tree and mix superseded cells into the live campaign.
# What replaced it: research/bondi_dipole/pack_runaway.sh, which packs the current tree.
# Kept for the reasoning in the header above.  Remove this block only if you
# deliberately want the old tree back.
printf '%s: superseded campaign; %s no longer exists.\n' "$(basename "$0")" "runs/bondi/rerun/" >&2
printf '  see research/bondi_dipole/docs/GPU_RUN_PAPER.md for the live campaign\n' >&2
exit 2
# ---------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/../.." && pwd)"
SIM_ROOT="$(cd -- "${ROOT}/.." && pwd)"
RUNS="${ROOT}/runs/bondi/rerun"
DEST="${ROOT}/results/bondi-dipole-runaway"

export ROOT SIM_ROOT
scrub() { python3 "${ROOT}/grteclyn-wrapper/src/grteclyn_wrapper/packaging/scrub_paths.py" "$@"; }

mkdir -p "${DEST}/campaign"

# Cell list: the convergence/big-box cells at the top level, plus the follow-up
# campaign one level down in next/.  The follow-up cells are packed under a
# "next_" prefix so a reader never mistakes them for the published convergence
# series -- they answer different questions and were run with different knobs
# (scaled solve tolerance, boundary sponge, a second phantom frequency).
cells=()
for celldir in "${RUNS}"/*/; do
  cell="$(basename "${celldir}")"
  case "${cell}" in published|experiments|logs|next) continue ;; esac
  cells+=("${cell}|${celldir}")
done
for celldir in "${RUNS}"/next/*/; do
  [[ -d "${celldir}" ]] || continue
  cell="$(basename "${celldir}")"
  [[ "${cell}" == _queue ]] && continue
  cells+=("next_${cell}|${celldir}")
done

for entry in "${cells[@]}"; do
  cell="${entry%%|*}"
  celldir="${entry#*|}"
  run="$(find "${celldir}" -maxdepth 1 -type d -name 'bondi_sg_*' | head -n 1)"
  if [[ -z "${run}" ]] || ! compgen -G "${run}/small_data/*.dat" > /dev/null; then
    echo "[pack-campaign] ${cell}: no time series yet -- skipping"
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
    [[ -f "${run}/${f}" ]] && cp "${run}/${f}" "${out}/"
  done

  # Evolution-side diagnostics arrive at every-step cadence; downsample to
  # dt = 0.5 exactly as pack_results.sh does for the published cells.
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

  found_txt=$(find "${out}" -maxdepth 1 \( -name '*.txt' -o -name '*.json' \))
  [[ -n "${found_txt}" ]] && scrub ${found_txt}

  # Curated frames: one every FRAME_STRIDE in simulation time (plus the final
  # frame), so a t = 60 cell keeps 7 per field and a t = 90 cell keeps 10.
  # Named by simulation time, matching pack_results.sh.  A plotfile's time is
  # step * dt with dt = dt_multiplier * L_full / N_full, all three read back
  # from the params file just packed; if any is missing the old
  # fraction-of-the-run selection is used instead, named by step.
  for field in scalar_activity_z chi_minus_1_z; do
    fdir="${run}/frames/${field}/frames"
    [[ -d "${fdir}" ]] || continue
    FIELD="${field}" OUT="${out}/frames" PARAMS="${out}/evolution_params.txt" \
      FRAME_STRIDE="${FRAME_STRIDE:-10}" python3 - "${fdir}" <<'PY'
import os, pathlib, re, shutil, sys

step_of = lambda p: int(re.search(r"(\d+)\s*$", p.stem).group(1))
frames = sorted(pathlib.Path(sys.argv[1]).glob("*.png"), key=step_of)
if not frames:
    raise SystemExit(0)

field = os.environ["FIELD"]
out = pathlib.Path(os.environ["OUT"])
stride = float(os.environ["FRAME_STRIDE"])

need = {"dt_multiplier", "L_full", "N_full"}
vals = {}
params = pathlib.Path(os.environ["PARAMS"])
if params.is_file():
    for line in params.read_text(encoding="utf-8").splitlines():
        m = re.match(r"\s*(dt_multiplier|L_full|N_full)\s*=\s*([0-9.eE+-]+)", line)
        if m:
            vals[m.group(1)] = float(m.group(2))

out.mkdir(parents=True, exist_ok=True)
if need <= set(vals) and vals["N_full"]:
    dt = vals["dt_multiplier"] * vals["L_full"] / vals["N_full"]
    times = [step_of(f) * dt for f in frames]
    targets = [k * stride for k in range(int(times[-1] // stride) + 1)]
    if times[-1] - targets[-1] >= stride / 2:   # keep the end of a run that
        targets.append(times[-1])               # stops between two decades
    picked = []
    for want in targets:
        i = min(range(len(times)), key=lambda j: abs(times[j] - want))
        if i not in picked:
            picked.append(i)
    for i in picked:
        shutil.copy2(frames[i], out / f"{field}_t{times[i]:04.0f}.png")
else:
    n = len(frames)
    for frac in (0.0, 1 / 3, 2 / 3, 1.0):
        i = min(n - 1, round(frac * (n - 1)))
        shutil.copy2(frames[i], out / f"{field}_step{step_of(frames[i]):05d}.png")
PY
  done

  echo "[pack-campaign] campaign/${cell}: $(find "${out}" -type f | wc -l) files"
done

# Derived tables: each analysis owns its own module and its own output files,
# and is rebuilt from whatever campaign cells are packed.
python3 "${DEST}/analysis/convergence_check.py" "${DEST}"
python3 "${DEST}/analysis/wave_check.py" "${DEST}"
python3 "${DEST}/analysis/constraint_check.py" "${DEST}"

echo "[pack-campaign] total campaign size: $(du -sh "${DEST}/campaign" | cut -f1)"
echo "[pack-campaign] done"
