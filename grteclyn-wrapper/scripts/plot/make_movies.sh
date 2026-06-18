#!/usr/bin/env bash
# Stitch pre-rendered PNG frames into mp4 movies for one or more episode dirs.
#
# Robust to gapped frame numbering (frames are named by sim step, e.g.
# frame_z_0000.png, frame_z_0015.png, ...), which the %04d-based make_movies.py
# cannot handle.  Uses ffmpeg's glob demuxer so any sorted PNG sequence works.
#
# Usage:
#   make_movies.sh EPISODE_DIR [EPISODE_DIR ...] [--framerate N] [--only chi_z K_z]
#
# For each EPISODE_DIR it looks under <EPISODE_DIR>/frames/<field>_<axis>/frames/
# and writes <EPISODE_DIR>/movies/movie_<field>_<axis>.mp4
set -euo pipefail

FRAMERATE=10
ONLY=()
DIRS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --framerate) FRAMERATE="$2"; shift 2;;
    --only) shift; while [[ $# -gt 0 && "$1" != --* ]]; do ONLY+=("$1"); shift; done;;
    *) DIRS+=("$1"); shift;;
  esac
done

if [[ ${#DIRS[@]} -eq 0 ]]; then
  echo "Usage: $0 EPISODE_DIR [EPISODE_DIR ...] [--framerate N] [--only chi_z K_z]" >&2
  exit 2
fi

command -v ffmpeg >/dev/null 2>&1 || { echo "ffmpeg not found on PATH" >&2; exit 1; }

want() {  # $1 = folder name; returns 0 if selected
  [[ ${#ONLY[@]} -eq 0 ]] && return 0
  local f; for f in "${ONLY[@]}"; do [[ "$f" == "$1" ]] && return 0; done; return 1
}

made=0
for ep in "${DIRS[@]}"; do
  ep="${ep%/}"
  froot="${ep}/frames"
  movies_dir="${ep}/movies"
  [[ -d "$froot" ]] || { echo "[skip] no frames dir: $froot" >&2; continue; }
  mkdir -p "$movies_dir"
  for fd in "$froot"/*/; do
    field_axis="$(basename "$fd")"
    frames_dir="${fd%/}/frames"
    [[ -d "$frames_dir" ]] || continue
    want "$field_axis" || continue
    shopt -s nullglob
    pngs=("$frames_dir"/*.png)
    shopt -u nullglob
    [[ ${#pngs[@]} -gt 0 ]] || continue
    out="${movies_dir}/movie_${field_axis}.mp4"
    echo "[movie] $ep :: $field_axis (${#pngs[@]} frames) -> movies/$(basename "$out")"
    ffmpeg -y -loglevel error -framerate "$FRAMERATE" \
      -pattern_type glob -i "${frames_dir}/*.png" \
      -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" \
      -c:v libx264 -pix_fmt yuv420p "$out"
    made=$((made+1))
  done
done
echo "[done] wrote $made movie(s)"
