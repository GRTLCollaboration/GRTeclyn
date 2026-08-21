# shellcheck shell=bash
# Bondi dipole -- frame colour scales, sourced by the bondi launchers.
#
# WHY THIS EXISTS
# The frame renderer's built-in colour limits (visualisation/process_wave/
# consume_plotfiles/config.py) were calibrated for the search campaign's matter
# rung -- scalar_lambda=640, mu=85333 -- whose lumps are several times taller
# than the stars this campaign evolves.  On the bondi rung (lambda=10240,
# mu=21845333) EVERY star peaks at phi_c ~ 0.0197 and is essentially frequency
# independent: the family scan (results/bondi-dipole-runaway/stars/
# star_radius.csv) gives 0.019695 at omega=0.55, 0.019912 at 0.615, 0.019720 at
# 0.75 and 0.019126 at 0.80, with the phantom sector inside 1% of the canonical
# at every rung -- a 4% spread across the whole family.  Against the built-in
# activity window of 0..0.20 that renders the star at 17% of the colour range,
# i.e. a near-black smudge on a black field, and it has been doing so since the
# first bondi cell ran -- the published omega=0.55 movies are washed out too.
#
# Two of the built-ins fail in the OTHER direction.  phi and Pi are the two
# quadratures of one rotating complex field, so whichever is near zero the
# other is near phi_c; measured live at t=22.4 in single_p_w080, phi peaked at
# 0.0035 while Pi peaked at 0.0151 against a Pi window of +-0.01 -- clipped.
# A single shared scale for all quadrature channels is the only correct choice:
# a phase sweep must not change how bright the star looks.
#
# The renderer already exposes GRTECLYN_FRAMES_ZLIM_<FIELD>=lo,hi, and an
# explicit override there beats both the built-in preset and any per-frame
# rescaling.  This module is only the bondi-specific policy for those knobs; it
# changes no shared default, so every other campaign renders exactly as before.
#
# WHY NOT PER-FRAME AUTO-SCALING
# GRTECLYN_FRAMES_AUTO_ZLIM=1 would make the lump visible in one line, and it
# is the wrong tool here.  This campaign's question is "does the star stay the
# same size for 120 time units".  A scale that re-fits itself to each frame
# renders a star that is quietly evaporating exactly like one that is holding
# still, because the colour range shrinks with the signal.  The scale has to be
# fixed for the comparison to mean anything.
#
# KNOBS
#   BONDI_FRAME_CONTRAST=0    leave every colour scale at the library default
#   BONDI_FRAME_PHI_MAX=A     saturation amplitude for the field quadratures
#                             (default 0.025 ~ 1.25 x phi_c, so a quadrature at
#                             full swing fills ~80% of the bar).  Tight on
#                             purpose: this screen is looking for a star that
#                             slowly changes size, and a generous scale is
#                             exactly what hides a slow change.  A collapse
#                             spike will clip, which is a legible signal in
#                             itself, not a loss.
#   BONDI_FRAME_RHO_MAX=R     saturation for rho_req (default 1.5e-3; measured
#                             peak on this rung is 8.1e-4 against a built-in
#                             window of 3.0e-3)
#   BONDI_FRAME_ZOOM=W        render a W-wide window instead of the whole box.
#                             Off by default: it makes the star bigger on
#                             screen but crops the domain, and the framing is
#                             then no longer comparable with the packed cells.
#                             W=32 is the natural choice for a lone star -- it
#                             just contains the r99 tail (~14) around a core
#                             parked at x=+4.
# A knob already set in the environment is left alone, so a caller can scale a
# single field by hand without disabling the rest.

bondi_frame_contrast_env() {
  [[ "${BONDI_FRAME_CONTRAST:-1}" != "0" ]] || return 0

  local a r act
  a="${BONDI_FRAME_PHI_MAX:-0.025}"
  r="${BONDI_FRAME_RHO_MAX:-1.5e-3}"
  # scalar_activity sums sqrt(phi^2 + Pi^2) over the bicomplex channels.  It is
  # phase independent (that is the point of it) and sits at ~1.75 phi_c:
  # measured 0.0344 live in single_p_w080.  1.6 x the quadrature scale puts
  # that peak at ~86% of the bar.
  act="$(awk -v a="${a}" 'BEGIN{printf "%.10g", 1.6*a}')"

  # Set one GRTECLYN_FRAMES_ZLIM_<FIELD>, unless the caller already set it.
  _bondi_zlim() {
    local var="GRTECLYN_FRAMES_ZLIM_$1"
    [[ -n "${!var:-}" ]] && return 0
    printf -v "${var}" '%s' "$2"
    export "${var?}"
  }

  local f
  for f in PHI PI PHI_LUMP0 PI_LUMP0 PHI_LUMP1 PI_LUMP1 PHI_LUMP2 PI_LUMP2 \
           PHI_LUMP_SUM PI_LUMP_SUM; do
    _bondi_zlim "${f}" "-${a},${a}"
  done
  for f in SCALAR_ACTIVITY LUMP_ACTIVITY; do
    _bondi_zlim "${f}" "0,${act}"
  done
  _bondi_zlim RHO_REQ "-${r},${r}"

  unset -f _bondi_zlim

  if [[ -n "${BONDI_FRAME_ZOOM:-}" ]]; then
    export GRTECLYN_FRAMES_ZOOM="${BONDI_FRAME_ZOOM}"
  fi

  echo "[bondi] frame contrast: quadratures +-${a}, activity 0..${act}, rho_req +-${r}, zoom=${BONDI_FRAME_ZOOM:-full box}"
}
