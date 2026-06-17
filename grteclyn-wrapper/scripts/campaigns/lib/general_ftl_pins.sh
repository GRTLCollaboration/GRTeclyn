#!/usr/bin/env bash
# v20 general_ftl matter-class pin bundles (MapElites.md). Sourced by general_ftl/run_all.sh.

ftl_general_ftl_static_inert() {
  printf '%s' "grtresna_shell_toroidal_velocity=0 grtresna_shell_poloidal_velocity=0 grtresna_shell_radial_velocity=0 grtresna_shell_omega=0"
}

ftl_general_ftl_wormhole_pins() {
  local inert
  inert="$(ftl_general_ftl_static_inert)"
  printf '%s' "grtresna_matter_layout=2 grtresna_shell_axis_theta=1.5708 grtresna_shell_axis_phi=0 grtresna_shell_static=1 ${inert}"
}

ftl_general_ftl_ring_pins() {
  local inert
  inert="$(ftl_general_ftl_static_inert)"
  printf '%s' "grtresna_matter_layout=3 grtresna_shell_static=1 ${inert}"
}

ftl_general_ftl_spin_pins() {
  printf '%s' "grtresna_matter_layout=0 grtresna_shell_static=0 grtresna_shell_toroidal_velocity=0 grtresna_shell_poloidal_velocity=0 grtresna_shell_radial_velocity=0 grtresna_shift_seed=0"
}
