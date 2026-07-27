#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4StateVariables.hpp"
#include "GRTresnaScalarLayout.hpp"

enum
{
    c_phi  = NUM_CCZ4_VARS,
    c_Pi,
    // Single complex scalar reuses the first lump slots (phi_lump0 / Pi_lump0):
    //   c_phi = phi_re+, c_Pi = Pi_re+, c_phi2 = phi_im+, c_Pi2 = Pi_im+.
    c_phi2 = NUM_CCZ4_VARS + 2,
    c_Pi2  = NUM_CCZ4_VARS + 3,

    // Two-complex-field (canonical + phantom) model: the PHANTOM (sign -1)
    // field reuses the phi_lump1 / Pi_lump1 / phi_lump2 / Pi_lump2 slots.  The
    // canonical (sign +1) field uses the four slots above.  The two matter
    // models are mutually exclusive at runtime, so the slot reuse is safe.
    c_phi_m  = NUM_CCZ4_VARS + 4, // phi_re-  (= phi_lump1 slot)
    c_Pi_m   = NUM_CCZ4_VARS + 5, // Pi_re-   (= Pi_lump1  slot)
    c_phi2_m = NUM_CCZ4_VARS + 6, // phi_im-  (= phi_lump2 slot)
    c_Pi2_m  = NUM_CCZ4_VARS + 7, // Pi_im-   (= Pi_lump2  slot)

    // Controller energy-momentum reservoir (ρ_c, j_{c i}). Always present in
    // the state vector; remain identically zero unless controller_reservoir_mode
    // >= 1. Appended after the lump slots so existing gridinit files (which
    // write fewer components) leave them at the pre-zero from initData.
    c_rho_ctrl = NUM_CCZ4_VARS + 2 + 2 * GRTRESNA_MAX_INDEPENDENT_SCALARS,
    c_jctrl1,
    c_jctrl2,
    c_jctrl3,

    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {
    "phi", "Pi",       "phi_lump0", "Pi_lump0", "phi_lump1", "Pi_lump1",
    "phi_lump2", "Pi_lump2", "phi_lump3", "Pi_lump3", "phi_lump4", "Pi_lump4",
    "rho_ctrl", "jctrl1", "jctrl2", "jctrl3"};
static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, 16> additional_parities = {
    BCParity::even, BCParity::even, BCParity::even, BCParity::even,
    BCParity::even, BCParity::even, BCParity::even, BCParity::even,
    BCParity::even, BCParity::even, BCParity::even, BCParity::even,
    BCParity::even,   // rho_ctrl
    BCParity::odd_x,  // jctrl1
    BCParity::odd_y,  // jctrl2
    BCParity::odd_z}; // jctrl3
static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4StateVariables::parities, additional_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */
