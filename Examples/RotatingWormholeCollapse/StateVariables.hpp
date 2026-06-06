#ifndef STATEVARIABLES_HPP
#define STATEVARIABLES_HPP

#include "ArrayTools.hpp"
#include "CCZ4StateVariables.hpp"

enum
{
    c_phi = NUM_CCZ4_VARS,
    c_Pi,
    c_dust_rho,
    c_dust_v1,
    c_dust_v2,
    c_dust_v3,
    c_teo_rho,
    c_teo_j1,
    c_teo_j2,
    c_teo_j3,
    c_teo_S11,
    c_teo_S12,
    c_teo_S13,
    c_teo_S22,
    c_teo_S23,
    c_teo_S33,
    
    NUM_VARS
};

namespace StateVariables
{
static const amrex::Vector<std::string> additional_names = {
    "phi",     "Pi",      "dust_rho", "dust_v1", "dust_v2", "dust_v3",
    "teo_rho", "teo_j1", "teo_j2", "teo_j3",
    "teo_S11", "teo_S12", "teo_S13", "teo_S22", "teo_S23", "teo_S33"};
static const amrex::Vector<std::string> names =
    ArrayTools::concatenate(CCZ4StateVariables::names, additional_names);

static const std::array<BCParity, 16> additional_parities = {
    BCParity::even,   // phi
    BCParity::even,   // Pi
    BCParity::even,   // dust_rho
    BCParity::odd_x,  // dust_v1
    BCParity::odd_y,  // dust_v2
    BCParity::odd_z,  // dust_v3
    BCParity::even,   // teo_rho
    BCParity::odd_x,  // teo_j1
    BCParity::odd_y,  // teo_j2
    BCParity::odd_z,  // teo_j3
    BCParity::even,   // teo_S11
    BCParity::odd_xy, // teo_S12
    BCParity::odd_xz, // teo_S13
    BCParity::even,   // teo_S22
    BCParity::odd_yz, // teo_S23
    BCParity::even};  // teo_S33
static const std::array<BCParity, NUM_VARS> parities =
    ArrayTools::concatenate(CCZ4StateVariables::parities, additional_parities);
} // namespace StateVariables

#endif /* STATEVARIABLES_HPP */