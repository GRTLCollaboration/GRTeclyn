#ifndef DERIVEDVARIABLES_HPP_
#define DERIVEDVARIABLES_HPP_

// AMReX includes
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

// GRTeclyn includes
#include "Coordinates.hpp"

// KleinGordon includes
#include "KleinGordonRHS.hpp"
#include "SineGordon.hpp"
#include "Wave.hpp"

// For unused variables in the function signature, you don't have to name them
// and it is best practice to comment them out. However, they must still be
// included in any derived variable function definition because AMReX expects a
// certain function signature

// NB: The analytic solution doesn't use the input MultiFab mf_in - if you want
// to use the current state variables, you will need to uncomment this.
AMREX_FORCE_INLINE void
calc_analytic_solution(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                       const amrex::MultiFab & /*mf_in*/,
                       const amrex::Geometry &geom, const amrex::Real time,
                       const int * /*bcomp*/, int /*scomp*/);

template <typename model_t>
AMREX_FORCE_INLINE void
calc_energy_density(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                    const amrex::MultiFab &mf_in,
                    const amrex::Geometry /*&geom*/, const amrex::Real /*time*/,
                    const int * /*bcomp*/, int /*scomp*/);

template <typename model_t>
AMREX_FORCE_INLINE void calc_analytic_mf_3d(amrex::MultiFab &mf_out, int dcomp,
                                            const amrex::Geometry &geom,
                                            const amrex::Real time);

template <typename model_t>
AMREX_FORCE_INLINE void calc_analytic_mf_1d(amrex::MultiFab &mf_out, int dcomp,
                                            const amrex::Geometry &geom,
                                            const amrex::Real time);

// This needs to be a .hpp file because there is a template definition inside
// The template arguments are expanded inline by the compiler in the .hpp file.
#include "DerivedVariables.impl.hpp"

#endif // DERIVEDVARIABLES_HPP_
