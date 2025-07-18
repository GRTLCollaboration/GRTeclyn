#ifndef DERIVEDVARIABLES_HPP_
#define DERIVEDVARIABLES_HPP_

// AMReX includes
#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

// GRTeclyn includes
#include "VarsTools.hpp"

// KleinGordon includes
#include "SineGordon.hpp"
#include "Wave.hpp"

void calc_analytic_solution(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                            const amrex::MultiFab & /*mf_in*/,
                            const amrex::Geometry &geom, const amrex::Real time,
                            const int * /*bcomp*/, int /*scomp*/);

#endif // DERIVEDVARIABLES_HPP_
