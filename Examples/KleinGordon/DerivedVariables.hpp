#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>

#include "VarsTools.hpp"

#include "InitialConditions.hpp"

void calc_derive_mf(amrex::MultiFab &mf_out, int dcomp, int /*numcomp*/,
                    const amrex::MultiFab &mf_in, const amrex::Geometry &geom,
                    const amrex::Real time, const int * /*bcomp*/,
                    int /*scomp*/);
