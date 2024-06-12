#ifndef DERIVE_H_
#define DERIVE_H_

#include <AMReX_BLFort.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <InitialConditions.H>

void derive_func_fab(const amrex::Box &bx, amrex::FArrayBox &derfab, int dcomp,
                     int /*numcomp*/, const amrex::FArrayBox &datfab,
                     const amrex::Geometry &geom, const amrex::Real time,
                     const int * /*bcomp*/, int /*scomp*/);

#endif /* DERIVE_H_ */
