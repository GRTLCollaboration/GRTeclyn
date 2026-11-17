/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef NULLBCFILL_HPP_
#define NULLBCFILL_HPP_

#include <AMReX_PhysBCFunct.H>

// AMReX requires us to set a boundary fill function even though we are not
// using Dirichlet or "User" BCs (see docs)

struct NullBCFill
{
    AMREX_GPU_DEVICE
    void operator()(const amrex::IntVect & /*iv*/,
                    const amrex::Array4<amrex::Real> & /*dest*/,
                    const int /*dcomp*/, const int /*numcomp*/,
                    const amrex::GeometryData & /*geom*/,
                    const amrex::Real /*time*/, const amrex::BCRec * /*bcr*/,
                    const int /*bcomp*/, const int /*orig_comp*/) const
    {
        // We don't have any external Dirichlet BC
    }
};

void null_bc_fill(const amrex::Box &box, amrex::FArrayBox &data,
                  const int dcomp, const int numcomp,
                  const amrex::Geometry &geom, const amrex::Real time,
                  const amrex::Vector<amrex::BCRec> &bcr, const int bcomp,
                  const int scomp)
{
    amrex::GpuBndryFuncFab<NullBCFill> bndry_func(NullBCFill{});
    bndry_func(box, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}

#endif /* NULLBCFILL_HPP_ */