/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#include "AHJacobianOp.hpp"

#include <AMReX_BLassert.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_RealBox.H>

#include <algorithm>

AHJacobianOp::AHJacobianOp(int n_rings, int ring_size, int max_grid_size)
    : m_n_rings(n_rings), m_ring_size(ring_size),
      m_num_particles(n_rings * ring_size)
{
    using namespace amrex;

    // Cell-centred ring grid: latitude along dim 0, longitude (phi) along
    // dim 1, degenerate dim 2.
    const Box domain(IntVect(0, 0, 0),
                     IntVect(m_n_rings - 1, m_ring_size - 1, 0));

    // Unit index-space extent; the operator works in grid index space, so the
    // physical coordinates are immaterial. Only phi (dim 1) is periodic.
    const RealBox real_box({AMREX_D_DECL(0.0, 0.0, 0.0)},
                           {AMREX_D_DECL(1.0, 1.0, 1.0)});
    const Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 1, 0)};
    m_geom = Geometry(domain, real_box, CoordSys::cartesian, is_periodic);

    // Split along theta only, keeping each ring (phi) contiguous, then let
    // AMReX distribute the boxes across ranks.
    m_ba.define(domain);
    m_ba.maxSize(IntVect(max_grid_size, m_ring_size, 1));
    m_dm.define(m_ba);

    // Single level, no coarsening: the antipodal/periodic topology lives in
    // the mat-vec, not in a geometric hierarchy.
    LPInfo info;
    info.setMaxCoarseningLevel(0);
    define({m_geom}, {m_ba}, {m_dm}, info);

    // Phi is periodic; theta (and the degenerate z) are registered as
    // homogeneous Neumann for bookkeeping only -- the mat-vec never reads the
    // ghost cells these fill.
    setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann, LinOpBCType::Periodic,
                              LinOpBCType::Neumann)},
                {AMREX_D_DECL(LinOpBCType::Neumann, LinOpBCType::Periodic,
                              LinOpBCType::Neumann)});
    setLevelBC(0, nullptr);
}

amrex::MultiFab AHJacobianOp::make_mf() const
{
    return amrex::MultiFab(m_ba, m_dm, 1, 1);
}

void AHJacobianOp::flat_to_mf(const std::vector<double> &flat,
                              amrex::MultiFab &mf) const
{
    using namespace amrex;
    AMREX_ALWAYS_ASSERT(static_cast<int>(flat.size()) == m_num_particles);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box &bx           = mfi.validbox();
        const Array4<Real> &arr = mf.array(mfi);
        const auto lo           = lbound(bx);
        const auto hi           = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                arr(i, j, 0) = flat[i * m_ring_size + j];
            }
        }
    }
}

void AHJacobianOp::mf_to_flat(const amrex::MultiFab &mf,
                              std::vector<double> &flat) const
{
    using namespace amrex;
    flat.assign(m_num_particles, 0.0);

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box &bx                 = mfi.validbox();
        const Array4<const Real> &arr = mf.const_array(mfi);
        const auto lo                 = lbound(bx);
        const auto hi                 = ubound(bx);
        for (int j = lo.y; j <= hi.y; ++j)
        {
            for (int i = lo.x; i <= hi.x; ++i)
            {
                flat[i * m_ring_size + j] = arr(i, j, 0);
            }
        }
    }

    // Each cell is owned by exactly one rank, so summing across ranks leaves
    // every rank holding the complete array.
    ParallelDescriptor::ReduceRealSum(flat.data(), m_num_particles);
}

void AHJacobianOp::applyBC(int /*amrlev*/, int /*mglev*/,
                           amrex::MultiFab & /*in*/, BCMode /*bc_mode*/,
                           StateMode /*s_mode*/,
                           const amrex::MLMGBndryT<amrex::MultiFab> * /*bndry*/,
                           bool /*skip_fillboundary*/) const
{
    // Intentionally empty: see the declaration in AHJacobianOp.hpp.
}

void AHJacobianOp::Fapply(int /*amrlev*/, int /*mglev*/, amrex::MultiFab &out,
                          const amrex::MultiFab &in) const
{
    AMREX_ALWAYS_ASSERT(m_matvec);

    std::vector<double> in_flat;
    mf_to_flat(in, in_flat);

    std::vector<double> out_flat(m_num_particles, 0.0);
    m_matvec(in_flat, out_flat);

    flat_to_mf(out_flat, out);
}

void AHJacobianOp::Fsmooth(int /*amrlev*/, int /*mglev*/,
                           amrex::MultiFab & /*sol*/,
                           const amrex::MultiFab & /*rhs*/,
                           int /*redblack*/) const
{
    amrex::Abort("AHJacobianOp::Fsmooth: relaxation is not implemented; the "
                 "single-level solve must use a bottom solver (BiCGStab) with "
                 "coarsening disabled");
}

void AHJacobianOp::FFlux(
    int /*amrlev*/, const amrex::MFIter & /*mfi*/,
    const amrex::Array<amrex::FArrayBox *, AMREX_SPACEDIM> & /*flux*/,
    const amrex::FArrayBox & /*sol*/, Location /*loc*/,
    int /*face_only*/) const
{
    // Single level: no inter-level fluxes are ever requested.
}
