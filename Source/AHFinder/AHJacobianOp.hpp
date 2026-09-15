/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHJACOBIANOP_HPP_
#define AHJACOBIANOP_HPP_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MultiFab.H>

#include <functional>
#include <vector>

// Matrix-free linear operator for one pseudo-transient-continuation (PTC) step
// of AHFinder. Each PTC step solves
//   (I/dt + J) delta_h = -Theta(h_n),   J = dTheta/dh
// over the theta/phi ring grid, with J the frozen-metric Jacobian applied
// through a finite-difference directional derivative. This class plugs that
// operator into AMReX's MLMG framework so the system can be solved by a
// matrix-free BiCGStab (see AHFinder::find()).
//
// The ring grid is laid out as a single-level, cell-centred domain of shape
// (n_rings, ring_size, 1): dimension 0 is latitude (theta), dimension 1 is
// longitude (phi) and is periodic, dimension 2 is degenerate. The antipodal
// coupling at the poles and the phi periodicity are handled inside the mat-vec
// (which reuses AHGeometry via the supplied callback on complete flat arrays),
// so the solver itself never has to represent that topology: coarsening is
// disabled and the theta boundary is registered as homogeneous Neumann purely
// for MLMG bookkeeping. Its ghost cells are never read for physics.
//
// The operator owns its own solver-side decomposition (a BoxArray split along
// theta with phi kept contiguous, plus an AMReX-built DistributionMapping),
// independent of AHFinder's particle decomposition. mf_to_flat()/flat_to_mf()
// bridge the two: because the mat-vec leaves the complete array on every rank
// (via a reduction), gather copies owned cells out and reduces to a full-length
// vector, and scatter fills only owned cells.
class AHJacobianOp : public amrex::MLCellLinOp
{
  public:
    // Signature of the mat-vec supplied by AHFinder. Given the complete
    // input vector "in" (length n_rings * ring_size, identical on every rank),
    // it must fill "out" with (I/dt + J) in, also complete on every rank.
    using MatVec = std::function<void(const std::vector<double> &in,
                                      std::vector<double> &out)>;

    // Build the operator and its single-level grid from the ring dimensions.
    // max_grid_size caps the theta extent of each box (phi is kept whole).
    AHJacobianOp(int n_rings, int ring_size, int max_grid_size = 8);

    ~AHJacobianOp() override = default;

    AHJacobianOp(const AHJacobianOp &)            = delete;
    AHJacobianOp(AHJacobianOp &&)                 = delete;
    AHJacobianOp &operator=(const AHJacobianOp &) = delete;
    AHJacobianOp &operator=(AHJacobianOp &&)      = delete;

    // Bind the mat-vec used by Fapply. AHFinder rebinds this each PTC step so
    // the callback closes over the current frozen metric and pseudo-timestep.
    void set_matvec(MatVec matvec) { m_matvec = std::move(matvec); }

    // Solver-side layout, so AHFinder can allocate matching MultiFabs.
    [[nodiscard]] const amrex::BoxArray &grids() const { return m_ba; }
    [[nodiscard]] const amrex::DistributionMapping &dmap() const
    {
        return m_dm;
    }
    [[nodiscard]] const amrex::Geometry &geom() const { return m_geom; }

    // A MultiFab on the operator's layout (1 component, 1 ghost cell).
    [[nodiscard]] amrex::MultiFab make_mf() const;

    // Scatter a complete flat vector into the owned cells of mf.
    void flat_to_mf(const std::vector<double> &flat, amrex::MultiFab &mf) const;

    // Gather the owned cells of mf into a complete flat vector: every rank
    // holds the full result after the reduction.
    void mf_to_flat(const amrex::MultiFab &mf,
                    std::vector<double> &flat) const;

  protected:
    // No-op boundary fill. Fapply reconstructs the complete input from valid
    // cells only (mf_to_flat) and handles the pole/phi topology itself, so it
    // never reads ghost cells; there is nothing for applyBC to do. Overriding
    // it also avoids MLCellLinOp's default, which calls a Fortran kernel that
    // is unavailable in this Fortran-free build.
    void applyBC(int amrlev, int mglev, amrex::MultiFab &in, BCMode bc_mode,
                 StateMode s_mode,
                 const amrex::MLMGBndryT<amrex::MultiFab> *bndry = nullptr,
                 bool skip_fillboundary = false) const override;

    // out = (I/dt + J) in. Reconstructs the complete input via mf_to_flat,
    // runs the bound mat-vec, and scatters the result back. The applyBC() that
    // MLCellLinOp runs before this only touches ghost cells, which are ignored.
    void Fapply(int amrlev, int mglev, amrex::MultiFab &out,
                const amrex::MultiFab &in) const override;

    // Never invoked: with coarsening disabled the single-level V-cycle calls
    // the bottom solver directly and does no relaxation sweeps.
    void Fsmooth(int amrlev, int mglev, amrex::MultiFab &sol,
                 const amrex::MultiFab &rhs, int redblack) const override;

    // No fluxes: single level, no reflux between AMR levels.
    void FFlux(int amrlev, const amrex::MFIter &mfi,
               const amrex::Array<amrex::FArrayBox *, AMREX_SPACEDIM> &flux,
               const amrex::FArrayBox &sol, Location loc,
               int face_only = 0) const override;

    // The Jacobian is dense within each row (Hessian stencil plus antipodal
    // coupling), so it is not a compact cross stencil.
    [[nodiscard]] bool isCrossStencil() const override { return false; }

    // I/dt makes the operator strictly non-singular.
    [[nodiscard]] bool isSingular(int /*amrlev*/) const override
    {
        return false;
    }
    [[nodiscard]] bool isBottomSingular() const override { return false; }

  private:
    int m_n_rings;
    int m_ring_size;
    int m_num_particles;

    amrex::Geometry m_geom;
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;

    MatVec m_matvec;
};

#endif /* AHJACOBIANOP_HPP_ */
