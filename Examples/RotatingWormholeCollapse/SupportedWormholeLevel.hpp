#ifndef SUPPORTEDWORMHOLELEVEL_HPP_
#define SUPPORTEDWORMHOLELEVEL_HPP_

#include "BHAMR.hpp"
#include "DefaultLevelFactory.hpp"
#include "GRAMRLevel.hpp"

class SupportedWormholeLevel : public GRAMRLevel
{
  public:
    // The wormhole has no punctures, but we base the AMR container on BHAMR
    // (a GRAMR child) purely to reuse its set-up ParticleInterpolator for
    // in-code Weyl4 / Psi4 spherical-harmonic extraction.  Puncture tracking
    // stays disabled (puncture_tracking.enabled = 0), so the tracker member is
    // constructed but never used.
    static constexpr int num_punctures = 1;

    static void variableSetUp();

    using GRAMRLevel::GRAMRLevel;

    //! Access the owning BHAMR (for its m_weyl_interpolator).
    BHAMR<num_punctures> *get_bhamr_ptr();

    void specificAdvance() override;
    void initData() override;
    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         const double a_time) override;
    void specificUpdateODE(amrex::MultiFab &a_soln) override;
    void pre_tag_cells() final;
    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
    void specificPostTimeStep() override;
};

#endif /* SUPPORTEDWORMHOLELEVEL_HPP_ */