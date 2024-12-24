/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#define PUNCTURETRACKER_HPP_

#include <AMReX_Particles.H>
#include <AMReX_RealVect.H>

#include "GRAMR.hpp"

//!  The class tracks the puncture locations by advecting them in the reverse
//!  direction to the shift. It is an amrex AoS ParticleContainer.
class PunctureTracker : public amrex::ParticleContainer<AMREX_SPACEDIM, 0>
{
  private:
    //! Params for puncture tracking
    int m_num_punctures{0};
    amrex::Vector<amrex::Real>
        m_puncture_coords; //!< the puncture location broadcast to all ranks

    std::string m_punctures_filename;
    std::string m_checkpoint_subdir;

    GRAMR *m_gr_amr{nullptr};

    bool m_initialized{false};
    bool m_started{false};

  public:
    //! The constructor
    using amrex::ParticleContainer<AMREX_SPACEDIM, 0>::ParticleContainer;

    //! Initialize the tracker. Note that this does not set up the underlying
    //! ParticleContainer
    void initialize(const amrex::Vector<amrex::Real> &initial_puncture_coords,
                    GRAMR *a_gr_amr,
                    const std::string &a_filename    = "punctures",
                    const std::string &a_output_path = "./");

    //! start the puncture tracker from the initial punctures
    void start_from_initial_punctures();

    //! restart the puncture tracker
    void restart(const std::string &a_restart_chk_dir);

    //! write punctures to the checkpoint directory
    void checkpoint(const std::string &a_chk_dir);

    //! Track the punctures and write out if requested
    void track(double a_time, double a_restart_time, double a_dt,
               const bool write_punctures = true);

#ifndef AMREX_USE_CUDA
  private: // CUDA doesn't allow lambdas in private functions
#endif

    //! set the initial punctures in the particle container
    void set_initial_punctures_pc();

  private:
    //! write the initial punctures to a file
    void write_initial_punctures() const;
};

#endif /* PUNCTURETRACKER_HPP_ */
