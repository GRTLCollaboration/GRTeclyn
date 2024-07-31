/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#define PUNCTURETRACKER_HPP_

#include <AMReX_Particles.H>
#include <AMReX_RealVect.H>

#include "GRAMR.hpp"

//!  The class tracks the puncture locations by integrating the shift at
//!  The puncture position
class PunctureTracker
{
  private:
    //! Params for puncture tracking
    int m_num_punctures{0};
    amrex::Vector<amrex::RealVect> m_puncture_coords;
    int m_update_level{}; //!< the level on which to update positions

    std::string m_punctures_filename;
    std::string m_checkpoint_subdir;

    // We will use the AMREX_SPACEDIM real attributes to store the
    // midpoint position/shift
    static std::unique_ptr<amrex::ParticleContainerPureSoA<AMREX_SPACEDIM, 0>>
        s_particle_container;

    using PunctureIter         = s_particle_container::ParIterType;
    using PunctureParticleType = s_particle_container::ParticleType;

    GRAMR *m_gr_amr{nullptr};

  public:
    //! The constructor
    PunctureTracker() = default;

    //! set puncture locations on start (or restart)
    //! this needs to be done before 'setupAMRObject'
    //! if the puncture locations are required for Tagging Criteria
    void
    initial_setup(const amrex::Vector<amrex::RealVect> &initial_puncture_coords,
                  GRAMR *a_gr_amr, const std::string &a_filename = "punctures",
                  const std::string &a_output_path = "./",
                  const int a_update_level         = 0);

    //! set puncture locations on start (or restart)
    void restart(int a_coarse_step);

    //! write punctures to the checkpoint directory
    void checkpoint(const std::string &a_chk_dir);

    //! Execute the tracking and write out
    void execute_tracking(double a_time, double a_restart_time, double a_dt,
                          const bool write_punctures = true);

    // function to get punctures
    [[nodiscard]] ALWAYS_INLINE const auto &get_puncture_coords() const
    {
        return m_puncture_coords;
    }

  private:
    //! set and write initial puncture locations
    void set_initial_punctures();

    //! Use the interpolator to get the value of the shift at
    //! given coords
    void interp_shift();

    //! Get a vector of the puncture coords - used for write out
    [[nodiscard]] std::vector<double> get_puncture_vector() const;
};

#endif /* PUNCTURETRACKER_HPP_ */
