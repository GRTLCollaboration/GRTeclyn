/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef PUNCTURETRACKER_HPP_
#define PUNCTURETRACKER_HPP_

#include <AMReX_Array.H>
#include <AMReX_Particles.H>

#include "GRAMR.hpp"

//!  The class tracks the puncture locations by advecting them in the reverse
//!  direction to the shift. It is an amrex AoS ParticleContainer.
template <unsigned int num_punctures>
class PunctureTracker : public amrex::ParticleContainer<AMREX_SPACEDIM, 0>
{
  public:
    static constexpr unsigned int num_puncture_coords =
        num_punctures * AMREX_SPACEDIM;

  private:
    amrex::Array<amrex::Real, num_puncture_coords> m_puncture_coords;

    std::string m_punctures_filename;
    std::string m_checkpoint_subdir;

    GRAMR *m_gr_amr{nullptr};

    bool m_initialized{false};
    bool m_started{false};

    double m_restart_time{0.0};

  public:
    //! The constructor
    using amrex::ParticleContainer<AMREX_SPACEDIM, 0>::ParticleContainer;

    //! Initialize the tracker. Note that this does not set up the underlying
    //! ParticleContainer
    void initialize(GRAMR *a_gr_amr);

    //! start the puncture tracker from the initial punctures
    void start_from_initial_punctures(
        const amrex::Array<amrex::Real, num_puncture_coords>
            &a_initial_puncture_coords);

    //! restart the puncture tracker
    void restart(const std::string &a_restart_chk_dir);

    //! write punctures to the checkpoint directory
    void checkpoint(const std::string &a_chk_dir);

    //! Track the punctures and write out if requested
    void track(double a_time, double a_dt, const bool a_write_punctures = true);

#ifndef AMREX_USE_CUDA
  private: // CUDA doesn't allow lambdas in private functions
#endif

    //! set the initial punctures in the particle container
    void set_initial_punctures_pc();

    //! update m_puncture_coords from the particle locations
    void update_puncture_coords();

    //! return the linear index of the coord in the idir direction for the
    //! ipuncture puncture in m_puncture_coords
    static AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
    linear_idx(int ipuncture, int idir)
    {
        return ipuncture * AMREX_SPACEDIM + idir;
    }

#ifdef AMREX_USE_CUDA
  private:
#endif
    //! write the initial punctures to a file
    void write_initial_punctures() const;

    //! SmallDataIO requires a std::vector to write the coords
    std::vector<amrex::Real> get_puncture_vector() const;
};

#include "PunctureTracker.impl.hpp"

#endif /* PUNCTURETRACKER_HPP_ */
