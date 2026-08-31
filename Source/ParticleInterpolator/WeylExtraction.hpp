/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef WEYLEXTRACTION_HPP_
#define WEYLEXTRACTION_HPP_

#include "SphericalExtraction.hpp"

/*!
   The class allows the user to extract data from the grid for the Weyl
   components over spherical shells at specified radii. The values may then be
   written to an output file, or integrated across the surfaces.
*/
class WeylExtraction : public SphericalExtraction<2>
{
  public:

    //! The constructor
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    WeylExtraction(const spherical_extraction_params_t &a_params,
                   amrex::Real a_dt, amrex::Real a_time, bool a_first_step,
                   amrex::Real a_restart_time = 0.0)
        : SphericalExtraction<2>(a_params, a_dt, a_time, a_first_step,
                                 a_restart_time)
    {
        amrex::Vector<BCParity> parities = {BCParity::even, BCParity::odd_xyz};
        this->add_derived_vars({0, 1}, parities, "Weyl4");
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    amrex::Vector<BCParity> parities = {BCParity::even, BCParity::odd_xyz};

    //! The old constructor which assumes it is called in specific_post_timestep
    //! so the first time step is when m_time == m_dt
    WeylExtraction(const spherical_extraction_params_t &a_params,
                   amrex::Real a_dt, amrex::Real a_time,
                   amrex::Real a_restart_time = 0.0)
        : WeylExtraction(a_params, a_dt, a_time, (a_dt == a_time),
                         a_restart_time)
    {
    }

    //! Execute the query
    void execute_query(ParticleInterpolator<2> *a_interpolator)
    {
        // extract the values of the Weyl scalars on the spheres
        this->extract(a_interpolator);

        if (this->m_params.write_extraction)
        {
            this->write_extraction(this->m_params.extraction_file_prefix);
        }

        // now calculate and write the requested spherical harmonic modes
        std::vector<std::pair<std::vector<amrex::ParticleReal>,
                              std::vector<amrex::ParticleReal>>>
            mode_integrals(m_num_modes);

        // note that this is normalised by multiplying by radius
        // NOLINTBEGIN(bugprone-easily-swappable-parameters)
        auto normalised_Weyl4_complex =
            [](std::vector<amrex::ParticleReal> Weyl4_reim_parts,
               amrex::ParticleReal r, amrex::ParticleReal, amrex::ParticleReal)
        {
            // here the std::vector<amrex::ParticleReal> passed will just have
            // the real and imaginary parts of the Weyl4 scalar as its
            // only components
            return std::make_pair(r * Weyl4_reim_parts[0],
                                  r * Weyl4_reim_parts[1]);
        };
        // NOLINTEND(bugprone-easily-swappable-parameters)

        // add the modes that will be integrated
        for (int imode = 0; imode < this->m_num_modes; ++imode)
        {
            const auto &mode                  = this->m_modes[imode];
            constexpr int spin_quantum_number = -2;
            this->add_mode_integrand(spin_quantum_number, mode.first,
                                     mode.second, normalised_Weyl4_complex,
                                     mode_integrals[imode]);
        }

        // do the integration over the surface
        this->integrate();

        // write the integrals
        for (int imode = 0; imode < this->m_num_modes; ++imode)
        {
            const auto &mode = this->m_modes[imode];
            std::string integrals_filename =
                this->m_params.integral_file_prefix +
                std::to_string(mode.first) + std::to_string(mode.second);
            std::vector<std::vector<amrex::ParticleReal>>
                integrals_for_writing = {
                    std::move(mode_integrals[imode].first),
                    std::move(mode_integrals[imode].second)};
            std::vector<std::string> labels = {"integral Re", "integral Im"};
            this->write_integrals(integrals_filename, integrals_for_writing,
                                  labels);
        }
    }
};

#endif /* WEYLEXTRACTION_HPP_ */