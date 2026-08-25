/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef LINEEXTRACTION_HPP_
#define LINEEXTRACTION_HPP_

#include "InterpolationQueryParticle.hpp"
#include "ParticleInterpolator.hpp"
#include "SmallDataIO.hpp"
#include <array>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

//! extract values of a (possibly multiple-component) variable along a line.
//! Line coords can be general and are provided by the user via start and end
//! coords (see code below on how the points of the line are defined).
template <int num_components> class LineExtraction
{
  public:
    struct derived_vars_t
    {
        std::string name;
        std::array<std::string, num_components> component_names;
        std::array<BCParity, num_components> parities;
    };

  private:
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    const int m_start_comp; // first component
    const int m_num_points; // number of points along the line
    const std::array<amrex::ParticleReal, AMREX_SPACEDIM>
        m_start_coords; // starting coords of the line
    const std::array<amrex::ParticleReal, AMREX_SPACEDIM>
        m_end_coords; // ending coords of the line
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    // for data write-out
    amrex::Real m_dt{};
    amrex::Real m_time{};
    amrex::Real m_restart_time{};
    bool m_first_step{};

    void
    interpolate_and_write(ParticleInterpolator<num_components> *interpolator,
                          const std::string &a_file_prefix,
                          VariableType a_variable_type,
                          const derived_vars_t &a_derived_vars) const
    {
        BL_PROFILE("LineExtraction::interpolate_and_write()");

        // build coordinates along x from start to end
        std::vector<amrex::ParticleReal> interp_x(m_num_points);
        std::vector<amrex::ParticleReal> interp_y(m_num_points);
        std::vector<amrex::ParticleReal> interp_z(m_num_points);

        for (int i = 0; i < m_num_points; ++i)
        {
            interp_x[i] =
                m_start_coords[0] + (amrex::ParticleReal(i) /
                                     amrex::ParticleReal(m_num_points - 1)) *
                                        (m_end_coords[0] - m_start_coords[0]);
            interp_y[i] =
                m_start_coords[1] + (amrex::ParticleReal(i) /
                                     amrex::ParticleReal(m_num_points - 1)) *
                                        (m_end_coords[1] - m_start_coords[1]);
            interp_z[i] =
                m_start_coords[2] + (amrex::ParticleReal(i) /
                                     amrex::ParticleReal(m_num_points - 1)) *
                                        (m_end_coords[2] - m_start_coords[2]);
        }

        std::vector<amrex::ParticleReal> interp_var_data(
            m_num_points * num_components, 0.0);

        InterpolationQueryParticle query(m_num_points);
        query.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data());

        for (int component = 0; component < num_components; ++component)
        {
            amrex::ParticleReal *out =
                interp_var_data.data() + component * m_num_points;
            query.addComp(m_start_comp + component, out, a_variable_type,
                          a_derived_vars.parities[component]);
        }

        interpolator->interp(query, false, a_derived_vars.name, m_time);

        SmallDataIO output_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW, m_first_step);

        const amrex::Vector<std::string> data_header_line(
            a_derived_vars.component_names.begin(),
            a_derived_vars.component_names.end());
        output_file.write_header_line(data_header_line, {"x", "y", "z"});

        for (int point = 0; point < m_num_points; ++point)
        {
            std::vector<amrex::Real> data_line(num_components);
            for (int component = 0; component < num_components; ++component)
            {
                data_line[component] = static_cast<amrex::Real>(
                    interp_var_data[component * m_num_points + point]);
            }
            const std::vector<amrex::Real> coordinates{
                static_cast<amrex::Real>(interp_x[point]),
                static_cast<amrex::Real>(interp_y[point]),
                static_cast<amrex::Real>(interp_z[point])};
            output_file.write_data_line(data_line, coordinates);
        }
    }

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    LineExtraction(
        int a_start_comp, int a_num_points,
        std::array<amrex::ParticleReal, AMREX_SPACEDIM> a_start_coords,
        std::array<amrex::ParticleReal, AMREX_SPACEDIM> a_end_coords,
        amrex::Real a_dt, amrex::Real a_time, amrex::Real a_restart_time,
        bool a_first_step)
        : m_start_comp(a_start_comp),
          m_num_points(
              (amrex::ParallelDescriptor::IOProcessor() ? a_num_points : 0)),
          m_start_coords(a_start_coords), m_end_coords(a_end_coords),
          m_dt(a_dt), m_time(a_time), m_restart_time(a_restart_time),
          m_first_step(a_first_step)
    {
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    ~LineExtraction() = default;

    //! Execute using particle-based interpolation
    void execute_query(ParticleInterpolator<num_components> *interpolator,
                       const std::string &a_file_prefix) const
    {
        derived_vars_t state_vars;
        state_vars.parities.fill(BCParity::undefined);
        for (int component = 0; component < num_components; ++component)
        {
            state_vars.component_names[component] =
                StateVariables::names[m_start_comp + component];
        }
        interpolate_and_write(interpolator, a_file_prefix, VariableType::state,
                              state_vars);
    }

    //! Execute a query for components of an AMReX derived record.
    void execute_query(ParticleInterpolator<num_components> *interpolator,
                       const std::string &a_file_prefix,
                       const derived_vars_t &a_derived_vars) const
    {
        interpolate_and_write(interpolator, a_file_prefix,
                              VariableType::derived, a_derived_vars);
    }
};

#endif /* LINEEXTRACTION_HPP_ */
