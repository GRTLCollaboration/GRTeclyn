/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef LINEEXTRACTION_HPP_
#define LINEEXTRACTION_HPP_

#include "FilesystemTools.hpp"
#include "InterpolationQueryParticle.hpp"
#include "LineExtractionParameters.hpp"
#include "ParticleInterpolator.hpp"
#include "SmallDataIO.hpp"
#include <array>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

//! extract values of a (possibly multiple-component) variable along a line.
//! Line coords are defined in a named parameter scope by a starting point,
//! direction and maximum radius.
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
    line_extraction_params_t m_params;
    int m_start_comp;   // first component
    int m_num_points{}; // number of points along the line on this rank
    // for data write-out
    amrex::Real m_dt{};
    amrex::Real m_time{};
    amrex::Real m_restart_time{};
    bool m_first_step{};

    void
    interpolate_and_write(ParticleInterpolator<num_components> *interpolator,
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
                m_params.start_coords[0] +
                (amrex::ParticleReal(i) /
                 amrex::ParticleReal(m_num_points - 1)) *
                    (m_params.end_coords[0] - m_params.start_coords[0]);
            interp_y[i] =
                m_params.start_coords[1] +
                (amrex::ParticleReal(i) /
                 amrex::ParticleReal(m_num_points - 1)) *
                    (m_params.end_coords[1] - m_params.start_coords[1]);
            interp_z[i] =
                m_params.start_coords[2] +
                (amrex::ParticleReal(i) /
                 amrex::ParticleReal(m_num_points - 1)) *
                    (m_params.end_coords[2] - m_params.start_coords[2]);
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

        SmallDataIO output_file(m_params.data_path + m_params.file_prefix, m_dt,
                                m_time, m_restart_time, SmallDataIO::NEW,
                                m_first_step);

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
    LineExtraction(const std::string &a_param_scope, int a_start_comp,
                   amrex::Real a_dt, amrex::Real a_time,
                   amrex::Real a_restart_time, bool a_first_step)
        : m_params(a_param_scope), m_start_comp(a_start_comp), m_dt(a_dt),
          m_time(a_time), m_restart_time(a_restart_time),
          m_first_step(a_first_step)
    {
        m_params.fill_params();
        if (m_params.enabled)
        {
            m_num_points = amrex::ParallelDescriptor::IOProcessor()
                               ? m_params.num_points
                               : 0;
            if (m_time < m_restart_time + 1.5 * m_dt || m_first_step)
            {
                FilesystemTools::ensure_directory_exists(m_params.data_path);
            }
        }
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    ~LineExtraction() = default;

    //! Execute using particle-based interpolation
    void execute_query(ParticleInterpolator<num_components> *interpolator) const
    {
        if (!m_params.enabled)
        {
            return;
        }
        derived_vars_t state_vars;
        state_vars.parities.fill(BCParity::undefined);
        for (int component = 0; component < num_components; ++component)
        {
            state_vars.component_names[component] =
                StateVariables::names[m_start_comp + component];
        }
        interpolate_and_write(interpolator, VariableType::state, state_vars);
    }

    //! Execute a query for components of an AMReX derived record.
    void execute_query(ParticleInterpolator<num_components> *interpolator,
                       const derived_vars_t &a_derived_vars) const
    {
        if (!m_params.enabled)
        {
            return;
        }
        interpolate_and_write(interpolator, VariableType::derived,
                              a_derived_vars);
    }
};

#endif /* LINEEXTRACTION_HPP_ */
