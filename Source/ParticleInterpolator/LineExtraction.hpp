#ifndef LINEEXTRACTION_HPP_
#define LINEEXTRACTION_HPP_

#include "InterpolationQueryParticle.hpp"
#include "ParticleInterpolator.hpp"
#include "SmallDataIO.hpp"
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

//! extract values of a (possibly multiple-component) variable along a line.
//! Line coords can be general and are provided by the user via start and end
//! coords (see code below on how the points of the line are defined).
template <int num_components> class LineExtraction
{
  private:
    const int m_start_comp; // first component
    const int m_num_points; // number of points along the line
    const std::array<double, AMREX_SPACEDIM>
        m_start_coords; // starting coords of the line
    const std::array<double, AMREX_SPACEDIM>
        m_end_coords; // ending coords of the line
    // for data write-out
    const double m_dt;
    const double m_time;
    const double m_restart_time;
    const bool m_first_step;

  public:
    LineExtraction(int a_start_comp, int a_num_points,
                   std::array<double, AMREX_SPACEDIM> a_start_coords,
                   std::array<double, AMREX_SPACEDIM> a_end_coords, double a_dt,
                   double a_time, bool a_restart_time, bool a_first_step)
        : m_start_comp(a_start_comp),
          m_num_points(
              (amrex::ParallelDescriptor::IOProcessor() ? a_num_points : 0)),
          m_start_coords(a_start_coords), m_end_coords(a_end_coords),
          m_dt(a_dt), m_time(a_time), m_restart_time(a_restart_time),
          m_first_step(a_first_step)
    {
    }

    ~LineExtraction() = default;

    //! Execute using particle-based interpolation
    void execute_query(ParticleInterpolator<num_components> *interpolator,
                       std::string a_file_prefix) const
    {
        BL_PROFILE("LineExtraction::execute_query()");

        // build coordinates along x from start to end
        std::vector<double> interp_x(m_num_points);
        std::vector<double> interp_y(m_num_points);
        std::vector<double> interp_z(m_num_points);

        for (int i = 0; i < m_num_points; ++i)
        {
            interp_x[i] =
                m_start_coords[0] + (double(i) / double(m_num_points - 1)) *
                                        (m_end_coords[0] - m_start_coords[0]);
            interp_y[i] =
                m_start_coords[1] + (double(i) / double(m_num_points - 1)) *
                                        (m_end_coords[1] - m_start_coords[1]);
            interp_z[i] =
                m_start_coords[2] + (double(i) / double(m_num_points - 1)) *
                                        (m_end_coords[2] - m_start_coords[2]);
        }

        // set up InterpolationQuery
        std::vector<double> interp_var_data(m_num_points * num_components, 0.0);

        InterpolationQueryParticle query(m_num_points);
        query.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data());

        // register components individually
        for (int k = 0; k < num_components; ++k)
        {
            // each component writes into its own stride
            double *out_k = interp_var_data.data() + k * m_num_points;
            query.addComp(
                m_start_comp + k, out_k,
                VariableType::state); // here state or derived? we will write
                                      // stuff into out_k array
        }

        interpolator->interp(query);

        const bool first_step =
            (std::abs(m_time) == m_dt); // random for now, needs a fix

        SmallDataIO output_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::NEW, first_step);

        std::vector<std::string> data_header_line(num_components);
        for (int icomp = 0; icomp < num_components; ++icomp)
        {
            data_header_line[icomp] =
                StateVariables::names[m_start_comp + icomp];
        }
        output_file.write_header_line(data_header_line, {"x", "y", "z"});

        for (int ipoint = 0; ipoint < m_num_points; ++ipoint)
        {
            std::vector<double> data_line(num_components);
            for (int icomp = 0; icomp < num_components; ++icomp)
            {
                data_line[icomp] =
                    interp_var_data[icomp * m_num_points + ipoint];
            }
            output_file.write_data_line(
                data_line,
                {interp_x[ipoint], interp_y[ipoint], interp_z[ipoint]});
        }
    }
};

#endif /* LINEEXTRACTION_HPP_ */