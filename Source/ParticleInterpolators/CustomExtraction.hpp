#ifndef CUSTOMEXTRACTION_HPP_
#define CUSTOMEXTRACTION_HPP_

#include "InterpolationQueryParticle.hpp"
#include "ParticleInterpolators.hpp"
#include "SmallDataIO.hpp"
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

// function to convert double to string with precision
inline std::string to_string_with_precision(double value, int precision)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << value;
    return out.str();
}

//! extract values of a (possibly multiple-component) variable along a line
//! (x-axis)
class CustomExtraction
{
  private:
    const int m_start_comp; // first component
    const int m_ncomp;      // number of components (<= AMREX_SPACEDIM)
    const int m_num_points; // number of points along the line
    const double m_L;       // length of the line
    const std::array<double, AMREX_SPACEDIM> m_origin; // line starts here
    // for data write-out
    const double m_dt;
    const double m_time;
    const double m_restart_time;
    const bool m_first_step;

  public:
    CustomExtraction(int a_start_comp, int a_ncomp, int a_num_points,
                     double a_L, std::array<double, AMREX_SPACEDIM> a_origin,
                     double a_dt, double a_time, bool a_restart_time,
                     bool a_first_step)
        : m_start_comp(a_start_comp), m_ncomp(a_ncomp),
          m_num_points(a_num_points), m_L(a_L), m_origin(a_origin), m_dt(a_dt),
          m_time(a_time), m_restart_time(a_restart_time),
          m_first_step(a_first_step)
    {
    }

    ~CustomExtraction() = default;

    //! Execute using particle-based interpolation
    void execute_query(ParticleInterpolators &particles,
                       std::string a_file_prefix) const
    {

        // build coordinates along x from origin to origin + L
        std::vector<double> interp_x(m_num_points);
        std::vector<double> interp_y(m_num_points, m_origin[1]);
        std::vector<double> interp_z(m_num_points, m_origin[2]);

        for (int i = 0; i < m_num_points; ++i)
        {
            interp_x[i] =
                m_origin[0] + (double(i) / double(m_num_points - 1)) * m_L;
        }

        // set up InterpolationQuery
        std::vector<double> interp_var_data(m_num_points * m_ncomp, 0.0);

        InterpolationQueryParticle query(m_num_points);
        query.setCoords(0, interp_x.data())
            .setCoords(1, interp_y.data())
            .setCoords(2, interp_z.data());

        // register components individually
        for (int k = 0; k < m_ncomp; ++k)
        {
            // each component writes into its own stride
            double *out_k = interp_var_data.data() + k * m_num_points;
            query.addComp(
                m_start_comp + k, out_k, Derivative::LOCAL,
                VariableType::state); // here state or derived? we will write
                                      // stuff into out_k array
        }

        // populate particles at the query points
        particles.populate_from_query(query);

        // interpolate from the AMR state into particle SOA
        particles.interpolate_to_particle();

        // fill the InterpolationQuery outputs from particle SOA with parity
        // rules applied
        particles.interp(query);

        const bool first_step =
            (std::abs(m_time) == m_dt); // random for now, needs a fix

        SmallDataIO output_file(a_file_prefix, m_dt, m_time, m_restart_time,
                                SmallDataIO::APPEND, first_step);

        if (first_step)
        {
            // Header: "x = <coord>" for each point
            std::vector<std::string> header_line(m_num_points);
            for (int i = 0; i < m_num_points; ++i)
            {
                header_line[i] =
                    "x = " + to_string_with_precision(interp_x[i], 2);
            }
            output_file.write_header_line(header_line);
        }

        if (m_ncomp == 1)
        {
            output_file.write_time_data_line(interp_var_data);
        }
        else
        {
            for (int k = 0; k < m_ncomp; ++k)
            {
                const double *row = interp_var_data.data() + k * m_num_points;
                output_file.write_time_data_line(
                    std::vector<double>(row, row + m_num_points));
            }
        }
    }
};

#endif /* CUSTOMEXTRACTION_HPP_ */
