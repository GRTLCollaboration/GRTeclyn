/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "PunctureTracker.hpp"
#include "AMReXParameters.hpp" // for writing data
#include "DimensionDefinitions.hpp"
#include "InterpolationQuery.hpp"
#include "SmallDataIO.hpp"   // for writing data
#include "UserVariables.hpp" // for writing data

//! Set punctures post restart
void PunctureTracker::initial_setup(
    const std::vector<std::array<double, AMREX_SPACEDIM>>
        &initial_puncture_coords,
    const std::string &a_filename, const std::string &a_output_path,
    const int a_min_level)
{
    if (!FilesystemTools::directory_exists(a_output_path))
    {
        FilesystemTools::mkdir_recursive(a_output_path);
    }

    m_punctures_filename = a_output_path + a_filename;

    // first set the puncture data
    // m_num_punctures is only set later
    m_puncture_coords = initial_puncture_coords;

    m_min_level = a_min_level;
}

void PunctureTracker::restart_punctures()
{
#if 0
//xxxxx
    int current_step = m_interpolator->getAMR().s_step;

    if (current_step == 0)
    {
        // if it is the first timestep, use the param values
        // rather than look for the output file, e.g. for when
        // restart from IC solver checkpoint
        set_initial_punctures();
    }
    else
    {
        // look for the current puncture location in the
        // puncture output file (it needs to exist!)
//xxxxx        read_in_punctures(current_step,
//                          m_interpolator->getAMR().getCurrentTime());
    }
#endif
}

//! set and write initial puncture locations
void PunctureTracker::set_initial_punctures()
{
    AMREX_ASSERT(m_puncture_coords.size() > 0); // sanity check

    m_num_punctures = m_puncture_coords.size();
    m_puncture_shift.resize(m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        // assume initial shift is always zero
        FOR (i)
        {
            m_puncture_shift[ipuncture][i] = 0.0;
        }
    }

    // now the write out to a new file
    bool first_step     = true;
    double dt           = 1.; // doesn't matter
    double time         = 0.;
    double restart_time = 0.;
    SmallDataIO punctures_file(m_punctures_filename, dt, time, restart_time,
                               SmallDataIO::APPEND, first_step);
    std::vector<std::string> header1_strings(AMREX_SPACEDIM * m_num_punctures);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        std::string idx = std::to_string(ipuncture + 1);
        header1_strings[AMREX_SPACEDIM * ipuncture + 0] = "x_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 1] = "y_" + idx;
        header1_strings[AMREX_SPACEDIM * ipuncture + 2] = "z_" + idx;
    }
    punctures_file.write_header_line(header1_strings);

    // use a vector for the write out
    punctures_file.write_time_data_line(get_puncture_vector());
}

//! Set punctures post restart
void PunctureTracker::read_in_punctures(int a_int_step, double a_current_time)
{
    // read them in from the Punctures file at current time m_time
    // NB opening in APPEND mode allows reading where m_restart_time
    // is greater than zero and m_time < m_restart_time + m_dt
    bool first_step = false;
    double dt       = (a_current_time / a_int_step);
    SmallDataIO punctures_file(m_punctures_filename, dt, a_current_time,
                               a_current_time, SmallDataIO::APPEND, first_step);

    // NB need to give the get function an empty vector to fill
    std::vector<double> puncture_vector;
    punctures_file.get_specific_data_line(puncture_vector, a_current_time);

    // check the data returned is the right size
    AMREX_ASSERT(puncture_vector.size() % AMREX_SPACEDIM == 0);

    m_num_punctures = puncture_vector.size() / AMREX_SPACEDIM;
    m_puncture_coords.resize(m_num_punctures);

    // remove any duplicate data from the file
    const bool keep_m_time_data = true;
    punctures_file.remove_duplicate_time_data(keep_m_time_data);

    // convert vector to list of coords
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        m_puncture_coords[ipuncture] = {
            puncture_vector[ipuncture * AMREX_SPACEDIM + 0],
            puncture_vector[ipuncture * AMREX_SPACEDIM + 1],
            puncture_vector[ipuncture * AMREX_SPACEDIM + 2]};
    }

    // set the coordinates and get the current shift
    interp_shift();

    // print out values into pout files
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        amrex::Print() << "Puncture " << ipuncture
                       << " restarted at : " << m_puncture_coords[ipuncture][0]
                       << " " << m_puncture_coords[ipuncture][1] << " "
                       << m_puncture_coords[ipuncture][2] << std::endl;
        amrex::Print() << " with shift vector : "
                       << m_puncture_shift[ipuncture][0] << " "
                       << m_puncture_shift[ipuncture][1] << " "
                       << m_puncture_shift[ipuncture][2] << std::endl;
        amrex::Print() << "at time = " << a_current_time << std::endl;
    }
}

//! Execute the tracking and write out
void PunctureTracker::execute_tracking(double a_time, double a_restart_time,
                                       double a_dt, const bool write_punctures)
{
    BL_PROFILE("PunctureTracker::execute_tracking");
    // leave if this is called at t=0, we don't want to move the puncture yet
    if (m_num_punctures == 0 || a_time == 0.)
    {
        return;
    }
    AMREX_ASSERT(m_interpolator != nullptr); // sanity check

    // get puncture coordinates and old shift value
    std::vector<std::array<double, AMREX_SPACEDIM>> old_shift =
        m_puncture_shift;
    AMREX_ASSERT(static_cast<int>(m_puncture_coords.size()) ==
                 m_num_punctures); // sanity check

    // new shift value
    interp_shift();

    // update puncture locations using second order update
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        FOR (i)
        {
            m_puncture_coords[ipuncture][i] +=
                -0.5 * a_dt *
                (m_puncture_shift[ipuncture][i] + old_shift[ipuncture][i]);
        }
    }

    // print them out
    if (write_punctures)
    {
        bool first_step = false;
        SmallDataIO punctures_file(m_punctures_filename, a_dt, a_time,
                                   a_restart_time, SmallDataIO::APPEND,
                                   first_step);

        // use a vector for the write out
        punctures_file.write_time_data_line(get_puncture_vector());
    }
}

//! Use the interpolator to get the value of the shift at
//! given coords
void PunctureTracker::interp_shift()
{
#if 0
//xxxxx
    BL_PROFILE("PunctureTracker::interp_shift");
    // resize the vector to the number of punctures
    m_puncture_shift.resize(m_num_punctures);

    // refresh interpolator
    bool fill_ghosts = false;
    m_interpolator->refresh(fill_ghosts);
    // only fill the ghosts we need
    m_interpolator->fill_multilevel_ghosts(
        VariableType::evolution, Interval(c_shift1, c_shift3), m_min_level);

    // set up shift and coordinate holders
    std::vector<double> interp_shift1(m_num_punctures);
    std::vector<double> interp_shift2(m_num_punctures);
    std::vector<double> interp_shift3(m_num_punctures);
    std::vector<double> interp_x(m_num_punctures);
    std::vector<double> interp_y(m_num_punctures);
    std::vector<double> interp_z(m_num_punctures);

    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        interp_x[ipuncture] = m_puncture_coords[ipuncture][0];
        interp_y[ipuncture] = m_puncture_coords[ipuncture][1];
        interp_z[ipuncture] = m_puncture_coords[ipuncture][2];
    }

    // setup query
    InterpolationQuery query(m_num_punctures);
    query.setCoords(0, interp_x.data())
        .setCoords(1, interp_y.data())
        .setCoords(2, interp_z.data())
        .addComp(c_shift1, interp_shift1.data())
        .addComp(c_shift2, interp_shift2.data())
        .addComp(c_shift3, interp_shift3.data());

    // engage!
    m_interpolator->interp(query);

    // put the shift values into the output array
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        m_puncture_shift[ipuncture] = {interp_shift1[ipuncture],
                                       interp_shift2[ipuncture],
                                       interp_shift3[ipuncture]};
    }
#endif
}

//! get a vector of the puncture coords - used for write out
std::vector<double> PunctureTracker::get_puncture_vector() const
{
    std::vector<double> puncture_vector;
    puncture_vector.resize(m_num_punctures * AMREX_SPACEDIM);
    for (int ipuncture = 0; ipuncture < m_num_punctures; ipuncture++)
    {
        puncture_vector[ipuncture * AMREX_SPACEDIM + 0] =
            m_puncture_coords[ipuncture][0];
        puncture_vector[ipuncture * AMREX_SPACEDIM + 1] =
            m_puncture_coords[ipuncture][1];
        puncture_vector[ipuncture * AMREX_SPACEDIM + 2] =
            m_puncture_coords[ipuncture][2];
    }
    return puncture_vector;
}
