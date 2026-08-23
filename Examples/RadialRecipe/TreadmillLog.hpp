/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TREADMILLLOG_HPP_
#define TREADMILLLOG_HPP_

#include "SmallDataIO.hpp"

#include <string>
#include <vector>

/// The column contract of `data/treadmill.dat`, and nothing else.
///
/// SINGLE RESPONSIBILITY.  Kept apart from the policy that decides when to
/// shift so that changing what is recorded cannot change what is done, and
/// vice versa.
///
/// This is a NEW file with its OWN columns.  No existing stream is widened:
/// `constraint_norms.dat`, `sector_dynamics.dat` and the rest are positional
/// contracts shared with live campaigns and with every archived cell, and a
/// column inserted into one of them silently reinterprets every reader.
///
/// It deliberately duplicates position information that the plotfile consumer
/// also derives.  That is the point: the consumer's numbers arrive at plotfile
/// cadence (hundreds of samples on a long run) and these arrive at the tracker's
/// cadence (tens of thousands).  The two must agree where they overlap, and
/// checking that they do is a validation gate rather than an assumption.
class TreadmillLog
{
  public:
    explicit TreadmillLog(std::string a_file_prefix)
        : m_file_prefix(std::move(a_file_prefix))
    {
    }

    /// One row.  `a_cells_shifted` is zero except on the steps where a shift
    /// actually happened, so a reader can find every seam by filtering it.
    /// `core0` / `core1` are the two sectors in the order they were configured,
    /// along the shift axis; which of them leads is a property of the run, not
    /// of this file, so the columns do not pretend to know.
    ///
    /// `a_row_spacing` is the time between ROWS OF THIS FILE, not the
    /// simulation timestep.  `SmallDataIO` uses it for exactly one thing --
    /// deciding whether this write is the first one after a restart, via
    /// `m_time < m_restart_time + m_dt` -- so a stream that writes every N
    /// steps must pass `N * dt`.  Passing the raw timestep makes that test
    /// false on every write, `remove_duplicate_time_data` never fires, and a
    /// restart leaves the rows it re-computes duplicated in the file.
    void write_row(double a_time, double a_row_spacing, double a_restart_time,
                   bool a_first_step, int a_step, double a_core0,
                   double a_core1, double a_midpoint_grid,
                   int a_cells_shifted, long a_odometer_cells,
                   double a_odometer_length) const
    {
        SmallDataIO file(m_file_prefix, a_row_spacing, a_time, a_restart_time,
                         SmallDataIO::APPEND, a_first_step);
        file.remove_duplicate_time_data();
        if (a_first_step)
        {
            file.write_header_line({"step", "core0_grid", "core1_grid",
                                    "midpoint_grid", "cells_shifted",
                                    "odometer_cells", "odometer_length",
                                    "midpoint_true"});
        }
        file.write_time_data_line({static_cast<double>(a_step), a_core0,
                                   a_core1, a_midpoint_grid,
                                   static_cast<double>(a_cells_shifted),
                                   static_cast<double>(a_odometer_cells),
                                   a_odometer_length,
                                   a_midpoint_grid + a_odometer_length});
    }

  private:
    std::string m_file_prefix;
};

#endif /* TREADMILLLOG_HPP_ */
