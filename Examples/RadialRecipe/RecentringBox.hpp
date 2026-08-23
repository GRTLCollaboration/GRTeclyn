/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef RECENTRINGBOX_HPP_
#define RECENTRINGBOX_HPP_

#include "GridTreadmill.hpp"
#include "SectorCoreTracker.hpp"
#include "TreadmillLog.hpp"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

/// The policy layer of the recentring box: WHEN to shift, and the bookkeeping
/// that makes the shift readable afterwards.
///
/// SINGLE RESPONSIBILITY.  It owns the decision and the odometer.  It does not
/// know how to translate a grid (`GridTreadmill.hpp`), how to find a star
/// (`SectorCoreTracker.hpp`) or what a column of output looks like
/// (`TreadmillLog.hpp`) -- it composes those three.
///
/// DEPENDENCY INVERSION.  It is handed the state MultiFabs, the geometry, the
/// sector field indices and the asymptotic values; it never includes the level
/// or the simulation parameters, so it cannot quietly grow a dependency on
/// either and can be exercised without an Amr hierarchy.
///
/// OPEN/CLOSED.  A different notion of "where the source is" is a different
/// `SectorFieldSet` list; a different way of inventing the leading shell is a
/// different fill mode; neither requires touching this class.
///
/// The whole thing is inert unless `treadmill_enabled` is set, so every
/// archived cell is bit-for-bit unaffected by its presence.
class RecentringBox
{
  public:
    /// What happened on one visit.  Returned rather than logged directly so a
    /// caller (or a test) can react without parsing a file.
    struct Step
    {
        bool tracked{false};      //!< a position measurement was taken
        bool shifted{false};      //!< a shift was performed on this step
        int cells_shifted{0};     //!< 0 unless shifted
        double midpoint_grid{0.}; //!< source midpoint, grid coordinates
        double core0{0.};
        double core1{0.};
    };

    /// One-time setup.  `a_sectors` lists the matter sectors whose cores define
    /// the source position; their midpoint along `a_params.axis` is what the
    /// threshold is compared against.
    ///
    /// Returns false, having printed why, if the configuration cannot be
    /// honoured exactly.  A caller that ignores the return value and runs
    /// anyway is choosing a silently wrong trajectory, so every failure here is
    /// escalated to an abort by `require_valid`.
    bool configure(const GridTreadmillParams &a_params, double a_dx,
                   const std::array<double, AMREX_SPACEDIM> &a_center,
                   const std::vector<SectorFieldSet> &a_sectors,
                   const std::vector<double> &a_asymptotic_values,
                   const std::string &a_output_prefix, int a_max_level,
                   const std::array<int, AMREX_SPACEDIM> &a_is_periodic,
                   int a_min_box_extent_on_axis, double a_sponge_inner_radius,
                   double a_source_reach)
    {
        m_params             = a_params;
        m_sectors            = a_sectors;
        m_asymptotic_values  = a_asymptotic_values;
        m_center             = a_center;
        m_odometer.cells     = 0;
        m_odometer.dx        = a_dx;
        m_log                = std::make_unique<TreadmillLog>(a_output_prefix);
        m_seeded             = false;
        m_configured         = false;
        m_wrote_header       = false;

        if (!m_params.enabled)
        {
            return true; // inert: nothing to validate, nothing to do
        }

        std::ostringstream why;

        if (m_params.axis < 0 || m_params.axis >= AMREX_SPACEDIM)
        {
            why << "treadmill_axis must be 0, 1 or 2\n";
        }
        else if (a_is_periodic[m_params.axis] != 0)
        {
            // The whole seam argument assumes data leaving the trailing face is
            // gone.  On a periodic axis it comes back round the other side.
            why << "the shift axis is periodic; the recentring box requires a "
                   "non-periodic axis\n";
        }

        // A shift that does not land on the mesh cannot be taken exactly, and
        // rounding it silently would reintroduce the sub-cell aliasing this
        // design exists to avoid.
        m_shift_cells =
            GridTreadmill::cells_for_length(m_params.threshold, a_dx);
        if (m_shift_cells < 1)
        {
            why << "treadmill_threshold (" << m_params.threshold
                << ") is not a whole number of cells of size " << a_dx << "\n";
        }

        if (a_max_level > 0)
        {
            // Refuse rather than do something plausible on a refined grid.
            why << "the recentring box is single-level only, but max_level = "
                << a_max_level << "\n";
        }

        // The leading sliver is filled from the outermost surviving layer,
        // which is only in the same box as the sliver when the shift is shorter
        // than a box.
        if (m_shift_cells >= 1 &&
            a_min_box_extent_on_axis < m_shift_cells + 1)
        {
            why << "the shift is " << m_shift_cells
                << " cells but the smallest box on that axis is only "
                << a_min_box_extent_on_axis
                << " cells; raise min_box_size or lower the threshold\n";
        }

        // If the source can reach the sponge while sitting at its furthest
        // allowed excursion, the run is being quietly damped and no amount of
        // care elsewhere recovers it.
        if (a_sponge_inner_radius > 0.0 &&
            m_params.threshold + a_source_reach > a_sponge_inner_radius)
        {
            why << "at the threshold the source reaches radius "
                << (m_params.threshold + a_source_reach)
                << ", which is inside the sponge at "
                << a_sponge_inner_radius << "\n";
        }

        if (m_sectors.empty())
        {
            why << "no matter sectors were configured to track\n";
        }

        const std::string message = why.str();
        if (!message.empty())
        {
            amrex::Print() << "[treadmill] refusing to start:\n" << message;
            return false;
        }

        m_configured = true;
        amrex::Print() << "[treadmill] recentring box armed: axis "
                       << m_params.axis << ", shift " << m_shift_cells
                       << " cells (" << m_params.threshold
                       << "), checked every " << m_params.check_interval
                       << " steps, fill mode " << m_params.fill_mode << "\n";
        return true;
    }

    /// Abort with the reason if `configure` refused.  Kept separate so a test
    /// can inspect the refusal instead of dying on it.
    static void require_valid(bool a_ok)
    {
        if (!a_ok)
        {
            amrex::Abort("[treadmill] invalid recentring-box configuration");
        }
    }

    [[nodiscard]] bool active() const
    {
        return m_params.enabled && m_configured;
    }

    [[nodiscard]] const GridTreadmill::Odometer &odometer() const
    {
        return m_odometer;
    }

    /// Track, decide, shift, and record.
    ///
    /// `a_states` are every time level that must move together.  The old level
    /// matters as much as the new one: `AmrLevel::RK` starts each step from the
    /// old data, so leaving it behind mixes two frames.
    Step advance(const std::vector<amrex::MultiFab *> &a_states,
                 const amrex::Geometry &a_geom, int a_step, double a_time,
                 double a_dt, double a_restart_time)
    {
        Step result;
        if (!active() || a_states.empty() || a_states[0] == nullptr)
        {
            return result;
        }
        const bool due = !m_wrote_header || (m_params.check_interval <= 1) ||
                         (a_step % m_params.check_interval == 0);
        if (!due)
        {
            return result;
        }

        result = measure(*a_states[0], a_geom);
        if (!result.tracked)
        {
            return result;
        }

        const double excursion =
            result.midpoint_grid - m_center[m_params.axis];
        if (std::abs(excursion) >= m_params.threshold)
        {
            // Always shift by exactly the threshold, never by the measured
            // excursion: the threshold is known to be a whole number of cells
            // and the excursion is not.  Checking often enough that the
            // overshoot is small is what keeps the source centred.
            //
            // The sign follows the source, so a mirror run travelling the other
            // way is handled by the same code rather than by a second path.
            const int cells = (excursion > 0.0) ? m_shift_cells : -m_shift_cells;
            for (amrex::MultiFab *state : a_states)
            {
                if (state == nullptr)
                {
                    continue;
                }
                GridTreadmill::shift(*state, m_params.axis, cells, a_geom,
                                     m_center, m_asymptotic_values,
                                     m_params.fill_mode);
            }
            m_odometer.add(cells);
            result.shifted       = true;
            result.cells_shifted = cells;

            // The source is now back at the centre; the tracker's search balls
            // must follow it, or the next measurement looks in the old place.
            const double moved = cells * m_odometer.dx;
            for (auto &seed : m_seeds)
            {
                seed[m_params.axis] -= moved;
            }
            result.core0 -= moved;
            result.core1 -= moved;
            result.midpoint_grid -= moved;

            amrex::Print() << "[treadmill] step " << a_step << " t=" << a_time
                           << ": shifted " << cells << " cells, odometer "
                           << m_odometer.cells << " cells ("
                           << m_odometer.length() << ")\n";
        }

        if (m_log)
        {
            // The header goes on the first row this module ever writes, not on
            // the first step of the run: this stream starts when the tracker
            // does, and on a restart the file already exists and must be
            // appended to rather than truncated.
            // The spacing between ROWS, not the timestep -- see the note on
            // write_row.  Rows land every check_interval steps, and passing
            // anything smaller stops the post-restart de-duplication firing.
            const double row_spacing =
                a_dt * static_cast<double>(std::max(1, m_params.check_interval));
            m_log->write_row(a_time, row_spacing, a_restart_time,
                             !m_wrote_header, a_step, result.core0,
                             result.core1, result.midpoint_grid,
                             result.cells_shifted, m_odometer.cells,
                             m_odometer.length());
            m_wrote_header = true;
        }
        return result;
    }

    /// Persist the odometer beside a checkpoint.  Written by the I/O rank only.
    void save(const std::string &a_dir) const
    {
        if (!active() || !amrex::ParallelDescriptor::IOProcessor())
        {
            return;
        }
        std::ofstream os(state_file(a_dir));
        os << m_odometer.to_string();
    }

    /// Restore the odometer after a restart.
    ///
    /// Returns false if the file is missing or unreadable.  The caller must
    /// treat that as fatal: resuming with the odometer at zero produces a
    /// trajectory that is wrong and looks fine, which is the worst failure
    /// available here.
    bool load(const std::string &a_dir)
    {
        if (!active())
        {
            return true;
        }
        std::ifstream is(state_file(a_dir));
        if (!is)
        {
            return false;
        }
        std::string key;
        long cells = 0;
        double dx  = 0.0;
        bool got_cells = false;
        while (is >> key)
        {
            if (key == "treadmill_odometer_cells")
            {
                is >> cells;
                got_cells = true;
            }
            else if (key == "treadmill_dx")
            {
                is >> dx;
            }
        }
        if (!got_cells)
        {
            return false;
        }
        if (dx > 0.0 && std::abs(dx - m_odometer.dx) > 1.0e-12 * dx)
        {
            amrex::Print() << "[treadmill] checkpoint was written at dx=" << dx
                           << " but this run has dx=" << m_odometer.dx << "\n";
            return false;
        }
        m_odometer.cells = cells;
        // The log file survived the restart, so its header did too: appending
        // is right and re-writing the header would truncate the history.
        m_wrote_header = true;
        // The cores have moved with the data, so the search must be re-seeded
        // from scratch rather than from stale positions.
        m_seeded = false;
        amrex::Print() << "[treadmill] restored odometer " << cells
                       << " cells (" << m_odometer.length() << ")\n";
        return true;
    }

  private:
    static std::string state_file(const std::string &a_dir)
    {
        return a_dir + "/Treadmill";
    }

    /// Locate every sector's core and reduce them to one number.
    Step measure(const amrex::MultiFab &a_state, const amrex::Geometry &a_geom)
    {
        Step out;
        if (!m_seeded)
        {
            m_seeds.assign(m_sectors.size(),
                           std::array<double, AMREX_SPACEDIM>{});
            for (auto &seed : m_seeds)
            {
                seed = m_center;
            }
        }

        double sum_axis = 0.0;
        int found       = 0;
        std::vector<double> axis_positions(m_sectors.size(), 0.0);
        for (size_t s = 0; s < m_sectors.size(); ++s)
        {
            // The first measurement has no previous position to search around,
            // so it searches the whole domain; later ones use a ball, which is
            // what keeps departed wake from pulling the answer backwards.
            const double ball = m_seeded ? m_params.ball_radius : -1.0;
            const SectorCore core = track_sector_core(
                a_state, a_geom, m_sectors[s], m_seeds[s][0], m_seeds[s][1],
                m_seeds[s][2], ball);
            if (!core.found)
            {
                continue;
            }
            m_seeds[s] = {AMREX_D_DECL(core.x, core.y, core.z)};
            axis_positions[s] = m_seeds[s][m_params.axis];
            sum_axis += axis_positions[s];
            ++found;
        }

        if (found != static_cast<int>(m_sectors.size()))
        {
            // A missing sector means the midpoint is not the midpoint any more.
            // Refuse to act on it rather than shift on half the pair.
            return out;
        }

        m_seeded          = true;
        out.tracked       = true;
        out.midpoint_grid = sum_axis / static_cast<double>(m_sectors.size());
        out.core0         = axis_positions.empty() ? 0.0 : axis_positions[0];
        out.core1 = (axis_positions.size() > 1) ? axis_positions[1] : out.core0;
        return out;
    }

    GridTreadmillParams m_params{};
    std::vector<SectorFieldSet> m_sectors{};
    std::vector<double> m_asymptotic_values{};
    std::array<double, AMREX_SPACEDIM> m_center{};
    std::vector<std::array<double, AMREX_SPACEDIM>> m_seeds{};
    GridTreadmill::Odometer m_odometer{};
    std::unique_ptr<TreadmillLog> m_log{};
    int m_shift_cells{0};
    bool m_seeded{false};
    bool m_configured{false};
    bool m_wrote_header{false};
};

/// Process-wide instance.  The odometer has to outlive any one call site --
/// it is written at checkpoint time, read at restart time and updated every
/// step -- and there is exactly one grid being recentred.
inline RecentringBox &recentring_box()
{
    static RecentringBox instance;
    return instance;
}

#endif /* RECENTRINGBOX_HPP_ */
