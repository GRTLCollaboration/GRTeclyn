/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPECTRALWISHES_HPP_
#define SPECTRALWISHES_HPP_

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include "AMRReductions.hpp"
#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "SmallDataIO.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>

//! The operations SpectralWishes can perform. Declared outside the class
//! template so the type does not depend on matter_t. Parameter values are
//! matched case-insensitively (e.g. "mean", "Mean" and "MEAN" all work).
//!   Mean:           coordinate average, int f dV / V
//!   Variance:       int (f - Mean)^2 dV / V
//!   VolumeMean:     average weighted by the proper volume,
//!                   int f sqrt(gamma) dV / int sqrt(gamma) dV,
//!                   sqrt(gamma) = chi^-3/2
//!   Integral:       int f dV
//!   VolumeIntegral: integral over the proper volume, int f sqrt(gamma) dV
//!   Min, Max:       minimum and maximum over the cells
//!   L1Norm:         int |f| dV / V
//!   L2Norm:         sqrt(int f^2 dV / V)
//! where V is the volume the wishes are computed over. On a single level every
//! cell has the same dV so e.g. Mean = sum(f) / N. For a magnitude variable
//! (first:last) f is sqrt(sum_c f_c^2) over the components first to last.
AMREX_ENUM(SpectralWishOperation, Mean, Variance, VolumeMean, Integral,
           VolumeIntegral, Min, Max, L1Norm, L2Norm);

//! How the levels in [min_level, max_level] are used (case-insensitive):
//!   Composite: one result per wish over all the levels, counting each point
//!              once with the finest data available there
//!   Each:      one result per wish and level, over that level alone
AMREX_ENUM(SpectralWishLevelMode, Composite, Each);

/*!
 * An AMRReductions reduction for the integrals SpectralWishes needs,
 *
 *     int w (f - shift)^power dV
 *
 * where w = 1, or w = sqrt(gamma) = chi^(-3/2) if chi_comp >= 0, and f is
 * component comp or, if num_comps > 1, the magnitude of the num_comps
 * components from comp. power = 0 gives the (proper) volume.
 */
struct SpectralWishMoment
{
    int comp          = 0;
    int num_comps     = 1;
    int power         = 1;
    amrex::Real shift = 0.0;
    //! Component holding chi to weight by the proper volume, or -1
    int chi_comp = -1;

    using reduce_op_t = amrex::ReduceOpSum;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    cell_value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
               int iz, amrex::Real a_cell_volume) const noexcept
    {
        amrex::Real result = a_cell_volume;
        if (chi_comp >= 0)
        {
            const amrex::Real chi  = a_cell(ix, iy, iz, chi_comp);
            result                /= chi * std::sqrt(chi);
        }
        const amrex::Real deviation =
            AMRReductions::detail::value(a_cell, ix, iy, iz, comp, num_comps) -
            shift;
        for (int ipower = 0; ipower < power; ++ipower)
        {
            result *= deviation;
        }
        return result;
    }

    [[nodiscard]] static amrex::Real finalize(amrex::Real a_integral,
                                              amrex::Real /*a_domain_volume*/)
    {
        return a_integral;
    }
};

/*!
 * Computes global quantities (e.g. the mean or variance) of state or derived
 * variables over one or more AMR levels with AMRReductions. The parameters
 * are read with a prefix (spectral_wishes by default, so several independent
 * sets of wishes can be used in one run). Each wish is an (operation,
 * variable) pair, e.g.
 *
 *     spectral_wishes.operations = Mean Variance Mean
 *     spectral_wishes.variables  = Ham  Ham      chi
 *
 * The two lists must have the same length. A variable is a state variable, a
 * component of a derived variable (e.g. Ham or Mom1) or the magnitude of
 * consecutive state variables or components of one derived variable, written
 * first:last (e.g. Mom1:Mom3). Each variable is derived only once per call to
 * compute(), however many wishes use it.
 *
 * By default the wishes are computed on the level SpectralWishes is called
 * from. Finer levels can be included with
 *
 *     spectral_wishes.min_level  = <calling level> # default
 *     spectral_wishes.max_level  = <calling level> # default, -1 = finest
 *     spectral_wishes.level_mode = Composite       # default, or Each
 *
 * Levels that do not exist (yet) give NaN. Coarser levels than the calling
 * level cannot be used as they are not at the same time.
 *
 * Layers of cells next to the faces of the domain can be left out, e.g. those
 * whose evolution is replaced by an outer boundary condition. The widths are
 * in level-0 cells for each direction (faces in periodic directions are never
 * excluded):
 *
 *     spectral_wishes.exclude_lo = 0 0 0 # default
 *     spectral_wishes.exclude_hi = 0 0 0 # default
 *
 * write() appends the wishes to <path><file_prefix>.dat every interval steps
 * of the calling level (-1 disables writing), with one column per wish (and
 * level in Each mode, e.g. Mean_Ham_lev1):
 *
 *     spectral_wishes.interval    = 1                  # default
 *     spectral_wishes.path        = <grteclyn.output_path>/data/ # default
 *     spectral_wishes.file_prefix = <prefix>           # default
 */
template <class matter_t> class SpectralWishes
{
  public:
    //! A single (operation, variable) pair
    struct wish_t
    {
        SpectralWishOperation operation;
        std::string variable;
        //! Key of the variable's source in m_sources: the name of its derive
        //! record for a derived variable or the variable itself otherwise
        std::string source;
        //! Component of the variable within what derive() gives for source
        int source_comp;
        //! Number of consecutive components (> 1 for a magnitude)
        int num_comps;
    };

    //! Constructor. The parameters are read from a_param_prefix.<name>.
    SpectralWishes(int state_index, const amrex::Real time, int lev, int ngrow,
                   std::string a_param_prefix = "spectral_wishes")
        : m_time(time), m_lev(lev), m_state_index(state_index), m_ngrow(ngrow),
          m_param_prefix(std::move(a_param_prefix))
    {
        fill_params();
    };

    void fill_params();

    //! Compute every wish, in the order given in the parameter file (for each
    //! level in turn in Each mode)
    amrex::Vector<amrex::Real> compute(GRAmr *gramr_ptr);

    //! The wishes, in the order given in the parameter file
    [[nodiscard]] const amrex::Vector<wish_t> &get_wishes() const
    {
        return m_wishes;
    }

    //! A label for each result of compute(), e.g. "Mean_Ham", or
    //! "Mean_Ham_lev1" in Each mode
    [[nodiscard]] amrex::Vector<std::string> get_labels() const;

    //! Whether write() writes on this step of the level
    [[nodiscard]] bool is_write_step(int a_step) const
    {
        return !m_wishes.empty() && m_interval > 0 && a_step % m_interval == 0;
    }

    //! Compute every wish and append them to the output file if this is a
    //! write step. Call after the level has been advanced (e.g. in
    //! specific_post_timestep) on every rank.
    void write(GRAmr *gramr_ptr);

    //! The output file name without the .dat extension
    [[nodiscard]] std::string get_file_prefix() const
    {
        return m_data_path + m_file_prefix;
    }

  protected:
    //! Data to derive: a derive record (state_comp < 0) or the state
    //! components [state_comp, state_comp + num_comps)
    struct source_t
    {
        std::string name;
        int num_comps;
        int state_comp;
    };

    //! Where a single (not first:last) variable lives
    struct variable_t
    {
        bool is_state;
        //! The variable for a state variable, else its derive record's name
        std::string source;
        //! The state component, or the component within the derive record
        int comp;
        //! Number of components of source
        int source_num_comps;
    };

    //! Convert a parameter value to an AMREX_ENUM (case-insensitively), with
    //! an error listing the valid options if there is no match
    template <class enum_t>
    static enum_t to_enum(const GRParmParse &a_pp, const char *a_name,
                          const std::string &a_value);

    //! Where a single variable lives, if it exists
    [[nodiscard]] std::optional<variable_t>
    find_variable(const std::string &a_variable) const;

    //! Make a wish and add what it needs to m_sources
    wish_t make_wish(SpectralWishOperation a_operation,
                     const std::string &a_variable);

    //! Whether the operation needs chi for the proper volume
    static bool uses_proper_volume(SpectralWishOperation a_operation)
    {
        return a_operation == SpectralWishOperation::VolumeMean ||
               a_operation == SpectralWishOperation::VolumeIntegral;
    }

    //! Add a_source to m_sources unless one with the same name is there
    void add_source(const source_t &a_source);

    //! First component of a_source in the data from derive_sources()
    [[nodiscard]] int source_start(const std::string &a_source) const;

    //! Every source, in the order of m_sources, derived into consecutive
    //! components of one MultiFab per level in [a_min_level, a_max_level]
    [[nodiscard]] amrex::Vector<amrex::MultiFab>
    derive_sources(GRAmr *gramr_ptr, int a_min_level, int a_max_level) const;

    //! Every wish computed over levels [a_min_level, a_max_level], counting
    //! each point once with the finest data available there
    [[nodiscard]] amrex::Vector<amrex::Real>
    compute_on_levels(GRAmr *gramr_ptr, int a_min_level, int a_max_level) const;

    //! State variable holding the conformal factor used for volume weights
    static constexpr const char *s_chi_name = "chi";

    matter_t m_matter;        //!< The matter object, e.g. a scalar field
    const amrex::Real m_time; // Time as measured by AmrLevel
    const int m_lev;          // AMR level
    int m_state_index;
    const std::string m_param_prefix; // Prefix of the parameters
    // Unused: the reductions only use valid cells so the variables are
    // derived without ghost cells. Kept so existing callers do not change.
    const int m_ngrow;
    amrex::Vector<wish_t> m_wishes;    // The wishes from the parameter file
    amrex::Vector<source_t> m_sources; // What to derive for the wishes
    int m_min_level{0};                // Coarsest level used
    int m_max_level{0};                // Finest level used
    SpectralWishLevelMode m_level_mode{SpectralWishLevelMode::Composite};
    // Level-0 cells to leave out next to the low and high faces
    AMRReductions::BoundaryExclusion m_exclusion;
    int m_interval{1};         // Steps between writes, -1 disables
    std::string m_data_path;   // Directory of the output file
    std::string m_file_prefix; // Output file name without extension
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
