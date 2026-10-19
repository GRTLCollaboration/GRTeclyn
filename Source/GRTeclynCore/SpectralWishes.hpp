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

#include <cmath>
#include <optional>
#include <string>

//! The operations SpectralWishes can perform. Declared outside the class
//! template so the type does not depend on matter_t. Parameter values are
//! matched case-insensitively (e.g. "mean", "Mean" and "MEAN" all work).
//!   Mean:       coordinate average, int f dV / V
//!   Variance:   int (f - Mean)^2 dV / V
//!   VolumeMean: average weighted by the proper volume,
//!               int f sqrt(gamma) dV / int sqrt(gamma) dV,
//!               sqrt(gamma) = chi^-3/2
//! where V is the volume the wishes are computed over. On a single level every
//! cell has the same dV so e.g. Mean = sum(f) / N.
AMREX_ENUM(SpectralWishOperation, Mean, Variance, VolumeMean);

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
 * variables with AMRReductions. Each wish is an (operation, variable) pair
 * read from the parameter file, e.g.
 *
 *     spectral_wishes.operations = Mean Variance Mean
 *     spectral_wishes.variables  = Ham  Ham      chi
 *
 * The two lists must have the same length. Each variable is derived only once
 * per call to compute(), however many wishes use it.
 *
 * write() appends the wishes to <path><file_prefix>.dat every interval steps
 * of the level (-1 disables writing), with one column per wish:
 *
 *     spectral_wishes.interval    = 1                  # default
 *     spectral_wishes.path        = <grteclyn.output_path>/data/ # default
 *     spectral_wishes.file_prefix = spectral_wishes    # default
 */
template <class matter_t> class SpectralWishes
{
  public:
    //! A single (operation, variable) pair
    struct wish_t
    {
        SpectralWishOperation operation;
        std::string variable;
        //! Name passed to derive(): the variable itself for a state
        //! variable or the name of its derive record for a derived variable
        std::string source;
        //! Component of the variable within what derive() gives for source
        int source_comp;
        //! Number of consecutive components (> 1 for a magnitude)
        int num_comps;
    };

    //! Constructor
    SpectralWishes(int state_index, const amrex::Real time, int lev, int ngrow)
        : m_time(time), m_lev(lev), m_state_index(state_index), m_ngrow(ngrow)
    {
        fill_params();
    };

    void fill_params();

    //! Compute every wish, in the order given in the parameter file
    amrex::Vector<amrex::Real> compute(GRAmr *gramr_ptr);

    //! The wishes, in the order given in the parameter file
    [[nodiscard]] const amrex::Vector<wish_t> &get_wishes() const
    {
        return m_wishes;
    }

    //! A label for each wish, e.g. "Mean_Ham"
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
    //! A name passed to derive() and the number of components it gives
    struct source_t
    {
        std::string name;
        int num_comps;
    };

    //! Find which source and component a variable lives in
    [[nodiscard]] wish_t make_wish(SpectralWishOperation a_operation,
                                   const std::string &a_variable) const;

    //! Add a_source to m_sources unless it is already there
    void add_source(const std::string &a_source);

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
    // Unused: the reductions only use valid cells so the variables are
    // derived without ghost cells. Kept so existing callers do not change.
    const int m_ngrow;
    amrex::Vector<wish_t> m_wishes;    // The wishes from the parameter file
    amrex::Vector<source_t> m_sources; // What to derive for the wishes
    int m_interval{1};                 // Steps between writes, -1 disables
    std::string m_data_path;           // Directory of the output file
    std::string m_file_prefix;         // Output file name without extension
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
