/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPECTRALWISHES_HPP_
#define SPECTRALWISHES_HPP_

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include "FilesystemTools.hpp"
#include "GRParmParse.hpp"
#include "SmallDataIO.hpp"

#include <array>
#include <map>
#include <memory>
#include <string>

//! The operations SpectralWishes can perform. Declared outside the class
//! template so the type does not depend on matter_t. Parameter values are
//! matched case-insensitively (e.g. "mean", "Mean" and "MEAN" all work).
//!   Mean:       coordinate average, sum(f) / N
//!   Variance:   sum((f - Mean)^2) / (N - 1)
//!   VolumeMean: average weighted by the proper volume,
//!               sum(f sqrt(gamma)) / sum(sqrt(gamma)), sqrt(gamma) = chi^-3/2
AMREX_ENUM(SpectralWishOperation, Mean, Variance, VolumeMean);

/*!
 * Computes global quantities (e.g. the mean or variance) of state or derived
 * variables on one AMR level. Each wish is an (operation, variable) pair read
 * from the parameter file, e.g.
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
        //! Component of the variable within the MultiFab derive() returns
        int comp;
    };

    //! Constructor
    SpectralWishes(int state_index, const amrex::Real time, int lev, int ngrow)
        : m_time(time), m_lev(lev), m_state_index(state_index), m_ngrow(ngrow)
    {
        fill_params();
    };

    void fill_params();

    //! Compute every wish, in the order given in the parameter file
    AMREX_FORCE_INLINE amrex::Vector<amrex::Real> compute(GRAmr *gramr_ptr);

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
    //! Find which MultiFab (source) and component a variable lives in
    [[nodiscard]] wish_t make_wish(SpectralWishOperation a_operation,
                                   const std::string &a_variable) const;

    //! Compute members
    static amrex::Real compute_mean(const amrex::MultiFab &a_mf, int a_comp,
                                    amrex::Long a_num_cells);

    static amrex::Real compute_variance(const amrex::MultiFab &a_mf, int a_comp,
                                        amrex::Long a_num_cells,
                                        amrex::Real a_mean);

    static amrex::Real compute_volume_mean(const amrex::MultiFab &a_mf,
                                           int a_comp,
                                           const amrex::MultiFab &a_chi_mf);

    //! State variable holding the conformal factor used for volume weights
    static constexpr const char *s_chi_name = "chi";

    matter_t m_matter;        //!< The matter object, e.g. a scalar field
    const amrex::Real m_time; // Time as measured by AmrLevel
    const int m_lev;          // AMR level
    int m_state_index;
    const int
        m_ngrow; // Number of ghost cells to use in Spectral Wish calculation,
                 // should match what is used in the state MF in most cases.
    amrex::Vector<wish_t> m_wishes; // The wishes from the parameter file
    int m_interval{1};              // Steps between writes, -1 disables
    std::string m_data_path;        // Directory of the output file
    std::string m_file_prefix;      // Output file name without extension
};

#include "SpectralWishes.impl.hpp"

#endif /* SPECTRALWISHES_HPP_ */
