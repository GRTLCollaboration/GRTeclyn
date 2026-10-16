/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Calculate the desired quantity

#ifndef SPECTRALWISHES_IMPL_HPP_
#define SPECTRALWISHES_IMPL_HPP_

#include "SpectralWishes.hpp"

#include <AMReX_ParReduce.H>

template <class matter_t> void SpectralWishes<matter_t>::fill_params()
{
    GRParmParse spectral_wishes_pp("spectral_wishes");

    std::vector<std::string> operations;
    std::vector<std::string> variables;
    spectral_wishes_pp.queryarr("operations", operations);
    spectral_wishes_pp.queryarr("variables", variables);

    if (operations.size() != variables.size())
    {
        // Report against a parameter that exists so its value is printed
        spectral_wishes_pp.error(
            variables.empty() ? "operations" : "variables",
            "spectral_wishes.operations and spectral_wishes.variables must "
            "have the same number of entries");
    }

    m_wishes.clear();
    m_wishes.reserve(static_cast<amrex::Long>(operations.size()));

    for (std::size_t iwish = 0; iwish < operations.size(); ++iwish)
    {
        SpectralWishOperation operation{};
        try
        {
            operation = amrex::getEnumCaseInsensitive<SpectralWishOperation>(
                operations[iwish]);
        }
        catch (const std::runtime_error &)
        {
            std::string valid_options;
            for (const auto &name :
                 amrex::getEnumNameStrings<SpectralWishOperation>())
            {
                valid_options += " " + name;
            }
            spectral_wishes_pp.error(
                "operations", "unknown operation \"" + operations[iwish] +
                                  "\". Valid options are:" + valid_options);
        }

        m_wishes.push_back(make_wish(operation, variables[iwish]));
    }

    // Output file parameters
    spectral_wishes_pp.query("interval", m_interval);

    GRParmParse grteclyn_pp("grteclyn");
    std::string output_path{"."};
    grteclyn_pp.query("output_path", output_path);
    m_data_path = output_path + "/data/";
    spectral_wishes_pp.query("path", m_data_path);

    m_file_prefix = "spectral_wishes";
    spectral_wishes_pp.query("file_prefix", m_file_prefix);
}

template <class matter_t> void SpectralWishes<matter_t>::write(GRAmr *gramr_ptr)
{
    AMREX_ASSERT(gramr_ptr != nullptr);

    // In post_timestep the step count has already been incremented, so the
    // first step of a run is 1
    const int step = gramr_ptr->levelSteps(m_lev);
    if (!is_write_step(step))
    {
        return;
    }

    const bool first_step          = (step == m_interval);
    const amrex::Real restart_time = gramr_ptr->get_restart_time();
    // SmallDataIO uses this to recognise the first write after a restart so
    // it must be the time between writes rather than the timestep
    const amrex::Real write_dt = m_interval * gramr_ptr->dtLevel(m_lev);

    if (first_step || m_time < restart_time + 1.5 * write_dt)
    {
        FilesystemTools::ensure_directory_exists(m_data_path);
    }

    // All ranks take part in the computation but only the IO processor
    // writes
    const auto results = compute(gramr_ptr);

    SmallDataIO output_file(get_file_prefix(), write_dt, m_time, restart_time,
                            SmallDataIO::APPEND, first_step);
    // Remove any data from a previous run beyond the restart time
    output_file.remove_duplicate_time_data();
    if (first_step)
    {
        output_file.write_header_line(get_labels());
    }
    output_file.write_time_data_line(results);
}

template <class matter_t>
typename SpectralWishes<matter_t>::wish_t
SpectralWishes<matter_t>::make_wish(SpectralWishOperation a_operation,
                                    const std::string &a_variable) const
{
    wish_t wish{a_operation, a_variable, a_variable, 0};

    int state_index{0};
    int state_comp{0};
    if (amrex::AmrLevel::isStateVariable(a_variable, state_index, state_comp))
    {
        // derive() returns a single component copy of a state variable
        AMREX_ALWAYS_ASSERT(state_index == m_state_index);
        return wish;
    }

    // Not a state variable so get the component number from the derive list
    const auto *derive_rec = amrex::AmrLevel::get_derive_lst().get(a_variable);
    if (derive_rec != nullptr)
    {
        // derive() returns every component of the derive record so find
        // this variable's position in it
        for (int icomp = 0; icomp < derive_rec->numDerive(); ++icomp)
        {
            if (derive_rec->variableName(icomp) == a_variable)
            {
                wish.source = derive_rec->name();
                wish.comp   = icomp;
                return wish;
            }
        }
    }

    GRParmParse spectral_wishes_pp("spectral_wishes");
    spectral_wishes_pp.error("variables",
                             "\"" + a_variable +
                                 "\" is not a state variable or a single "
                                 "component of a derived variable");
    return wish;
}

template <class matter_t>
amrex::Vector<std::string> SpectralWishes<matter_t>::get_labels() const
{
    amrex::Vector<std::string> labels;
    labels.reserve(m_wishes.size());
    for (const auto &wish : m_wishes)
    {
        labels.push_back(amrex::getEnumNameString(wish.operation) + "_" +
                         wish.variable);
    }
    return labels;
}

template <class matter_t>
AMREX_FORCE_INLINE amrex::Vector<amrex::Real>
SpectralWishes<matter_t>::compute(GRAmr *gramr_ptr)
{
    AMREX_ASSERT(gramr_ptr != nullptr);

    amrex::Vector<amrex::Real> results;
    results.reserve(m_wishes.size());

    const amrex::Long num_cells = gramr_ptr->CountCells(m_lev);

    // Derive each source once, however many wishes use it (e.g. Ham and Mom1
    // both come from the same constraints calculation)
    std::map<std::string, std::unique_ptr<amrex::MultiFab>> derived;

    for (const auto &wish : m_wishes)
    {
        auto &source_mf = derived[wish.source];
        if (source_mf == nullptr)
        {
            source_mf = gramr_ptr->derive(wish.source, m_time, m_lev, m_ngrow);
        }

        const amrex::Real mean = compute_mean(*source_mf, wish.comp, num_cells);

        amrex::Real result{0.};

        // No default case so the compiler warns (-Wswitch) if a new
        // SpectralWishOperation is added without being handled here. Invalid
        // parameter values are already rejected in fill_params().
        switch (wish.operation)
        {
        case SpectralWishOperation::Mean:
            result = mean;
            break;
        case SpectralWishOperation::Variance:
            result = compute_variance(*source_mf, wish.comp, num_cells, mean);
            break;
        }

        results.push_back(result);
    }

    return results;
};

template <class matter_t>
amrex::Real SpectralWishes<matter_t>::compute_mean(const amrex::MultiFab &a_mf,
                                                   int a_comp,
                                                   amrex::Long a_num_cells)
{
    // MultiFab::sum only sums valid cells and reduces over MPI ranks
    return a_mf.sum(a_comp) / static_cast<amrex::Real>(a_num_cells);
}

template <class matter_t>
amrex::Real
SpectralWishes<matter_t>::compute_variance(const amrex::MultiFab &a_mf,
                                           int a_comp, amrex::Long a_num_cells,
                                           amrex::Real a_mean)
{
    const auto &mf_arrays = a_mf.const_arrays();

    // Sum of squared deviations over valid cells without modifying a_mf
    amrex::Real sum_sq = amrex::ParReduce(
        amrex::TypeList<amrex::ReduceOpSum>{}, amrex::TypeList<amrex::Real>{},
        a_mf, amrex::IntVect(0),
        [=] AMREX_GPU_DEVICE(int ibox, int i, int j,
                             int k) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            const amrex::Real deviation =
                mf_arrays[ibox](i, j, k, a_comp) - a_mean;
            return {deviation * deviation};
        });

    // ParReduce only reduces over this rank's boxes
    amrex::ParallelAllReduce::Sum(sum_sq,
                                  amrex::ParallelDescriptor::Communicator());

    return sum_sq / static_cast<amrex::Real>(a_num_cells - 1);
}

#endif /* SPECTRALWISHES_IMPL_HPP_ */
