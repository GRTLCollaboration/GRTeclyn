/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Calculate the desired quantity

#ifndef SPECTRALWISHES_IMPL_HPP_
#define SPECTRALWISHES_IMPL_HPP_

#include "SpectralWishes.hpp"

template <class matter_t>
template <class enum_t>
enum_t SpectralWishes<matter_t>::to_enum(const GRParmParse &a_pp,
                                         const char *a_name,
                                         const std::string &a_value)
{
    try
    {
        return amrex::getEnumCaseInsensitive<enum_t>(a_value);
    }
    catch (const std::runtime_error &)
    {
        std::string valid_options;
        for (const auto &name : amrex::getEnumNameStrings<enum_t>())
        {
            valid_options += " " + name;
        }
        a_pp.error(a_name, "unknown option \"" + a_value +
                               "\". Valid options are:" + valid_options);
    }
    return enum_t{};
}

template <class matter_t> void SpectralWishes<matter_t>::fill_params()
{
    GRParmParse spectral_wishes_pp(m_param_prefix);

    std::vector<std::string> operations;
    std::vector<std::string> variables;
    spectral_wishes_pp.queryarr("operations", operations);
    spectral_wishes_pp.queryarr("variables", variables);

    if (operations.size() != variables.size())
    {
        // Report against a parameter that exists so its value is printed
        spectral_wishes_pp.error(
            variables.empty() ? "operations" : "variables",
            "operations and variables must have the same number of "
            "entries");
    }

    m_wishes.clear();
    m_wishes.reserve(static_cast<amrex::Long>(operations.size()));

    for (std::size_t iwish = 0; iwish < operations.size(); ++iwish)
    {
        const auto operation = to_enum<SpectralWishOperation>(
            spectral_wishes_pp, "operations", operations[iwish]);

        int chi_state_index{0};
        int chi_comp{0};
        if (operation == SpectralWishOperation::VolumeMean &&
            !amrex::AmrLevel::isStateVariable(s_chi_name, chi_state_index,
                                              chi_comp))
        {
            spectral_wishes_pp.error("operations",
                                     "VolumeMean needs the state variable chi");
        }

        m_wishes.push_back(make_wish(operation, variables[iwish]));
    }

    // Everything to derive: each wish's source and chi for volume weights
    m_sources.clear();
    for (const auto &wish : m_wishes)
    {
        add_source(wish.source);
        if (wish.operation == SpectralWishOperation::VolumeMean)
        {
            add_source(s_chi_name);
        }
    }

    // Levels: by default only the calling level. Coarser levels are not at
    // the same time as the calling level after it has been advanced.
    int amr_max_level{0};
    GRParmParse amr_pp("amr");
    amr_pp.query("max_level", amr_max_level);

    m_min_level = m_lev;
    spectral_wishes_pp.query("min_level", m_min_level);
    if (m_min_level < m_lev || m_min_level > amr_max_level)
    {
        spectral_wishes_pp.error(
            "min_level", "must be between the level the wishes are computed "
                         "from (" +
                             std::to_string(m_lev) + ") and amr.max_level");
    }

    m_max_level = m_lev;
    spectral_wishes_pp.query("max_level", m_max_level);
    if (m_max_level < 0)
    {
        m_max_level = amr_max_level;
    }
    if (m_max_level < m_min_level || m_max_level > amr_max_level)
    {
        spectral_wishes_pp.error("max_level",
                                 "must be -1 (the finest level) or between "
                                 "min_level and amr.max_level");
    }

    std::string level_mode{"Composite"};
    spectral_wishes_pp.query("level_mode", level_mode);
    m_level_mode = to_enum<SpectralWishLevelMode>(spectral_wishes_pp,
                                                  "level_mode", level_mode);

    // Output file parameters
    spectral_wishes_pp.query("interval", m_interval);

    GRParmParse grteclyn_pp("grteclyn");
    std::string output_path{"."};
    grteclyn_pp.query("output_path", output_path);
    m_data_path = output_path + "/data/";
    spectral_wishes_pp.query("path", m_data_path);

    m_file_prefix = m_param_prefix;
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
    wish_t wish{.operation   = a_operation,
                .variable    = a_variable,
                .source      = a_variable,
                .source_comp = 0,
                .num_comps   = 1};

    int state_index{0};
    int state_comp{0};
    if (amrex::AmrLevel::isStateVariable(a_variable, state_index, state_comp))
    {
        // derive() gives a single component copy of a state variable
        AMREX_ALWAYS_ASSERT(state_index == m_state_index);
        return wish;
    }

    // Not a state variable so get the component number from the derive list
    const auto *derive_rec = amrex::AmrLevel::get_derive_lst().get(a_variable);
    if (derive_rec != nullptr)
    {
        // derive() gives every component of the derive record so find
        // this variable's position in it
        for (int icomp = 0; icomp < derive_rec->numDerive(); ++icomp)
        {
            if (derive_rec->variableName(icomp) == a_variable)
            {
                wish.source      = derive_rec->name();
                wish.source_comp = icomp;
                return wish;
            }
        }
    }

    GRParmParse spectral_wishes_pp(m_param_prefix);
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
    const auto add_labels = [&](const std::string &a_suffix)
    {
        for (const auto &wish : m_wishes)
        {
            labels.push_back(amrex::getEnumNameString(wish.operation) + "_" +
                             wish.variable + a_suffix);
        }
    };

    switch (m_level_mode)
    {
    case SpectralWishLevelMode::Composite:
        add_labels("");
        break;
    case SpectralWishLevelMode::Each:
        for (int level = m_min_level; level <= m_max_level; ++level)
        {
            add_labels("_lev" + std::to_string(level));
        }
        break;
    }
    return labels;
}

template <class matter_t>
void SpectralWishes<matter_t>::add_source(const std::string &a_source)
{
    for (const auto &source : m_sources)
    {
        if (source.name == a_source)
        {
            return;
        }
    }

    int num_comps{1};
    int state_index{0};
    int state_comp{0};
    if (!amrex::AmrLevel::isStateVariable(a_source, state_index, state_comp))
    {
        const auto *derive_rec =
            amrex::AmrLevel::get_derive_lst().get(a_source);
        AMREX_ALWAYS_ASSERT(derive_rec != nullptr);
        num_comps = derive_rec->numDerive();
    }
    m_sources.push_back({.name = a_source, .num_comps = num_comps});
}

template <class matter_t>
int SpectralWishes<matter_t>::source_start(const std::string &a_source) const
{
    int start{0};
    for (const auto &source : m_sources)
    {
        if (source.name == a_source)
        {
            return start;
        }
        start += source.num_comps;
    }
    amrex::Abort("SpectralWishes: " + a_source + " has not been derived");
    return -1;
}

template <class matter_t>
amrex::Vector<amrex::MultiFab>
SpectralWishes<matter_t>::derive_sources(GRAmr *gramr_ptr, int a_min_level,
                                         int a_max_level) const
{
    int num_comps{0};
    for (const auto &source : m_sources)
    {
        num_comps += source.num_comps;
    }

    amrex::Vector<amrex::MultiFab> hierarchy;
    hierarchy.reserve(a_max_level - a_min_level + 1);
    for (int level = a_min_level; level <= a_max_level; ++level)
    {
        // No ghost cells as the reductions only use valid cells
        auto &level_data = hierarchy.emplace_back(
            gramr_ptr->boxArray(level), gramr_ptr->DistributionMap(level),
            num_comps, 0);

        int dcomp{0};
        for (const auto &source : m_sources)
        {
            gramr_ptr->getLevel(level).derive(source.name, m_time, level_data,
                                              dcomp);
            dcomp += source.num_comps;
        }
    }
    return hierarchy;
}

template <class matter_t>
amrex::Vector<amrex::Real> SpectralWishes<matter_t>::compute(GRAmr *gramr_ptr)
{
    AMREX_ASSERT(gramr_ptr != nullptr);

    // Levels finer than the finest level present have no data (yet), but
    // still give a result so the output has the same columns at every step
    const int finest_level        = gramr_ptr->finestLevel();
    const auto compute_if_present = [&](int a_min_level, int a_max_level)
    {
        if (a_min_level > finest_level)
        {
            return amrex::Vector<amrex::Real>(
                m_wishes.size(), std::numeric_limits<amrex::Real>::quiet_NaN());
        }
        return compute_on_levels(gramr_ptr, a_min_level,
                                 std::min(a_max_level, finest_level));
    };

    amrex::Vector<amrex::Real> results;
    switch (m_level_mode)
    {
    case SpectralWishLevelMode::Composite:
        results = compute_if_present(m_min_level, m_max_level);
        break;
    case SpectralWishLevelMode::Each:
        for (int level = m_min_level; level <= m_max_level; ++level)
        {
            const auto level_results = compute_if_present(level, level);
            results.insert(results.end(), level_results.begin(),
                           level_results.end());
        }
        break;
    }
    return results;
}

template <class matter_t>
amrex::Vector<amrex::Real>
SpectralWishes<matter_t>::compute_on_levels(GRAmr *gramr_ptr, int a_min_level,
                                            int a_max_level) const
{
    amrex::Vector<amrex::Real> results;
    if (m_wishes.empty())
    {
        return results;
    }
    results.reserve(m_wishes.size());

    const auto hierarchy = derive_sources(gramr_ptr, a_min_level, a_max_level);
    const auto data      = amrex::GetVecOfConstPtrs(hierarchy);

    amrex::Vector<amrex::Geometry> geometry;
    amrex::Vector<amrex::IntVect> refinement_ratios;
    for (int level = a_min_level; level <= a_max_level; ++level)
    {
        geometry.push_back(gramr_ptr->Geom(level));
        if (level < a_max_level)
        {
            refinement_ratios.push_back(gramr_ptr->refRatio(level));
        }
    }
    const AMRReductions::AMRReducer reducer(std::move(geometry),
                                            std::move(refinement_ratios));

    // AMRReducer::compute takes a fixed list of reductions while the wishes
    // are only known at run time, so integrate one moment at a time
    const auto integrate = [&](const SpectralWishMoment &a_moment)
    { return reducer.compute(data, a_moment)[0]; };

    const amrex::Real volume = integrate({.power = 0});
    std::optional<amrex::Real> proper_volume;

    for (const auto &wish : m_wishes)
    {
        const SpectralWishMoment first_moment{
            .comp      = source_start(wish.source) + wish.source_comp,
            .num_comps = wish.num_comps};

        amrex::Real result{0.};

        // No default case so the compiler warns (-Wswitch) if a new
        // SpectralWishOperation is added without being handled here. Invalid
        // parameter values are already rejected in fill_params().
        switch (wish.operation)
        {
        case SpectralWishOperation::Mean:
            result = integrate(first_moment) / volume;
            break;
        case SpectralWishOperation::Variance:
        {
            // Two passes, as the mean can be much larger than the deviations
            auto second_moment  = first_moment;
            second_moment.power = 2;
            second_moment.shift = integrate(first_moment) / volume;
            result              = integrate(second_moment) / volume;
            break;
        }
        case SpectralWishOperation::VolumeMean:
        {
            const int chi_comp = source_start(s_chi_name);
            if (!proper_volume)
            {
                proper_volume = integrate({.power = 0, .chi_comp = chi_comp});
            }
            auto weighted_moment     = first_moment;
            weighted_moment.chi_comp = chi_comp;
            result = integrate(weighted_moment) / *proper_volume;
            break;
        }
        }

        results.push_back(result);
    }

    return results;
}

#endif /* SPECTRALWISHES_IMPL_HPP_ */
