/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(AMRREDUCTIONS_HPP)
#error "This file should only be included through AMRReductions.hpp"
#endif

#ifndef AMRREDUCTIONS_IMPL_HPP
#define AMRREDUCTIONS_IMPL_HPP

template <VariableType var_t>
AMRReductions<var_t>::AMRReductions(const GRAMR &a_gramr,
                                    const int a_base_level)
    : m_base_level(a_base_level) // xxxxx,
// xxxxx      m_coarsest_dx(a_gramr.get_gramrlevels()[0]->get_dx())
{
    set_level_data_vect(a_gramr);
    set_ref_ratios_vect(a_gramr);
    set_domain_volume();
}

template <VariableType var_t>
void AMRReductions<var_t>::set_level_data_vect(const GRAMR &a_gramr)
{
#if 0
//xxxxx
    // this function already checks if the level pointers are null
    const auto gramrlevel_ptrs = a_gramr.get_gramrlevels();
    int num_levels = gramrlevel_ptrs.size();
    m_level_data_ptrs.resize(num_levels);

    for (int ilev = 0; ilev < num_levels; ++ilev)
    {
        m_level_data_ptrs[ilev] = const_cast<GRLevelData *>(
            &gramrlevel_ptrs[ilev]->getLevelData(var_t));
    }
#endif
}

template <VariableType var_t>
void AMRReductions<var_t>::set_ref_ratios_vect(const GRAMR &a_gramr)
{
#if 0
//xxxxx
    // this function already checks if the level pointers are null
    const auto gramrlevel_ptrs = a_gramr.get_gramrlevels();
    int num_levels = gramrlevel_ptrs.size();
    m_ref_ratios.resize(num_levels);

    for (int ilev = 0; ilev < num_levels; ++ilev)
    {
        m_ref_ratios[ilev] = gramrlevel_ptrs[ilev]->refRatio();
    }
#endif
}

template <VariableType var_t> void AMRReductions<var_t>::set_domain_volume()
{
    // first check if m_level_data_ptrs has been set
    AMREX_ASSERT((m_level_data_ptrs.size() > 0) &&
                 (m_level_data_ptrs[0] != nullptr));

    // first calculate the volume assuming each cell on the coarsest level has
    // unit length
    int cell_volume = 0;
    // xxxxx
    // m_level_data_ptrs[0]->disjointBoxLayout().physDomain().size().product();

    // multiply by dx_coarsest to get real volume
    m_domain_volume =
        pow(m_coarsest_dx, AMREX_SPACEDIM) * static_cast<double>(cell_volume);
}

#if 0
//xxxxx
template <VariableType var_t>
amrex::Real AMRReductions<var_t>::min(const Interval &a_vars) const
{
    AMREX_ASSERT(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    BL_PROFILE("AMRReductions::min");
    return computeMin(m_level_data_ptrs, m_ref_ratios, a_vars, m_base_level);
}
#endif

template <VariableType var_t>
amrex::Real AMRReductions<var_t>::min(const int a_var) const
{
    return 0; // xxxxxmin(Interval(a_var, a_var));
}

#if 0
//xxxxx
template <VariableType var_t>
amrex::Real AMRReductions<var_t>::max(const Interval &a_vars) const
{
    AMREX_ASSERT(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    BL_PROFILE("AMRReductions::max");
    return computeMax(m_level_data_ptrs, m_ref_ratios, a_vars, m_base_level);
}
#endif

template <VariableType var_t>
amrex::Real AMRReductions<var_t>::max(const int a_var) const
{
    return 0; // xxxxxmax(Interval(a_var, a_var));
}

#if 0
//xxxxx
template <VariableType var_t>
amrex::Real AMRReductions<var_t>::norm(const Interval &a_vars,
                                const int a_norm_exponent,
                                const bool a_normalize_by_volume) const
{
    AMREX_ASSERT(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    BL_PROFILE("AMRReductions::norm");
    amrex::Real norm = computeNorm(m_level_data_ptrs, m_ref_ratios, m_coarsest_dx,
                            a_vars, a_norm_exponent, m_base_level);
    if (a_normalize_by_volume)
    {
        norm /=
            pow(m_domain_volume, 1.0 / static_cast<double>(a_norm_exponent));
    }

    return norm;
}
#endif

template <VariableType var_t>
amrex::Real AMRReductions<var_t>::norm(const int a_var,
                                       const int a_norm_exponent,
                                       const bool a_normalize_by_volume) const
{
    return 0; // xxxxxnorm(Interval(a_var, a_var), a_norm_exponent,
              // a_normalize_by_volume);
}

#if 0
//xxxxx
template <VariableType var_t>
amrex::Real AMRReductions<var_t>::sum(const Interval &a_vars) const
{
    AMREX_ASSERT(a_vars.begin() >= 0 && a_vars.end() < m_num_vars);
    BL_PROFILE("AMRReductions::sum");
    return computeSum(m_level_data_ptrs, m_ref_ratios, m_coarsest_dx, a_vars,
                      m_base_level);
}
#endif

template <VariableType var_t>
amrex::Real AMRReductions<var_t>::sum(const int a_var) const
{
    return 0; // sum(Interval(a_var, a_var));
}

template <VariableType var_t>
amrex::Real AMRReductions<var_t>::get_domain_volume() const
{
    return m_domain_volume;
}

#endif
