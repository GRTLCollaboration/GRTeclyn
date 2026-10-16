/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

// Calculate the desired quantity

#ifndef SPECTRALWISHES_IMPL_HPP_
#define SPECTRALWISHES_IMPL_HPP_

#include "SpectralWishes.hpp"
#include "StateVariables.hpp"

template <class matter_t>
AMREX_FORCE_INLINE amrex::Real SpectralWishes<matter_t>::get_ncomp()
{
    int ncomp{0};

    // isStateVariable also fills in ncomp as well if it is a state variable
    if (amrex::AmrLevel::isStateVariable(m_var_name, m_state_index, ncomp))
    {
        AMREX_ASSERT(ncomp < NUM_VARS);
    }
    else
    {
        // Not a state variable so get the component number from the derive list
        auto &derive_lst = amrex::AmrLevel::get_derive_lst();
        const auto d     = derive_lst.get(m_var_name);

        for (int i = 0; i < d->numDerive(); ++i)
        {
            if (d->variableName(i) == m_var_name)
            {
                ncomp = i;
                break;
            }
        }

        AMREX_ASSERT(ncomp < d->numDerive());
    }

    return ncomp;
};

template <class matter_t>
AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute(GRAmr *gramr_ptr)
{
    amrex::Real result{0.};

    // No default case so the compiler warns (-Wswitch) if a new
    // SpectralWishOperation is added without being handled here. Invalid
    // parameter values are already rejected in fill_params().
    switch (m_operation)
    {
    case SpectralWishOperation::Mean:
        result = compute_mean(gramr_ptr);
        break;
    case SpectralWishOperation::Variance:
        result = compute_variance(gramr_ptr);
        break;
    }

    return result;
};

template <class matter_t>
AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute_mean(GRAmr *gramr_ptr)

{
    AMREX_ASSERT(gramr_ptr != nullptr);

    amrex::Real mean{0.};

    auto out_mf = gramr_ptr->derive(m_var_name, m_time, m_lev, m_ngrow);
    auto ncomp  = get_ncomp();
    mean        = out_mf->sum(ncomp);

    // Calculate total number of cells on this level

    auto n_cells = gramr_ptr->CountCells(m_lev);

    mean /= static_cast<amrex::Real>(n_cells);

    return mean;
}

template <class matter_t>
AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute_variance(GRAmr *gramr_ptr)
{

    AMREX_ASSERT(gramr_ptr != nullptr);

    amrex::Real var{0.};
    amrex::Real mean = compute_mean(gramr_ptr);

    amrex::Real mean_sq{0.};

    int ncomp = get_ncomp();

    auto out_mf = gramr_ptr->derive(m_var_name, m_time, m_lev, m_ngrow);

    out_mf->plus(-1.0 * mean, m_ngrow);
    amrex::MultiFab::Multiply(*out_mf, *out_mf, ncomp, ncomp, 1, m_ngrow);
    auto n_cells = gramr_ptr->CountCells(m_lev);

    var = (out_mf->sum(ncomp)) / static_cast<amrex::Real>(n_cells - 1);

    return (var);
}

#endif /* SPECTRALWISHES_IMPL_HPP_ */
