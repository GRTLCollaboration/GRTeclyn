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
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::get_ncomp()
{
    int ncomp{0};

    // isStateVariable also fills in ncomp as well if it is a state variable
    if (amrex::AmrLevel::isStateVariable(m_var_name, m_state_index, ncomp))
    {
        AMREX_ASSERT(ncomp < NUM_VARS);
    }
    else
    {
        std::vector<std::string> derive_names;
        auto &derive_lst = amrex::AmrLevel::get_derive_lst();
        const auto d     = derive_lst.get(m_var_name);
        // const std::list<amrex::DeriveRec> &dlist = derive_lst.dlist();
        // for (auto const &d : dlist)
        // {
        //       if (amrex::Amr::isDerivePlotVar(d.name()))

        for (int i = 0; i < d->numDerive(); ++i)
        {
            if (d->variableName(i) == m_var_name)
            {
                ncomp = i;
                amrex::Print()
                    << d->variableName(ncomp) << " " << ncomp << "\n";

                break;
                //            derive_names.push_back(d.name());
                //            num_derive += d.numDerive();
            }
            else
            {
                amrex::Print() << d->variableName(i) << "\n";
                //                    ncomp_derive++;
            }
        }

        AMREX_ASSERT(ncomp < d->numDerive());
    }

    // mean = derive_mf->sum(
    //     ncomp_derive); // sum needs a component so check this;

    return ncomp;
};

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute(GRAmr *gramr_ptr)
{
    amrex::Real result;

    if (m_operation == "Mean")
    {
        result = compute_mean(gramr_ptr);
    }
    else if (m_operation == "Variance")
    {
        result = compute_variance(gramr_ptr);
    }
    else
    {
        amrex::Abort("SpectralWishes: Option not found\n");
    }

    return result;
};

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute_mean(GRAmr *gramr_ptr)

{
    AMREX_ASSERT(gramr_ptr != nullptr);

    //  const std::string var_name = m_var_name;
    amrex::Real mean{0.};

    //    int ncomp{0}; // this will be set in isStateVariable or get_derived_mf

    // if (amrex::AmrLevel::isStateVariable(m_var_name, m_state_index, ncomp))
    // {
    //     // Use the state multifab
    //     // Sum all the values in ncomp

    //     //        auto out_mf = get_state_data(state_index);
    //     mean = src_mf.sum(ncomp);
    // }
    // else
    // {
    //     auto out_mf = get_derived_mf(gramr_ptr, ncomp, ngrow);
    //     mean        = out_mf->sum(ncomp);
    // }
    auto out_mf = gramr_ptr->derive(m_var_name, m_time, m_lev, m_ngrow);
    auto ncomp  = get_ncomp();
    mean        = out_mf->sum(ncomp);

    //    AMREX_ASSERT(ncomp_state < NUM_VARS);

    // Calculate total number of cells on this level

    // auto n_cells = amrex::AmrLevel::countCells();

    auto n_cells = gramr_ptr->CountCells(m_lev);

    mean /= static_cast<amrex::Real>(n_cells);

    amrex::AllPrint() << "Mean of " << m_var_name << " on Rank "
                      << amrex::ParallelDescriptor::MyProc() << ": " << mean
                      << "\n";
    amrex::AllPrint() << "Number of cells on Rank "
                      << amrex::ParallelDescriptor::MyProc() << ": " << n_cells
                      << "\n";

    return mean;
}

template <class matter_t>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
SpectralWishes<matter_t>::compute_variance(GRAmr *gramr_ptr)
{

    AMREX_ASSERT(gramr_ptr != nullptr);

    amrex::Real var{0.};
    amrex::Real mean = compute_mean(gramr_ptr);

    amrex::Real mean_sq{0.};

    int ncomp = get_ncomp();

    auto out_mf = gramr_ptr->derive(m_var_name, m_time, m_lev, m_ngrow);

    // if (amrex::AmrLevel::isStateVariable(m_var_name, m_state_index, ncomp))
    // {
    //     // Use the state multifab
    //     // Sum all the values in ncomp

    //     //        auto out_mf = get_state_data(state_index);
    //     //        mean_sq = amrex::MultiFab::Dot(src_mf, ncomp, 1, ngrow);

    //   //      mf = std::make_unique<MultiFab>(state[index].boxArray(), dmap,
    //   1, ngrow, MFInfo(), *m_factory);
    //   //      FillPatch(*this,*mf,ngrow,time,index,scomp,1,0);

    //     amrex::MultiFab::Multiply(src_mf, src_mf, ncomp, ncomp, 1, ngrow);
    // }
    // else
    // {
    //     auto out_mf = get_derived_mf(gramr_ptr, ncomp, ngrow);

    //     //    mean_sq = amrex::MultiFab::Dot(*out_mf, ncomp, 1, ngrow);
    //     amrex::MultiFab::Multiply(*out_mf, *out_mf, ncomp, ncomp, 1, ngrow);
    // }

    out_mf->plus(-1.0 * mean, m_ngrow);
    amrex::MultiFab::Multiply(*out_mf, *out_mf, ncomp, ncomp, 1, m_ngrow);
    auto n_cells = gramr_ptr->CountCells(m_lev);

    var = (out_mf->sum(ncomp)) / static_cast<amrex::Real>(n_cells - 1);

    //    var = mean_sq - mean * mean;

    amrex::AllPrint() << "Variance of " << m_var_name << " on Rank "
                      << amrex::ParallelDescriptor::MyProc() << ": " << var
                      << "\n";

    return (var);
}

#endif /* SPECTRALWISHES_IMPL_HPP_ */
