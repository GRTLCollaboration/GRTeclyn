/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CELL_HPP_
#define CELL_HPP_

#include <AMReX_Array4.H>

using namespace amrex::disabled;

template <template <typename> class vars_t, class data_t>
AMREX_GPU_HOST_DEVICE void store_vars(const amrex::CellData<data_t> &cell,
                                      vars_t<data_t> &vars)
{
    vars.enum_mapping(
        [&](const int &ivar, const data_t &var)
        {
            // NOLINTNEXTLINE(clang-analyzer-core.uninitialized.Assign)
            cell[ivar] = var;
        });
}

template <template <typename> class vars_t, class data_t>
AMREX_GPU_HOST_DEVICE void load_vars(const amrex::CellData<data_t> &cell,
                                     vars_t<std::remove_const_t<data_t>> &vars)
{
    vars.enum_mapping([&](const int &ivar, std::remove_const_t<data_t> &var)
                      { var = cell[ivar]; });
}

template <template <typename> class vars_t, class data_t>
AMREX_GPU_HOST_DEVICE auto load_vars(const amrex::CellData<data_t> &cell)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-member-init)
    vars_t<std::remove_const_t<data_t>> vars;
    load_vars(cell, vars);
    return vars;
}

#endif /* CELL_HPP_ */
