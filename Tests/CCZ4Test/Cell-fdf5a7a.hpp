/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef CELL_FDF5A7A_HPP_
#define CELL_FDF5A7A_HPP_

#include <AMReX_Array4.H>

using namespace amrex::disabled;

// Namespace to avoid conflicts with current code
namespace Old
{
template <template <typename> class vars_t, class data_t>
AMREX_GPU_DEVICE void store_vars(amrex::CellData<data_t> const &cell,
                                 vars_t<data_t> &vars)
{
    vars.enum_mapping([&](const int &ivar, data_t const &var)
                      { cell[ivar] = var; });
}

template <template <typename> class vars_t, class data_t>
AMREX_GPU_DEVICE void load_vars(amrex::CellData<data_t> const &cell,
                                vars_t<std::remove_const_t<data_t>> &vars)
{
    vars.enum_mapping([&](const int &ivar, std::remove_const_t<data_t> &var)
                      { var = cell[ivar]; });
}

template <template <typename> class vars_t, class data_t>
AMREX_GPU_DEVICE auto load_vars(amrex::CellData<data_t> const &cell)
{
    vars_t<std::remove_const_t<data_t>> vars;
    load_vars(cell, vars);
    return vars;
}
} // namespace Old

#endif /* CELL_FDF5A7A_HPP_ */
