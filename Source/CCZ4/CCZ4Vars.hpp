/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4VARS_HPP_
#define CCZ4VARS_HPP_

#include "StateVariables.hpp"
#include "Tensor.hpp"
#include "AMReX_Array4.H"

class CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    CCZ4Vars(const amrex::CellData<const amrex::Real> &input_cell_data)
        : cell_data(input_cell_data)
    {
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &chi() const
    {
        return cell_data[c_chi];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &h(int i, int j) const
    {
        return cell_data[sym_var_idx(c_h11, i, j)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &K() const
    {
        return cell_data[c_K];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &A(int i, int j) const
    {
        return cell_data[sym_var_idx(c_A11, i, j)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Theta() const
    {
        return cell_data[c_Theta];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Gamma(int i) const
    {
        return cell_data[c_Gamma1 + i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &lapse() const
    {
        return cell_data[c_lapse];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &shift(int i) const
    {
        return cell_data[c_shift1 + i];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &B(int i) const
    {
        return cell_data[c_B1 + i];
    }

    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const amrex::CellData<const amrex::Real> &cell_data;
};

#endif /* CCZ4VARS_HPP */
