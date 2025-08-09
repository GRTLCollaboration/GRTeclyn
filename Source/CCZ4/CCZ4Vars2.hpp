/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4VARS2_HPP_
#define CCZ4VARS2_HPP_

#include "StateVariables.hpp"
#include "Tensor.hpp"

class CCZ4Vars2
{
  public:
    AMREX_GPU_DEVICE inline CCZ4Vars2(
        const amrex::CellData<amrex::Real> &input_cell_data)
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
        return cell_data[var_idx(c_h11, i, j)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &K() const
    {
        return cell_data[c_K];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &A(int i, int j) const
    {
        return cell_data[var_idx(c_A11, i, j)];
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

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_var(amrex::Real var,
                                                       int ivar)
    {
        cell_data[ivar] = var;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_chi(amrex::Real chi)
    {
        cell_data[c_chi] = chi;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    store_h(Tensor<2, amrex::Real> h_LL)
    {
        FOR2 (i, j)
        {
            cell_data[var_idx(c_h11, i, j)] = h_LL[i][j];
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_K(amrex::Real K)
    {
        cell_data[c_K] = K;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    store_A(Tensor<2, amrex::Real> A_LL)
    {
        FOR2 (i, j)
        {
            cell_data[var_idx(c_A11, i, j)] = A_LL[i][j];
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_Theta(amrex::Real Theta)
    {
        cell_data[c_Theta] = Theta;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    store_Gamma(Tensor<1, amrex::Real> Gamma_U)
    {
        FOR (i)
        {
            cell_data[c_Gamma1 + i] = Gamma_U[i];
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_lapse(amrex::Real lapse)
    {
        cell_data[c_lapse] = lapse;
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    store_shift(Tensor<1, amrex::Real> shift_U)
    {
        FOR (i)
        {
            cell_data[c_shift1 + i] = shift_U[i];
        }
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void store_B(Tensor<1, amrex::Real> B_U)
    {
        FOR (i)
        {
            cell_data[c_B1 + i] = B_U[i];
        }
    }

    const amrex::CellData<amrex::Real> &cell_data;
};

#endif /* CCZ4VARS2_HPP */