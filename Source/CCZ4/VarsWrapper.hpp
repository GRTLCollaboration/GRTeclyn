/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef VARSWRAPPER_HPP_
#define VARSWRAPPER_HPP_

// A function to return the right index for the tensor
template <std::integral T> constexpr int var_idx(T ivar, T i, T j)
{
    return ivar + i + j + ((i * j != 0) ? 1 : 0);
}

class VarsWrapper
{
  public:
    AMREX_GPU_DEVICE inline VarsWrapper(
        const amrex::CellData<amrex::Real const> &input_cell_data)
        : cell_data(input_cell_data)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &chi()
    {
        return cell_data[c_chi];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &h(int i, int j)
    {
        return cell_data[var_idx(c_h11, i, j)];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &K()
    {
        return cell_data[c_K];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &A(int i, int j)
    {
        return cell_data[var_idx(c_A11, i, j)];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Theta()
    {
        return cell_data[c_Theta];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Gamma(int i)
    {
        return cell_data[c_Gamma1 + i];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &lapse()
    {
        return cell_data[c_lapse];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &shift(int i)
    {
        return cell_data[c_shift1 + i];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &B(int i)
    {
        return cell_data[c_B1 + i];
    }

    const amrex::CellData<amrex::Real const> &cell_data;
};

#endif /* VARSWRAPPER_HPP */