/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef CCZ4VARS2_HPP_
#define CCZ4VARS2_HPP_

// A function to return the right index for the tensor
[[nodiscard]]
int var_idx(int ivar, int i, int j)
{
    return ivar + i + j + ((i * j != 0) ? 1 : 0);
}

class CCZ4Vars2
{
  public:
    AMREX_GPU_DEVICE inline CCZ4Vars2(
        const amrex::CellData<amrex::Real const> &input_cell_data)
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

    const amrex::CellData<amrex::Real const> &cell_data;
};

#endif /* CCZ4VARS2_HPP */