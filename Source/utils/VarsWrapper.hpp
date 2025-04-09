/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef VARSWRAPPER_HPP_
#define VARSWRAPPER_HPP_

template <typename data_t> class VarsWrapper
{
  public:
    AMREX_GPU_DEVICE inline VarsWrapper(
        const amrex::CellData<data_t const> &input_cell_data)
        : cell_data(input_cell_data)
    {
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &chi() const
    {
        return cell_data[c_chi];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &h(int i, int j) const
    {
        return cell_data[SYMM_INDEX(c_h11, i, j)];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &K() const
    {
        return cell_data[c_K];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &A(int i, int j) const
    {
        return cell_data[SYMM_INDEX(c_A11, i, j)];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &Theta() const
    {
        return cell_data[c_Theta];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &Gamma(int i) const
    {
        return cell_data[c_Gamma1 + i];
    }

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &lapse() const
    {
        return cell_data[c_lapse];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &shift(int i) const
    {
        return cell_data[c_shift1 + i];
    }
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t &B(int i) const
    {
        return cell_data[c_B1 + i];
    }

    const amrex::CellData<data_t const> &cell_data;
};

#endif /* VARSWRAPPER_HPP */