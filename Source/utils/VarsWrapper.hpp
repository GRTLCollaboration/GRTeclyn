/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

 #ifndef VARSWRAPPER_HPP_
 #define VARSWRAPPER_HPP_

template <typename data_t>
class VarsWrapper
{
    public:
        AMREX_GPU_DEVICE inline VarsWrapper(const amrex::CellData<data_t const> &input_cell_data)
            : cell_data(input_cell_data)
        {
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &chi() { return cell_data[c_chi]; }
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &h(int i, int j) 
        {
            return cell_data[SYMM_INDEX(c_h11, i, j)];
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &K() { return cell_data[c_K]; }
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &A(int i, int j) 
        {
            return cell_data[SYMM_INDEX(c_A11, i, j)];
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &Theta() { return cell_data[c_Theta]; }
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &Gamma(int i) 
        {
            return cell_data[c_Gamma1 +i];
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &lapse() { return cell_data[c_lapse]; }
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &shift(int i) { return cell_data[c_shift1 + i]; }
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE const data_t &B(int i) { return cell_data[c_B1 + i]; }

        const amrex::CellData<data_t const> &cell_data;
};

#endif /* VARSWRAPPER_HPP */