#ifndef GRTRESNA_INDEPENDENT_SCALARS_VARS_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_VARS_HPP_

#include "CCZ4Vars.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "Tensor.hpp"

class GRTresnaIndependentScalarsVars : public CCZ4Vars
{
  public:
    AMREX_GPU_DEVICE
    GRTresnaIndependentScalarsVars(
        const amrex::CellData<const amrex::Real> &input_cell_data)
        : CCZ4Vars(input_cell_data)
    {
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &phi(int k) const
    {
        return cell_data[c_phi_lump_index(k)];
    }

    [[nodiscard]]
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE const amrex::Real &Pi(int k) const
    {
        return cell_data[c_Pi_lump_index(k)];
    }
};

#endif /* GRTRESNA_INDEPENDENT_SCALARS_VARS_HPP_ */
