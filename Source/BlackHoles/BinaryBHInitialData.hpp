/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BINARYBHINITIALDATA_HPP_
#define BINARYBHINITIALDATA_HPP_

#include "BoostedBHInitialData.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "Tensor.hpp"
#include <array>

#include <AMReX_Geometry.H>

enum Lapse
{
    ONE,
    PRE_COLLAPSED,
    CHI
};

class BinaryBHInitialData
{
  protected:
    amrex::Geometry m_geom;
    BoostedBHInitialData bh1;
    BoostedBHInitialData bh2;
    int m_initial_lapse;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    AMREX_FORCE_INLINE
    BinaryBHInitialData(BoostedBHInitialData::params_t a_bh1_params,
                        BoostedBHInitialData::params_t a_bh2_params,
                        const amrex::Geometry &a_geom,
                        int a_initial_lapse = Lapse::PRE_COLLAPSED);
    // NOLINTEND(bugprone-easily-swappable-parameters)

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE void
    init_data(int i, int j, int k,
              const amrex::CellData<amrex::Real> &cell) const;

  protected:
    [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
    compute_chi(Coordinates coords) const;

    [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE Tensor<2, amrex::Real>
    compute_A(amrex::Real chi, Coordinates coords) const;
};

#include "BinaryBHInitialData.impl.hpp"

#endif /* BINARYBHINITIALDATA_HPP_ */
