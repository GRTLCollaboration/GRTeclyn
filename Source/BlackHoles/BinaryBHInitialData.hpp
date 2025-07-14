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

enum Lapse
{
    ONE,
    PRE_COLLAPSED,
    CHI
};

class BinaryBHInitialData
{
  protected:
    double m_dx;
    BoostedBHInitialData bh1;
    BoostedBHInitialData bh2;
    int m_initial_lapse;

  public:
    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    BinaryBHInitialData(BoostedBHInitialData::params_t a_bh1_params,
                        BoostedBHInitialData::params_t a_bh2_params,
                        double a_dx,
                        int a_initial_lapse = Lapse::PRE_COLLAPSED);
    // NOLINTEND(bugprone-easily-swappable-parameters)

    AMREX_GPU_DEVICE void
    init_data(int i, int j, int k,
              const amrex::CellData<amrex::Real> &cell) const;

  protected:
    AMREX_GPU_DEVICE [[nodiscard]] amrex::Real
    compute_chi(Coordinates coords) const;

    AMREX_GPU_DEVICE [[nodiscard]] Tensor<2, amrex::Real>
    compute_A(amrex::Real chi, Coordinates coords) const;
};

#endif /* BINARYBHINITIALDATA_HPP_ */
