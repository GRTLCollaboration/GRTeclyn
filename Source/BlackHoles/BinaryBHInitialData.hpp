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
    AMREX_FORCE_INLINE
    BinaryBHInitialData(BoostedBHInitialData::params_t a_bh1_params,
                        BoostedBHInitialData::params_t a_bh2_params,
                        double a_dx,
                        int a_initial_lapse = Lapse::PRE_COLLAPSED);
    // NOLINTEND(bugprone-easily-swappable-parameters)

    AMREX_GPU_DEVICE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &state) const;

  protected:
    [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE amrex::Real
    compute_chi(Coordinates coords) const;

    [[nodiscard]] AMREX_FORCE_INLINE AMREX_GPU_DEVICE Tensor::Rank2
    compute_A(amrex::Real chi, Coordinates coords) const;
};

#include "BinaryBHInitialData.impl.hpp"

#endif /* BINARYBHINITIALDATA_HPP_ */
