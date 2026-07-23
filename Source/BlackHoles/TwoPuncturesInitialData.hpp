/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifdef USE_TWOPUNCTURES

#ifndef TWOPUNCTURESINITIALDATA_HPP_
#define TWOPUNCTURESINITIALDATA_HPP_

#include "CCZ4Vars.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total number of components
#include "TensorAlgebra.hpp"
#include "TwoPunctures.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"
#include <array>

//! This compute class sets the initial data computed by TwoPunctures on the
//! grid
class TwoPuncturesInitialData
{
  protected:
    amrex::Real m_dx;
    std::array<amrex::Real, AMREX_SPACEDIM> m_center;
    const TP::TwoPunctures &m_two_punctures;

  public:
    template <class data_t> using Vars = CCZ4Vars::VarsWithGauge<data_t>;

    TwoPuncturesInitialData(
        const amrex::Real a_dx,
        const std::array<amrex::Real, AMREX_SPACEDIM> a_center,
        const TP::TwoPunctures &a_two_punctures)
        : m_dx(a_dx), m_center(a_center), m_two_punctures(a_two_punctures)
    {
    }

    void compute(Cell<amrex::Real> current_cell) const;

  protected:
    void interpolate_tp_vars(const Coordinates<amrex::Real> &coords,
                             Tensor::Rank2 &out_h_phys,
                             Tensor::Rank2 &out_extrinsic_K,
                             amrex::Real &out_lapse, Tensor::Rank1 &out_shift,
                             amrex::Real &out_Theta,
                             Tensor::Rank1 &out_Z3) const;
};

#include "TwoPuncturesInitialData.impl.hpp"

#endif /* TWOPUNCTURESINITIALDATA_HPP_ */
#endif /* USE_TWOPUNCTURES */
