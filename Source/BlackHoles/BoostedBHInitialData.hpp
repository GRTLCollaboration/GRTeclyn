/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BOOSTEDBHINITIALDATA_HPP_
#define BOOSTEDBHINITIALDATA_HPP_
/**
 * BOOSTED SCHWARZSCHILD BLACK HOLE
 * Baumgarte & Shapiro, pp. 73-74
 * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
 */

#include "Coordinates.hpp"
#include "Tensor.hpp"
#include <array>

class BoostedBHInitialData
{

  public:
    struct params_t
    {
        double mass;
        std::array<double, AMREX_SPACEDIM> center;
        std::array<double, AMREX_SPACEDIM> momentum;
    };

    BoostedBHInitialData(params_t a_params);

    // conformal factor
    AMREX_GPU_DEVICE amrex::Real psi_minus_one(Coordinates a_coords) const;

    // extrinsic curvature
    AMREX_GPU_DEVICE Tensor<2, amrex::Real> Aij(Coordinates a_coords) const;

  private:
    params_t m_params;

    AMREX_GPU_DEVICE amrex::Real center_dist(Coordinates a_coords) const;

    AMREX_GPU_DEVICE amrex::Real psi0(amrex::Real a_r) const;

    AMREX_GPU_DEVICE amrex::Real psi2(amrex::Real a_r,
                                      amrex::Real a_cos_theta) const;

    AMREX_GPU_DEVICE amrex::Real psi2_0(amrex::Real a_r) const;

    AMREX_GPU_DEVICE amrex::Real psi2_2(amrex::Real a_r) const;
};

#endif /*BOOSTEDBHINITIALDATA_HPP_*/
