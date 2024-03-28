/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef BOOSTEDBH_HPP_
#define BOOSTEDBH_HPP_
/**
 * BOOSTED SCHWARZSCHILD BLACK HOLE
 * Baumgarte & Shapiro, pp. 73-74
 * NB: \bar{A} as defined in the book is psi^{-6} * \bar{A}_{BSSN}
 */

#include "Coordinates.hpp"
#include "Tensor.hpp"
#include <array>

class BoostedBH
{

  public:
    struct params_t
    {
        double mass;
        std::array<double, AMREX_SPACEDIM> center;
        std::array<double, AMREX_SPACEDIM> momentum;
    };

    BoostedBH(params_t a_params);

    // conformal factor
    template <class data_t>
    AMREX_GPU_DEVICE data_t psi_minus_one(Coordinates<data_t> a_coords) const;

    // extrinsic curvature
    template <class data_t>
    AMREX_GPU_DEVICE Tensor<2, data_t> Aij(Coordinates<data_t> a_coords) const;

  private:
    params_t m_params;

    template <class data_t>
    AMREX_GPU_DEVICE data_t center_dist(Coordinates<data_t> a_coords) const;

    template <class data_t> AMREX_GPU_DEVICE data_t psi0(data_t a_r) const;

    template <class data_t>
    AMREX_GPU_DEVICE data_t psi2(data_t a_r, data_t a_cos_theta) const;

    template <class data_t> AMREX_GPU_DEVICE data_t psi2_0(data_t a_r) const;

    template <class data_t> AMREX_GPU_DEVICE data_t psi2_2(data_t a_r) const;
};

#include "BoostedBH.impl.hpp"

#endif /*BOOSTEDBH_HPP_*/
