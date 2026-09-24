/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(EPPLEYPACKET_HPP_)
#error "This file should only be included through EppleyPacket.hpp"
#endif

#ifndef EPPLEYPACKET_IMPL_HPP_
#define EPPLEYPACKET_IMPL_HPP_

#include "EppleyPacket.hpp"
#include <AMReX_REAL.H>

// Get F and its derivatives where x = r \pm t
AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketDerivs
EppleyPacket::get_F_derivs(amrex::Real x) const
{
    amrex::Real A = this->m_params.amplitude, sigma = this->m_params.sigma,
                r0 = this->m_params.radial_offset;

    // --- temporary variables ---
    amrex::Real sigma2    = sigma * sigma;
    amrex::Real xp        = r0 + x;
    amrex::Real xm        = r0 - x;
    amrex::Real exp_plus  = exp(-(xp * xp) / sigma2);
    amrex::Real exp_minus = exp(-(xm * xm) / sigma2);
    amrex::Real sigma4    = sigma2 * sigma2;
    amrex::Real sigma6    = sigma4 * sigma2;
    amrex::Real sigma8    = sigma4 * sigma4;

    // --- F and its derivatives ---
    amrex::Real F0, F1, F2, F3, F4;

    F0 = (exp_plus + exp_minus) * x;

    F1 = exp_minus + exp_plus + (2 * exp_minus * xm * x) / sigma2 -
         (2 * exp_plus * x * xp) / sigma2;

    F2 = (4 * (exp_minus * xm - exp_plus * xp)) / sigma2 -
         (2 * x *
          (exp_minus * sigma2 + exp_plus * sigma2 - 2 * exp_minus * xm * xm -
           2 * exp_plus * xp * xp)) /
             sigma4;

    F3 = (-2 * (3 * sigma2 *
                    (exp_minus * sigma2 + exp_plus * sigma2 -
                     2 * exp_minus * xm * xm - 2 * exp_plus * xp * xp) +
                2 * x *
                    (3 * exp_minus * sigma2 * xm - 2 * exp_minus * pow(xm, 3) -
                     3 * exp_plus * sigma2 * xp + 2 * exp_plus * pow(xp, 3)))) /
         sigma6;

    F4 =
        (4 *
         (-4 * sigma2 *
              (3 * exp_minus * sigma2 * xm - 2 * exp_minus * pow(xm, 3) -
               3 * exp_plus * sigma2 * xp + 2 * exp_plus * pow(xp, 3)) +
          x * (3 * exp_minus * sigma4 + 3 * exp_plus * sigma4 -
               12 * exp_minus * sigma2 * xm * xm + 4 * exp_minus * pow(xm, 4) -
               12 * exp_plus * sigma2 * xp * xp + 4 * exp_plus * pow(xp, 4)))) /
        sigma8;

    // Add amplitude and factor 1/2 here and return the result
    return EppleyPacketDerivs{A * F0 / 2., A * F1 / 2., A * F2 / 2.,
                              A * F3 / 2., A * F4 / 2.};
}

// Auxiliary functions. In the end we want the superposition, so the returned
// coefficients are the difference between the out- and ingoing coefficients.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE EvenEppleyPacketCoefficients
EppleyPacket::get_ABC(amrex::Real r) const
{
    // t = 0
    amrex::Real x_out               = -r;
    amrex::Real x_in                = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in  = get_F_derivs(x_in);

    // Compute inverse powers of r
    amrex::Real r2_inv = 1.0 / (r * r);
    amrex::Real r3_inv = r2_inv / r;
    amrex::Real r4_inv = r3_inv / r;
    amrex::Real r5_inv = r4_inv / r;

    // Compute out coefficients
    amrex::Real A_out = 3 * F_derivs_out.F2 * r3_inv +
                        9. * F_derivs_out.F1 * r4_inv +
                        9. * F_derivs_out.F0 * r5_inv;
    amrex::Real B_out =
        -1. * F_derivs_out.F3 * r2_inv - 3. * F_derivs_out.F2 * r3_inv -
        6. * F_derivs_out.F1 * r4_inv - 6. * F_derivs_out.F0 * r5_inv;
    amrex::Real C_out =
        0.25 * F_derivs_out.F4 / r + 0.5 * F_derivs_out.F3 * r2_inv +
        2.25 * F_derivs_out.F2 * r3_inv + 5.25 * F_derivs_out.F1 * r4_inv +
        5.25 * F_derivs_out.F0 * r5_inv;

    // Compute in coefficients
    amrex::Real A_in = 3 * F_derivs_in.F2 * r3_inv -
                       9. * F_derivs_in.F1 * r4_inv +
                       9. * F_derivs_in.F0 * r5_inv;
    amrex::Real B_in =
        1. * F_derivs_in.F3 * r2_inv - 3. * F_derivs_in.F2 * r3_inv +
        6. * F_derivs_in.F1 * r4_inv - 6. * F_derivs_in.F0 * r5_inv;
    amrex::Real C_in =
        0.25 * F_derivs_in.F4 / r - 0.5 * F_derivs_in.F3 * r2_inv +
        2.25 * F_derivs_in.F2 * r3_inv - 5.25 * F_derivs_in.F1 * r4_inv +
        5.25 * F_derivs_in.F0 * r5_inv;

    return EvenEppleyPacketCoefficients{A_out - A_in, B_out - B_in,
                                        C_out - C_in};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE OddEppleyPacketCoefficients
EppleyPacket::get_KL(amrex::Real r) const
{
    // t = 0
    amrex::Real x_out               = -r;
    amrex::Real x_in                = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in  = get_F_derivs(x_in);

    // Compute inverse powers of r
    amrex::Real r2_inv = 1.0 / (r * r);
    amrex::Real r3_inv = r2_inv / r;
    amrex::Real r4_inv = r3_inv / r;

    // Compute out coefficients
    amrex::Real K_out = F_derivs_out.F2 * r2_inv +
                        3. * F_derivs_out.F1 * r3_inv +
                        3. * F_derivs_out.F0 * r4_inv;
    amrex::Real L_out = F_derivs_out.F3 / r + 2. * F_derivs_out.F2 * r2_inv +
                        3. * F_derivs_out.F1 * r3_inv +
                        3. * F_derivs_out.F0 * r4_inv;

    // Compute in coefficients
    amrex::Real K_in = F_derivs_in.F2 * r2_inv - 3. * F_derivs_in.F1 * r3_inv +
                       3. * F_derivs_in.F0 * r4_inv;
    amrex::Real L_in = -1. * F_derivs_in.F3 / r + 2. * F_derivs_in.F2 * r2_inv -
                       3. * F_derivs_in.F1 * r3_inv +
                       3. * F_derivs_in.F0 * r4_inv;

    return OddEppleyPacketCoefficients{K_out - K_in, L_out - L_in};
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketMetricComponents
EppleyPacket::get_metric_components(amrex::Real x, amrex::Real y,
                                    amrex::Real z) const
{
    if (m_type == EppleyPacketType::even_m0)
    {
        return get_metric_components_even_m0(x, y, z);
    }
    if (m_type == EppleyPacketType::even_m2)
    {
        return get_metric_components_even_m2(x, y, z);
    }
    // m_type == EppleyPacketType::odd_m2
    return get_metric_components_odd_m2(x, y, z);
}

// ------------- m = 0 Even parity EppleyPacket -----------------

AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketMetricComponents
EppleyPacket::get_metric_components_even_m0(amrex::Real x, amrex::Real y,
                                            amrex::Real z) const
{
    amrex::Real x2 = x * x, y2 = y * y, z2 = z * z;
    amrex::Real r2 = x2 + y2 + z2;
    amrex::Real r  = sqrt(r2);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r2 / (m_params.regularize_r * m_params.regularize_r));
    // Propagate the regularization to r2
    r2                 = r * r;
    amrex::Real r2_inv = 1. / r2;
    amrex::Real xy2    = x2 + y2;
    amrex::Real r4_inv = r2_inv * r2_inv;
    // Coefficients for the even parity case
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    amrex::Real A_tot                   = coeffs.A;
    amrex::Real B_tot                   = coeffs.B;
    amrex::Real C_tot                   = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx = 1. +
                     (-1. + 3. * y2 * r2_inv + 3. * x2 * z2 * r4_inv) * A_tot -
                     6. * z2 * x2 * B_tot * r4_inv +
                     3. * (-y2 * r2_inv + x2 * z2 * r4_inv) * C_tot;
    components.gxy = 3. * x * y *
                     (-1. * A_tot * xy2 - 2 * z2 * B_tot + (r2 + z2) * C_tot) *
                     r4_inv;
    components.gxz =
        3. * x * z * (z2 * A_tot + (xy2 - z2) * B_tot - xy2 * C_tot) * r4_inv;
    components.gyy = 1. +
                     (-1. + 3. * x2 * r2_inv + 3. * y2 * z2 * r4_inv) * A_tot -
                     6. * z2 * y2 * B_tot * r4_inv +
                     3. * (-x2 * r2_inv + y2 * z2 * r4_inv) * C_tot;
    components.gyz =
        3. * y * z * (z2 * A_tot + (xy2 - z2) * B_tot - xy2 * C_tot) * r4_inv;
    components.gzz = 1. + (-1. + 3. * z2 * z2 * r4_inv) * A_tot +
                     6. * z2 * xy2 * B_tot * r4_inv +
                     3. * xy2 * xy2 * C_tot * r4_inv;
    return components;
}

// ------------- m = 2 Even parity EppleyPacket -----------------

AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketMetricComponents
EppleyPacket::get_metric_components_even_m2(amrex::Real x, amrex::Real y,
                                            amrex::Real z) const
{
    amrex::Real x2 = x * x, y2 = y * y, z2 = z * z;
    amrex::Real r2 = x2 + y2 + z2;
    amrex::Real r  = sqrt(r2);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r2 / (m_params.regularize_r * m_params.regularize_r));
    // Propagate the regularization to r2
    r2                                  = r * r;
    amrex::Real r2_inv                  = 1. / r2;
    amrex::Real r4_inv                  = r2_inv * r2_inv;
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    amrex::Real A_tot                   = coeffs.A;
    amrex::Real B_tot                   = coeffs.B;
    amrex::Real C_tot                   = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx =
        1. + ((x2 - z2) * r2_inv - x2 * (z2 + 2 * y2) * r4_inv) * A_tot +
        2 * x2 * (z2 + 2 * y2) * B_tot * r4_inv +
        ((y2 + 2 * z2) * r2_inv - x2 * (z2 + 2 * y2) * r4_inv) * C_tot;
    components.gxy = x * y * (x2 - y2) * (A_tot - 2 * B_tot + C_tot) * r4_inv;
    components.gxz = x * z *
                     ((2 * x2 + z2) * A_tot + (z2 + 3 * y2 - x2) * B_tot -
                      (x2 + 2 * z2 + 3 * y2) * C_tot) *
                     r4_inv;
    components.gyy =
        1. + ((z2 - y2) * r2_inv + y2 * (z2 + 2 * x2) * r4_inv) * A_tot -
        2 * y2 * (z2 + 2 * x2) * B_tot * r4_inv +
        (-(x2 + 2 * z2) * r2_inv + y2 * (z2 + 2 * x2) * r4_inv) * C_tot;
    components.gyz = y * z *
                     (-(2 * y2 + z2) * A_tot - (z2 + 3 * x2 - y2) * B_tot +
                      (y2 + 2 * z2 + 3 * x2) * C_tot) *
                     r4_inv;
    components.gzz =
        1. + ((y2 * y2 - x2 * x2) * A_tot - 2. * z2 * (x2 - y2) * B_tot +
              (x2 - y2) * (r2 + z2) * C_tot) *
                 r4_inv;
    return components;
}

// -------------- m = 2 Odd parity EppleyPacket -----------------

AMREX_GPU_DEVICE AMREX_FORCE_INLINE EppleyPacketMetricComponents
EppleyPacket::get_metric_components_odd_m2(amrex::Real x, amrex::Real y,
                                           amrex::Real z) const
{
    amrex::Real x2 = x * x, y2 = y * y, z2 = z * z;
    amrex::Real r2 = x2 + y2 + z2;
    amrex::Real r  = sqrt(r2);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r2 / (m_params.regularize_r * m_params.regularize_r));
    // Propagate the regularization to r2
    r2                                 = r * r;
    OddEppleyPacketCoefficients coeffs = get_KL(r);
    amrex::Real K_tot                  = coeffs.K;
    amrex::Real L_tot                  = coeffs.L;
    amrex::Real r2_inv                 = 1. / r2;
    amrex::Real r3_inv                 = r2_inv / r;

    EppleyPacketMetricComponents components;
    components.gxx =
        1. + 8. * x2 * z * r3_inv * K_tot - 2. * (y2 + z2) * z * r3_inv * L_tot;
    components.gxy = 0.;
    components.gxz = 4. * (r2 - 2 * x2) * x * r3_inv * K_tot +
                     2. * (y2 + z2) * x * r3_inv * L_tot;
    components.gyy =
        1. - 8. * y2 * z * r3_inv * K_tot + 2. * (x2 + z2) * z * r3_inv * L_tot;
    components.gyz = -4. * (r2 - 2 * y2) * y * r3_inv * K_tot -
                     2. * (x2 + z2) * y * r3_inv * L_tot;
    components.gzz = 1. - 2. * (x2 - y2) * z * r3_inv * (4. * K_tot + L_tot);
    return components;
}

#endif /* EPPLEYPACKET_IMPL_HPP_ */
