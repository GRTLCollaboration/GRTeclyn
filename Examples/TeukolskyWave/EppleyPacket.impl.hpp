/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#if !defined(EPPLEYPACKET_HPP_)
#error "This file should only be included through EppleyPacket.hpp"
#endif

#ifndef EPPLEYPACKET_IMPL_HPP_
#define EPPLEYPACKET_IMPL_HPP_

#include "EppleyPacket.hpp"
#include <AMReX_REAL.H>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// F function and its derivatives where x = r \pm t
EppleyPacketDerivs EppleyPacket::get_F_derivs(amrex::Real x) const
{
    amrex::Real A = this->m_params.amplitude, sigma = this->m_params.sigma,
                r0 = this->m_params.radial_offset;

    // --- paste generated temporaries here ---
    amrex::Real exp_plus  = exp(-pow(x + r0, 2) / pow(sigma, 2));
    amrex::Real exp_minus = exp(-pow(x - r0, 2) / pow(sigma, 2));
    amrex::Real sigma2    = sigma * sigma;
    amrex::Real sigma4    = sigma2 * sigma2;
    amrex::Real sigma6    = sigma4 * sigma2;
    amrex::Real sigma8    = sigma4 * sigma4;

    // --- paste generated F0..F4 assignments here ---
    amrex::Real F0, F1, F2, F3, F4;

    F0 = (exp_plus + exp_minus) * x;

    F1 = exp_minus + exp_plus + (2 * exp_minus * (r0 - x) * x) / sigma2 -
         (2 * exp_plus * x * (r0 + x)) / sigma2;

    F2 = (4 * (exp_minus * (r0 - x) - exp_plus * (r0 + x))) / sigma2 -
         (2 * x *
          (exp_minus * sigma2 + exp_plus * sigma2 -
           2 * exp_minus * pow(r0 - x, 2) - 2 * exp_plus * pow(r0 + x, 2))) /
             sigma4;

    F3 = (-2 * (3 * sigma2 *
                    (exp_minus * sigma2 + exp_plus * sigma2 -
                     2 * exp_minus * pow(r0 - x, 2) -
                     2 * exp_plus * pow(r0 + x, 2)) +
                2 * x *
                    (3 * exp_minus * sigma2 * (r0 - x) -
                     2 * exp_minus * pow(r0 - x, 3) -
                     3 * exp_plus * sigma2 * (r0 + x) +
                     2 * exp_plus * pow(r0 + x, 3)))) /
         sigma6;

    F4 = (4 * (-4 * sigma2 *
                   (3 * exp_minus * sigma2 * (r0 - x) -
                    2 * exp_minus * pow(r0 - x, 3) -
                    3 * exp_plus * sigma2 * (r0 + x) +
                    2 * exp_plus * pow(r0 + x, 3)) +
               x * (3 * exp_minus * sigma4 + 3 * exp_plus * sigma4 -
                    12 * exp_minus * sigma2 * pow(r0 - x, 2) +
                    4 * exp_minus * pow(r0 - x, 4) -
                    12 * exp_plus * sigma2 * pow(r0 + x, 2) +
                    4 * exp_plus * pow(r0 + x, 4)))) /
         sigma8;

    return EppleyPacketDerivs{A * F0 / 2., A * F1 / 2., A * F2 / 2.,
                              A * F3 / 2., A * F4 / 2.};
}

// Auxiliary functions. In the end we want the superposition, so these are also
// implemented as get_X_tot
EvenEppleyPacketCoefficients EvenEppleyPacket::get_ABC(amrex::Real r) const
{
    amrex::Real x_out               = -r;
    amrex::Real x_in                = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in  = get_F_derivs(x_in);

    // Compute out coefficients
    amrex::Real A_out = 3 * F_derivs_out.F2 / pow(r, 3) +
                        9. * F_derivs_out.F1 / pow(r, 4) +
                        9. * F_derivs_out.F0 / pow(r, 5);
    amrex::Real B_out =
        -1. * F_derivs_out.F3 / (r * r) - 3. * F_derivs_out.F2 / pow(r, 3) -
        6. * F_derivs_out.F1 / pow(r, 4) - 6. * F_derivs_out.F0 / pow(r, 5);
    amrex::Real C_out =
        0.25 * F_derivs_out.F4 / r + 0.5 * F_derivs_out.F3 / (r * r) +
        2.25 * F_derivs_out.F2 / pow(r, 3) +
        5.25 * F_derivs_out.F1 / pow(r, 4) + 5.25 * F_derivs_out.F0 / pow(r, 5);

    // Compute in coefficients
    amrex::Real A_in = 3 * F_derivs_in.F2 / pow(r, 3) -
                       9. * F_derivs_in.F1 / pow(r, 4) +
                       9. * F_derivs_in.F0 / pow(r, 5);
    amrex::Real B_in =
        1. * F_derivs_in.F3 / (r * r) - 3. * F_derivs_in.F2 / pow(r, 3) +
        6. * F_derivs_in.F1 / pow(r, 4) - 6. * F_derivs_in.F0 / pow(r, 5);
    amrex::Real C_in =
        0.25 * F_derivs_in.F4 / r - 0.5 * F_derivs_in.F3 / (r * r) +
        2.25 * F_derivs_in.F2 / pow(r, 3) - 5.25 * F_derivs_in.F1 / pow(r, 4) +
        5.25 * F_derivs_in.F0 / pow(r, 5);

    return EvenEppleyPacketCoefficients{A_out - A_in, B_out - B_in,
                                        C_out - C_in};
}

OddEppleyPacketCoefficients OddEppleyPacket::get_KL(amrex::Real r) const
{
    amrex::Real x_out               = -r;
    amrex::Real x_in                = r;
    EppleyPacketDerivs F_derivs_out = get_F_derivs(x_out);
    EppleyPacketDerivs F_derivs_in  = get_F_derivs(x_in);

    // Compute out coefficients
    amrex::Real K_out = F_derivs_out.F2 / pow(r, 2) +
                        3. * F_derivs_out.F1 / pow(r, 3) +
                        3. * F_derivs_out.F0 / pow(r, 4);
    amrex::Real L_out = F_derivs_out.F3 / r + 2. * F_derivs_out.F2 / pow(r, 2) +
                        3. * F_derivs_out.F1 / pow(r, 3) +
                        3. * F_derivs_out.F0 / pow(r, 4);

    // Compute in coefficients
    amrex::Real K_in = F_derivs_in.F2 / pow(r, 2) -
                       3. * F_derivs_in.F1 / pow(r, 3) +
                       3. * F_derivs_in.F0 / pow(r, 4);
    amrex::Real L_in =
        -1. * F_derivs_in.F3 / r + 2. * F_derivs_in.F2 / pow(r, 2) -
        3. * F_derivs_in.F1 / pow(r, 3) + 3. * F_derivs_in.F0 / pow(r, 4);

    return OddEppleyPacketCoefficients{K_out - K_in, L_out - L_in};
}

// ------------- m = 0 EppleyPacket -----------------

EppleyPacketMetricComponents
EvenEppleyPacketM0::get_metric_components(amrex::Real x, amrex::Real y,
                                          amrex::Real z) const
{
    amrex::Real r = sqrt(x * x + y * y + z * z);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r * r / (m_params.regularize_r * m_params.regularize_r));
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    amrex::Real A_tot                   = coeffs.A;
    amrex::Real B_tot                   = coeffs.B;
    amrex::Real C_tot                   = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx =
        1. +
        (-1. + 3. * y * y / (r * r) + 3. * x * x * z * z / pow(r, 4)) * A_tot -
        6. * z * z * x * x * B_tot / pow(r, 4) +
        3. * (-y * y / (r * r) + x * x * z * z / pow(r, 4)) * C_tot;
    components.gxy = 3. * x * y *
                     (-1. * A_tot * (x * x + y * y) - 2 * z * z * B_tot +
                      (r * r + z * z) * C_tot) /
                     pow(r, 4);
    components.gxz = 3. * x * z *
                     (z * z * A_tot + (x * x + y * y - z * z) * B_tot -
                      (x * x + y * y) * C_tot) /
                     pow(r, 4);
    components.gyy =
        1. +
        (-1. + 3. * x * x / (r * r) + 3. * y * y * z * z / pow(r, 4)) * A_tot -
        6. * z * z * y * y * B_tot / pow(r, 4) +
        3. * (-x * x / (r * r) + y * y * z * z / pow(r, 4)) * C_tot;
    components.gyz = 3. * y * z *
                     (z * z * A_tot + (x * x + y * y - z * z) * B_tot -
                      (x * x + y * y) * C_tot) /
                     pow(r, 4);
    components.gzz = 1. + (-1. + 3. * pow(z, 4) / pow(r, 4)) * A_tot +
                     6. * z * z * (x * x + y * y) * B_tot / pow(r, 4) +
                     3. * (x * x + y * y) * (x * x + y * y) * C_tot / pow(r, 4);
    return components;
}

// ------------- m = 2 EppleyPacket -----------------

EppleyPacketMetricComponents
EvenEppleyPacketM2::get_metric_components(amrex::Real x, amrex::Real y,
                                          amrex::Real z) const
{
    amrex::Real r = sqrt(x * x + y * y + z * z);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r * r / (m_params.regularize_r * m_params.regularize_r));
    EvenEppleyPacketCoefficients coeffs = get_ABC(r);
    amrex::Real A_tot                   = coeffs.A;
    amrex::Real B_tot                   = coeffs.B;
    amrex::Real C_tot                   = coeffs.C;
    EppleyPacketMetricComponents components;
    components.gxx =
        1. +
        ((x * x - z * z) / (r * r) - x * x * (z * z + 2 * y * y) / pow(r, 4)) *
            A_tot +
        2 * x * x * (z * z + 2 * y * y) * B_tot / pow(r, 4) +
        ((y * y + 2 * z * z) / (r * r) -
         x * x * (z * z + 2 * y * y) / pow(r, 4)) *
            C_tot;
    components.gxy =
        x * y * (x * x - y * y) * (A_tot - 2 * B_tot + C_tot) / pow(r, 4);
    components.gxz =
        x * z *
        ((2 * x * x + z * z) * A_tot + (z * z + 3 * y * y - x * x) * B_tot -
         (x * x + 2 * z * z + 3 * y * y) * C_tot) /
        pow(r, 4);
    components.gyy =
        1. +
        ((z * z - y * y) / (r * r) + y * y * (z * z + 2 * x * x) / pow(r, 4)) *
            A_tot -
        2 * y * y * (z * z + 2 * x * x) * B_tot / pow(r, 4) +
        (-(x * x + 2 * z * z) / (r * r) +
         y * y * (z * z + 2 * x * x) / pow(r, 4)) *
            C_tot;
    components.gyz =
        y * z *
        (-(2 * y * y + z * z) * A_tot - (z * z + 3 * x * x - y * y) * B_tot +
         (y * y + 2 * z * z + 3 * x * x) * C_tot) /
        pow(r, 4);
    components.gzz = 1. + ((pow(y, 4) - pow(x, 4)) * A_tot -
                           2. * z * z * (x * x - y * y) * B_tot +
                           (x * x - y * y) * (r * r + z * z) * C_tot) /
                              pow(r, 4);
    return components;
}

// -------------- m = 2 Odd parity EppleyPacket -----------------

EppleyPacketMetricComponents
OddEppleyPacketM2::get_metric_components(amrex::Real x, amrex::Real y,
                                         amrex::Real z) const
{
    amrex::Real r = sqrt(x * x + y * y + z * z);
    // regularize at the origin
    r = r + m_params.regularize_r *
                exp(-r * r / (m_params.regularize_r * m_params.regularize_r));
    OddEppleyPacketCoefficients coeffs = get_KL(r);
    amrex::Real K_tot                  = coeffs.K;
    amrex::Real L_tot                  = coeffs.L;
    amrex::Real r3_inv                 = 1. / (r * r * r);

    EppleyPacketMetricComponents components;
    components.gxx = 1. + 8. * x * x * z * r3_inv * K_tot -
                     2. * (y * y + z * z) * z * r3_inv * L_tot;
    components.gxy = 0.;
    components.gxz = 4. * (r * r - 2. * x * x) * x * r3_inv * K_tot +
                     2. * (y * y + z * z) * x * r3_inv * L_tot;
    components.gyy = 1. - 8. * y * y * z * r3_inv * K_tot +
                     2. * (x * x + z * z) * z * r3_inv * L_tot;
    components.gyz = -4. * (r * r - 2. * y * y) * y * r3_inv * K_tot -
                     2. * (x * x + z * z) * y * r3_inv * L_tot;
    components.gzz =
        1. - 2. * (x * x - y * y) * z * r3_inv * (4. * K_tot + L_tot);
    return components;
}

#endif /* EPPLEYPACKET_IMPL_HPP_ */
