/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
// Our headers
#include "DimensionDefinitions.hpp"
#include "StateVariables.hpp"

// AMReX headers
#include "AMReX_Array4.H"
#include "AMReX_IntVect.H"
#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
random_ccz4_initial_data(const amrex::IntVect &a_iv,
                         const amrex::Array4<amrex::Real> &a_array,
                         const amrex::RealVect &a_coords)
{
    amrex::Real x = a_coords[0];
    amrex::Real y = a_coords[1];
    amrex::Real z = a_coords[2];

    // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    // NOLINTNEXTLINE(readability-identifier-length)
    amrex::Real g[3][3];
    amrex::Real g_UU[3][3];
    amrex::Real chi; // NOLINT(cppcoreguidelines-init-variables)
    {
        g[0][0] = 1.36778 + 2.39731 * x + 4.53541 * x * x +
                  19.9771 * x * y * y * y + 6.13801 * y * z + 5.65185 * z * z +
                  9.35842 * z * z * z * z;
        g[0][1] = -0.07646 - 0.48786 * x - 0.75098 * x * x -
                  1.73683 * x * y * y * y + 1.71676 * y * z + 1.03662 * z * z +
                  0.35630 * z * z * z * z;
        g[0][2] = -0.10083 + 0.12445 * x - 1.26649 * x * x -
                  1.95052 * x * y * y * y + 0.73091 * y * z - 1.49835 * z * z -
                  2.39024 * z * z * z * z;
        g[1][1] = 0.84072 + 2.31163 * x + 3.32275 * x * x +
                  15.1662 * x * y * y * y + 8.48730 * y * z + 3.05098 * z * z +
                  17.8448 * z * z * z * z;
        g[1][2] = -0.42495 - 0.33464 * x - 0.47012 * x * x -
                  7.38477 * x * y * y * y + 0.41896 * y * z - 1.36394 * z * z +
                  5.25894 * z * z * z * z;
        g[2][2] = 0.60995 + 1.30428 * x + 3.86237 * x * x +
                  22.7614 * x * y * y * y + 6.93818 * y * z + 4.39250 * z * z +
                  19.0244 * z * z * z * z;
        g[1][0] = g[0][1];
        g[2][0] = g[0][2];
        g[2][1] = g[1][2];

        const amrex::Real detg =
            g[0][0] * g[1][1] * g[2][2] + 2 * g[0][1] * g[0][2] * g[1][2] -
            g[0][0] * g[1][2] * g[1][2] - g[1][1] * g[0][2] * g[0][2] -
            g[2][2] * g[0][1] * g[0][1];
        g_UU[0][0] = (g[1][1] * g[2][2] - g[1][2] * g[1][2]) / detg;
        g_UU[0][1] = (g[0][2] * g[1][2] - g[0][1] * g[2][2]) / detg;
        g_UU[0][2] = (g[0][1] * g[1][2] - g[0][2] * g[1][1]) / detg;
        g_UU[1][1] = (g[0][0] * g[2][2] - g[0][2] * g[0][2]) / detg;
        g_UU[1][2] = (g[0][1] * g[0][2] - g[0][0] * g[1][2]) / detg;
        g_UU[2][2] = (g[0][0] * g[1][1] - g[0][1] * g[0][1]) / detg;
        g_UU[1][0] = g_UU[0][1];
        g_UU[2][0] = g_UU[0][2];
        g_UU[2][1] = g_UU[1][2];

        chi = std::pow(std::fabs(detg), -1.0 / GR_SPACEDIM);

        a_array(a_iv, c_chi) = chi;
        a_array(a_iv, c_h11) = chi * g[0][0];
        a_array(a_iv, c_h12) = chi * g[0][1];
        a_array(a_iv, c_h13) = chi * g[0][2];
        a_array(a_iv, c_h22) = chi * g[1][1];
        a_array(a_iv, c_h23) = chi * g[1][2];
        a_array(a_iv, c_h33) = chi * g[2][2];
    }

    {
        // NOLINTNEXTLINE(readability-identifier-length)
        amrex::Real K[3][3];
        K[0][0] = -0.16238 - 0.74295 * x + 0.51595 * x * x -
                  6.60239 * x * y * y * y - 0.76401 * y * z - 1.81131 * z * z -
                  3.88228 * z * z * z * z;
        K[0][1] = 0.15054 - 0.60088 * x - 0.15428 * x * x +
                  3.16779 * x * y * y * y - 2.00687 * y * z - 1.35442 * z * z -
                  0.67601 * z * z * z * z;
        K[0][2] = -0.02174 - 0.36243 * x + 0.81531 * x * x +
                  4.34918 * x * y * y * y + 0.90419 * y * z - 0.85088 * z * z -
                  6.45097 * z * z * z * z;
        K[1][1] = -0.47653 - 0.43889 * x + 0.87342 * x * x +
                  4.24684 * x * y * y * y + 0.26290 * y * z + 1.90095 * z * z +
                  3.69515 * z * z * z * z;
        K[1][2] = 0.37472 + 0.03657 * x - 0.10327 * x * x -
                  0.95744 * x * y * y * y - 1.20800 * y * z - 0.43064 * z * z -
                  0.25419 * z * z * z * z;
        K[2][2] = 0.34184 + 0.21495 * x - 0.73195 * x * x +
                  7.81626 * x * y * y * y + 2.48359 * y * z + 1.89657 * z * z -
                  4.10980 * z * z * z * z;
        K[1][0] = K[0][1];
        K[2][0] = K[0][2];
        K[2][1] = K[1][2];

        amrex::Real trK = 0;
        FOR (i, j)
        {
            trK += g_UU[i][j] * K[i][j];
        }

        a_array(a_iv, c_K)   = trK;
        a_array(a_iv, c_A11) = chi * (K[0][0] - trK * g[0][0] / GR_SPACEDIM);
        a_array(a_iv, c_A12) = chi * (K[0][1] - trK * g[0][1] / GR_SPACEDIM);
        a_array(a_iv, c_A13) = chi * (K[0][2] - trK * g[0][2] / GR_SPACEDIM);
        a_array(a_iv, c_A22) = chi * (K[1][1] - trK * g[1][1] / GR_SPACEDIM);
        a_array(a_iv, c_A23) = chi * (K[1][2] - trK * g[1][2] / GR_SPACEDIM);
        a_array(a_iv, c_A33) = chi * (K[2][2] - trK * g[2][2] / GR_SPACEDIM);
    }

    a_array(a_iv, c_Theta) = 0.27579 + 0.25791 * x + 1.40488 * x * x +
                             5.68276 * x * y * y * y + 3.04325 * y * z +
                             1.81250 * z * z + 1.01832 * z * z * z * z;
    a_array(a_iv, c_Gamma1) = -0.49482 + 0.89227 * x + 0.05571 * x * x -
                              5.38570 * x * y * y * y + 0.13979 * y * z -
                              0.68588 * z * z - 4.39964 * z * z * z * z;
    a_array(a_iv, c_Gamma2) = -0.09082 - 0.31017 * x + 1.06980 * x * x +
                              7.81524 * x * y * y * y - 1.65016 * y * z -
                              0.53352 * z * z - 3.20997 * z * z * z * z;
    a_array(a_iv, c_Gamma3) = -0.42367 + 0.03891 * x - 0.87898 * x * x +
                              6.67657 * x * y * y * y - 3.44662 * y * z -
                              0.19655 * z * z + 2.97524 * z * z * z * z;

    a_array(a_iv, c_lapse) = 0.73578 + 0.36898 * x + 0.64348 * x * x +
                             9.33487 * x * y * y * y + 0.99469 * y * z +
                             0.20515 * z * z + 8.88385 * z * z * z * z;
    a_array(a_iv, c_shift1) = 0.00000 + 0.18795 * x - 0.52389 * x * x -
                              4.14079 * x * y * y * y + 0.73135 * y * z -
                              0.27057 * z * z + 3.24187 * z * z * z * z;
    a_array(a_iv, c_shift2) = 0.00000 - 0.30316 * x - 0.15184 * x * x -
                              0.48815 * x * y * y * y + 2.45991 * y * z -
                              0.79248 * z * z + 7.14007 * z * z * z * z;
    a_array(a_iv, c_shift3) = 0.00000 + 0.68835 * x - 0.52219 * x * x -
                              7.50449 * x * y * y * y - 2.35372 * y * z -
                              0.21476 * z * z + 4.36363 * z * z * z * z;
    a_array(a_iv, c_B1) = -0.26928 + 0.35045 * x - 0.48884 * x * x +
                          2.72465 * x * y * y * y - 2.59022 * y * z -
                          0.27384 * z * z + 0.38748 * z * z * z * z;
    a_array(a_iv, c_B2) = 0.40234 + 0.26741 * x + 1.94822 * x * x -
                          0.78276 * x * y * y * y + 2.12346 * y * z +
                          0.69086 * z * z - 4.47639 * z * z * z * z;
    a_array(a_iv, c_B3) = 0.40313 + 0.00569 * x - 1.12452 * x * x -
                          5.49255 * x * y * y * y - 2.21932 * y * z +
                          0.49523 * z * z + 1.29460 * z * z * z * z;
    // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
}
