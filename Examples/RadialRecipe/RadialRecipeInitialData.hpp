#ifndef RADIALRECIPEINITIALDATA_HPP_
#define RADIALRECIPEINITIALDATA_HPP_

#include "CCZ4StateVariables.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "SphericalHarmonics.hpp"
#include "StateVariables.hpp"
#include <AMReX_REAL.H>

#include <array>
#include <cmath>

class RadialRecipeInitialData
{
  public:
    static constexpr int MAX_BASES = 16;
    static constexpr int MAX_ANGULAR_MODES = 8;
    static constexpr int MAX_YLM_MODES = 8;

    struct AngularMode
    {
        int ell{0};
        double amplitude{0.0};
        double radial_center{0.0};
        double radial_width{1.0};
    };

    struct YlmMode
    {
        int ell{0};
        int em{0};
        double amplitude{0.0};
        double radial_center{0.0};
        double radial_width{1.0};
    };

    struct params_t
    {
        int num_bases{4};
        std::array<double, AMREX_SPACEDIM> grid_center{};
        double basis_width{1.0};
        double basis_radius_max{8.0};

        double chi_asymptotic{1.0};
        double alpha_asymptotic{1.0};
        double beta_asymptotic{0.0};
        double K_asymptotic{0.0};
        double phi_asymptotic{0.0};
        double Pi_asymptotic{0.0};

        std::array<double, MAX_BASES> chi_coeffs{};
        std::array<double, MAX_BASES> alpha_coeffs{};
        std::array<double, MAX_BASES> beta_coeffs{};
        std::array<double, MAX_BASES> K_coeffs{};
        std::array<double, MAX_BASES> phi_coeffs{};
        std::array<double, MAX_BASES> Pi_coeffs{};

        int num_chi_angular_modes{0};
        std::array<AngularMode, MAX_ANGULAR_MODES> chi_angular_modes{};

        int num_chi_Ylm_modes{0};
        std::array<YlmMode, MAX_YLM_MODES> chi_Ylm_modes{};

        // Axisymmetric (Legendre) angular deformations of the gauge/extrinsic
        // fields, enabling non-spherical warp/wormhole candidates beyond the
        // purely radial recipe.
        int num_lapse_angular_modes{0};
        std::array<AngularMode, MAX_ANGULAR_MODES> lapse_angular_modes{};

        int num_beta_angular_modes{0};
        std::array<AngularMode, MAX_ANGULAR_MODES> beta_angular_modes{};

        int num_K_angular_modes{0};
        std::array<AngularMode, MAX_ANGULAR_MODES> K_angular_modes{};
    };

    RadialRecipeInitialData(params_t a_params, double a_dx)
        : m_params(a_params), m_dx(a_dx)
    {
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, amrex::Array4<data_t> cell) const
    {
        amrex::IntVect grid_index(i, j, k);
        Coordinates coords(grid_index, m_dx, m_params.grid_center);

        const data_t dx = coords.x;
        const data_t dy = coords.y;
        const data_t dz = coords.z;
        const data_t r = std::sqrt(dx * dx + dy * dy + dz * dz);

        data_t chi = radial_profile((data_t)m_params.chi_asymptotic,
                                    m_params.chi_coeffs, r);
        chi += chi_angular_contribution(r, dx, dy, dz);
        data_t lapse = radial_profile((data_t)m_params.alpha_asymptotic,
                                      m_params.alpha_coeffs, r);
        lapse += legendre_angular_sum(m_params.num_lapse_angular_modes,
                                      m_params.lapse_angular_modes, r, dz);
        data_t beta_x = radial_profile((data_t)m_params.beta_asymptotic,
                                       m_params.beta_coeffs, r);
        beta_x += legendre_angular_sum(m_params.num_beta_angular_modes,
                                       m_params.beta_angular_modes, r, dz);
        data_t K = radial_profile((data_t)m_params.K_asymptotic,
                                  m_params.K_coeffs, r);
        K += legendre_angular_sum(m_params.num_K_angular_modes,
                                  m_params.K_angular_modes, r, dz);
        const data_t phi = radial_profile((data_t)m_params.phi_asymptotic,
                                          m_params.phi_coeffs, r);
        const data_t Pi = radial_profile((data_t)m_params.Pi_asymptotic,
                                         m_params.Pi_coeffs, r);

        if (chi < (data_t)1.0e-10)
        {
            chi = (data_t)1.0e-10;
        }
        if (lapse < (data_t)1.0e-10)
        {
            lapse = (data_t)1.0e-10;
        }

        cell(i, j, k, c_chi) = chi;
        cell(i, j, k, c_h11) = (data_t)1.0;
        cell(i, j, k, c_h12) = (data_t)0.0;
        cell(i, j, k, c_h13) = (data_t)0.0;
        cell(i, j, k, c_h22) = (data_t)1.0;
        cell(i, j, k, c_h23) = (data_t)0.0;
        cell(i, j, k, c_h33) = (data_t)1.0;
        cell(i, j, k, c_A11) = (data_t)0.0;
        cell(i, j, k, c_A12) = (data_t)0.0;
        cell(i, j, k, c_A13) = (data_t)0.0;
        cell(i, j, k, c_A22) = (data_t)0.0;
        cell(i, j, k, c_A23) = (data_t)0.0;
        cell(i, j, k, c_A33) = (data_t)0.0;
        cell(i, j, k, c_Theta) = (data_t)0.0;
        cell(i, j, k, c_Gamma1) = (data_t)0.0;
        cell(i, j, k, c_Gamma2) = (data_t)0.0;
        cell(i, j, k, c_Gamma3) = (data_t)0.0;
        cell(i, j, k, c_lapse) = lapse;
        cell(i, j, k, c_shift1) = beta_x;
        cell(i, j, k, c_shift2) = (data_t)0.0;
        cell(i, j, k, c_shift3) = (data_t)0.0;
        cell(i, j, k, c_B1) = (data_t)0.0;
        cell(i, j, k, c_B2) = (data_t)0.0;
        cell(i, j, k, c_B3) = (data_t)0.0;
        cell(i, j, k, c_K) = K;
        cell(i, j, k, c_phi) = phi;
        cell(i, j, k, c_Pi) = Pi;
    }

  private:
    params_t m_params;
    double m_dx;

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    legendre_p(int ell, data_t mu) const
    {
        if (ell == 0)
        {
            return (data_t)1.0;
        }
        if (ell == 1)
        {
            return mu;
        }
        if (ell == 2)
        {
            return (data_t)0.5 * ((data_t)3.0 * mu * mu - (data_t)1.0);
        }
        if (ell == 3)
        {
            return (data_t)0.5 *
                   ((data_t)5.0 * mu * mu * mu - (data_t)3.0 * mu);
        }
        return (data_t)0.0;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    radial_gaussian_bump(data_t r, data_t center, data_t width) const
    {
        const data_t width_sq = width * width;
        const data_t dr = r - center;
        return std::exp(-dr * dr / ((data_t)2.0 * width_sq));
    }

    // Axisymmetric (Legendre) angular sum shared by chi, lapse, shift and K:
    // sum_n amp_n * exp(-(r-rc_n)^2/(2 rw_n^2)) * P_{ell_n}(cos theta).
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    legendre_angular_sum(int num_modes,
                         const std::array<AngularMode, MAX_ANGULAR_MODES> &modes,
                         data_t r, data_t dz) const
    {
        data_t delta        = (data_t)0.0;
        const data_t r_safe = std::max(r, (data_t)1.0e-12);
        const data_t mu     = dz / r_safe;
        for (int n = 0; n < num_modes; ++n)
        {
            const auto &mode = modes[n];
            const data_t radial = radial_gaussian_bump(
                r, (data_t)mode.radial_center, (data_t)mode.radial_width);
            delta += (data_t)mode.amplitude * radial * legendre_p(mode.ell, mu);
        }
        return delta;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    chi_angular_contribution(data_t r, data_t dx, data_t dy, data_t dz) const
    {
        data_t delta = legendre_angular_sum(m_params.num_chi_angular_modes,
                                            m_params.chi_angular_modes, r, dz);

        for (int n = 0; n < m_params.num_chi_Ylm_modes; ++n)
        {
            const auto &mode = m_params.chi_Ylm_modes[n];
            const data_t radial = radial_gaussian_bump(
                r, (data_t)mode.radial_center, (data_t)mode.radial_width);
            const auto ylm = SphericalHarmonics::spin_Y_lm(
                static_cast<amrex::Real>(dx), static_cast<double>(dy),
                static_cast<double>(dz), 0, mode.ell, mode.em);
            delta += (data_t)mode.amplitude * radial * (data_t)ylm.Real;
        }

        return delta;
    }

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE data_t
    radial_profile(data_t asymptotic, const std::array<double, MAX_BASES> &coeffs,
                   data_t r) const
    {
        const int num_bases = m_params.num_bases;
        const data_t width_sq =
            (data_t)(m_params.basis_width * m_params.basis_width);
        const data_t inv_two_width_sq = (data_t)0.5 / width_sq;

        data_t value = asymptotic;
        if (num_bases <= 1)
        {
            const data_t dr = r - (data_t)0.0;
            value += (data_t)coeffs[0] *
                     std::exp(-dr * dr * inv_two_width_sq);
            return value;
        }

        const data_t radius_max = (data_t)m_params.basis_radius_max;
        const data_t inv_denom = (data_t)1.0 / (data_t)(num_bases - 1);
        for (int n = 0; n < num_bases; ++n)
        {
            const data_t node = radius_max * (data_t)n * inv_denom;
            const data_t dr = r - node;
            value += (data_t)coeffs[n] *
                     std::exp(-dr * dr * inv_two_width_sq);
        }
        return value;
    }
};

#endif /* RADIALRECIPEINITIALDATA_HPP_ */
