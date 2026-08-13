/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SPHERICALEXTRACTION_HPP_
#define SPHERICALEXTRACTION_HPP_

#include "SphericalGeometry.hpp"
#include "SphericalHarmonics.hpp"
#include "SurfaceExtraction.hpp"

// Parameters

#include "SphericalExtractionParameters.hpp"

//! A child class of SurfaceExtraction for extraction on spherical shells
template <int num_components>
class SphericalExtraction
    : public SurfaceExtraction<SphericalGeometry, num_components>
{
  public:
    using Base   = SurfaceExtraction<SphericalGeometry, num_components>;
    using vars_t = typename Base::vars_t;
    using mode_integrals_t =
        std::pair<std::vector<double>, std::vector<double>>;
    using params_t = spherical_extraction_params_t;

  protected:
    std::array<double, AMREX_SPACEDIM> m_center;
    int m_num_modes;
    std::vector<std::pair<int, int>> m_modes{};

  public:
    SphericalExtraction(const params_t &a_params, double a_dt, double a_time,
                        bool a_first_step, double a_restart_time = 0.0)
        : Base(a_params.get_surface_extraction_params(), a_dt, a_time,
               a_first_step, a_restart_time),
          m_center(a_params.center), m_num_modes(a_params.num_modes),
          m_modes(a_params.modes)
    {
    }

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    SphericalExtraction(const params_t &a_params,
                        const std::vector<vars_t> &a_vars, double a_dt,
                        double a_time, bool a_first_step,
                        double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        this->add_vars(a_vars);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    // NOLINTBEGIN(bugprone-easily-swappable-parameters)
    SphericalExtraction(const params_t &a_params,
                        const std::vector<int> &a_state_vars, double a_dt,
                        double a_time, bool a_first_step,
                        double a_restart_time = 0.0)
        : SphericalExtraction(a_params, a_dt, a_time, a_first_step,
                              a_restart_time)
    {
        this->add_state_vars(a_state_vars);
    }
    // NOLINTEND(bugprone-easily-swappable-parameters)

    // alias this long type used for complex functions defined on the surface
    // and dependent on the interpolated data
    using complex_function_t = std::function<std::pair<double, double>(
        std::vector<double> &a_data_here, double r, double theta, double phi)>;

    //! Add the integrand corresponding to the spin-weighted spherical harmonic
    //! decomposition of a complex-valued function, a_function
    //! (normalised by 1/r^2), over each spherical shell
    // NOLINTBEGIN(readability-identifier-length)
    void add_mode_integrand(
        int es, int el, int em, const complex_function_t &a_function,
        mode_integrals_t &out_integrals,
        const IntegrationMethod &a_method_theta = IntegrationMethod::simpson,
        const IntegrationMethod &a_method_phi   = IntegrationMethod::trapezium,
        const bool a_broadcast_integral         = false)
    {
        auto integrand_re = [center = m_center, &geom = this->m_geom, es, el,
                             em, a_function](std::vector<double> &a_data_here,
                                             double r, double theta, double phi)
        {
            // note that spin_Y_lm requires the coordinates with the center
            // at the origin
            double x = geom.get_grid_coord(0, r, theta, phi) - center[0];
            double y = geom.get_grid_coord(1, r, theta, phi) - center[1];
            double z = geom.get_grid_coord(2, r, theta, phi) - center[2];
            SphericalHarmonics::Y_lm_t Y_lm =
                SphericalHarmonics::spin_Y_lm(x, y, z, es, el, em);
            auto function_here = a_function(a_data_here, r, theta, phi);
            return (function_here.first * Y_lm.Real +
                    function_here.second * Y_lm.Im) /
                   (r * r);
        };
        this->add_integrand(integrand_re, out_integrals.first, a_method_theta,
                            a_method_phi, a_broadcast_integral);

        auto integrand_im = [center = m_center, &geom = this->m_geom, es, el,
                             em, a_function](std::vector<double> &a_data_here,
                                             double r, double theta, double phi)
        {
            // note that spin_Y_lm requires the coordinates with the center
            // at the origin
            double x = geom.get_grid_coord(0, r, theta, phi) - center[0];
            double y = geom.get_grid_coord(1, r, theta, phi) - center[1];
            double z = geom.get_grid_coord(2, r, theta, phi) - center[2];
            SphericalHarmonics::Y_lm_t Y_lm =
                SphericalHarmonics::spin_Y_lm(x, y, z, es, el, em);
            auto function_here = a_function(a_data_here, r, theta, phi);
            return (function_here.second * Y_lm.Real -
                    function_here.first * Y_lm.Im) /
                   (r * r);
        };
        this->add_integrand(integrand_im, out_integrals.second, a_method_theta,
                            a_method_phi, a_broadcast_integral);
    }

    //! If you only want to extract one mode, you can use this function which
    //! calls add_mode_integrand, then integrate and returns the integrals
    std::pair<std::vector<double>, std::vector<double>> integrate_mode(
        int es, int el, int em, const complex_function_t &a_function,
        const IntegrationMethod &a_method_theta = IntegrationMethod::simpson,
        const IntegrationMethod &a_method_phi   = IntegrationMethod::trapezium)
    {
        this->m_integrands.clear();
        this->m_integration_methods.clear();
        this->m_integrals.clear();

        std::pair<std::vector<double>, std::vector<double>> integrals;
        add_mode_integrand(es, el, em, a_function, integrals, a_method_theta,
                           a_method_phi);
        this->integrate();
        return integrals;
    }
    // NOLINTEND(readability-identifier-length)
};

#endif /* SPHERICALEXTRACTION_HPP_ */