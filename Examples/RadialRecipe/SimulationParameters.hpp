#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

#include "ExternalGridInitialData.hpp"
#include "GRParmParse.hpp"
#include "RadialRecipeInitialData.hpp"
#include "SimulationParametersBase.hpp"

#include <sstream>
#include <string>

class SimulationParameters : public SimulationParametersBase
{
  public:
    SimulationParameters(GRParmParse &pp) : SimulationParametersBase(pp)
    {
        read_shared_params(pp);
        read_recipe_params(pp);
        check_params();
    }

    void read_shared_params(GRParmParse &pp)
    {
        pp.load("calculate_constraint_norms", calculate_constraint_norms, false);
        pp.load("calculate_energy_conditions", calculate_energy_conditions,
                false);
        pp.load("calculate_curvature_invariants",
                calculate_curvature_invariants, false);

        // Matter sector selection.  The constrained initial-data recipe may
        // require exotic (phantom, rho <= 0) matter to source a wormhole/warp
        // geometry.  When recipe_exotic_matter is set the level evolves an
        // ExoticScalarField whose stress-energy is the canonical one scaled by
        // -recipe_support_strength, so the matter that is actually evolved
        // matches the matter the geometry was reconstructed for.
        pp.load("recipe_exotic_matter", recipe_exotic_matter, false);
        pp.load("recipe_support_strength", recipe_support_strength, 1.0);

        pp.load("recipe_initial_data_file", recipe_initial_data_file,
                std::string(""));
        if (!recipe_initial_data_file.empty())
        {
            external_grid_params.gridinit_file = recipe_initial_data_file;
            external_grid_params.grid_center = center;
        }
    }

    void read_recipe_params(GRParmParse &pp)
    {
        pp.load("recipe_num_bases", recipe_params.num_bases, 4);
        pp.load("center", recipe_params.grid_center, center);
        pp.load("recipe_basis_width", recipe_params.basis_width, 1.0);
        pp.load("recipe_basis_radius_max", recipe_params.basis_radius_max, 8.0);

        pp.load("recipe_chi_asymptotic", recipe_params.chi_asymptotic, 1.0);
        pp.load("recipe_alpha_asymptotic", recipe_params.alpha_asymptotic, 1.0);
        pp.load("recipe_beta_asymptotic", recipe_params.beta_asymptotic, 0.0);
        pp.load("recipe_K_asymptotic", recipe_params.K_asymptotic, 0.0);
        pp.load("recipe_phi_asymptotic", recipe_params.phi_asymptotic, 0.0);
        pp.load("recipe_Pi_asymptotic", recipe_params.Pi_asymptotic, 0.0);

        load_coeff_array(pp, "recipe_chi_coeff", recipe_params.chi_coeffs);
        load_coeff_array(pp, "recipe_alpha_coeff", recipe_params.alpha_coeffs);
        load_coeff_array(pp, "recipe_beta_coeff", recipe_params.beta_coeffs);
        load_coeff_array(pp, "recipe_K_coeff", recipe_params.K_coeffs);
        load_coeff_array(pp, "recipe_phi_coeff", recipe_params.phi_coeffs);
        load_coeff_array(pp, "recipe_Pi_coeff", recipe_params.Pi_coeffs);

        pp.load("recipe_num_chi_angular_modes",
                recipe_params.num_chi_angular_modes, 0);
        load_angular_modes(pp, "recipe_chi_mode",
                           recipe_params.num_chi_angular_modes,
                           recipe_params.chi_angular_modes);

        pp.load("recipe_num_chi_Ylm_modes", recipe_params.num_chi_Ylm_modes, 0);
        load_Ylm_modes(pp);

        pp.load("recipe_num_lapse_angular_modes",
                recipe_params.num_lapse_angular_modes, 0);
        load_angular_modes(pp, "recipe_lapse_mode",
                           recipe_params.num_lapse_angular_modes,
                           recipe_params.lapse_angular_modes);

        pp.load("recipe_num_beta_angular_modes",
                recipe_params.num_beta_angular_modes, 0);
        load_angular_modes(pp, "recipe_beta_mode",
                           recipe_params.num_beta_angular_modes,
                           recipe_params.beta_angular_modes);

        pp.load("recipe_num_K_angular_modes",
                recipe_params.num_K_angular_modes, 0);
        load_angular_modes(pp, "recipe_K_mode",
                           recipe_params.num_K_angular_modes,
                           recipe_params.K_angular_modes);
    }

    void check_params()
    {
        check_parameter("recipe_num_bases", recipe_params.num_bases,
                        recipe_params.num_bases > 0 &&
                            recipe_params.num_bases <=
                                RadialRecipeInitialData::MAX_BASES,
                        "must be between 1 and MAX_BASES");

        check_parameter("recipe_basis_width", recipe_params.basis_width,
                        recipe_params.basis_width > 0.0, "must be positive");

        check_parameter("recipe_basis_radius_max",
                        recipe_params.basis_radius_max,
                        recipe_params.basis_radius_max > 0.0,
                        "must be positive");

        check_parameter("recipe_num_chi_angular_modes",
                        recipe_params.num_chi_angular_modes,
                        recipe_params.num_chi_angular_modes >= 0 &&
                            recipe_params.num_chi_angular_modes <=
                                RadialRecipeInitialData::MAX_ANGULAR_MODES,
                        "must be between 0 and MAX_ANGULAR_MODES");

        check_parameter("recipe_num_chi_Ylm_modes",
                        recipe_params.num_chi_Ylm_modes,
                        recipe_params.num_chi_Ylm_modes >= 0 &&
                            recipe_params.num_chi_Ylm_modes <=
                                RadialRecipeInitialData::MAX_YLM_MODES,
                        "must be between 0 and MAX_YLM_MODES");

        check_angular_mode_count("recipe_num_lapse_angular_modes",
                                 recipe_params.num_lapse_angular_modes);
        check_angular_mode_count("recipe_num_beta_angular_modes",
                                 recipe_params.num_beta_angular_modes);
        check_angular_mode_count("recipe_num_K_angular_modes",
                                 recipe_params.num_K_angular_modes);
    }

    void check_angular_mode_count(const char *name, int count)
    {
        check_parameter(name, count,
                        count >= 0 &&
                            count <= RadialRecipeInitialData::MAX_ANGULAR_MODES,
                        "must be between 0 and MAX_ANGULAR_MODES");
    }

    bool calculate_constraint_norms{};
    bool calculate_energy_conditions{};
    bool calculate_curvature_invariants{};
    bool recipe_exotic_matter{};
    double recipe_support_strength{1.0};

    std::string recipe_initial_data_file;
    ExternalGridInitialData::params_t external_grid_params{};

    RadialRecipeInitialData::params_t recipe_params{};

  private:
    void load_coeff_array(GRParmParse &pp, const char *prefix,
                          std::array<double, RadialRecipeInitialData::MAX_BASES>
                              &coeffs)
    {
        for (int n = 0; n < recipe_params.num_bases; ++n)
        {
            std::ostringstream key;
            key << prefix << "_" << n;
            pp.load(key.str().c_str(), coeffs[n], 0.0);
        }
    }

    void load_angular_modes(
        GRParmParse &pp, const char *prefix, int num_modes,
        std::array<RadialRecipeInitialData::AngularMode,
                   RadialRecipeInitialData::MAX_ANGULAR_MODES> &modes)
    {
        for (int n = 0; n < num_modes; ++n)
        {
            auto &mode = modes[n];
            std::ostringstream ell_key, amp_key, rc_key, rw_key;
            ell_key << prefix << "_ell_" << n;
            amp_key << prefix << "_amp_" << n;
            rc_key << prefix << "_rc_" << n;
            rw_key << prefix << "_rw_" << n;
            pp.load(ell_key.str().c_str(), mode.ell, 0);
            pp.load(amp_key.str().c_str(), mode.amplitude, 0.0);
            pp.load(rc_key.str().c_str(), mode.radial_center, 0.0);
            pp.load(rw_key.str().c_str(), mode.radial_width, 1.0);
        }
    }

    void load_Ylm_modes(GRParmParse &pp)
    {
        for (int n = 0; n < recipe_params.num_chi_Ylm_modes; ++n)
        {
            auto &mode = recipe_params.chi_Ylm_modes[n];
            std::ostringstream l_key, m_key, amp_key, rc_key, rw_key;
            l_key << "recipe_chi_Ylm_l_" << n;
            m_key << "recipe_chi_Ylm_m_" << n;
            amp_key << "recipe_chi_Ylm_amp_" << n;
            rc_key << "recipe_chi_Ylm_rc_" << n;
            rw_key << "recipe_chi_Ylm_rw_" << n;
            pp.load(l_key.str().c_str(), mode.ell, 0);
            pp.load(m_key.str().c_str(), mode.em, 0);
            pp.load(amp_key.str().c_str(), mode.amplitude, 0.0);
            pp.load(rc_key.str().c_str(), mode.radial_center, 0.0);
            pp.load(rw_key.str().c_str(), mode.radial_width, 1.0);
        }
    }
};

#endif /* SIMULATIONPARAMETERS_HPP */
