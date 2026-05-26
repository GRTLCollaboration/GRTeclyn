#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

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
    }

    void read_recipe_params(GRParmParse &pp)
    {
        pp.load("recipe_num_bases", recipe_params.num_bases, 4);
        pp.load("center", recipe_params.grid_center, center);
        pp.load("recipe_basis_width", recipe_params.basis_width, 1.0);
        pp.load("recipe_basis_radius_max", recipe_params.basis_radius_max, 8.0);

        pp.load("recipe_chi_asymptotic", recipe_params.chi_asymptotic, 1.0);
        pp.load("recipe_alpha_asymptotic", recipe_params.alpha_asymptotic, 1.0);
        pp.load("recipe_K_asymptotic", recipe_params.K_asymptotic, 0.0);
        pp.load("recipe_phi_asymptotic", recipe_params.phi_asymptotic, 0.0);
        pp.load("recipe_Pi_asymptotic", recipe_params.Pi_asymptotic, 0.0);

        load_coeff_array(pp, "recipe_chi_coeff", recipe_params.chi_coeffs);
        load_coeff_array(pp, "recipe_alpha_coeff", recipe_params.alpha_coeffs);
        load_coeff_array(pp, "recipe_K_coeff", recipe_params.K_coeffs);
        load_coeff_array(pp, "recipe_phi_coeff", recipe_params.phi_coeffs);
        load_coeff_array(pp, "recipe_Pi_coeff", recipe_params.Pi_coeffs);
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
    }

    bool calculate_constraint_norms{};

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
};

#endif /* SIMULATIONPARAMETERS_HPP */
