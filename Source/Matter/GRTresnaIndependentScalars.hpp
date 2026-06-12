#ifndef GRTRESNA_INDEPENDENT_SCALARS_HPP_
#define GRTRESNA_INDEPENDENT_SCALARS_HPP_

#include "CCZ4Geometry.hpp"
#include "Coordinates.hpp"
#include "DimensionDefinitions.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRTresnaIndependentScalarsAdvecVars.hpp"
#include "GRTresnaIndependentScalarsD1Vars.hpp"
#include "GRTresnaIndependentScalarsD2Vars.hpp"
#include "GRTresnaIndependentScalarsVars.hpp"
#include "GRTresnaScalarLayout.hpp"
#include "GRTresnaScalarPotential.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"

#include <array>
#include <sstream>
#include <string>

//! Independent signed scalar lumps matching the GRTresna constraint solve.
class GRTresnaIndependentScalars
{
  public:
    GRTresnaIndependentScalars() { load_from_inputs(); }

    GRTresnaIndependentScalars(
        int a_num_fields,
        const std::array<int, GRTRESNA_MAX_INDEPENDENT_SCALARS> &a_signs,
        double a_mass, double a_lambda = 0.0)
        : m_num_fields(a_num_fields), m_signs(a_signs), m_mass_param(a_mass),
          m_lambda_param(a_lambda),
          m_potential(GRTresnaScalarPotential(a_mass, a_lambda))
    {
    }

    void load_from_inputs()
    {
        GRParmParse pp;
        std::string model;
        pp.load("recipe_matter_model", model, std::string(""));
        if (model != "grtresna_independent_scalars")
        {
            return;
        }
        pp.load("recipe_num_scalar_fields", m_num_fields, 0);
        pp.load("recipe_scalar_mass", m_mass_param, 0.0);
        pp.load("recipe_scalar_lambda", m_lambda_param, 0.0);
        m_potential = GRTresnaScalarPotential(m_mass_param, m_lambda_param);

        std::string signs_line;
        pp.load("recipe_scalar_field_signs", signs_line, std::string(""));
        if (!signs_line.empty())
        {
            std::istringstream iss(signs_line);
            int sign = 0;
            int idx  = 0;
            while (iss >> sign && idx < GRTRESNA_MAX_INDEPENDENT_SCALARS)
            {
                m_signs[idx++] = sign;
            }
            return;
        }
        for (int k = 0; k < GRTRESNA_MAX_INDEPENDENT_SCALARS; ++k)
        {
            std::ostringstream key;
            key << "recipe_scalar_field_sign_" << k;
            pp.load(key.str().c_str(), m_signs[k], 1);
        }
    }

    [[nodiscard]] int num_fields() const { return m_num_fields; }

    using Vars      = GRTresnaIndependentScalarsVars;
    using D1Vars    = GRTresnaIndependentScalarsD1Vars;
    using D2Vars    = GRTresnaIndependentScalarsD2Vars;
    using AdvecVars = GRTresnaIndependentScalarsAdvecVars;

    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL) const;

    [[nodiscard]]
    AMREX_GPU_DEVICE emtensor_t
    compute_emtensor(const Vars &vars, const D1Vars &d1,
                     const Tensor<2, amrex::Real> &h_UU,
                     const Tensor<3, amrex::Real> &chris_ULL,
                     const Coordinates &coords, amrex::Real time) const;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    add_matter_rhs(const amrex::CellData<amrex::Real> &rhs, const Vars &vars,
                   const D1Vars &d1, const D2Vars &d2,
                   const AdvecVars &advec) const;

  private:
    int m_num_fields{};
    std::array<int, GRTRESNA_MAX_INDEPENDENT_SCALARS> m_signs{};
    double m_mass_param{};
    double m_lambda_param{};
    GRTresnaScalarPotential m_potential;

    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE amrex::Real
    phi_sum(const Vars &vars) const
    {
        amrex::Real sum = 0.0;
        for (int k = 0; k < m_num_fields; ++k)
        {
            sum += vars.phi(k);
        }
        return sum;
    }
};

#include "GRTresnaIndependentScalars.impl.hpp"

#endif /* GRTRESNA_INDEPENDENT_SCALARS_HPP_ */
