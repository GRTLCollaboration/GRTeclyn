/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef KLEINGORDONRHS_HPP_
#define KLEINGORDONRHS_HPP_

#include "Cell.hpp"
#include "FourthOrderDerivatives.hpp"
#include "TensorAlgebra.hpp"
#include "VarsTools.hpp"

// Problem specific includes
#include "Potential.hpp"
#include "StateVariables.hpp"

template <class deriv_t = FourthOrderDerivatives, class potential_t = Potential>
class KleinGordonRHS
{
  public:

    // NOLINTNEXTLINE(bugprone-easily-swappable-parameters)
    KleinGordonRHS(amrex::Real a_sigma, amrex::Real a_dx,
                   const potential_t a_potential)
        : m_sigma(a_sigma), m_deriv(a_dx), m_potential(a_potential) {};

    template <class data_t> struct Vars
    {
        data_t phi, Pi;

        template <typename mapping_function_t>
        AMREX_GPU_DEVICE void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
            VarsTools::define_enum_mapping(mapping_function, c_Pi, Pi);
        }
    };

    template <class data_t> struct Diff2Vars
    {
        data_t phi;

        template <typename mapping_function_t>
        AMREX_GPU_DEVICE void enum_mapping(mapping_function_t mapping_function)
        {
            VarsTools::define_enum_mapping(mapping_function, c_phi, phi);
        }
    };

    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t const> &input,
            const amrex::Array4<data_t> &output) const;


  private:
    amrex::Real m_sigma;
    deriv_t m_deriv;
    potential_t m_potential;

    template <class data_t, template <typename> class vars_t,
              template <typename> class diff2_vars_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    rhs_equation(vars_t<data_t> &rhs, const vars_t<data_t> &vars,
                 const vars_t<Tensor<1, data_t>> &d1,
                 const diff2_vars_t<Tensor<2, data_t>> &d2) const;
};

#include "KleinGordonRHS.impl.hpp"

#endif
