/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SCALARFIELDLEVEL_HPP_
#define SCALARFIELDLEVEL_HPP_

#include "DefaultLevelBld.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRAmrLevel.hpp"
#include "Potential.hpp"
#include "ScalarField.hpp"
#include "ScalarFieldAmr.hpp"
#include "SixthOrderDerivatives.hpp"

/// Evolution level for a real scalar field minimally coupled to gravity.
class ScalarFieldLevel : public GRAmrLevel
{
  public:
    using GRAmrLevel::GRAmrLevel;

    template <class deriv_t = FourthOrderDerivatives>
    using ScalarFieldWithPotential = ScalarField<Potential, deriv_t>;

    ScalarFieldAmr *get_scalar_field_amr_ptr();

    static void variableSetUp();

    void specific_advance() override;

    void initData() override;

    void specificEvalRHS(amrex::MultiFab &a_soln, amrex::MultiFab &a_rhs,
                         amrex::Real a_time) override;

    void specificUpdateODE(amrex::MultiFab &a_soln) override;

    void specific_post_timestep() override;

    void tag_cells(amrex::TagBoxArray &a_tag_box_array,
                   amrex::Real a_regrid_threshold) final;
};

#endif /* SCALARFIELDLEVEL_HPP_ */
