/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INITIALSCALARDATA_HPP_
#define INITIALSCALARDATA_HPP_

#include "Cell.hpp"
#include "Coordinates.hpp"
#include "MatterCCZ4RHS.hpp"
#include "ScalarField.hpp"
#include "StateVariables.hpp" //This files needs NUM_VARS - total no. components
#include "Tensor.hpp"
#include "VarsTools.hpp"
#include "simd.hpp"

//! Class which sets the initial scalar field matter config
class InitialScalarData
{
  public:
    //! A structure for the input params for scalar field properties and initial
    //! conditions
    struct params_t
    {
        double amplitude; //!< Amplitude of bump in initial SF bubble
        std::array<double, AMREX_SPACEDIM>
            center;   //!< Centre of perturbation in initial SF bubble
        double width; //!< Width of bump in initial SF bubble
    };

    //! The constructor
    InitialScalarData(params_t a_params, double a_dx)
        : m_dx(a_dx), m_params(a_params)
    {
    }

    //! Function to compute the value of all the initial vars on the grid
    template <class data_t>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    compute(int i, int j, int k, const amrex::Array4<data_t> &cell) const
    {
        MatterCCZ4RHS<ScalarField<>>::Vars<data_t> vars;
        VarsTools::assign(vars, 0.); // Set only the non-zero components below

        // start with unit lapse and flat metric (must be relaxed for chi)
        vars.lapse = 1.0;
        vars.chi   = 1.0;

        // conformal metric is flat
        FOR (index)
            vars.h[index][index] = 1.;

        // where am i?
        amrex::IntVect pos(i, j, k);
        Coordinates<data_t> coords(pos, m_dx, m_params.center);
        data_t rr  = coords.get_radius();
        data_t rr2 = rr * rr;

        // calculate the field value
        vars.phi = m_params.amplitude *
                   (1.0 + 0.01 * rr2 * exp(-pow(rr / m_params.width, 2.0)));
        vars.Pi = 0;

        // store the vars
        //        cell(i, j, k, c_phi) = phi;
        //        cell(i, j, k, c_Pi)  = 0.0;

        // Store the initial values of the variables
        store_vars(cell.cellData(i, j, k), vars);
    }

  protected:
    double m_dx;
    const params_t m_params; //!< The matter initial condition params
};

#endif /* INITIALSCALARDATA_HPP_ */
