/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHFINDERPARAMETERS_HPP_
#define AHFINDERPARAMETERS_HPP_

#include "GRParmParse.hpp"

#include <AMReX_REAL.H>

// The knobs of AHFinder's implicit pseudo-transient-continuation (PTC) solve,
// read from the "ah_finder" scope of the input file. Defaults are registered
// (and validated) by check_params(), which runs from
// SimulationParameters::check_params() at amrex::Initialize() time, so
// fill_params() can just get() them.
//
// Each PTC iteration solves the linear system
//   (I/dt + J) delta_h = -Theta(h_n),   h_{n+1} = h_n + delta_h
// with J = dTheta/dh (frozen-metric Jacobian) via a matrix-free BiCGStab, and
// grows the pseudo-timestep dt with an SER rule (dt *= theta_old/theta_new) so
// the method approaches Newton as the residual falls. The SER bounds
// (m_min_dt, m_dt_shrink, m_dt_grow, m_theta_floor) remain hard-coded in
// AHFinder::init().
struct ah_finder_params_t
{
    // Convergence threshold on the inf-norm of the expansion Theta.
    amrex::Real tolerance{};

    // Target growth factor for the adaptive pseudo-timestep dt: each step
    // scales dt by r * theta_old / theta_new, clamped to
    // [m_dt_shrink, m_dt_grow].
    amrex::Real r{};

    // Maximum number of PTC iterations before find() gives up.
    int max_iter{};

    // Relative and absolute tolerances passed to the per-step BiCGStab solve.
    amrex::Real linear_rel_tol{};
    amrex::Real linear_abs_tol{};

    static void check_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        amrex::Real tolerance = 1e-4;
        ah_finder_pp.queryAdd("tolerance", tolerance);
        if (tolerance <= 0.0)
        {
            ah_finder_pp.error("tolerance", "must be > 0 or find() cannot "
                                            "terminate");
        }

        amrex::Real r = 1.15;
        ah_finder_pp.queryAdd("r", r);
        if (r <= 0.0)
        {
            ah_finder_pp.error("r", "must be > 0");
        }

        int max_iter = 200;
        ah_finder_pp.queryAdd("max_iter", max_iter);
        if (max_iter <= 0)
        {
            ah_finder_pp.error("max_iter", "must be > 0");
        }

        amrex::Real linear_rel_tol = 1e-6;
        ah_finder_pp.queryAdd("linear_rel_tol", linear_rel_tol);
        if (linear_rel_tol <= 0.0)
        {
            ah_finder_pp.error("linear_rel_tol", "must be > 0");
        }

        amrex::Real linear_abs_tol = 0.0;
        ah_finder_pp.queryAdd("linear_abs_tol", linear_abs_tol);
        if (linear_abs_tol < 0.0)
        {
            ah_finder_pp.error("linear_abs_tol", "must be >= 0");
        }
    }

    void fill_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        ah_finder_pp.get("tolerance", tolerance);
        ah_finder_pp.get("r", r);
        ah_finder_pp.get("max_iter", max_iter);
        ah_finder_pp.get("linear_rel_tol", linear_rel_tol);
        ah_finder_pp.get("linear_abs_tol", linear_abs_tol);
    }
};

#endif /* AHFINDERPARAMETERS_HPP_ */
