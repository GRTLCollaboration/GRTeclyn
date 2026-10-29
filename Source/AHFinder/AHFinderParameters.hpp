/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AHFINDERPARAMETERS_HPP_
#define AHFINDERPARAMETERS_HPP_

#include "GRParmParse.hpp"

#include <AMReX_REAL.H>

// The knobs of AHFinder's pseudo-time evolution that are worth varying between
// runs, read from the "ah_finder" scope of the input file. Defaults are
// registered (and validated) by check_params(), which runs from
// SimulationParameters::check_params() at amrex::Initialize() time, so
// fill_params() can just get() them.
//
// The surface radius h is driven towards the horizon by pseudo-transient
// continuation of h_dot = -Theta, with Theta the expansion: each step solves
//   (J + I/dt) delta_h = -Theta,   J = dTheta/dh,
// with a banded LU solve, and dt grows as Theta falls so the step tends to a
// Newton step. This continues until the inf-norm of Theta drops below
// `tolerance`, `max_iter` iterations have been taken, or the pseudo-timestep
// has sat at its minimum for `max_stall_iter` consecutive iterations. The
// remaining timestepping bounds (m_min_dt, m_dt_shrink)
// and the finite-difference scale m_fd_sigma are still
// hard-coded in AHFinder::init().
struct ah_finder_params_t
{
    // Convergence threshold on the inf-norm of the expansion Theta.
    amrex::Real tolerance{};

    // Target growth factor for the adaptive pseudo-timestep: each step scales
    // dt by r * theta_old / theta_new, shrinking by no more than
    // m_dt_shrink. Growth is unbounded.
    amrex::Real r{};

    // Pseudo-timestep used for the first step of each find(). Subsequent
    // steps adapt it as above.
    amrex::Real initial_dt{};

    // Maximum number of pseudo-time iterations before find() gives up
    // without convergence.
    int max_iter{};

    // Maximum number of consecutive iterations the pseudo-timestep may sit at
    // its minimum (m_min_dt) before find() treats the solve as stalled and
    // gives up without convergence.
    int max_stall_iter{};

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

        amrex::Real initial_dt = 1e-2;
        ah_finder_pp.queryAdd("initial_dt", initial_dt);
        if (initial_dt <= 0.0)
        {
            ah_finder_pp.error("initial_dt", "must be > 0");
        }

        int max_iter = 200;
        ah_finder_pp.queryAdd("max_iter", max_iter);
        if (max_iter <= 0)
        {
            ah_finder_pp.error("max_iter", "must be > 0");
        }

        int max_stall_iter = 50;
        ah_finder_pp.queryAdd("max_stall_iter", max_stall_iter);
        if (max_stall_iter <= 0)
        {
            ah_finder_pp.error("max_stall_iter", "must be > 0");
        }
    }

    void fill_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        ah_finder_pp.get("tolerance", tolerance);
        ah_finder_pp.get("r", r);
        ah_finder_pp.get("initial_dt", initial_dt);
        ah_finder_pp.get("max_iter", max_iter);
        ah_finder_pp.get("max_stall_iter", max_stall_iter);
    }
};

#endif /* AHFINDERPARAMETERS_HPP_ */
