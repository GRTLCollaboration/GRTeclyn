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
// Newton step.
//
// If a step does not reduce the inf-norm of Theta, it is backtracked along the
// same direction, trying alpha * delta_h with alpha = backtrack_factor^k for
// k = 1..max_backtracks, without re-solving. If none of these reduce Theta, h
// is left unchanged for that iteration and dt is shrunk.
//
// This continues until the inf-norm of Theta drops below `tolerance`,
// `max_iter` iterations have been taken, or backtracking fails with dt
// already at its minimum, so no smaller step can be tried. The
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

    // Factor by which the step length alpha is reduced after each rejected
    // trial step.
    amrex::Real backtrack_factor{};

    // Maximum number of shortened trial steps after a rejected full step
    // before the iteration is abandoned (h reverted, dt shrunk). 0 disables
    // backtracking, so a step that does not reduce Theta is simply rejected.
    int max_backtracks{};

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

        amrex::Real backtrack_factor = 0.5;
        ah_finder_pp.queryAdd("backtrack_factor", backtrack_factor);
        if (backtrack_factor <= 0.0 || backtrack_factor >= 1.0)
        {
            ah_finder_pp.error("backtrack_factor", "must be in (0, 1)");
        }

        int max_backtracks = 5;
        ah_finder_pp.queryAdd("max_backtracks", max_backtracks);
        if (max_backtracks < 0)
        {
            ah_finder_pp.error("max_backtracks", "must be >= 0");
        }
    }

    void fill_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        ah_finder_pp.get("tolerance", tolerance);
        ah_finder_pp.get("r", r);
        ah_finder_pp.get("initial_dt", initial_dt);
        ah_finder_pp.get("max_iter", max_iter);
        ah_finder_pp.get("backtrack_factor", backtrack_factor);
        ah_finder_pp.get("max_backtracks", max_backtracks);
    }
};

#endif /* AHFINDERPARAMETERS_HPP_ */
