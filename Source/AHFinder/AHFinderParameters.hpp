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
// The surface radius h is evolved as a damped wave towards the horizon:
//   h_dot = v - eta * h
//   v_dot = -c^2 * Theta
// with Theta the expansion, until the inf-norm of Theta drops below
// `tolerance`. The remaining timestepping bounds (m_min_dt, m_dt_shrink,
// m_dt_grow, m_theta_floor) are still hard-coded in AHFinder::init().
struct ah_finder_params_t
{
    // Damping coefficient eta. Too small and the surface rings; too large and
    // it crawls.
    amrex::Real eta{};

    // Wave speed c driving v towards the zero-expansion surface.
    amrex::Real c{};

    // Convergence threshold on the inf-norm of the expansion Theta.
    amrex::Real tolerance{};

    // Target growth factor for the adaptive pseudo-timestep: each step scales
    // dt by r * theta_old / theta_new, clamped to [m_dt_shrink, m_dt_grow].
    amrex::Real r{};

    // Safety factor on the CFL cap: dt is never allowed above
    // cfl_factor * (smallest coordinate spacing between neighbouring points
    // on the ring grid). Larger values let dt grow further before the cap
    // bites, at the cost of stability.
    amrex::Real cfl_factor{};

    static void check_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        amrex::Real eta = 3.0;
        ah_finder_pp.queryAdd("eta", eta);
        if (eta <= 0.0)
        {
            ah_finder_pp.error("eta", "must be > 0 for the surface to damp");
        }

        amrex::Real c = 1.0;
        ah_finder_pp.queryAdd("c", c);
        if (c <= 0.0)
        {
            ah_finder_pp.error("c", "must be > 0");
        }

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

        amrex::Real cfl_factor = 7.0;
        ah_finder_pp.queryAdd("cfl_factor", cfl_factor);
        if (cfl_factor <= 0.0)
        {
            ah_finder_pp.error("cfl_factor", "must be > 0");
        }
    }

    void fill_params()
    {
        GRParmParse ah_finder_pp("ah_finder");

        ah_finder_pp.get("eta", eta);
        ah_finder_pp.get("c", c);
        ah_finder_pp.get("tolerance", tolerance);
        ah_finder_pp.get("r", r);
        ah_finder_pp.get("cfl_factor", cfl_factor);
    }
};

#endif /* AHFINDERPARAMETERS_HPP_ */
