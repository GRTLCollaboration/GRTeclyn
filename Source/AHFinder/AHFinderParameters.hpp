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
// with J = dTheta/dh (frozen-metric Jacobian) via a matrix-free GMRES, and
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

    // Relative and absolute convergence tolerances of the per-step GMRES
    // solve.
    amrex::Real linear_rel_tol{};
    amrex::Real linear_abs_tol{};

    // Maximum GMRES iterations per PTC step. Each one costs a full Theta
    // evaluation over the ring grid, which is cheap next to the metric
    // interpolation a rejected step forces, so the budget wants to be generous:
    // measured Krylov counts run from ~75 at small dt to ~300 as dt grows and
    // (I/dt + J) becomes ill-conditioned. Capping below that truncates the
    // solve, and the resulting poor directions cause step rejections that cost
    // far more than the Krylov iterations saved -- see
    // docs/ah_finder_implicit_ptc.md.
    int linear_max_iter{};

    // GMRES restart length (Krylov vectors kept between restarts). Longer
    // restarts converge in fewer iterations but cost memory
    // (linear_restart_length * num_particles doubles) and orthogonalisation
    // work quadratic in the length. The default is large enough that the solve
    // never restarts at the resolutions tested; restarting part-way through
    // these solves was measured to stagnate rather than converge.
    int linear_restart_length{};

    // Multiplier on the Brown-Saad finite-difference step used by the Jacobian
    // mat-vec, eps = fd_eps_scale * sqrt(macheps) * (1 + ||h||) / ||v||. The
    // Brown-Saad value (scale 1) balances FD truncation error against roundoff
    // in Theta, which is the right trade-off for the frozen Jacobian: there the
    // mat-vec only re-evaluates theta_from_metric(), a smooth function of h.
    // With the metric unfrozen the mat-vec also re-runs interpolate_metric(),
    // which is only piecewise smooth in the particle positions (stencil and AMR
    // level selection change across cell boundaries), so a larger eps may be
    // needed to step over that noise.
    amrex::Real fd_eps_scale{};

    // Whether the Jacobian mat-vec re-interpolates the metric at the perturbed
    // surface. 0 (default) freezes the metric at h_n, so J omits the
    // d(metric)/dh term and is approximate but cheap. 1 gives the exact Newton
    // Jacobian at roughly the cost of a metric interpolation per Krylov
    // iteration instead of per PTC step -- around 50x more expensive.
    int unfreeze_jacobian{};

    // Hybrid policies: pay for the exact Jacobian only on selected PTC steps
    // and use the cheap frozen one otherwise. Ignored when unfreeze_jacobian
    // is 1 (which unfreezes unconditionally).
    //
    // The choice is necessarily per *PTC step*, not per mat-vec: GMRES builds
    // its Arnoldi basis assuming one fixed linear operator, so mixing frozen
    // and unfrozen applies inside a single solve would break the Krylov
    // relation. Each solve therefore uses one operator throughout; only which
    // operator varies from step to step.

    // Use the exact Jacobian on every Nth PTC step. 0 (default) disables.
    int unfreeze_every{};

    // Use the exact Jacobian on the step immediately after a rejected one. A
    // rejection means the frozen direction was an ascent direction, which is
    // exactly the failure the missing d(metric)/dh term causes, so this spends
    // the expensive Jacobian only where the cheap one has demonstrably failed.
    int unfreeze_on_reject{};

    // Newton-at-frozen-cost: the constant c in J_exact ~= J_frozen + c I. When
    // positive, dt is capped at 1/c, which makes (I/dt + J_frozen) the exact
    // Newton Jacobian at the top of the SER ramp and a Levenberg-Marquardt
    // damped one below it -- all using the cheap frozen mat-vec, with no metric
    // re-interpolation inside the solve. 0 (default) disables the cap.
    //
    // The cap also keeps dt clear of the pole of (I/dt + J_frozen) at
    // dt = -1/lambda_min ~ 0.125, which sits above 1/c ~ 0.084; walking into
    // that pole is what produces the accept/reject cycle.
    //
    // c is a property of the surface and the data (~10.6 at the initial guess,
    // ~12.0 at the horizon for the equal-mass binary test), so a hand-set value
    // does not transfer between problems. Prefer newton_shift_auto.
    amrex::Real newton_shift{};

    // Apply the measured shift as a diagonal, diag(c(x)), rather than as the
    // scalar mean c. Requires newton_shift_auto, which is what measures the
    // diagonal; costs nothing extra, since the scalar is only its mean.
    //
    // The scalar path caps dt at 1/c, which is the same as using a per-step
    // shift of max(1/dt_SER, c). The diagonal path drops the cap and applies
    // max(1/dt_SER, c(x_ip)) particle by particle instead, so it reduces
    // exactly to the scalar behaviour when c(x) is constant. Taking the max
    // with 1/dt also keeps the shift positive where c(x) is not, which happens
    // on strongly distorted surfaces.
    //
    // 0 (default) uses the scalar. Worth enabling when the surface is far from
    // spherical: on the individual horizon of the test binary c(x) varies by
    // 0.6% and the diagonal is not worth having, but on the common horizon of a
    // close binary it varies by 86% and the scalar costs a 5x slowdown -- see
    // docs/ah_finder_implicit_ptc.md.
    int newton_shift_diagonal{};
    int newton_shift_analytic{};

    // Measure c at every PTC step instead of using the fixed newton_shift.
    // Costs two extra metric interpolations per step (one unfrozen mat-vec plus
    // the interpolation that restores the frozen metric), which is roughly the
    // cost of one rejected step and vastly less than unfreezing the solve.
    int newton_shift_auto{};

    // Print a Rayleigh-quotient probe of the frozen and unfrozen Jacobians on a
    // few low-order modes at the top of every PTC step. Diagnostic only: it
    // does not change the trajectory, and costs a handful of extra metric
    // interpolations per step.
    int jacobian_diagnostic{};

    // How many times the line search may halve the step length alpha before
    // the step is rejected outright. 0 disables backtracking, recovering the
    // plain accept-or-reject behaviour. Each backtrack costs one
    // interpolate_metric + theta_from_metric, far less than the linear solve
    // that a rejection discards.
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

        int linear_max_iter = 1000;
        ah_finder_pp.queryAdd("linear_max_iter", linear_max_iter);
        if (linear_max_iter <= 0)
        {
            ah_finder_pp.error("linear_max_iter", "must be > 0");
        }

        int linear_restart_length = 1000;
        ah_finder_pp.queryAdd("linear_restart_length", linear_restart_length);
        if (linear_restart_length <= 0)
        {
            ah_finder_pp.error("linear_restart_length", "must be > 0");
        }

        amrex::Real fd_eps_scale = 1.0;
        ah_finder_pp.queryAdd("fd_eps_scale", fd_eps_scale);
        if (fd_eps_scale <= 0.0)
        {
            ah_finder_pp.error("fd_eps_scale", "must be > 0");
        }

        int unfreeze_jacobian = 0;
        ah_finder_pp.queryAdd("unfreeze_jacobian", unfreeze_jacobian);
        if (unfreeze_jacobian != 0 && unfreeze_jacobian != 1)
        {
            ah_finder_pp.error("unfreeze_jacobian", "must be 0 or 1");
        }

        int unfreeze_every = 0;
        ah_finder_pp.queryAdd("unfreeze_every", unfreeze_every);
        if (unfreeze_every < 0)
        {
            ah_finder_pp.error("unfreeze_every", "must be >= 0");
        }

        int unfreeze_on_reject = 0;
        ah_finder_pp.queryAdd("unfreeze_on_reject", unfreeze_on_reject);
        if (unfreeze_on_reject != 0 && unfreeze_on_reject != 1)
        {
            ah_finder_pp.error("unfreeze_on_reject", "must be 0 or 1");
        }

        amrex::Real newton_shift = 0.0;
        ah_finder_pp.queryAdd("newton_shift", newton_shift);
        if (newton_shift < 0.0)
        {
            ah_finder_pp.error("newton_shift", "must be >= 0");
        }

        int newton_shift_auto = 0;
        ah_finder_pp.queryAdd("newton_shift_auto", newton_shift_auto);
        if (newton_shift_auto != 0 && newton_shift_auto != 1)
        {
            ah_finder_pp.error("newton_shift_auto", "must be 0 or 1");
        }

        int newton_shift_diagonal = 0;
        ah_finder_pp.queryAdd("newton_shift_diagonal", newton_shift_diagonal);
        if (newton_shift_diagonal != 0 && newton_shift_diagonal != 1)
        {
            ah_finder_pp.error("newton_shift_diagonal", "must be 0 or 1");
        }
        // Compute c by transporting the frozen fields to the displaced radius
        // with their own derivatives, rather than re-interpolating the metric
        // at a perturbed surface and differencing. Cheaper (three partial
        // derivative queries and two Theta evaluations, against two full metric
        // interpolations) and free of the double finite difference the
        // interpolating route carries. Costs the second-derivative queries,
        // which are only registered when this is on.
        int newton_shift_analytic = 0;
        ah_finder_pp.queryAdd("newton_shift_analytic", newton_shift_analytic);
        if (newton_shift_analytic != 0 && newton_shift_analytic != 1)
        {
            ah_finder_pp.error("newton_shift_analytic", "must be 0 or 1");
        }
        if (newton_shift_analytic != 0 && newton_shift_auto == 0)
        {
            ah_finder_pp.error("newton_shift_analytic",
                               "requires newton_shift_auto = 1, which is what "
                               "triggers the measurement");
        }

        if (newton_shift_diagonal != 0 && newton_shift_auto == 0)
        {
            ah_finder_pp.error("newton_shift_diagonal",
                               "requires newton_shift_auto = 1, which is what "
                               "measures the diagonal");
        }

        int jacobian_diagnostic = 0;
        ah_finder_pp.queryAdd("jacobian_diagnostic", jacobian_diagnostic);
        if (jacobian_diagnostic != 0 && jacobian_diagnostic != 1)
        {
            ah_finder_pp.error("jacobian_diagnostic", "must be 0 or 1");
        }

        int max_backtracks = 4;
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
        ah_finder_pp.get("max_iter", max_iter);
        ah_finder_pp.get("linear_rel_tol", linear_rel_tol);
        ah_finder_pp.get("linear_abs_tol", linear_abs_tol);
        ah_finder_pp.get("linear_max_iter", linear_max_iter);
        ah_finder_pp.get("linear_restart_length", linear_restart_length);
        ah_finder_pp.get("fd_eps_scale", fd_eps_scale);
        ah_finder_pp.get("unfreeze_jacobian", unfreeze_jacobian);
        ah_finder_pp.get("unfreeze_every", unfreeze_every);
        ah_finder_pp.get("unfreeze_on_reject", unfreeze_on_reject);
        ah_finder_pp.get("newton_shift", newton_shift);
        ah_finder_pp.get("newton_shift_auto", newton_shift_auto);
        ah_finder_pp.get("newton_shift_diagonal", newton_shift_diagonal);
        ah_finder_pp.get("newton_shift_analytic", newton_shift_analytic);
        ah_finder_pp.get("jacobian_diagnostic", jacobian_diagnostic);
        ah_finder_pp.get("max_backtracks", max_backtracks);
    }
};

#endif /* AHFINDERPARAMETERS_HPP_ */
