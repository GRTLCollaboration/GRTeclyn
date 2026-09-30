#if !defined(AHFINDER_HPP_)
#error "This file should only be included through AHFinder.hpp"
#endif

#ifndef AHFINDER_IMPL_HPP_
#define AHFINDER_IMPL_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>

#include "CCZ4StateVariables.hpp"
#include "DefaultLevelBld.hpp"
#include "Derivative.hpp"
#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

template <int num_components>
void AHFinder<num_components>::init(GRAmr *gramr_ptr)
{
    // tolerance, r, max_iter and the linear-solve tolerances come
    // from the "ah_finder" scope of the input file; see AHFinderParameters.hpp
    // for defaults and meaning.
    m_params.fill_params();

    m_min_dt            = 1e-4;
    m_dt_shrink         = 0.5;
    m_dt_grow           = 1.25;
    m_theta_floor       = 1e-12;
    m_backtrack_factor  = 0.5;
    m_shift_gate_margin = 2.0;

    this->setup_metric_query();

    // Set up interpolator
    this->setup(gramr_ptr);

    m_geometry.set_surface_data(&m_state.h, &m_gamma_LL);

    // Matrix-free operator for the per-step linear solve, laid out on the same
    // ring grid the surface is discretised on.
    m_jac_op = std::make_unique<AHJacobianOp>(m_geometry.n_rings(),
                                              m_geometry.ring_size());

    // Initialise h for the particles, then place them
    this->init_particle_vals();
    this->set_particle_positions(m_state.h);
}

template <int num_components> void AHFinder<num_components>::find()
{
    // The Jacobian mat-vec closes over the current pseudo-timestep m_dt and the
    // frozen residual m_theta_n, both refreshed at the top of each PTC step.
    // For a direction "in" it returns (I/dt + J) in, with J dTheta/dh applied
    // as a finite-difference directional derivative of theta_from_metric about
    // the frozen surface m_state.h:
    //   J in ~= (Theta(h_n + eps*in) - Theta(h_n)) / eps
    // and eps set by the Brown-Saad rule. theta_from_metric() leaves the full
    // array on every rank, so "in"/"out" are complete flat vectors throughout.
    //
    // The shift multiplying "in" is 1/dt in the usual case. With
    // newton_shift_diagonal it is instead max(1/dt, c(x_ip)) particle by
    // particle: where the SER ramp has grown past the local c the operator
    // becomes J_frozen + diag(c) = J_exact there, and elsewhere it keeps the
    // heavier PTC damping. That is the pointwise version of capping dt at 1/c,
    // and coincides with it when c(x) is constant.
    m_jac_op->set_matvec(
        [this](const std::vector<double> &in, std::vector<double> &out)
        {
            std::vector<double> jv;
            this->jacobian_apply(in, m_unfreeze_now, jv);

            if (m_params.newton_shift_diagonal != 0 &&
                !m_newton_shift_diag.empty())
            {
                const double inv_dt = 1.0 / m_dt;
                for (int ip = 0; ip < m_num_particles; ++ip)
                {
                    out[ip] =
                        in[ip] * std::max(inv_dt, m_newton_shift_diag[ip]) +
                        jv[ip];
                }
            }
            else
            {
                for (int ip = 0; ip < m_num_particles; ++ip)
                    out[ip] = in[ip] / m_dt + jv[ip];
            }
        });

    int n_iter = 0;

    // Freeze the metric and residual at the initial surface. The loop keeps the
    // invariant that on entry m_gamma_LL and m_theta_n are frozen at m_state.h
    // and theta_old == inf_norm(m_theta_n).
    this->interpolate_metric(m_state.h);
    this->theta_from_metric(m_state.h, m_theta_n);

    double theta_old = inf_norm(m_theta_n);

    // Global pseudo-timestep
    m_dt = 2;

    const bool io_proc = amrex::ParallelDescriptor::IOProcessor();

    amrex::Print() << "\n AHFinder expansion Theta inf norm = " << theta_old
                   << "\n";

    std::ofstream theta_log;
    std::ofstream dt_log;

    if (io_proc)
    {
        theta_log.open("theta_vs_iter.csv");
        dt_log.open("dt_vs_iter.csv");
        std::filesystem::create_directory("particles");

        theta_log << n_iter << "," << theta_old << std::endl;
        dt_log << n_iter << "," << m_dt << std::endl;
    }

    auto write_particles = [&](int iter)
    {
        if (!io_proc)
            return;

        std::ofstream pfile("particles/particles_" + std::to_string(iter) +
                            ".csv");
        pfile << "x,y,z\n";
        for (int ip = 0; ip < m_num_particles; ++ip)
            pfile << interp_coords_x[ip] << "," << interp_coords_y[ip] << ","
                  << interp_coords_z[ip] << "\n";
    };
    write_particles(n_iter);

    // MultiFabs on the operator's ring-grid layout for the RHS -Theta(h_n) and
    // the solution increment delta_h.
    amrex::MultiFab rhs   = m_jac_op->make_mf();
    amrex::MultiFab delta = m_jac_op->make_mf();

    const double t_solve_start = amrex::second();

    // Whether the previous PTC step was rejected, for the unfreeze-on-reject
    // policy.
    bool last_rejected = false;

    while (theta_old > m_params.tolerance && n_iter < m_params.max_iter)
    {
        // Pick this step's Jacobian, once, before the solve: the exact
        // (unfrozen) one either unconditionally, on a fixed cadence, or only
        // after the cheap one has demonstrably failed. Held fixed for the
        // whole GMRES solve below so the Krylov iteration sees a consistent
        // linear operator.
        m_unfreeze_now = (m_params.unfreeze_jacobian != 0) ||
                         (m_params.unfreeze_on_reject != 0 && last_rejected) ||
                         (m_params.unfreeze_every > 0 &&
                          n_iter % m_params.unfreeze_every == 0);

        if (m_unfreeze_now)
            m_n_unfrozen_steps++;

        if (m_params.jacobian_diagnostic != 0)
            this->jacobian_diagnostic(n_iter);

        // Newton at frozen cost. Since J_exact ~= J_frozen + c I, the operator
        // (I/dt + J_frozen) *is* the exact Newton Jacobian when 1/dt = c. So
        // capping dt at 1/c turns the top of the SER ramp into a true Newton
        // step, computed entirely with the cheap frozen mat-vec and no metric
        // re-interpolation inside the solve; below the cap the step is a
        // Levenberg-Marquardt-damped Newton step, which is the globalisation
        // one wants anyway.
        //
        // The cap also keeps dt off the pole of (I/dt + J_frozen) at
        // dt = -1/lambda_min(J_frozen) ~ 0.125, which lies *above* 1/c ~ 0.084.
        // Uncapped, SER walks dt into that pole every few iterations and the
        // step is rejected; that is the accept/reject cycle documented in
        // docs/ah_finder_implicit_ptc.md.
        //
        // Measuring c costs two extra metric interpolations, so it is gated on
        // whether c can actually affect this step. Both uses of it -- the
        // scalar cap dt = min(dt, 1/c) and the diagonal max(1/dt, c(x)) -- only
        // bite once 1/dt has fallen to c; while 1/dt is above it the measured
        // value is computed and discarded. That masked phase is also the one in
        // which c drifts (13.5% over a typical solve, monotonically upwards),
        // whereas by the time it binds it has converged to ~1e-5, so reusing a
        // stale value costs nothing where it is used. m_shift_gate_margin is
        // the safety factor against c rising further while the gate is shut.
        if (m_params.newton_shift_auto != 0)
        {
            const bool no_valid_shift = (m_newton_shift <= 0.0);
            const bool shift_binds =
                (1.0 / m_dt) < m_shift_gate_margin * m_newton_shift;

            if (no_valid_shift || shift_binds)
            {
                m_newton_shift =
                    (m_params.newton_shift_analytic != 0)
                        ? this->measure_newton_shift_analytic(n_iter)
                        : this->measure_newton_shift(n_iter);
                m_n_shift_measured++;
            }
        }
        else
        {
            m_newton_shift = m_params.newton_shift;
        }

        // The scalar cap and the diagonal are alternatives, not cumulative: the
        // diagonal mat-vec already applies max(1/dt, c(x)) pointwise, so
        // capping dt as well would only raise the floor of that max and
        // re-impose the mean on every particle.
        if (m_newton_shift > 0.0 && m_params.newton_shift_diagonal == 0)
        {
            m_dt = std::min(m_dt, 1.0 / m_newton_shift);
        }

        // RHS = -Theta(h_n) (frozen at the current surface).
        std::vector<double> minus_theta(m_num_particles);
        for (int ip = 0; ip < m_num_particles; ++ip)
            minus_theta[ip] = -m_theta_n[ip];
        m_jac_op->flat_to_mf(minus_theta, rhs);

        // Solve (I/dt + J) delta_h = -Theta(h_n) with a single-level,
        // matrix-free GMRES (coarsening disabled in the operator). MLMG is
        // used only as the operator host; GMRESMLMG drives the Krylov
        // iteration through it.
        delta.setVal(0.0);
        amrex::MLMG mlmg(*m_jac_op);
        mlmg.setVerbose(0);
        mlmg.setBottomVerbose(0);

        amrex::GMRESMLMG gmres(mlmg);
        // Run unpreconditioned: there are no geometric multigrid levels to
        // precondition with, and Fsmooth (relaxation) is not implemented for
        // this matrix-free operator. With the preconditioner off, GMRES only
        // ever calls MLMG::applyPrecond (a plain mat-vec) -- the V-cycle and
        // the bottom solver are never entered.
        gmres.usePrecond(false);
        gmres.setVerbose(0);
        // As dt grows the frozen-metric operator (I/dt + J) becomes
        // ill-conditioned and GMRES stagnates. An inexact linear solve is
        // acceptable for PTC, so cap the Krylov work per step -- each GMRES
        // iteration costs a full Theta evaluation over the ring grid -- and
        // keep whatever increment the solver reached. A step that makes no
        // progress is caught by the rejection test below and answered by
        // shrinking dt, so no exception handling is needed here: unlike
        // MLMG::solve(), GMRES simply returns with a failure status.
        gmres.setMaxIters(m_params.linear_max_iter);
        gmres.getGMRES().setRestartLength(m_params.linear_restart_length);
        gmres.solve(delta, rhs, m_params.linear_rel_tol,
                    m_params.linear_abs_tol);

        std::vector<double> delta_h(m_num_particles);
        m_jac_op->mf_to_flat(delta, delta_h);

        // Step amplification: ||delta_h|| / (dt ||Theta||). If the I/dt term
        // dominated the operator this would be 1; for a definite elliptic J it
        // stays below 1. A value far above 1 means (I/dt + J) is close to
        // singular along some direction, which is the signature of an
        // indefinite J whose negative eigenvalues have been cancelled by the
        // 1/dt shift rather than regularised by it.
        double dh_norm = 0.0;
        double tn_norm = 0.0;
        for (int ip = 0; ip < m_num_particles; ++ip)
        {
            dh_norm += delta_h[ip] * delta_h[ip];
            tn_norm += m_theta_n[ip] * m_theta_n[ip];
        }
        dh_norm = std::sqrt(dh_norm);
        tn_norm = std::sqrt(tn_norm);

        const double amplification =
            (tn_norm > 0.0) ? dh_norm / (m_dt * tn_norm) : 0.0;

        // Backtracking line search along delta_h. The linear solve is by far
        // the most expensive part of the step, so rather than discarding it
        // when the full increment overshoots, retry the same direction at
        // shorter step lengths h_n + alpha * delta_h, halving alpha until the
        // residual falls. delta_h solves an *approximate* (frozen-metric)
        // Jacobian system, so near convergence the full step is regularly too
        // long even though the direction is good; without this the step is
        // thrown away and the whole solve repeated at a smaller dt.
        const std::vector<double> h_backup = m_state.h;

        double alpha     = 1.0;
        double theta_new = theta_old;
        bool accepted    = false;

        for (int ls = 0; ls <= m_params.max_backtracks; ++ls)
        {
            for (int ip = 0; ip < m_num_particles; ++ip)
                m_state.h[ip] = h_backup[ip] + alpha * delta_h[ip];

            // Evaluate Theta at the trial state, freezing its metric.
            // interpolate_metric() places the particles itself.
            this->interpolate_metric(m_state.h);
            this->theta_from_metric(m_state.h, m_theta_vals);

            theta_new = inf_norm(m_theta_vals);

            if (theta_new < theta_old)
            {
                accepted = true;
                break;
            }

            alpha *= m_backtrack_factor;
        }

        n_iter++;

        if (accepted)
        {
            // Accept: the frozen metric/residual now hold at the trial surface.
            m_theta_n = m_theta_vals;

            if (alpha < 1.0)
            {
                // The full step overshot, so the linear model is optimistic at
                // this dt: damp rather than grow. alpha is the measured
                // over-prediction factor, so reuse it as the shrink factor.
                m_dt = std::max(m_dt * alpha, m_min_dt);
            }
            else
            {
                // SER: grow dt as the residual falls, approaching Newton.
                m_dt = update_dt(m_dt, theta_old, theta_new);
            }

            theta_old     = theta_new;
            last_rejected = false;
        }
        else
        {
            // Reject: the residual rose for every alpha tried, and rose
            // proportionally to alpha, so delta_h is an *ascent* direction --
            // not merely an overshoot. Shortening the step cannot fix that;
            // only recomputing the direction from a more strongly regularised
            // (I/dt dominant) operator can. Restore h_n with its frozen
            // metric/residual and cut dt hard: shrinking by the gentle
            // m_dt_shrink instead makes the solver ratchet down over several
            // iterations, each costing a full GMRES solve, before it reaches a
            // dt whose direction is usable again.
            m_state.h = h_backup;
            this->interpolate_metric(m_state.h);
            this->theta_from_metric(m_state.h, m_theta_n);
            m_dt          = std::max(m_dt * m_dt_shrink, m_min_dt);
            last_rejected = true;
        }

        amrex::Print() << " AHFinder iter " << n_iter
                       << ": theta = " << theta_old << ", dt = " << m_dt
                       << ", t = " << amrex::second() - t_solve_start << " s"
                       << ", interp = " << m_n_interp
                       << ", theta_evals = " << m_n_theta
                       << ", |dh| = " << dh_norm
                       << ", ampl = " << amplification;

        if (m_newton_shift > 0.0)
        {
            amrex::Print() << ", shift = " << m_newton_shift;

            // The spread of c(x) about that mean is what decides whether the
            // scalar is an adequate summary: a few per cent on a near-spherical
            // surface, but comparable to the mean itself once the surface is
            // strongly distorted, at which point newton_shift_diagonal is
            // needed. Reported alongside the mean so the choice is visible.
            if (m_params.newton_shift_diagonal != 0 &&
                !m_newton_shift_diag.empty())
            {
                double sumsq = 0.0;
                for (int ip = 0; ip < m_num_particles; ++ip)
                {
                    const double d  = m_newton_shift_diag[ip] - m_newton_shift;
                    sumsq          += d * d;
                }

                amrex::Print() << " +- " << std::sqrt(sumsq / m_num_particles)
                               << " (diag)";
            }
        }

        if (m_params.jacobian_diagnostic != 0)
        {
            // The convergence test uses the inf-norm, which is set by a single
            // particle and so says nothing about whether the rest of the
            // surface is converging at the same rate. Print the rms alongside
            // it, and where the inf-norm currently lives: if the two rates
            // differ, the observed contraction is a property of the merit
            // function rather than of the step.
            double sumsq = 0.0;
            int arg_ip   = 0;
            for (int ip = 0; ip < m_num_particles; ++ip)
            {
                sumsq += m_theta_n[ip] * m_theta_n[ip];
                if (std::abs(m_theta_n[ip]) > std::abs(m_theta_n[arg_ip]))
                    arg_ip = ip;
            }

            amrex::Print() << ", rms = " << std::sqrt(sumsq / m_num_particles)
                           << ", argmax = " << arg_ip << " (ring "
                           << arg_ip / m_geometry.ring_size() << " of "
                           << m_geometry.n_rings() << ")";
        }

        amrex::Print() << "\n";

        if (io_proc)
        {
            theta_log << n_iter << "," << theta_old << std::endl;
            dt_log << n_iter << "," << m_dt << std::endl;
        }
        write_particles(n_iter);
    }

    if (io_proc)
    {
        theta_log.close();
        dt_log.close();
    }

    const double t_solve = amrex::second() - t_solve_start;

    amrex::AllPrint() << "\n AHFinder converged with inf norm of theta = "
                      << theta_old << " in " << n_iter << " iterations\n";

    // Name the Jacobian policy that was in force, so a log can be attributed to
    // a run configuration without going back to the input file.
    std::string jac_policy = "frozen";
    if (m_params.unfreeze_jacobian != 0)
    {
        jac_policy = "unfrozen";
    }
    else if (m_params.unfreeze_every > 0 || m_params.unfreeze_on_reject != 0)
    {
        jac_policy = "hybrid";
        if (m_params.unfreeze_every > 0)
        {
            jac_policy += ", every " + std::to_string(m_params.unfreeze_every);
        }
        if (m_params.unfreeze_on_reject != 0)
        {
            jac_policy += ", on reject";
        }
    }

    amrex::Print() << " AHFinder timing: solve = " << t_solve << " s"
                   << " (jacobian " << jac_policy << "; " << m_n_unfrozen_steps
                   << " of " << n_iter << " steps exact; shift measured on "
                   << m_n_shift_measured << " of " << n_iter << " steps)\n"
                   << "   interpolate_metric: " << m_n_interp << " calls, "
                   << m_t_interp << " s (" << 100.0 * m_t_interp / t_solve
                   << "%)\n"
                   << "   theta_from_metric:  " << m_n_theta << " calls, "
                   << m_t_theta << " s (" << 100.0 * m_t_theta / t_solve
                   << "%)\n";

    // Report the converged surface's area and irreducible mass
    // (Christodoulou formula: M = sqrt(A / 16 pi)).
    const amrex::Real area = m_geometry.area();
    amrex::AllPrint() << " AHFinder surface area = " << area << "\n";

    const amrex::Real mass = std::sqrt(area / (16.0 * M_PI));
    amrex::AllPrint() << " AHFinder irreducible mass = " << mass << "\n";
}

template <int num_components>
void AHFinder<num_components>::jacobian_apply(const std::vector<double> &v,
                                              bool unfreeze,
                                              std::vector<double> &jv)
{
    const double macheps = std::numeric_limits<double>::epsilon();

    double h_norm = 0.0;
    double v_norm = 0.0;
    for (int ip = 0; ip < m_num_particles; ++ip)
    {
        h_norm += m_state.h[ip] * m_state.h[ip];
        v_norm += v[ip] * v[ip];
    }
    h_norm = std::sqrt(h_norm);
    v_norm = std::sqrt(v_norm);

    jv.assign(m_num_particles, 0.0);

    // If the direction is (numerically) zero, J v = 0.
    if (v_norm == 0.0)
    {
        return;
    }

    const double eps =
        m_params.fd_eps_scale * std::sqrt(macheps) * (1.0 + h_norm) / v_norm;

    std::vector<double> h_pert(m_num_particles);
    for (int ip = 0; ip < m_num_particles; ++ip)
        h_pert[ip] = m_state.h[ip] + eps * v[ip];

    std::vector<double> theta_pert(m_num_particles);

    if (unfreeze)
    {
        // Exact Newton Jacobian: re-interpolate the metric at the perturbed
        // surface so the directional derivative picks up the d(metric)/dh term
        // too. This costs a full particle re-query per apply rather than per
        // PTC step, and it leaves m_gamma_LL interpolated at h_pert -- every
        // consumer outside the mat-vec (the line search, the reject branch and
        // jacobian_diagnostic()) calls interpolate_metric() itself before using
        // it, so the frozen invariant is re-established there rather than here.
        this->interpolate_metric(h_pert);
    }

    this->theta_from_metric(h_pert, theta_pert);

    for (int ip = 0; ip < m_num_particles; ++ip)
        jv[ip] = (theta_pert[ip] - m_theta_n[ip]) / eps;
}

template <int num_components>
void AHFinder<num_components>::jacobian_diagnostic(int n_iter)
{
    // Probe both Jacobians at the *same* surface with the same test modes, so
    // frozen and unfrozen are compared as operators rather than through the
    // trajectories they produce.
    //
    // Unfreezing adds a zeroth-order (multiplicative) term
    // (dTheta/dgamma) (dir^k d_k gamma) to what is otherwise the elliptic
    // principal part in h. What this probe established (see
    // docs/ah_finder_implicit_ptc.md) is that the term is, to within 0.2%
    // across modes whose own eigenvalues span a factor of 200, a constant
    // multiple of the identity: J_exact ~= J_frozen + c I with c ~= 12. It is
    // the *frozen* operator that is indefinite -- its l=0 eigenvalue is about
    // -8, so (I/dt + J_frozen) has a pole at dt ~ 0.125, which is what limits
    // dt and what the step-amplification column in the iteration log tracks.
    //
    // The Rayleigh quotient q = <v, J v> / <v, v> is a one-mat-vec estimate of
    // the eigenvalue J sees along v (J is not symmetric, so q is indicative,
    // not an eigenvalue). Modes are chosen to span the low end of the spectrum,
    // where the elliptic part is weakest and the potential can dominate: the
    // uniform l=0 breathing mode, and l=1/l=2 modes aligned with the binary
    // axis (x), which is where the surface comes closest to the punctures.
    const char *names[3] = {"l0    ", "l1_x  ", "l2_x  "};
    std::array<std::vector<double>, 3> modes;

    for (auto &m : modes)
        m.resize(m_num_particles);

    for (int ip = 0; ip < m_num_particles; ++ip)
    {
        const Tensor::Rank1 dir = m_geometry.direction(ip);
        modes[0][ip]            = 1.0;
        modes[1][ip]            = dir(0);
        modes[2][ip]            = 3.0 * dir(0) * dir(0) - 1.0;
    }

    amrex::Print() << " AHFinder jacobian probe, iter " << n_iter + 1
                   << ", dt = " << m_dt << "\n";

    std::vector<double> jv;
    for (int im = 0; im < 3; ++im)
    {
        const std::vector<double> &v = modes[im];

        double vv = 0.0;
        for (int ip = 0; ip < m_num_particles; ++ip)
            vv += v[ip] * v[ip];

        double q[2];
        for (int variant = 0; variant < 2; ++variant)
        {
            this->jacobian_apply(v, variant == 1, jv);

            double vjv = 0.0;
            for (int ip = 0; ip < m_num_particles; ++ip)
                vjv += v[ip] * jv[ip];

            q[variant] = vjv / vv;
        }

        // dt at which (I/dt + J) is singular along this mode, if q < 0.
        const double dt_sing = (q[1] < 0.0) ? -1.0 / q[1] : -1.0;

        amrex::Print() << "   mode " << names[im] << " q_frozen = " << q[0]
                       << ", q_exact = " << q[1]
                       << ", dt_singular(exact) = " << dt_sing << "\n";
    }

    // jacobian_apply() with unfreeze left the metric at the last perturbed
    // surface; restore the loop invariant.
    this->interpolate_metric(m_state.h);
}

template <int num_components>
double AHFinder<num_components>::measure_newton_shift(int n_iter)
{
    // The exact Jacobian differs from the frozen one by a term that is, to
    // within 0.2% across modes whose own eigenvalues span a factor of 200, a
    // constant multiple of the identity:
    //
    //   J_exact ~= J_frozen + c I
    //
    // (see docs/ah_finder_implicit_ptc.md). One mode therefore suffices to
    // estimate c; the uniform l=0 breathing mode is used because it is the mode
    // the PTC step is overwhelmingly built from. Two Rayleigh quotients, so two
    // metric interpolations: one inside the unfrozen apply, one to restore the
    // frozen metric at m_state.h afterwards.
    const std::vector<double> v(m_num_particles, 1.0);
    std::vector<double> jv[2];

    double q[2];
    for (int variant = 0; variant < 2; ++variant)
    {
        this->jacobian_apply(v, variant == 1, jv[variant]);

        double vjv = 0.0;
        for (int ip = 0; ip < m_num_particles; ++ip)
            vjv += jv[variant][ip];

        q[variant] = vjv / static_cast<double>(m_num_particles);
    }

    // The difference of the two applies is the dropped term evaluated on the
    // uniform mode, particle by particle: since that term is zeroth order in
    // delta_h it is a multiplication operator, so
    //   c(x_ip) = (J_exact - J_frozen)_ip,ip = jv_unfrozen[ip] - jv_frozen[ip].
    // The scalar shift returned below is its mean; newton_shift_diagonal makes
    // the mat-vec use this array directly instead, which is worth doing exactly
    // when the surface is distorted enough for the two to differ.
    m_newton_shift_diag.resize(m_num_particles);
    for (int ip = 0; ip < m_num_particles; ++ip)
        m_newton_shift_diag[ip] = jv[1][ip] - jv[0][ip];

    // Whether the diagonal is genuinely constant or only looks constant to the
    // low-order probe modes is what this dump answers -- see
    // docs/ah_finder_implicit_ptc.md.
    if (m_params.jacobian_diagnostic != 0 &&
        amrex::ParallelDescriptor::IOProcessor())
    {
        std::ofstream pfile("shift_profile_" + std::to_string(n_iter) + ".csv");
        pfile << "ip,ring,dir_x,dir_y,dir_z,h,c\n";
        for (int ip = 0; ip < m_num_particles; ++ip)
        {
            const Tensor::Rank1 dir = m_geometry.direction(ip);
            pfile << ip << "," << ip / m_geometry.ring_size() << "," << dir(0)
                  << "," << dir(1) << "," << dir(2) << "," << m_state.h[ip]
                  << "," << jv[1][ip] - jv[0][ip] << "\n";
        }
    }

    this->interpolate_metric(m_state.h);

    return q[1] - q[0];
}

template <int num_components>
void AHFinder<num_components>::interpolate_shift_data()
{
    // The extra derivative data the analytic shift needs, at the particle
    // positions the last interpolate_metric() already placed: no
    // set_particle_positions(), and no refresh flag on the queries, because the
    // surface has not moved. find() guarantees that -- its loop invariant is
    // that the metric is frozen at m_state.h on entry to each step, and
    // jacobian_diagnostic() restores it if it ran.
    //
    // Counted as an interpolation in the cost report even though it is three
    // partial queries rather than a full call, so the reported interpolation
    // time stays honest.
    const double t_start = amrex::second();
    m_n_interp++;

    this->interp(m_metric_query_deriv2, false);
    this->interp(m_metric_query_d2a, false);
    this->interp(m_metric_query_d2b, false);

    m_t_interp += amrex::second() - t_start;
}

template <int num_components>
double AHFinder<num_components>::measure_newton_shift_analytic(int n_iter)
{
    // The same c(x) that measure_newton_shift() finite-differences, but reached
    // without moving the surface.
    //
    // Displacing particle ip by delta_h moves its sample point along
    // dir(ip), so every field Phi that theta_from_metric() reads changes by
    // delta_h * dir^k d_k Phi. Rather than re-interpolating at the displaced
    // surface, transport the frozen arrays by that amount and re-evaluate. The
    // chain rule is applied exactly -- there is no finite difference in h, and
    // so none of the double-FD error the interpolating version carries. What
    // remains is a central difference in *field* space, where Theta is a smooth
    // algebraic function, taken at the cbrt(macheps) step that is optimal for
    // it.
    //
    // theta_from_metric() reads 14 state values and 21 first-derivative values,
    // so both have to be transported: the values need d_k Phi and the
    // derivatives need d_j d_k Phi.
    this->interpolate_shift_data();

    const double macheps = std::numeric_limits<double>::epsilon();

    std::array<std::vector<double>, 14> dr_state;
    for (auto &v : dr_state)
        v.assign(m_num_particles, 0.0);

    std::array<std::array<std::vector<double>, 7>, 3> dr_d1;
    for (auto &a : dr_d1)
        for (auto &v : a)
            v.assign(m_num_particles, 0.0);

    double phi_sq = 0.0;
    double dr_sq  = 0.0;

    for (int ip = m_start; ip < m_start + m_n_local; ++ip)
    {
        const Tensor::Rank1 n_L = m_geometry.direction(ip);

        for (int c = 0; c < 14; ++c)
        {
            const double d = n_L(0) * m_metric_dx[c][ip] +
                             n_L(1) * m_metric_dy[c][ip] +
                             n_L(2) * m_metric_dz[c][ip];
            dr_state[c][ip] = d;

            phi_sq += m_metric_state[c][ip] * m_metric_state[c][ip];
            dr_sq  += d * d;
        }

        for (int c = 0; c < 7; ++c)
        {
            FOR (j)
            {
                double d = 0.0;
                FOR (k)
                {
                    d += n_L(k) * m_metric_d2[sym2_idx(j, k)][c][ip];
                }
                dr_d1[j][c][ip]  = d;
                dr_sq           += d * d;
            }

            phi_sq += m_metric_dx[c][ip] * m_metric_dx[c][ip] +
                      m_metric_dy[c][ip] * m_metric_dy[c][ip] +
                      m_metric_dz[c][ip] * m_metric_dz[c][ip];
        }
    }

    // The step has to be identical on every rank, so the norms are reduced
    // before it is formed.
    amrex::ParallelDescriptor::ReduceRealSum(phi_sq);
    amrex::ParallelDescriptor::ReduceRealSum(dr_sq);

    if (!(dr_sq > 0.0))
    {
        // A field with no radial variation contributes no shift.
        m_newton_shift_diag.assign(m_num_particles, 0.0);
        return 0.0;
    }

    // cbrt rather than sqrt: this is a central difference, whose truncation
    // error is O(eps^2) against a roundoff floor of O(macheps/eps).
    const double eps =
        std::cbrt(macheps) * (1.0 + std::sqrt(phi_sq)) / std::sqrt(dr_sq);

    const auto base_state = m_metric_state;
    const auto base_dx    = m_metric_dx;
    const auto base_dy    = m_metric_dy;
    const auto base_dz    = m_metric_dz;

    // Set every field to base + s * (its radial derivative). s = 0 restores the
    // frozen state exactly, so the surface is left as it was found.
    auto transport = [&](double s)
    {
        for (int ip = m_start; ip < m_start + m_n_local; ++ip)
        {
            for (int c = 0; c < 14; ++c)
            {
                m_metric_state[c][ip] = base_state[c][ip] + s * dr_state[c][ip];
            }
            for (int c = 0; c < 7; ++c)
            {
                m_metric_dx[c][ip] = base_dx[c][ip] + s * dr_d1[0][c][ip];
                m_metric_dy[c][ip] = base_dy[c][ip] + s * dr_d1[1][c][ip];
                m_metric_dz[c][ip] = base_dz[c][ip] + s * dr_d1[2][c][ip];
            }

            // gamma_ij = h_ij / chi is derived, not interpolated, so it has to
            // be rebuilt from the transported values rather than transported
            // itself. Only this rank's slice is touched; theta_from_metric()
            // reads no other.
            const double chi = m_metric_state[c_chi][ip];
            FOR (i, j)
            {
                m_gamma_LL[ip](i, j) =
                    m_metric_state[sym_var_idx(c_h11, i, j)][ip] / chi;
            }
        }
    };

    std::vector<double> theta_p;
    std::vector<double> theta_m;

    transport(+eps);
    this->theta_from_metric(m_state.h, theta_p);
    transport(-eps);
    this->theta_from_metric(m_state.h, theta_m);
    transport(0.0);

    m_newton_shift_diag.resize(m_num_particles);

    double sum = 0.0;
    for (int ip = 0; ip < m_num_particles; ++ip)
    {
        m_newton_shift_diag[ip]  = (theta_p[ip] - theta_m[ip]) / (2.0 * eps);
        sum                     += m_newton_shift_diag[ip];
    }

    if (m_params.jacobian_diagnostic != 0 &&
        amrex::ParallelDescriptor::IOProcessor())
    {
        std::ofstream pfile("shift_profile_" + std::to_string(n_iter) + ".csv");
        pfile << "ip,ring,dir_x,dir_y,dir_z,h,c\n";
        for (int ip = 0; ip < m_num_particles; ++ip)
        {
            const Tensor::Rank1 dir = m_geometry.direction(ip);
            pfile << ip << "," << ip / m_geometry.ring_size() << "," << dir(0)
                  << "," << dir(1) << "," << dir(2) << "," << m_state.h[ip]
                  << "," << m_newton_shift_diag[ip] << "\n";
        }
    }

    return sum / static_cast<double>(m_num_particles);
}

template <int num_components>
void AHFinder<num_components>::init_particle_vals()
{
    const double r0 = m_geometry.guess_radius();

    // The initial surface is the sphere r = guess_radius
    m_state.h.assign(m_num_particles, r0);
}

template <int num_components>
void AHFinder<num_components>::set_particle_positions(
    const std::vector<double> &h)
{
    const std::array<double, AMREX_SPACEDIM> &center = m_geometry.center();

    // Set particles' positions based on their radius from the centre
    for (int id = 0; id < m_num_particles; ++id)
    {
        const double r          = h[id];
        const Tensor::Rank1 dir = m_geometry.direction(id);

        interp_coords_x[id] = center[0] + r * dir(0);
        interp_coords_y[id] = center[1] + r * dir(1);
        interp_coords_z[id] = center[2] + r * dir(2);

        amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> coords = {
            interp_coords_x[id], interp_coords_y[id], interp_coords_z[id]};
        this->check_domain(coords);
    }
}

template <int num_components>
amrex::Real AHFinder<num_components>::update_dt(amrex::Real dt,
                                                double theta_old,
                                                double theta_new) const
{
    // SER: grow dt as the residual falls so the implicit step approaches
    // Newton. The implicit solve (I/dt + J) is unconditionally stable, so dt
    // is not capped; robustness against an over-large step is provided by the
    // step-rejection safeguard in find(), which shrinks dt on any step that
    // increases the residual.

    // Scale dt by the ratio of old to new theta, clamped so it can't grow or
    // shrink too fast in a single step.
    double ratio = (std::abs(theta_new) > m_theta_floor)
                       ? m_params.r * theta_old / theta_new
                       : m_params.r;
    ratio        = std::max(ratio, m_dt_shrink);
    ratio        = std::min(ratio, m_dt_grow);

    dt *= ratio;

    // Keep dt above the floor.
    dt = std::max(dt, m_min_dt);

    return dt;
}

template <int num_components>
double AHFinder<num_components>::inf_norm(std::vector<double> arr)
{
    double max_el = std::abs(arr[0]);
    for (auto &&i : arr)
    {
        if (std::abs(i) > max_el)
            max_el = std::abs(i);
    }

    return max_el;
}

template <int num_components>
void AHFinder<num_components>::setup_metric_query()
{
    for (auto &v : m_metric_state)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dx)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dy)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dz)
        v.resize(m_num_particles);
    for (auto &pair : m_metric_d2)
        for (auto &v : pair)
            v.resize(m_num_particles);

    m_metric_query_state.setCoords(0, interp_coords_x.data() + m_start)
        .setCoords(1, interp_coords_y.data() + m_start)
        .setCoords(2, interp_coords_z.data() + m_start);
    m_metric_query_deriv.setCoords(0, interp_coords_x.data() + m_start)
        .setCoords(1, interp_coords_y.data() + m_start)
        .setCoords(2, interp_coords_z.data() + m_start);
    for (auto *q :
         {&m_metric_query_deriv2, &m_metric_query_d2a, &m_metric_query_d2b})
    {
        q->setCoords(0, interp_coords_x.data() + m_start)
            .setCoords(1, interp_coords_y.data() + m_start)
            .setCoords(2, interp_coords_z.data() + m_start);
    }

    // chi, h_ij, K, A_ij (values only).
    m_metric_query_state.addComp(c_chi, m_metric_state[c_chi].data() + m_start,
                                 VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query_state.addComp(
            comp, m_metric_state[comp].data() + m_start, VariableType::state);
    }
    m_metric_query_state.addComp(c_K, m_metric_state[c_K].data() + m_start,
                                 VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_A11, i, j);
        m_metric_query_state.addComp(
            comp, m_metric_state[comp].data() + m_start, VariableType::state);
    }

    m_metric_query_deriv.addComp(c_chi, m_metric_dx[c_chi].data() + m_start,
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dx);
    m_metric_query_deriv.addComp(c_chi, m_metric_dy[c_chi].data() + m_start,
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dy);
    m_metric_query_deriv.addComp(c_chi, m_metric_dz[c_chi].data() + m_start,
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dz);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query_deriv.addComp(comp, m_metric_dx[comp].data() + m_start,
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dx);
        m_metric_query_deriv.addComp(comp, m_metric_dy[comp].data() + m_start,
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dy);
        m_metric_query_deriv.addComp(comp, m_metric_dz[comp].data() + m_start,
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dz);
    }

    if (m_params.newton_shift_analytic == 0)
    {
        return;
    }

    // d_k of K and A_ij. theta_from_metric() does not read these, but the
    // analytic shift does: transporting K and A_ij to the displaced radius
    // needs their gradients just as much as chi and h_ij need theirs.
    const std::array<Derivative, 3> d1{Derivative::dx, Derivative::dy,
                                       Derivative::dz};
    std::array<std::vector<double> *, 3> d1_out{
        m_metric_dx.data(), m_metric_dy.data(), m_metric_dz.data()};

    for (int k = 0; k < 3; ++k)
    {
        m_metric_query_deriv2.addComp(c_K, d1_out[k][c_K].data() + m_start,
                                      VariableType::state, BCParity::undefined,
                                      d1[k]);
        FOR2_SYM(i, j)
        {
            int comp = sym_var_idx(c_A11, i, j);
            m_metric_query_deriv2.addComp(
                comp, d1_out[k][comp].data() + m_start, VariableType::state,
                BCParity::undefined, d1[k]);
        }
    }

    // d_j d_k of chi and h_ij, six symmetric pairs over seven components = 42
    // entries, split into two queries because a single query's entry count
    // cannot exceed num_components.
    const std::array<Derivative, 6> d2{Derivative::dxdx, Derivative::dydy,
                                       Derivative::dzdz, Derivative::dxdy,
                                       Derivative::dxdz, Derivative::dydz};

    for (int p = 0; p < 6; ++p)
    {
        InterpolationQueryParticle &q =
            (p < 3) ? m_metric_query_d2a : m_metric_query_d2b;

        q.addComp(c_chi, m_metric_d2[p][c_chi].data() + m_start,
                  VariableType::state, BCParity::undefined, d2[p]);
        FOR2_SYM(i, j)
        {
            int comp = sym_var_idx(c_h11, i, j);
            q.addComp(comp, m_metric_d2[p][comp].data() + m_start,
                      VariableType::state, BCParity::undefined, d2[p]);
        }
    }
}

template <int num_components>
void AHFinder<num_components>::interpolate_metric(const std::vector<double> &h)
{
    // Place particles on the surface r = h and interpolate the CCZ4 metric
    // (chi, h_ij, K, A_ij) and its first derivatives onto them. Fills
    // m_metric_state / m_metric_dx/dy/dz (valid on this rank's local slice)
    // and the physical 3-metric m_gamma_LL (reduced to the full grid). These
    // are held fixed while theta_from_metric() varies h, so the Jacobian used
    // by the PTC solve is evaluated with the metric frozen at this surface.
    const double t_start = amrex::second();
    m_n_interp++;

    this->set_particle_positions(h);

    m_gamma_LL.assign(m_num_particles, Tensor::Rank2{0.0});

    this->interp(m_metric_query_state, true);
    this->interp(m_metric_query_deriv, false);

    for (int ip = m_start; ip < m_start + m_n_local; ++ip)
    {
        const double chi = m_metric_state[c_chi][ip];

        // gamma_ij = h_ij / chi.
        FOR (i, j)
        {
            const int h_comp     = sym_var_idx(c_h11, i, j);
            const double h_ij    = m_metric_state[h_comp][ip];
            m_gamma_LL[ip](i, j) = h_ij / chi;
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(
        reinterpret_cast<amrex::Real *>(m_gamma_LL.data()),
        m_num_particles * AMREX_SPACEDIM * AMREX_SPACEDIM);

    m_t_interp += amrex::second() - t_start;
}

template <int num_components>
void AHFinder<num_components>::theta_from_metric(const std::vector<double> &h,
                                                 std::vector<double> &theta_out)
{
    // Evaluate the expansion Theta at every grid point for the surface radius
    // h, using the metric frozen by the last interpolate_metric(). Only h and
    // its ring-grid derivatives enter, so this is cheap (local work plus one
    // reduce) and is what the matrix-free Jacobian differentiates.
    const double t_start = amrex::second();
    m_n_theta++;

    m_geometry.set_h_derivatives(h);

    theta_out.assign(m_num_particles, 0.0);

    amrex::GpuArray<const double *, 14> state_ptr;
    for (int c = 0; c < 14; ++c)
        state_ptr[c] = m_metric_state[c].data();
    amrex::GpuArray<const double *, 7> dx_ptr, dy_ptr, dz_ptr;
    for (int c = 0; c < 7; ++c)
    {
        dx_ptr[c] = m_metric_dx[c].data();
        dy_ptr[c] = m_metric_dy[c].data();
        dz_ptr[c] = m_metric_dz[c].data();
    }
    amrex::GpuArray<amrex::GpuArray<const double *, 7>, 3> d1_metric_ptr{
        dx_ptr, dy_ptr, dz_ptr};

    const double *h_ptr = h.data();

    double *theta_ptr = theta_out.data();

    for (int ip = m_start; ip < m_start + m_n_local; ++ip)
    {
        using namespace TensorAlgebra;

        double r   = h_ptr[ip];
        double chi = state_ptr[c_chi][ip];
        double K   = state_ptr[c_K][ip];

        // Physical 3-metric gamma_ij is frozen (from interpolate_metric); the
        // extrinsic curvature K_ij = (A_ij + (1/3) h_ij K)/chi is rebuilt here
        // from the frozen interpolated values.
        const Tensor::Rank2 &gamma_LL = m_gamma_LL[ip];
        Tensor::Rank2 K_LL;
        FOR (i, j)
        {
            int h_comp  = sym_var_idx(c_h11, i, j);
            int A_comp  = sym_var_idx(c_A11, i, j);
            double h_ij = state_ptr[h_comp][ip];
            double A_ij = state_ptr[A_comp][ip];
            K_LL(i, j)  = (A_ij + (1.0 / 3.0) * h_ij * K) / chi;
        }
        Tensor::Rank2 gamma_UU = compute_inverse(gamma_LL);

        // d_k(gamma_ij) from d1(chi), d1(h_ij) (product rule on
        // gamma_ij = h_ij/chi).
        Tensor::Rank1 d1_chi;
        FOR (k)
        {
            d1_chi(k) = d1_metric_ptr[k][c_chi][ip];
        }

        Tensor::Rank3 d1_gamma_LL; // (k, i, j) = d_k gamma_ij
        FOR (k, i, j)
        {
            int h_comp     = sym_var_idx(c_h11, i, j);
            double d1h_kij = d1_metric_ptr[k][h_comp][ip];
            d1_gamma_LL(k, i, j) =
                d1h_kij / chi - gamma_LL(i, j) * d1_chi(k) / chi;
        }

        // d_k(gamma^ij) = -gamma^im gamma^jn d_k(gamma_mn).
        Tensor::Rank3 d1_gamma_UU;
        FOR (k, i, j)
        {
            d1_gamma_UU(k, i, j) = 0.0;
            FOR (m, n)
            {
                d1_gamma_UU(k, i, j) -=
                    gamma_UU(i, m) * gamma_UU(j, n) * d1_gamma_LL(k, m, n);
            }
        }

        // V^j = sum_k d_k(gamma^kj) (divergence of the inverse
        // metric).
        Tensor::Rank1 V_U;
        FOR (j)
        {
            V_U(j) = 0.0;
            FOR (k)
            {
                V_U(j) += d1_gamma_UU(k, k, j);
            }
        }

        // Level-set gradient F_i = n_i - (grad h)_i, with n_i the flat
        // radial covector (x_i - center_i)/r -- which is exactly the
        // grid point's unit direction, since set_particle_positions()
        // built x_i as center_i + r * dir_i from this same h.
        const Tensor::Rank1 n_L      = m_geometry.direction(ip);
        const Tensor::Rank1 grad_h_L = m_geometry.grad_h(ip);
        Tensor::Rank1 F_L;
        FOR (i)
        {
            F_L(i) = n_L(i) - grad_h_L(i);
        }

        // Hess(F) = Hess(r) - Hess(h), with
        // Hess(r)_ij = (delta_ij - n_i n_j)/r (flat identity).
        const Tensor::Rank2 hess_h_LL = m_geometry.hess_h(ip);

        Tensor::Rank2 hess_F_LL;
        FOR (i, j)
        {
            hess_F_LL(i, j) =
                (delta(i, j) - n_L(i) * n_L(j)) / r - hess_h_LL(i, j);
        }

        // lambda = sqrt(gamma^ij F_i F_j); s_i = F_i/lambda;
        // s^i = gamma^ij s_j.
        amrex::Real lambda_sq = compute_dot_product(F_L, F_L, gamma_UU);
        amrex::Real lambda    = std::sqrt(lambda_sq);

        Tensor::Rank1 s_L;
        FOR (i)
        {
            s_L(i) = F_L(i) / lambda;
        }
        Tensor::Rank1 s_U = raise_all(s_L, gamma_UU);

        // d_k(lambda), from differentiating
        // lambda^2 = gamma^mn F_m F_n.
        Tensor::Rank1 d_lambda;
        FOR (k)
        {
            amrex::Real term1 = 0.0;
            FOR (m, n)
            {
                term1 += d1_gamma_UU(k, m, n) * F_L(m) * F_L(n);
            }
            amrex::Real term2 = 0.0;
            FOR (m, n)
            {
                term2 += gamma_UU(m, n) * hess_F_LL(m, k) * F_L(n);
            }
            d_lambda(k) = (term1 + 2.0 * term2) / (2.0 * lambda);
        }

        // d_i(ln sqrt(gamma)) = (1/2) gamma^jk d_i(gamma_jk)
        // (Jacobi's formula).
        Tensor::Rank1 d_ln_sqrt_gamma;
        FOR (k)
        {
            d_ln_sqrt_gamma(k) = 0.0;
            FOR (i, j)
            {
                d_ln_sqrt_gamma(k) +=
                    0.5 * gamma_UU(i, j) * d1_gamma_LL(k, i, j);
            }
        }

        // d_i s^i = (1/lambda)[V^j F_j + gamma^ij Hess(F)_ij
        //           - (d_i lambda) s^i]
        amrex::Real V_dot_F      = compute_dot_product(V_U, F_L);
        amrex::Real trace_hess_F = compute_trace(hess_F_LL, gamma_UU);
        amrex::Real lambda_dot_s = compute_dot_product(s_U, d_lambda);
        amrex::Real div_s = (V_dot_F + trace_hess_F - lambda_dot_s) / lambda;

        amrex::Real s_dot_dlnsqrtgamma =
            compute_dot_product(s_U, d_ln_sqrt_gamma);

        amrex::Real s_K_s = 0.0;
        FOR (i, j)
        {
            s_K_s += s_U(i) * s_U(j) * K_LL(i, j);
        }

        // Theta = D_i s^i - K + s^i s^j K_ij
        theta_ptr[ip] = div_s + s_dot_dlnsqrtgamma - K + s_K_s;
    }

    amrex::ParallelDescriptor::ReduceRealSum(theta_out.data(), m_num_particles);

    m_t_theta += amrex::second() - t_start;
}

#endif /* AHFINDER_IMPL_HPP_ */
