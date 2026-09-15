#if !defined(AHFINDER_HPP_)
#error "This file should only be included through AHFinder.hpp"
#endif

#ifndef AHFINDER_IMPL_HPP_
#define AHFINDER_IMPL_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

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

template <int num_components>
void AHFinder<num_components>::init(GRAmr *gramr_ptr)
{
    // tolerance, r, max_iter and the linear-solve tolerances come
    // from the "ah_finder" scope of the input file; see AHFinderParameters.hpp
    // for defaults and meaning.
    m_params.fill_params();

    m_min_dt      = 1e-4;
    m_dt_shrink   = 0.8;
    m_dt_grow     = 1.25;
    m_theta_floor = 1e-12;

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
    m_jac_op->set_matvec(
        [this](const std::vector<double> &in, std::vector<double> &out)
        {
            const double macheps = std::numeric_limits<double>::epsilon();

            double h_norm = 0.0;
            double v_norm = 0.0;
            for (int ip = 0; ip < m_num_particles; ++ip)
            {
                h_norm += m_state.h[ip] * m_state.h[ip];
                v_norm += in[ip] * in[ip];
            }
            h_norm = std::sqrt(h_norm);
            v_norm = std::sqrt(v_norm);

            // If the direction is (numerically) zero, J in = 0 and the apply
            // reduces to the I/dt term.
            if (v_norm == 0.0)
            {
                for (int ip = 0; ip < m_num_particles; ++ip)
                    out[ip] = in[ip] / m_dt;
                return;
            }

            const double eps =
                std::sqrt(macheps) * (1.0 + h_norm) / v_norm;

            std::vector<double> h_pert(m_num_particles);
            for (int ip = 0; ip < m_num_particles; ++ip)
                h_pert[ip] = m_state.h[ip] + eps * in[ip];

            std::vector<double> theta_pert(m_num_particles);
            this->theta_from_metric(h_pert, theta_pert);

            for (int ip = 0; ip < m_num_particles; ++ip)
                out[ip] = in[ip] / m_dt +
                          (theta_pert[ip] - m_theta_n[ip]) / eps;
        });

    int n_iter = 0;

    // Freeze the metric and residual at the initial surface. The loop keeps the
    // invariant that on entry m_gamma_LL and m_theta_n are frozen at m_state.h
    // and theta_old == inf_norm(m_theta_n).
    this->set_particle_positions(m_state.h);
    this->interpolate_metric(m_state.h);
    this->theta_from_metric(m_state.h, m_theta_n);

    double theta_old = inf_norm(m_theta_n);

    // Global pseudo-timestep
    m_dt = 1e-2;

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

    while (theta_old > m_params.tolerance && n_iter < m_params.max_iter)
    {
        // RHS = -Theta(h_n) (frozen at the current surface).
        std::vector<double> minus_theta(m_num_particles);
        for (int ip = 0; ip < m_num_particles; ++ip)
            minus_theta[ip] = -m_theta_n[ip];
        m_jac_op->flat_to_mf(minus_theta, rhs);

        // Solve (I/dt + J) delta_h = -Theta(h_n) with a single-level,
        // matrix-free BiCGStab (coarsening disabled in the operator).
        delta.setVal(0.0);
        amrex::MLMG mlmg(*m_jac_op);
        mlmg.setBottomSolver(amrex::BottomSolver::bicgstab);
        // No geometric multigrid levels exist, so the pre/post V-cycle sweeps
        // never run; disable the post-bottom smoothing too, since Fsmooth
        // (relaxation) is not implemented for this matrix-free operator.
        mlmg.setBottomSmooth(0);
        mlmg.setFinalSmooth(0);
        mlmg.setVerbose(0);
        mlmg.setBottomVerbose(0);
        // As dt grows the frozen-metric operator (I/dt + J) becomes
        // ill-conditioned; an inexact linear solve is acceptable for PTC, so
        // take whatever increment the solver reaches instead of aborting.
        mlmg.setThrowException(true);
        try
        {
            mlmg.solve({&delta}, {&rhs}, m_params.linear_rel_tol,
                       m_params.linear_abs_tol);
        }
        catch (const std::exception &)
        {
            // Keep the partially-converged increment in delta.
        }

        std::vector<double> delta_h(m_num_particles);
        m_jac_op->mf_to_flat(delta, delta_h);

        // Trial update h_{n+1} = h_n + delta_h, keeping h_n so the step can be
        // rejected if it does not reduce the residual.
        const std::vector<double> h_backup = m_state.h;
        for (int ip = 0; ip < m_num_particles; ++ip)
            m_state.h[ip] += delta_h[ip];

        // Evaluate Theta at the trial state, freezing its metric.
        this->set_particle_positions(m_state.h);
        this->interpolate_metric(m_state.h);
        this->theta_from_metric(m_state.h, m_theta_vals);

        const double theta_new = inf_norm(m_theta_vals);

        n_iter++;

        if (theta_new < theta_old)
        {
            // Accept: the frozen metric/residual now hold at the trial surface.
            m_theta_n = m_theta_vals;
            // SER: grow dt as the residual falls, approaching Newton.
            m_dt      = update_dt(m_dt, theta_old, theta_new);
            theta_old = theta_new;
        }
        else
        {
            // Reject: restore h_n and its frozen metric/residual, and shrink dt
            // so the next step is more strongly regularised (I/dt dominant). A
            // step that fails to reduce the residual is rejected here, including
            // the zero increment BiCGStab returns when it breaks down on the
            // ill-conditioned (I/dt + J) at large dt -- shrinking dt then pulls
            // the operator back into the regime the matrix-free solve handles,
            // so dt self-limits at the largest value the solver supports.
            m_state.h = h_backup;
            this->set_particle_positions(m_state.h);
            this->interpolate_metric(m_state.h);
            this->theta_from_metric(m_state.h, m_theta_n);
            m_dt = std::max(m_dt * m_dt_shrink, m_min_dt);
        }

        amrex::Print() << " AHFinder iter " << n_iter << ": theta = "
                       << theta_old << ", dt = " << m_dt << "\n";

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

    amrex::AllPrint() << "\n AHFinder converged with inf norm of theta = "
                      << theta_old << " in " << n_iter << " iterations\n";

    // Report the converged surface's area and irreducible mass
    // (Christodoulou formula: M = sqrt(A / 16 pi)).
    const amrex::Real area = m_geometry.area();
    amrex::AllPrint() << " AHFinder surface area = " << area << "\n";

    const amrex::Real mass = std::sqrt(area / (16.0 * M_PI));
    amrex::AllPrint() << " AHFinder irreducible mass = " << mass << "\n";
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
amrex::Real
AHFinder<num_components>::update_dt(amrex::Real dt, double theta_old,
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

    m_metric_query_state.setCoords(0, interp_coords_x.data() + m_start)
        .setCoords(1, interp_coords_y.data() + m_start)
        .setCoords(2, interp_coords_z.data() + m_start);
    m_metric_query_deriv.setCoords(0, interp_coords_x.data() + m_start)
        .setCoords(1, interp_coords_y.data() + m_start)
        .setCoords(2, interp_coords_z.data() + m_start);

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
}

template <int num_components>
void AHFinder<num_components>::theta_from_metric(const std::vector<double> &h,
                                                 std::vector<double> &theta_out)
{
    // Evaluate the expansion Theta at every grid point for the surface radius
    // h, using the metric frozen by the last interpolate_metric(). Only h and
    // its ring-grid derivatives enter, so this is cheap (local work plus one
    // reduce) and is what the matrix-free Jacobian differentiates.
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
}

#endif /* AHFINDER_IMPL_HPP_ */
