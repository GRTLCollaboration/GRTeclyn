#if !defined(AHFINDER_HPP_)
#error "This file should only be included through AHFinder.hpp"
#endif

#ifndef AHFINDER_IMPL_HPP_
#define AHFINDER_IMPL_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>

#include "CCZ4StateVariables.hpp"
#include "DefaultLevelBld.hpp"
#include "Derivative.hpp"
#include "GRAmr.hpp"
#include "ParticleInterpolator.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include <filesystem>
#include <fstream>

template <int num_components>
void AHFinder<num_components>::init(GRAmr *gramr_ptr)
{
    m_tol    = 1e-4;
    m_c      = 1.0;
    m_min_dt = 1e-4;
    m_r      = 1.15;
    m_eta    = 3;

    m_dt_shrink   = 0.8;
    m_dt_grow     = 1.25;
    m_theta_floor = 1e-12;

    this->setup_metric_query();

    // Set up interpolator
    this->setup(gramr_ptr);

    m_geometry.set_surface_data(&m_state.h, &m_gamma_LL);

    // Initialise h, v and dt values for particles, then place the particles
    this->init_particle_vals();
    this->set_particle_positions(m_state.h);
}

template <int num_components> void AHFinder<num_components>::find()
{
    // Create amrex time integrator
    // Allows for different (explicit) time stepping methods
    m_integrator = std::make_unique<amrex::TimeIntegrator<AHState>>(m_state);
    m_integrator->set_rhs([this](AHState &rhs, AHState &state, amrex::Real time)
                          { this->compute_rhs(rhs, state, time); });

    int n_iter = 0;

    this->set_particle_positions(m_state.h);
    compute_theta(m_state.h);

    double theta_old = inf_norm(m_theta_vals);

    // Global pseudo-timeste
    amrex::Real dt = 1e-2;

    const bool io_proc = amrex::ParallelDescriptor::IOProcessor();

    // Temp logging
    amrex::Print() << "\n AHFinder expansion Theta inf "
                      "norm = "
                   << theta_old << "\n";

    std::ofstream theta_log;
    std::ofstream dt_log;

    if (io_proc)
    {
        theta_log.open("theta_vs_iter.csv");
        dt_log.open("dt_vs_iter.csv");
        std::filesystem::create_directory("particles");

        theta_log << n_iter << "," << theta_old << std::endl;
        dt_log << n_iter << "," << dt << std::endl;
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

    AHState new_state = m_state;

    while (theta_old > m_tol)
    {
        // Advance one pseudo-time step with the AMReX integrator.
        m_integrator->set_time_step(dt);
        m_integrator->advance(m_state, new_state, n_iter * dt, dt);
        std::swap(m_state, new_state);

        // Evaluate Theta at the new state
        this->set_particle_positions(m_state.h);
        compute_theta(m_state.h);

        double theta_new = inf_norm(m_theta_vals);

        amrex::Print() << "\n theta_old = " << theta_old << "\n";

        amrex::Print() << "\n theta_new = " << theta_new << "\n";

        amrex::Print() << "-------------------------\n";

        // Adapt the global timestep for the next step.
        dt = update_dt(dt, theta_old, theta_new, m_state.h);
        amrex::Print() << "dt = " << dt << "\n";

        theta_old = theta_new;

        n_iter++;

        if (io_proc)
        {
            theta_log << n_iter << "," << theta_old << std::endl;
            dt_log << n_iter << "," << dt << std::endl;
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

    // Since h = v - eta * h, start velocity at eta * h so we don't
    // immediately collapse inwards
    m_state.v.assign(m_num_particles, m_eta * r0);
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
void AHFinder<num_components>::compute_rhs(AHState &rhs, AHState &state,
                                           amrex::Real /* time */)
{
    this->set_particle_positions(state.h);

    // Calculate Theta
    compute_theta(state.h);

    //   h_dot = v - eta * h
    //   v_dot = -c^2 * Theta
    rhs.h.assign(m_num_particles, 0.0);
    rhs.v.assign(m_num_particles, 0.0);
    for (int id = 0; id < m_num_particles; ++id)
    {
        rhs.h[id] = state.v[id] - m_eta * state.h[id];
        rhs.v[id] = -std::pow(m_c, 2) * m_theta_vals[id];
    }
}

template <int num_components>
amrex::Real
AHFinder<num_components>::update_dt(amrex::Real dt, double theta_old,
                                    double theta_new,
                                    const std::vector<double> &h) const
{
    // Update dt based on ratio of improvement of theta

    // Limit timestep by CFL condition of closest pair of particles
    const double cfl_factor = 7.0;
    const double max_dt     = cfl_factor * m_geometry.min_ring_spacing(h);

    // Update time step based on ratio of old to new theta
    // Ensure it doesn't grow or shrink too fast.
    double ratio = (std::abs(theta_new) > m_theta_floor)
                       ? m_r * theta_old / theta_new
                       : m_r;
    ratio        = std::max(ratio, m_dt_shrink);
    ratio        = std::min(ratio, m_dt_grow);

    dt *= ratio;

    // Ensure the timestep doesn't grow too large or small.
    dt = std::max(dt, m_min_dt);
    dt = std::min(dt, max_dt);

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
void AHFinder<num_components>::compute_theta(const std::vector<double> &h)
{

    m_geometry.set_h_derivatives(h);

    m_theta_vals.assign(m_num_particles, 0.0);

    this->interp(m_metric_query_state, true);
    this->interp(m_metric_query_deriv, false);

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

    double *theta_ptr = m_theta_vals.data();

    for (int ip = m_start; ip < m_start + m_n_local; ++ip)
    {
        using namespace TensorAlgebra;

        double r   = h_ptr[ip];
        double chi = state_ptr[c_chi][ip];
        double K   = state_ptr[c_K][ip];

        // Physical metric and extrinsic curvature from CCZ4
        // variables: gamma_ij = h_ij/chi,
        // K_ij = (A_ij + (1/3) h_ij K)/chi. gamma_ij is written directly
        // into the persistent m_gamma_LL[ip], shared with AHGeometry.
        Tensor::Rank2 K_LL;
        FOR (i, j)
        {
            int h_comp           = sym_var_idx(c_h11, i, j);
            int A_comp           = sym_var_idx(c_A11, i, j);
            double h_ij          = state_ptr[h_comp][ip];
            double A_ij          = state_ptr[A_comp][ip];
            m_gamma_LL[ip](i, j) = h_ij / chi;
            K_LL(i, j)           = (A_ij + (1.0 / 3.0) * h_ij * K) / chi;
        }
        Tensor::Rank2 gamma_UU = compute_inverse(m_gamma_LL[ip]);

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
                d1h_kij / chi - m_gamma_LL[ip](i, j) * d1_chi(k) / chi;
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

    amrex::ParallelDescriptor::ReduceRealSum(m_theta_vals.data(),
                                             m_num_particles);
}

#endif /* AHFINDER_IMPL_HPP_ */
