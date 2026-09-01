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

    this->generate_spherical_query();
    this->setup_metric_query();

    // Set up interpolator
    this->setup(gramr_ptr);

    // Populate so we can access the particle data
    this->populate_from_query(query);

    // Initialise h, v and dt values for particles
    this->init_particle_vals();
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

    // Temp logging
    amrex::AllPrint() << "\n AHFinder expansion Theta inf "
                         "norm = "
                      << theta_old << "\n";

    std::ofstream theta_log("theta_vs_iter.csv");
    theta_log << n_iter << "," << theta_old << std::endl;

    std::ofstream dt_log("dt_vs_iter.csv");
    dt_log << n_iter << "," << dt << std::endl;

    std::filesystem::create_directory("particles");
    auto write_particles = [&](int iter)
    {
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

        amrex::AllPrint() << "\n theta_old = " << theta_old << "\n";

        amrex::AllPrint() << "\n theta_new = " << theta_new << "\n";

        amrex::AllPrint() << "-------------------------\n";

        // Adapt the global timestep for the next step.
        dt = update_dt(dt, theta_old, theta_new, m_state.h);
        amrex::Print() << "dt = " << dt << "\n";

        theta_old = theta_new;

        n_iter++;

        theta_log << n_iter << "," << theta_old << std::endl;
        dt_log << n_iter << "," << dt << std::endl;
        write_particles(n_iter);
    }

    theta_log.close();
    dt_log.close();

    amrex::AllPrint() << "\n AHFinder converged with inf norm of theta = "
                      << theta_old << " in " << n_iter << " iterations\n";
}

template <int num_components>
void AHFinder<num_components>::generate_spherical_query()
{
    // Generate particle coordinates such that they are laid out
    // in rings on a sphere

    // Aim for ~twice as many longitudes as latitudes
    m_n_rings  = 1;
    int target = std::max(1, (int)std::round(std::sqrt(m_num_particles / 2.0)));
    for (int n = 1; n <= m_num_particles; ++n)
    {
        if (m_num_particles % n == 0 &&
            std::abs(n - target) < std::abs(m_n_rings - target))
        {
            m_n_rings = n;
        }
    }
    m_ring_size = m_num_particles / m_n_rings;

    for (int i = 0; i < m_n_rings; ++i)
    {
        // Offset theta so we don't get points on the poles
        double theta = (i + 0.5) * M_PI / m_n_rings;
        for (int j = 0; j < m_ring_size; ++j)
        {
            double phi = j * 2. * M_PI / m_ring_size;
            int idx    = i * m_ring_size + j;

            // Calculate normalised vector pointing to centre
            m_dir_x[idx] = cos(phi) * sin(theta);
            m_dir_y[idx] = sin(phi) * sin(theta);
            m_dir_z[idx] = cos(theta);

            interp_coords_x[idx] = m_center[0] + m_guess_radius * m_dir_x[idx];
            interp_coords_y[idx] = m_center[1] + m_guess_radius * m_dir_y[idx];
            interp_coords_z[idx] = m_center[2] + m_guess_radius * m_dir_z[idx];
        }
    }

    query.setCoords(0, interp_coords_x.data())
        .setCoords(1, interp_coords_y.data())
        .setCoords(2, interp_coords_z.data());
}

template <int num_components>
amrex::GpuArray<int, 4> AHFinder<num_components>::neighbours(int j) const
{
    // Get the 4 neighbours of a particle (north, south, east and west).
    // If we are at the top/bottom ring, we can take the north/south
    // neighbour respectively as the opposite point on the same ring
    int ring_num = j / m_ring_size;
    int ring_pos = j % m_ring_size;

    int opposite_point =
        ring_num * m_ring_size + (ring_pos + m_ring_size / 2) % m_ring_size;

    int north, south, east, west;

    // North/south neighbours are next ring above/below.
    if (ring_num > 0)
    {
        north = (ring_num - 1) * m_ring_size + ring_pos;
    }
    else
    {
        north = opposite_point;
    }

    if (ring_num < m_n_rings - 1)
    {
        south = (ring_num + 1) * m_ring_size + ring_pos;
    }
    else
    {
        south = opposite_point;
    }

    // East/west neighbours are adjacent in the same ring
    // Modulo to wrap around
    east = ring_num * m_ring_size + (ring_pos + 1) % m_ring_size;
    west = ring_num * m_ring_size + (ring_pos - 1 + m_ring_size) % m_ring_size;

    return {north, south, east, west};
}

template <int num_components>
void AHFinder<num_components>::init_particle_vals()
{
    for (int id = 0; id < m_num_particles; ++id)
    {
        // Height from centre
        m_state.h[id] =
            std::sqrt(std::pow(interp_coords_x[id] - m_center[0], 2) +
                      std::pow(interp_coords_y[id] - m_center[1], 2) +
                      std::pow(interp_coords_z[id] - m_center[2], 2));

        // Since h = v - eta * h, start velocity at eta * h so we don't
        // immediately collapse inwards
        m_state.v[id] = m_eta * m_state.h[id];
    }
}

template <int num_components>
void AHFinder<num_components>::set_particle_positions(
    const std::vector<double> &h)
{
    // Set particles' positions based on their radius from the centre
    for (int id = 0; id < m_num_particles; ++id)
    {
        const double r = h[id];

        interp_coords_x[id] = m_center[0] + r * m_dir_x[id];
        interp_coords_y[id] = m_center[1] + r * m_dir_y[id];
        interp_coords_z[id] = m_center[2] + r * m_dir_z[id];

        amrex::GpuArray<amrex::ParticleReal, AMREX_SPACEDIM> coords = {
            interp_coords_x[id], interp_coords_y[id], interp_coords_z[id]};
        this->check_domain(coords);
    }

    // Push positions to particles
    this->update_particle_positions(query);
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
    double cfl_factor = 7.0;
    double max_dt     = std::numeric_limits<double>::max();

    // Limit timestep by CFL condition of closest pair of particles
    for (int id = 0; id < m_num_particles; ++id)
    {
        const int ring_num = id / m_ring_size;
        const double theta = (ring_num + 0.5) * M_PI / m_n_rings;
        const double r     = h[id];

        const double d_theta = r * M_PI / m_n_rings;
        const double d_phi   = r * sin(theta) * 2.0 * M_PI / m_ring_size;

        max_dt = std::min(max_dt, cfl_factor * std::min(d_theta, d_phi));
    }

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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::GpuArray<amrex::Real, 3>
AHFinder<num_components>::ring_gradient(const double *field_ptr, int north_idx,
                                        int south_idx, int east_idx,
                                        int west_idx, double d_theta,
                                        double d_phi, double theta, double phi,
                                        double r)
{
    // Perform finite differences with a particle and its
    // neighbours in spherical (theta/phi) coordinates and convert
    // them to cartesian (x/y/z)
    amrex::Real d_dtheta =
        (field_ptr[south_idx] - field_ptr[north_idx]) / (2.0 * d_theta);
    amrex::Real d_dphi =
        (field_ptr[east_idx] - field_ptr[west_idx]) / (2.0 * d_phi);

    // Convert to Cartesian via the spherical Jacobian
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double sin_phi   = sin(phi);
    double cos_phi   = cos(phi);

    amrex::Real dx =
        d_dtheta * cos_theta * cos_phi / r - d_dphi * sin_phi / (r * sin_theta);
    amrex::Real dy =
        d_dtheta * cos_theta * sin_phi / r + d_dphi * cos_phi / (r * sin_theta);
    amrex::Real dz = -d_dtheta * sin_theta / r;

    return {dx, dy, dz};
}

template <int num_components>
void AHFinder<num_components>::h_derivs(const std::vector<double> &h)
{
    // Calculate dh/dx, dh/dy, dh/dz for all particles.
    const double *h_ptr = h.data();
    double *dhdx_ptr    = m_dhdx.data();
    double *dhdy_ptr    = m_dhdy.data();
    double *dhdz_ptr    = m_dhdz.data();

    const double d_theta = M_PI / m_n_rings;
    const double d_phi   = 2.0 * M_PI / m_ring_size;

    for (int id = 0; id < m_num_particles; ++id)
    {
        const int self_ring_num = id / m_ring_size;
        const int self_ring_pos = id % m_ring_size;
        const double theta      = (self_ring_num + 0.5) * M_PI / m_n_rings;
        const double phi        = self_ring_pos * 2.0 * M_PI / m_ring_size;
        const double r          = h_ptr[id];

        const amrex::GpuArray<int, 4> nb = this->neighbours(id);

        const amrex::GpuArray<amrex::Real, 3> grad_h = ring_gradient(
            h_ptr, nb[0], nb[1], nb[2], nb[3], d_theta, d_phi, theta, phi, r);

        dhdx_ptr[id] = grad_h[0];
        dhdy_ptr[id] = grad_h[1];
        dhdz_ptr[id] = grad_h[2];
    }
}

template <int num_components>
void AHFinder<num_components>::h_hessian(const std::vector<double> &h)
{
    // Calculate the Hessian of h
    // Must be called after h_derivs() has populated m_dhdx/m_dhdy/m_dhdz.
    const double *h_ptr    = h.data();
    const double *dhdx_ptr = m_dhdx.data();
    const double *dhdy_ptr = m_dhdy.data();
    const double *dhdz_ptr = m_dhdz.data();

    double *d2h_xx_ptr = m_d2h_xx.data();
    double *d2h_yy_ptr = m_d2h_yy.data();
    double *d2h_zz_ptr = m_d2h_zz.data();
    double *d2h_xy_ptr = m_d2h_xy.data();
    double *d2h_xz_ptr = m_d2h_xz.data();
    double *d2h_yz_ptr = m_d2h_yz.data();

    const double d_theta = M_PI / m_n_rings;
    const double d_phi   = 2.0 * M_PI / m_ring_size;

    for (int id = 0; id < m_num_particles; ++id)
    {
        const int self_ring_num = id / m_ring_size;
        const int self_ring_pos = id % m_ring_size;
        const double theta      = (self_ring_num + 0.5) * M_PI / m_n_rings;
        const double phi        = self_ring_pos * 2.0 * M_PI / m_ring_size;
        const double r          = h_ptr[id];

        const amrex::GpuArray<int, 4> nb = this->neighbours(id);

        // d/dx_j(dh/dx), d/dx_j(dh/dy), d/dx_j(dh/dz)
        const amrex::GpuArray<amrex::Real, 3> grad_x =
            ring_gradient(dhdx_ptr, nb[0], nb[1], nb[2], nb[3], d_theta, d_phi,
                          theta, phi, r);
        const amrex::GpuArray<amrex::Real, 3> grad_y =
            ring_gradient(dhdy_ptr, nb[0], nb[1], nb[2], nb[3], d_theta, d_phi,
                          theta, phi, r);
        const amrex::GpuArray<amrex::Real, 3> grad_z =
            ring_gradient(dhdz_ptr, nb[0], nb[1], nb[2], nb[3], d_theta, d_phi,
                          theta, phi, r);

        d2h_xx_ptr[id] = grad_x[0];
        d2h_yy_ptr[id] = grad_y[1];
        d2h_zz_ptr[id] = grad_z[2];

        d2h_xy_ptr[id] = 0.5 * (grad_x[1] + grad_y[0]);
        d2h_xz_ptr[id] = 0.5 * (grad_x[2] + grad_z[0]);
        d2h_yz_ptr[id] = 0.5 * (grad_y[2] + grad_z[1]);
    }
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

    m_metric_query_state.setCoords(0, interp_coords_x.data())
        .setCoords(1, interp_coords_y.data())
        .setCoords(2, interp_coords_z.data());
    m_metric_query_deriv.setCoords(0, interp_coords_x.data())
        .setCoords(1, interp_coords_y.data())
        .setCoords(2, interp_coords_z.data());

    // chi, h_ij, K, A_ij (values only).
    m_metric_query_state.addComp(c_chi, m_metric_state[c_chi].data(),
                                 VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query_state.addComp(comp, m_metric_state[comp].data(),
                                     VariableType::state);
    }
    m_metric_query_state.addComp(c_K, m_metric_state[c_K].data(),
                                 VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_A11, i, j);
        m_metric_query_state.addComp(comp, m_metric_state[comp].data(),
                                     VariableType::state);
    }

    m_metric_query_deriv.addComp(c_chi, m_metric_dx[c_chi].data(),
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dx);
    m_metric_query_deriv.addComp(c_chi, m_metric_dy[c_chi].data(),
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dy);
    m_metric_query_deriv.addComp(c_chi, m_metric_dz[c_chi].data(),
                                 VariableType::state, BCParity::undefined,
                                 Derivative::dz);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query_deriv.addComp(comp, m_metric_dx[comp].data(),
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dx);
        m_metric_query_deriv.addComp(comp, m_metric_dy[comp].data(),
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dy);
        m_metric_query_deriv.addComp(comp, m_metric_dz[comp].data(),
                                     VariableType::state, BCParity::undefined,
                                     Derivative::dz);
    }
}

template <int num_components>
void AHFinder<num_components>::compute_theta(const std::vector<double> &h)
{

    h_derivs(h);
    h_hessian(h);

    this->interp(m_metric_query_state, false);
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

    const double *pos_ptr[AMREX_SPACEDIM] = {
        interp_coords_x.data(), interp_coords_y.data(), interp_coords_z.data()};

    const double *dhdx_ptr = m_dhdx.data();
    const double *dhdy_ptr = m_dhdy.data();
    const double *dhdz_ptr = m_dhdz.data();

    const double *d2h_xx_ptr = m_d2h_xx.data();
    const double *d2h_yy_ptr = m_d2h_yy.data();
    const double *d2h_zz_ptr = m_d2h_zz.data();
    const double *d2h_xy_ptr = m_d2h_xy.data();
    const double *d2h_xz_ptr = m_d2h_xz.data();
    const double *d2h_yz_ptr = m_d2h_yz.data();

    double *theta_ptr = m_theta_vals.data();

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> center = {
        m_center[0], m_center[1], m_center[2]};

    for (int ip = 0; ip < m_num_particles; ++ip)
    {
        using namespace TensorAlgebra;

        double r   = h_ptr[ip];
        double chi = state_ptr[c_chi][ip];
        double K   = state_ptr[c_K][ip];

        // Physical metric and extrinsic curvature from CCZ4
        // variables: gamma_ij = h_ij/chi,
        // K_ij = (A_ij + (1/3) h_ij K)/chi.
        Tensor::Rank2 gamma_LL;
        Tensor::Rank2 K_LL;
        FOR (i, j)
        {
            int h_comp     = sym_var_idx(c_h11, i, j);
            int A_comp     = sym_var_idx(c_A11, i, j);
            double h_ij    = state_ptr[h_comp][ip];
            double A_ij    = state_ptr[A_comp][ip];
            gamma_LL(i, j) = h_ij / chi;
            K_LL(i, j)     = (A_ij + (1.0 / 3.0) * h_ij * K) / chi;
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

        // Level-set gradient F_i = n_i - (grad h)_i, with
        // n_i = (x_i - center_i)/r the flat radial covector
        Tensor::Rank1 n_L;
        FOR (i)
        {
            n_L(i) = (pos_ptr[i][ip] - center[i]) / r;
        }
        Tensor::Rank1 grad_h_L = {dhdx_ptr[ip], dhdy_ptr[ip], dhdz_ptr[ip]};
        Tensor::Rank1 F_L;
        FOR (i)
        {
            F_L(i) = n_L(i) - grad_h_L(i);
        }

        // Hess(F) = Hess(r) - Hess(h), with
        // Hess(r)_ij = (delta_ij - n_i n_j)/r (flat identity).
        Tensor::Rank2 hess_h_LL;
        hess_h_LL(0, 0) = d2h_xx_ptr[ip];
        hess_h_LL(1, 1) = d2h_yy_ptr[ip];
        hess_h_LL(2, 2) = d2h_zz_ptr[ip];
        hess_h_LL(0, 1) = hess_h_LL(1, 0) = d2h_xy_ptr[ip];
        hess_h_LL(0, 2) = hess_h_LL(2, 0) = d2h_xz_ptr[ip];
        hess_h_LL(1, 2) = hess_h_LL(2, 1) = d2h_yz_ptr[ip];

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
}

#endif /* AHFINDER_IMPL_HPP_ */
