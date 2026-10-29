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
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>

inline void AHFinder::init(GRAmr *gramr_ptr)
{
    // tolerance, r, initial_dt, max_iter and max_stall_iter come
    // from the "ah_finder" scope of the input file; see
    // AHFinderParameters.hpp for defaults and meaning.
    m_params.fill_params();

    m_min_dt    = 1e-4;
    m_dt_shrink = 0.5;
    m_fd_sigma  = 1.0;

    this->setup_metric_query();

    // Set up interpolator
    this->setup(gramr_ptr);

    m_geometry.set_surface_data(&m_h, &m_gamma_LL);

    this->setup_jacobian();

    // Initialise h for particles, then place the particles
    this->init_particle_vals();
    this->set_particle_positions(m_h);
}

inline void AHFinder::find()
{
    int n_iter = 0;

    // The grid may have evolved since the last find(), so refill the ghost
    // cells on the first interpolation
    m_fill_ghosts = true;

    this->set_particle_positions(m_h);
    interp_metric();
    theta_from_metric(m_h, m_theta_vals);

    double theta_old = inf_norm(m_theta_vals);

    // Global pseudo-timestep
    amrex::Real dt = m_params.initial_dt;

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

    std::vector<double> delta_h(m_num_particles);
    std::vector<double> jacobian(static_cast<std::size_t>(m_ldab) *
                                 m_num_particles);

    // Number of consecutive iterations with dt at its minimum
    int n_stall = 0;

    while (theta_old > m_params.tolerance && n_iter < m_params.max_iter)
    {
        // One pseudo-transient continuation step of h_dot = -Theta:
        // solve (J + I/dt) delta_h = -Theta.
        assemble_jacobian(jacobian);
        if (!lapack_solve(jacobian, dt, delta_h))
        {
            amrex::Print() << "\n AHFinder: banded LU factorisation of "
                              "J + I/dt failed\n";
            break;
        }

        for (int id = 0; id < m_num_particles; ++id)
        {
            m_h[id] += delta_h[id];
        }

        // Evaluate Theta at the new state
        this->set_particle_positions(m_h);
        interp_metric();
        theta_from_metric(m_h, m_theta_vals);

        double theta_new = inf_norm(m_theta_vals);

        amrex::Print() << "\n theta_old = " << theta_old << "\n";

        amrex::Print() << "\n theta_new = " << theta_new << "\n";

        amrex::Print() << "-------------------------\n";

        // Adapt the global timestep for the next step.
        dt = update_dt(dt, theta_old, theta_new);
        amrex::Print() << "dt = " << dt << "\n";

        n_stall = (dt <= m_min_dt) ? n_stall + 1 : 0;

        theta_old = theta_new;

        n_iter++;

        if (io_proc)
        {
            theta_log << n_iter << "," << theta_old << std::endl;
            dt_log << n_iter << "," << dt << std::endl;
        }
        write_particles(n_iter);

        if (n_stall >= m_params.max_stall_iter)
        {
            amrex::Print()
                << "\n AHFinder stalled: dt has been at its minimum ("
                << m_min_dt << ") for " << n_stall
                << " consecutive iterations\n";
            break;
        }
    }

    if (io_proc)
    {
        theta_log.close();
        dt_log.close();
    }

    if (theta_old > m_params.tolerance)
    {
        amrex::Print() << "\n AHFinder failed to converge in " << n_iter
                       << " iterations, inf norm of theta = " << theta_old
                       << "\n";
        return;
    }

    amrex::AllPrint() << "\n AHFinder converged with inf norm of theta = "
                      << theta_old << " in " << n_iter << " iterations\n";

    // Report the converged surface's area and irreducible mass
    // (Christodoulou formula: M = sqrt(A / 16 pi)).
    const amrex::Real area = m_geometry.area(m_start, m_n_local);
    amrex::Print() << " AHFinder surface area = " << area << "\n";

    const amrex::Real mass = std::sqrt(area / (16.0 * M_PI));
    amrex::Print() << " AHFinder irreducible mass = " << mass << "\n";
}

inline void AHFinder::init_particle_vals()
{
    const double r0 = m_geometry.guess_radius();

    // The initial surface is the sphere r = guess_radius
    m_h.assign(m_num_particles, r0);
}

inline void AHFinder::set_particle_positions(const std::vector<double> &h)
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

inline amrex::Real AHFinder::update_dt(amrex::Real dt, double theta_old,
                                       double theta_new) const
{
    // Switched evolution relaxation: grow dt as Theta improves, so the
    // implicit step tends to a Newton step near convergence.
    double ratio = m_params.r * theta_old / theta_new;
    // Limit shrinking of timestep
    ratio = std::max(ratio, m_dt_shrink);

    dt *= ratio;

    // Ensure the timestep doesn't get too small.
    dt = std::max(dt, m_min_dt);

    return dt;
}

inline double AHFinder::inf_norm(std::vector<double> arr)
{
    double max_el = std::abs(arr[0]);
    for (auto &&i : arr)
    {
        if (std::abs(i) > max_el)
            max_el = std::abs(i);
    }

    return max_el;
}

inline void AHFinder::setup_metric_query()
{
    for (auto &v : m_metric_state)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dx)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dy)
        v.resize(m_num_particles);
    for (auto &v : m_metric_dz)
        v.resize(m_num_particles);

    m_metric_query.setCoords(0, interp_coords_x.data() + m_start)
        .setCoords(1, interp_coords_y.data() + m_start)
        .setCoords(2, interp_coords_z.data() + m_start);

    // chi, h_ij, K, A_ij (values only).
    m_metric_query.addComp(c_chi, m_metric_state[c_chi].data() + m_start,
                           VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query.addComp(comp, m_metric_state[comp].data() + m_start,
                               VariableType::state);
    }
    m_metric_query.addComp(c_K, m_metric_state[c_K].data() + m_start,
                           VariableType::state);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_A11, i, j);
        m_metric_query.addComp(comp, m_metric_state[comp].data() + m_start,
                               VariableType::state);
    }

    m_metric_query.addComp(c_chi, m_metric_dx[c_chi].data() + m_start,
                           VariableType::state, BCParity::undefined,
                           Derivative::dx);
    m_metric_query.addComp(c_chi, m_metric_dy[c_chi].data() + m_start,
                           VariableType::state, BCParity::undefined,
                           Derivative::dy);
    m_metric_query.addComp(c_chi, m_metric_dz[c_chi].data() + m_start,
                           VariableType::state, BCParity::undefined,
                           Derivative::dz);
    FOR2_SYM(i, j)
    {
        int comp = sym_var_idx(c_h11, i, j);
        m_metric_query.addComp(comp, m_metric_dx[comp].data() + m_start,
                               VariableType::state, BCParity::undefined,
                               Derivative::dx);
        m_metric_query.addComp(comp, m_metric_dy[comp].data() + m_start,
                               VariableType::state, BCParity::undefined,
                               Derivative::dy);
        m_metric_query.addComp(comp, m_metric_dz[comp].data() + m_start,
                               VariableType::state, BCParity::undefined,
                               Derivative::dz);
    }
}

inline void AHFinder::interp_metric()
{
    // Particle positions move with h, so refresh them on every call. The
    // grid is frozen during find(), so only fill the ghost cells once.
    this->interp(m_metric_query, true, "", 0.0, m_fill_ghosts);
    m_fill_ghosts = false;
}

inline void AHFinder::theta_from_metric(const std::vector<double> &h,
                                        std::vector<double> &theta)
{
    m_geometry.set_h_derivatives(h);

    theta.assign(m_num_particles, 0.0);

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

    double *theta_ptr = theta.data();

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
        // grid point's unit direction, since the surface point is
        // center_i + r * dir_i. Taken from the direction rather than the
        // particle position so it stays consistent with h when the metric
        // is frozen at a different h (assemble_jacobian()).
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

    amrex::ParallelDescriptor::ReduceRealSum(theta.data(), m_num_particles);
}

inline void AHFinder::setup_jacobian()
{
    const int n         = m_num_particles;
    const int ring_size = m_geometry.ring_size();
    const int n_rings   = m_geometry.n_rings();

    // Rows each column appears in: the transpose of derivative_stencil()
    std::vector<std::vector<int>> rows_of(n);
    for (int i = 0; i < n; ++i)
    {
        const amrex::GpuArray<int, 9> stencil =
            m_geometry.derivative_stencil(i);
        for (int s = 0; s < 9; ++s)
        {
            rows_of[stencil[s]].push_back(i);
        }
    }

    // Greedy colouring: give each column the lowest colour not already used
    // by a column it shares a row with. used[c] == j marks colour c as taken
    // while colouring column j.
    m_colour.assign(n, -1);
    m_n_colours = 0;
    std::vector<int> used;
    for (int j = 0; j < n; ++j)
    {
        for (const int i : rows_of[j])
        {
            const amrex::GpuArray<int, 9> stencil =
                m_geometry.derivative_stencil(i);
            for (int s = 0; s < 9; ++s)
            {
                const int c = m_colour[stencil[s]];
                if (c >= 0)
                {
                    used[c] = j;
                }
            }
        }

        int c = 0;
        while (c < m_n_colours && used[c] == j)
        {
            ++c;
        }
        if (c == m_n_colours)
        {
            ++m_n_colours;
            used.push_back(-1);
        }
        m_colour[j] = c;
    }

    // Lookup table: the one column of each colour that each row depends on
    m_colour_column.assign(static_cast<std::size_t>(n) * m_n_colours, -1);
    for (int i = 0; i < n; ++i)
    {
        const amrex::GpuArray<int, 9> stencil =
            m_geometry.derivative_stencil(i);
        for (int s = 0; s < 9; ++s)
        {
            const int j                                    = stencil[s];
            m_colour_column[i * m_n_colours + m_colour[j]] = j;
        }
    }

    // Band ordering: fold each ring (0, n-1, 1, n-2, ...) so neighbours
    // across the periodic seam sit close together
    m_band_order.resize(n);
    m_band_index.resize(n);
    for (int ring = 0; ring < n_rings; ++ring)
    {
        for (int p = 0; p < ring_size; ++p)
        {
            const int fold = (p % 2 == 0) ? p / 2 : ring_size - 1 - p / 2;
            const int k    = ring * ring_size + p;
            const int idx  = ring * ring_size + fold;

            m_band_order[k]   = idx;
            m_band_index[idx] = k;
        }
    }

    // Bandwidths of the sparsity pattern in the band ordering
    m_kl = 0;
    m_ku = 0;
    for (int i = 0; i < n; ++i)
    {
        const amrex::GpuArray<int, 9> stencil =
            m_geometry.derivative_stencil(i);
        for (int s = 0; s < 9; ++s)
        {
            const int d = m_band_index[i] - m_band_index[stencil[s]];
            m_kl        = std::max(m_kl, d);
            m_ku        = std::max(m_ku, -d);
        }
    }
    m_ldab = 2 * m_kl + m_ku + 1;

    amrex::Print() << "AHFinder Jacobian: " << m_n_colours
                   << " colours, bandwidths kl = " << m_kl << ", ku = " << m_ku
                   << " (ring size " << ring_size << ")\n";
}

inline void AHFinder::assemble_jacobian(std::vector<double> &jacobian)
{
    const int n = m_num_particles;

    // Finite-difference step at each point, scaled to h
    const double sqrt_eps = std::sqrt(std::numeric_limits<double>::epsilon());
    std::vector<double> fd_eps(n);
    for (int j = 0; j < n; ++j)
    {
        fd_eps[j] = m_fd_sigma * sqrt_eps * (1.0 + std::abs(m_h[j]));
    }

    // Perturbed h and the Theta evaluated there
    std::vector<double> h_pert(n);
    std::vector<double> theta_pert(n);

    std::fill(jacobian.begin(), jacobian.end(), 0.0);

    // Perturb all points of one colour at once and update the
    // points in the jacobian they contribute to
    for (int c = 0; c < m_n_colours; ++c)
    {
        for (int j = 0; j < n; ++j)
        {
            h_pert[j] = (m_colour[j] == c) ? m_h[j] + fd_eps[j] : m_h[j];
        }

        theta_from_metric(h_pert, theta_pert);

        for (int i = 0; i < n; ++i)
        {
            const int j = m_colour_column[i * m_n_colours + c];
            if (j >= 0)
            {
                jacobian[band_offset(i, j)] =
                    (theta_pert[i] - m_theta_vals[i]) / fd_eps[j];
            }
        }
    }

    // Perturb all points, interpolate and calculate new value of
    // Theta, calculate delta Theta with finite difference and apply
    // resulting update to the diagonal of the jacobian
    for (int j = 0; j < n; ++j)
    {
        h_pert[j] = m_h[j] + fd_eps[j];
    }
    this->set_particle_positions(h_pert);
    interp_metric();

    theta_from_metric(m_h, theta_pert);

    for (int i = 0; i < n; ++i)
    {
        jacobian[band_offset(i, i)] +=
            (theta_pert[i] - m_theta_vals[i]) / fd_eps[i];
    }
}

inline bool AHFinder::lapack_solve(std::vector<double> &jacobian,
                                   amrex::Real dt, std::vector<double> &delta_h)
{
    const int n    = m_num_particles;
    const int root = amrex::ParallelDescriptor::IOProcessorNumber();
    int info       = 0;

    // Every rank holds the full J and Theta, but only one does the solve;
    // the result is then broadcast to the others
    if (amrex::ParallelDescriptor::MyProc() == root)
    {
        const int nrhs = 1;
        std::vector<int> ipiv(n);
        std::vector<double> rhs(n);

        // Apply pseudo timestep to diagonal
        for (int k = 0; k < n; ++k)
        {
            jacobian[m_kl + m_ku + k * m_ldab] += 1.0 / dt;
        }

        for (int k = 0; k < n; ++k)
        {
            rhs[k] = -m_theta_vals[m_band_order[k]];
        }

        dgbsv_(&n, &m_kl, &m_ku, &nrhs, jacobian.data(), &m_ldab, ipiv.data(),
               rhs.data(), &n, &info);

        if (info == 0)
        {
            for (int k = 0; k < n; ++k)
            {
                delta_h[m_band_order[k]] = rhs[k];
            }
        }
    }

    amrex::ParallelDescriptor::Bcast(&info, 1, root);
    if (info != 0)
    {
        return false;
    }

    amrex::ParallelDescriptor::Bcast(delta_h.data(), n, root);

    return true;
}

#endif /* AHFINDER_IMPL_HPP_ */
