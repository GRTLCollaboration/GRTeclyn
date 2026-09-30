#ifndef AHFINDER_HPP_
#define AHFINDER_HPP_

#include <AMReX_Array.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particles.H>
#include <algorithm>
#include <array>
#include <memory>

#include "AHFinderParameters.hpp"
#include "AHFinderState.hpp"
#include "AHGeometry.hpp"
#include "AHJacobianOp.hpp"
#include "ParticleInterpolator.hpp"
#include "Tensor.hpp"

template <int num_components>
class AHFinder : public ParticleInterpolator<num_components>
{
  private:
    int m_num_particles;
    int m_n_local;
    int m_start;

    // PTC/solver parameters read from the "ah_finder" scope of the input file
    // (tolerance, r, max_iter, linear_rel_tol, linear_abs_tol).
    // Filled by init().
    ah_finder_params_t m_params{};

    // Smallest permitted pseudo-timestep, bounds on the per-iteration change
    // in dt, and the magnitude of theta below which the SER ratio is not
    // trusted. Not input parameters: these are guard rails on the adaptive
    // timestep rather than knobs to tune per run.
    //
    // m_dt_shrink is what dt is multiplied by when the line search fails
    // outright. That failure means no step length along delta_h helps, i.e.
    // the direction itself is ascent, so dt has to fall far enough to change
    // the direction in one go; shrinking gently just makes the solver ratchet
    // down over several iterations, each costing a full GMRES solve. It is
    // also the lower clamp on the SER ratio in update_dt(), but that clamp
    // cannot bind for r >= 1: update_dt() is only reached on an accepted full
    // step, where theta_new < theta_old and so the ratio exceeds r.
    amrex::Real m_min_dt;
    amrex::Real m_dt_shrink;
    amrex::Real m_dt_grow;
    amrex::Real m_theta_floor;

    // Factor by which the line search shortens the step length alpha on each
    // backtrack. Also a guard rail rather than a tuning knob; how many times
    // it may be applied is the max_backtracks input parameter.
    amrex::Real m_backtrack_factor;

    // Safety factor in the gate on measure_newton_shift(): c is re-measured
    // once 1/dt falls below m_shift_gate_margin * c, i.e. one factor of two
    // before the shift can actually influence the step. Also a guard rail
    // rather than a knob -- c climbs monotonically in every case measured, so
    // the margin only has to cover one step's worth of that drift.
    amrex::Real m_shift_gate_margin;

    // Coords for particleinterpolator query
    std::vector<double> interp_coords_x{};
    std::vector<double> interp_coords_y{};
    std::vector<double> interp_coords_z{};

    // State storing the surface radius h for all particles. Stored off
    // particles since we need h from other particles to compute its
    // derivative, and this cannot be accessed from another particle if they
    // are not on the same tile
    AHState m_state{};

    // Owns the ring (latitude x longitude) grid: the per-particle
    // directions, the finite-difference stencil, the derivatives of h, and
    // the surface diagnostics (area)
    AHGeometry m_geometry;

    // Physical 3-metric gamma_ij at each particle (flat-indexed as
    // i * m_geometry.ring_size() + j), interpolated once per PTC step in
    // interpolate_metric() and then held frozen while theta_from_metric()
    // varies h. AHGeometry is given a pointer to this in init(), so it always
    // reads the latest values without a separate copy.
    std::vector<Tensor::Rank2> m_gamma_LL{};

    // Matrix-free Jacobian operator (I/dt + J) for the per-step GMRES
    // solve, and the state the operator's mat-vec closes over: the residual
    // Theta(h_n) frozen at the current surface and the current pseudo-timestep.
    std::unique_ptr<AHJacobianOp> m_jac_op;
    std::vector<double> m_theta_n{};
    amrex::Real m_dt{};

    // Output arrays for interpolation queries
    std::array<std::vector<double>, 14> m_metric_state{};
    // First derivatives. theta_from_metric() only reads the first 7 (chi and
    // h_ij, which build the Christoffels); the analytic shift additionally
    // needs d_k of K and A_ij, hence 14 rather than 7.
    std::array<std::vector<double>, 14> m_metric_dx{};
    std::array<std::vector<double>, 14> m_metric_dy{};
    std::array<std::vector<double>, 14> m_metric_dz{};

    // Second derivatives d_j d_k of chi and h_ij, symmetric in (j, k) and
    // indexed by sym2_idx() in the order xx, yy, zz, xy, xz, yz. Only filled
    // when ah_finder.newton_shift_analytic is set: they are what lets the
    // radial transport below move the *first* derivatives that
    // theta_from_metric() reads, not just the values.
    std::array<std::array<std::vector<double>, 7>, 6> m_metric_d2{};

    // Flat index into m_metric_d2 for the symmetric derivative pair (j, k).
    static constexpr int sym2_idx(int j, int k)
    {
        return (j == k) ? j : (2 + j + k);
    }

    // Split up queries num_components doubles as both the
    // query's flat scratch-array size and the number of contiguous grid
    // comps FillPatch fetches, so a query's total (comp, derivative) entry
    // count can't exceed the simulation's total number of state variables.
    InterpolationQueryParticle m_metric_query_state;
    InterpolationQueryParticle m_metric_query_deriv;
    // Only issued by interpolate_shift_data(), for the analytic shift: d_k of
    // K and A_ij, and the six second derivatives of chi and h_ij split across
    // two queries to stay inside the num_components entry limit.
    InterpolationQueryParticle m_metric_query_deriv2;
    InterpolationQueryParticle m_metric_query_d2a;
    InterpolationQueryParticle m_metric_query_d2b;

    std::vector<double> m_theta_vals{};

    // Whether the mat-vec re-interpolates the metric for the PTC step now in
    // progress. Set once per step by the unfreeze policy (see
    // AHFinderParameters) and held constant for the whole GMRES solve, so the
    // Krylov method always sees a single consistent linear operator.
    bool m_unfreeze_now{};

    // Number of PTC steps that used the exact (unfrozen) Jacobian, and the
    // number that actually re-measured c rather than reusing the last value.
    long m_n_unfrozen_steps{};
    long m_n_shift_measured{};

    // The shift c currently capping dt at 1/c, either the fixed
    // ah_finder.newton_shift or the value measured this step. 0 when the cap
    // is disabled.
    double m_newton_shift{};

    // The per-particle diagonal c(x) measured this step, of which
    // m_newton_shift is the mean. Filled by measure_newton_shift(); only read
    // by the mat-vec when ah_finder.newton_shift_diagonal is set, in which case
    // the scalar dt cap is dropped in its favour.
    std::vector<double> m_newton_shift_diag{};

    // Cost counters for the two kernels the PTC step is built from, reported
    // at the end of find(). interpolate_metric() is the expensive one (particle
    // placement, two interpolator queries, an MPI reduce); theta_from_metric()
    // is local tensor algebra plus one reduce. Which of them dominates depends
    // entirely on ah_finder.unfreeze_jacobian: frozen, the metric is
    // interpolated once per PTC step; unfrozen, once per Krylov iteration.
    long m_n_interp{};
    long m_n_theta{};
    double m_t_interp{};
    double m_t_theta{};

    static int local_count(int num_particles)
    {
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        return num_particles / nprocs +
               (myproc < num_particles % nprocs ? 1 : 0);
    }

    static int local_start(int num_particles)
    {
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int myproc = amrex::ParallelDescriptor::MyProc();

        return myproc * (num_particles / nprocs) +
               std::min(myproc, num_particles % nprocs);
    }

    void init_particle_vals();

    // Set particles' coordinates according to their distance from the centre
    void set_particle_positions(const std::vector<double> &h);

    // Grow the pseudo-timestep as the residual falls (SER rule), so the
    // implicit step approaches Newton. The implicit solve is unconditionally
    // stable, so dt is not capped.
    amrex::Real update_dt(amrex::Real dt, double theta_old,
                          double theta_new) const;

    void setup_metric_query();

    // Interpolate and freeze the CCZ4 metric at the surface r = h (expensive:
    // particle placement plus the interpolator queries). Run once per PTC
    // step, at h_n.
    void interpolate_metric(const std::vector<double> &h);

    // Evaluate the expansion Theta at every grid point for surface radius h,
    // using the metric frozen by the last interpolate_metric(). Cheap and
    // purely local (plus one reduce); this is what the Jacobian mat-vec
    // differentiates.
    void theta_from_metric(const std::vector<double> &h,
                           std::vector<double> &theta_out);

    // Apply the finite-difference directional derivative J v about the current
    // surface m_state.h, with the metric either frozen at h_n or
    // re-interpolated at the perturbed surface. Shared by the GMRES mat-vec and
    // the spectral probe; the I/dt term is *not* included.
    void jacobian_apply(const std::vector<double> &v, bool unfreeze,
                        std::vector<double> &jv);

    // Rayleigh-quotient probe of the frozen and unfrozen Jacobians on a few
    // low-order modes, printed when ah_finder.jacobian_diagnostic is set.
    void jacobian_diagnostic(int n_iter);

    // Estimate the constant shift c in J_exact ~= J_frozen + c I at the current
    // surface, from the l=0 mode. Costs two metric interpolations.
    double measure_newton_shift(int n_iter);

    // Interpolate the extra derivative data the analytic shift needs, at the
    // particle positions the last interpolate_metric() already placed.
    void interpolate_shift_data();

    // Same quantity as measure_newton_shift(), but obtained by transporting the
    // frozen fields to the displaced radius with their own derivatives instead
    // of re-interpolating at a perturbed surface. Costs three partial queries
    // and two theta evaluations rather than two full metric interpolations.
    double measure_newton_shift_analytic(int n_iter);

    double inf_norm(std::vector<double>);

  public:

    using Base         = ParticleInterpolator<num_components>;
    using ParIterType  = typename Base::ParIterType;
    using ParticleType = typename Base::ParticleType;
    using Base::Base;

    AHFinder(int num_particles,
             const std::array<double, AMREX_SPACEDIM> &center,
             double guess_radius = 1.0)
        : m_num_particles(num_particles), m_n_local(local_count(num_particles)),
          m_start(local_start(num_particles)), interp_coords_x(num_particles),
          interp_coords_y(num_particles), interp_coords_z(num_particles),
          m_state(std::vector<double>(num_particles)),
          m_geometry(num_particles, center, guess_radius),
          m_gamma_LL(num_particles), m_theta_n(num_particles),
          m_metric_query_state(m_n_local), m_metric_query_deriv(m_n_local),
          m_metric_query_deriv2(m_n_local), m_metric_query_d2a(m_n_local),
          m_metric_query_d2b(m_n_local), m_theta_vals(num_particles)
    {
    }

    void init(GRAmr *gramr_ptr);

    void find();
};

#include "AHFinder.impl.hpp"

#endif /* AHFINDER_HPP_ */
