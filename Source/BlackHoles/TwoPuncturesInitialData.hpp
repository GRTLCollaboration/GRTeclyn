/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifdef USE_TWOPUNCTURES

#ifndef TWOPUNCTURESINITIALDATA_HPP_
#define TWOPUNCTURESINITIALDATA_HPP_

#include "ArrayTools.hpp"
#include "BoundaryConditions.hpp"
#include "CCZ4Vars.hpp"
#include "Coordinates.hpp"
#include "GRParmParse.hpp"
#include "Tensor.hpp"
#include "TensorAlgebra.hpp"
#include "TwoPunctures.hpp"
#include <array>

//! This compute class stores, manages and runs the TwoPunctures object which it
//! uses to set the initial data.
class TwoPuncturesInitialData
{
  private:
    double m_dx;
    std::array<double, AMREX_SPACEDIM> m_center;
    // Let's assume we only need one TwoPunctures object for the whole
    // simulation
    static inline TP::TwoPunctures s_two_punctures;
    static inline bool s_parameters_read      = false;
    static inline bool s_two_punctures_solved = false;

  public:
    TwoPuncturesInitialData(const double a_dx) : m_dx(a_dx)
    {
        GRParmParse pp;
        pp.get("geometry.center", m_center);
        read_parameters();
    }

    void solve()
    {
        if (s_two_punctures_solved)
        {
            return;
        }
        s_two_punctures.Run();
        s_two_punctures_solved = true;
    }

    void set_bh_params(BoostedBHInitialData::params_t &bh1_params,
                       BoostedBHInitialData::params_t &bh2_params)
    {
        if (!s_two_punctures_solved)
        {
            // Need to solve to get bare masses if using target masses
            // We could use target masses in that case but that would be
            // inconsistent.
            solve();
        }

        bh1_params.mass = s_two_punctures.par_m_plus;
        bh2_params.mass = s_two_punctures.par_m_minus;

        bh1_params.center = m_center;
        bh2_params.center = m_center;
        int offset_dir    = (s_two_punctures.swap_xz) ? 2 : 0;
        bh1_params.center[offset_dir] +=
            s_two_punctures.center_offset[offset_dir] + s_two_punctures.par_b;
        bh2_params.center[offset_dir] +=
            s_two_punctures.center_offset[offset_dir] - s_two_punctures.par_b;
    }

    AMREX_FORCE_INLINE void
    operator()(int ix, int iy, int iz,
               const amrex::Array4<amrex::Real> &a_state) const
    {
        AMREX_ASSERT(s_two_punctures_solved);

        const amrex::CellData<amrex::Real> &state_cell_data =
            a_state.cellData(ix, iy, iz);

        for (int i = 0; i < NUM_CCZ4_VARS; ++i)
        {
            // Set only the non-zero components explicitly below
            state_cell_data[i] = 0.0;
        }

        Coordinates coords(amrex::IntVect(ix, iy, iz), m_dx, m_center);
        Tensor::Sym12Rank2 h_phys, K_tensor;
        Tensor::Rank1 shift, Z3;
        amrex::Real lapse, Theta;

        interpolate_tp_vars(coords, h_phys, K_tensor, state_cell_data[c_lapse],
                            shift, state_cell_data[c_Theta], Z3);

        using namespace TensorAlgebra;
        // analytically set Bowen-York properties below (e.g. conformal
        // flatness, tracefree K)

        // metric variables
        amrex::Real chi = std::pow(compute_determinant_sym(h_phys), -1. / 3.);
        state_cell_data[c_chi] = chi;
        FOR (i)
        {
            // Bowen-York data is conformally flat
            state_cell_data[sym_var_idx(c_h11, i, i)] = 1.0;
        }

        amrex::Real trace_A = 0.0;

        // extrinsic curvature
        FOR (i, j)
        {
            state_cell_data[sym_var_idx(c_A11, i, j)] = chi * K_tensor(i, j);
            // conformal flatness
            trace_A += state_cell_data[sym_var_idx(c_A11, i, j)] *
                       TensorAlgebra::delta(i, j);
        }

        // We choose not to use TensorAlgebra::make_trace_free so that we don't
        // need to instantiate a temporary tensor object.
        amrex::Real one_over_gr_spacedim =
            1. / static_cast<amrex::Real>(GR_SPACEDIM);
        FOR (i)
        {
            // conformal flatness
            state_cell_data[sym_var_idx(c_A11, i, i)] -=
                one_over_gr_spacedim * trace_A;
        }

        // gauge
        FOR (i)
        {
            state_cell_data[c_shift1 + i] = shift(i);
        }

        // Z4 variables
    }

    static void read_parameters()
    {
        if (s_parameters_read)
        {
            return;
        }
        GRParmParse tp_pp("two_punctures");

        s_two_punctures.verbose = false;
        tp_pp.queryAdd("verbose", s_two_punctures.verbose);

        // default to using bare masses rather than solving for target masses
        bool calculate_target_masses = false;
        tp_pp.queryAdd("calculate_target_masses", calculate_target_masses);
        s_two_punctures.give_bare_mass = !calculate_target_masses;

        if (calculate_target_masses)
        {
            tp_pp.get("target_mass_plus", s_two_punctures.target_M_plus);
            tp_pp.get("target_mass_minus", s_two_punctures.target_M_minus);

            s_two_punctures.adm_tol = 1.0e-10;
            tp_pp.queryAdd("adm_tol", s_two_punctures.adm_tol);
            if (s_two_punctures.verbose)
            {
                amrex::Print()
                    << "TwoPunctures will solve for bare masses to "
                       "achieve target ADM masses of "
                    << s_two_punctures.target_M_plus << " and "
                    << s_two_punctures.target_M_minus << " with tolerance "
                    << s_two_punctures.adm_tol << "\n";
            }
        }
        else
        {
            tp_pp.get("mass_plus", s_two_punctures.par_m_plus);
            tp_pp.get("mass_minus", s_two_punctures.par_m_minus);
            if (s_two_punctures.verbose)
            {
                amrex::Print()
                    << "TwoPunctures will use bare masses of "
                    << std::setprecision(16) << s_two_punctures.par_m_plus
                    << " and " << s_two_punctures.par_m_minus << "\n";
                // reset precision
                amrex::Print() << std::setprecision(6);
            }
        }

        std::array<amrex::Real, AMREX_SPACEDIM> momentum_plus{};
        std::array<amrex::Real, AMREX_SPACEDIM> momentum_minus{};
        tp_pp.get("momentum_plus", momentum_plus);
        tp_pp.get("momentum_minus", momentum_minus);

        std::array<amrex::Real, AMREX_SPACEDIM> spin_plus{};
        std::array<amrex::Real, AMREX_SPACEDIM> spin_minus{};
        tp_pp.get("spin_plus", spin_plus);
        tp_pp.get("spin_minus", spin_minus);

        // TwoPunctures implements swap_xz as an exchange of the x and z axes,
        // not as a proper rotation. These momentum and spin components use
        // that internal convention. In particular, with swap_xz enabled, a
        // physical spin in the positive z direction is specified by a
        // negative x component here.
        FOR (i)
        {
            s_two_punctures.par_P_plus[i]  = momentum_plus[i];
            s_two_punctures.par_P_minus[i] = momentum_minus[i];
            s_two_punctures.par_S_plus[i]  = spin_plus[i];
            s_two_punctures.par_S_minus[i] = spin_minus[i];
        }

        if (s_two_punctures.verbose)
        {
            amrex::Print() << "The corresponding momenta are"
                           << "\nP_plus = ";
            FOR (i)
            {
                amrex::Print() << s_two_punctures.par_P_plus[i] << " ";
            }
            amrex::Print() << "\nP_minus = ";
            FOR (i)
            {
                amrex::Print() << s_two_punctures.par_P_minus[i] << " ";
            }

            amrex::Print() << "\nThe corresponding spins are"
                           << "\nS_plus = ";
            FOR (i)
            {
                amrex::Print() << s_two_punctures.par_S_plus[i] << " ";
            }
            amrex::Print() << "\nS_minus = ";
            FOR (i)
            {
                amrex::Print() << s_two_punctures.par_S_minus[i] << " ";
            }
        }

        // default to Taylor expansion interpolation as it is much faster
        bool use_spectral_interpolation = false;
        tp_pp.queryAdd("use_spectral_interpolation",
                       use_spectral_interpolation);
        s_two_punctures.grid_setup_method =
            (use_spectral_interpolation) ? "evaluation" : "Taylor expansion";

        s_two_punctures.initial_lapse = "psi^n";
        tp_pp.queryAdd("initial_lapse", s_two_punctures.initial_lapse);

        if (s_two_punctures.initial_lapse != "twopunctures-antisymmetric" &&
            s_two_punctures.initial_lapse != "twopunctures-averaged" &&
            s_two_punctures.initial_lapse != "psi^n" &&
            s_two_punctures.initial_lapse != "brownsville")
        {
            tp_pp.error("initial_lapse",
                        "must be one of 'twopunctures-antisymmetric', "
                        "'twopunctures-averaged', 'psi^n', or 'brownsville'");
        }
        if (s_two_punctures.initial_lapse == "psi^n")
        {
            s_two_punctures.initial_lapse_psi_exponent = -2.0;
            tp_pp.queryAdd("initial_lapse_psi_exponent",
                           s_two_punctures.initial_lapse_psi_exponent);
        }

        // spectral grid parameters
        s_two_punctures.npoints_A = 30;
        tp_pp.queryAdd("num_points_A", s_two_punctures.npoints_A);
        s_two_punctures.npoints_B = 30;
        tp_pp.queryAdd("num_points_B", s_two_punctures.npoints_B);
        s_two_punctures.npoints_phi = 16;
        tp_pp.queryAdd("num_points_phi", s_two_punctures.npoints_phi);
        if (s_two_punctures.npoints_phi < 4 ||
            s_two_punctures.npoints_phi % 4 != 0)
        {
            tp_pp.error("num_points_phi", "must be >= 4 and divisible by 4");
        }

        // solver parameters
        s_two_punctures.Newton_tol = 1.0e-10;
        tp_pp.queryAdd("solver_tolerance", s_two_punctures.Newton_tol);
        s_two_punctures.Newton_maxit = 5;
        tp_pp.queryAdd("solver_max_iterations", s_two_punctures.Newton_maxit);
        s_two_punctures.TP_epsilon = 1.0e-6;
        tp_pp.queryAdd("epsilon", s_two_punctures.TP_epsilon);
        s_two_punctures.TP_Tiny = 0.0;
        tp_pp.queryAdd("tiny", s_two_punctures.TP_Tiny);
        s_two_punctures.TP_Extend_Radius = 0.0;
        tp_pp.queryAdd("extend_radius", s_two_punctures.TP_Extend_Radius);

        // BH positions
        amrex::Real offset_plus{};
        amrex::Real offset_minus{};
        tp_pp.get("offset_plus", offset_plus);
        tp_pp.get("offset_minus", offset_minus);

        s_two_punctures.swap_xz = false;
        tp_pp.queryAdd("swap_xz", s_two_punctures.swap_xz);

        double center_offset_xz = 0.5 * (offset_plus + offset_minus);
        int offset_dir          = (s_two_punctures.swap_xz) ? 2 : 0;
        s_two_punctures.center_offset[offset_dir] = center_offset_xz;
        s_two_punctures.par_b = 0.5 * (offset_plus - offset_minus);

        // debug output
        s_two_punctures.do_residuum_debug_output = false;
        tp_pp.queryAdd("do_residuum_debug_output",
                       s_two_punctures.do_residuum_debug_output);
        s_two_punctures.do_initial_debug_output = false;
        tp_pp.queryAdd("do_initial_debug_output",
                       s_two_punctures.do_initial_debug_output);

        // Irrelevant parameters set to default value
        s_two_punctures.keep_u_around                   = false;
        s_two_punctures.use_sources                     = false;
        s_two_punctures.rescale_sources                 = true;
        s_two_punctures.use_external_initial_guess      = false;
        s_two_punctures.multiply_old_lapse              = false;
        s_two_punctures.schedule_in_ADMBase_InitialData = true;
        s_two_punctures.solve_momentum_constraint       = false;
        s_two_punctures.metric_type                     = "something else";
        s_two_punctures.conformal_storage = "not conformal at all";
        s_two_punctures.conformal_state   = 0;
        s_two_punctures.mp                = 0;
        s_two_punctures.mm                = 0;
        s_two_punctures.mp_adm            = 0;
        s_two_punctures.mm_adm            = 0;

        s_parameters_read = true;
    }

    static void check_params()
    {
        // These checks are mostly taken from the Einstein Toolkit thorn
        // documentation:
        // https://einsteintoolkit.org/thornguide/EinsteinInitialData/TwoPunctures/documentation.html

        GRParmParse tp_pp("two_punctures");

        bool verbose = false;
        tp_pp.queryAdd("verbose", verbose);

        bool calculate_target_masses = false;
        tp_pp.queryAdd("calculate_target_masses", calculate_target_masses);

        double mass_plus{};
        double mass_minus{};
        if (calculate_target_masses)
        {
            tp_pp.get("target_mass_plus", mass_plus);
            tp_pp.get("target_mass_minus", mass_minus);

            double adm_tol = 1.0e-10;
            tp_pp.queryAdd("adm_tol", adm_tol);
            if (adm_tol <= 0.0)
            {
                tp_pp.error("adm_tol", "must be > 0.0");
            }
        }
        else
        {
            tp_pp.get("mass_plus", mass_plus);
            tp_pp.get("mass_minus", mass_minus);
        }

        const char *mass_plus_name =
            calculate_target_masses ? "target_mass_plus" : "mass_plus";
        const char *mass_minus_name =
            calculate_target_masses ? "target_mass_minus" : "mass_minus";
        if (mass_plus < 0.0)
        {
            tp_pp.error(mass_plus_name, "must be >= 0.0");
        }
        if (mass_minus < 0.0)
        {
            tp_pp.error(mass_minus_name, "must be >= 0.0");
        }

        std::array<amrex::Real, AMREX_SPACEDIM> momentum_plus{};
        std::array<amrex::Real, AMREX_SPACEDIM> momentum_minus{};
        std::array<amrex::Real, AMREX_SPACEDIM> spin_plus{};
        std::array<amrex::Real, AMREX_SPACEDIM> spin_minus{};
        tp_pp.get("momentum_plus", momentum_plus);
        tp_pp.get("momentum_minus", momentum_minus);
        tp_pp.get("spin_plus", spin_plus);
        tp_pp.get("spin_minus", spin_minus);

        // The target masses are puncture ADM masses, so |P| / M must remain
        // below one. A bare puncture mass is not the physical ADM mass, so it
        // can only support a warning for the same nominal ratio.
        auto check_momentum = [&tp_pp, calculate_target_masses](
                                  const char *momentum_name,
                                  const auto &momentum, const char *mass_name,
                                  double mass)
        {
            const double momentum_magnitude =
                std::sqrt(ArrayTools::norm2(momentum));
            if (momentum_magnitude == 0.0)
            {
                return;
            }
            if (mass == 0.0)
            {
                tp_pp.error(momentum_name,
                            "must be zero when the mass is zero");
            }
            if (momentum_magnitude >= mass)
            {
                const std::string message =
                    "gives |P| / " + std::string(mass_name) + " >= 1";
                if (calculate_target_masses)
                {
                    tp_pp.error(momentum_name, message);
                }
                else
                {
                    tp_pp.warning(momentum_name,
                                  message +
                                      ", but this uses a bare mass rather than "
                                      "the physical ADM mass");
                }
            }
        };
        check_momentum("momentum_plus", momentum_plus, mass_plus_name,
                       mass_plus);
        check_momentum("momentum_minus", momentum_minus, mass_minus_name,
                       mass_minus);

        auto check_spin = [&tp_pp, calculate_target_masses](
                              const char *spin_name, const auto &spin,
                              const char *mass_name, double mass)
        {
            const double spin_magnitude = std::sqrt(ArrayTools::norm2(spin));
            if (spin_magnitude == 0.0)
            {
                return;
            }
            if (mass == 0.0)
            {
                tp_pp.error(spin_name, "must be zero when the mass is zero");
            }
            if (spin_magnitude >= mass * mass)
            {
                const std::string message =
                    "gives |J| / " + std::string(mass_name) + "^2 >= 1";
                if (calculate_target_masses)
                {
                    tp_pp.error(spin_name, message);
                }
                else
                {
                    tp_pp.warning(spin_name,
                                  message +
                                      ", but this uses a bare mass rather than "
                                      "the physical black-hole mass");
                }
            }
        };
        check_spin("spin_plus", spin_plus, mass_plus_name, mass_plus);
        check_spin("spin_minus", spin_minus, mass_minus_name, mass_minus);

        bool use_spectral_interpolation = false;
        tp_pp.queryAdd("use_spectral_interpolation",
                       use_spectral_interpolation);

        std::string initial_lapse = "psi^n";
        tp_pp.queryAdd("initial_lapse", initial_lapse);
        if (initial_lapse != "twopunctures-antisymmetric" &&
            initial_lapse != "twopunctures-averaged" &&
            initial_lapse != "psi^n" && initial_lapse != "brownsville")
        {
            tp_pp.error("initial_lapse",
                        "must be one of 'twopunctures-antisymmetric', "
                        "'twopunctures-averaged', 'psi^n', or 'brownsville'");
        }
        if (initial_lapse == "psi^n")
        {
            double initial_lapse_psi_exponent = -2.0;
            tp_pp.queryAdd("initial_lapse_psi_exponent",
                           initial_lapse_psi_exponent);
        }

        int num_points_a = 30;
        tp_pp.queryAdd("num_points_A", num_points_a);
        if (num_points_a < 4)
        {
            tp_pp.error("num_points_A", "must be >= 4");
        }

        int num_points_b = 30;
        tp_pp.queryAdd("num_points_B", num_points_b);
        if (num_points_b < 4)
        {
            tp_pp.error("num_points_B", "must be >= 4");
        }

        int num_points_phi = 16;
        tp_pp.queryAdd("num_points_phi", num_points_phi);
        if (num_points_phi < 4 || num_points_phi % 4 != 0)
        {
            tp_pp.error("num_points_phi", "must be >= 4 and divisible by 4");
        }

        int solver_max_iterations = 5;
        tp_pp.queryAdd("solver_max_iterations", solver_max_iterations);
        if (solver_max_iterations < 0)
        {
            tp_pp.error("solver_max_iterations", "must be >= 0");
        }

        double solver_tolerance = 1.0e-10;
        tp_pp.queryAdd("solver_tolerance", solver_tolerance);
        if (solver_tolerance < 0.0)
        {
            tp_pp.error("solver_tolerance", "must be >= 0.0");
        }

        double epsilon = 1.0e-6;
        tp_pp.queryAdd("epsilon", epsilon);
        if (epsilon < 0.0)
        {
            tp_pp.error("epsilon", "must be >= 0.0");
        }

        double tiny = 0.0;
        tp_pp.queryAdd("tiny", tiny);
        if (tiny < 0.0)
        {
            tp_pp.error("tiny", "must be >= 0.0");
        }

        double extend_radius = 0.0;
        tp_pp.queryAdd("extend_radius", extend_radius);
        if (extend_radius < 0.0)
        {
            tp_pp.error("extend_radius", "must be >= 0.0");
        }

        double offset_plus{};
        double offset_minus{};
        tp_pp.get("offset_plus", offset_plus);
        tp_pp.get("offset_minus", offset_minus);

        bool swap_xz = false;
        tp_pp.queryAdd("swap_xz", swap_xz);

        // TwoPunctures always solves with the punctures on its internal x
        // axis. When swap_xz is enabled, it exchanges x and z during
        // interpolation, so the physical separation and center offset belong
        // in the z component here.
        GRParmParse geom_pp("geometry");
        std::array<double, AMREX_SPACEDIM> center{};
        std::array<double, AMREX_SPACEDIM> prob_extent{};
        geom_pp.get("center", center);
        geom_pp.get("prob_extent", prob_extent);

        GRParmParse boundary_pp("boundary");
        const auto lo_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "lo_condition");
        const auto hi_condition = BoundaryConditions::params_t::read_conditions(
            boundary_pp, "hi_condition");

        const int offset_dir = swap_xz ? 2 : 0;
        const double domain_lo =
            lo_condition[offset_dir] == BoundaryConditions::REFLECTIVE_BC
                ? -prob_extent[offset_dir]
                : 0.0;
        const double domain_hi =
            hi_condition[offset_dir] == BoundaryConditions::REFLECTIVE_BC
                ? 2.0 * prob_extent[offset_dir]
                : prob_extent[offset_dir];
        if (center[offset_dir] + offset_plus < domain_lo ||
            center[offset_dir] + offset_plus > domain_hi)
        {
            tp_pp.warning("offset_plus",
                          "places the puncture outside the computational "
                          "domain after applying reflective symmetry");
        }
        if (center[offset_dir] + offset_minus < domain_lo ||
            center[offset_dir] + offset_minus > domain_hi)
        {
            tp_pp.warning("offset_minus",
                          "places the puncture outside the computational "
                          "domain after applying reflective symmetry");
        }

        bool do_residuum_debug_output = false;
        tp_pp.queryAdd("do_residuum_debug_output", do_residuum_debug_output);
        bool do_initial_debug_output = false;
        tp_pp.queryAdd("do_initial_debug_output", do_initial_debug_output);
    }

    static const TP::TwoPunctures &get_two_punctures()
    {
        return s_two_punctures;
    }

  private:
    void interpolate_tp_vars(const Coordinates &coords,
                             Tensor::Sym12Rank2 &out_h_phys,
                             Tensor::Sym12Rank2 &out_K_tensor,
                             amrex::Real &out_lapse, Tensor::Rank1 &out_shift,
                             amrex::Real &out_Theta,
                             Tensor::Rank1 &out_Z3) const
    {
        amrex::Real coords_array[AMREX_SPACEDIM];
        coords_array[0] = coords.x;
        coords_array[1] = coords.y;
        coords_array[2] = coords.z;

        using namespace TP::Z4VectorShortcuts;
        amrex::Real TP_state[Qlen];
        s_two_punctures.Interpolate(coords_array, TP_state);

        // metric
        out_h_phys(0, 0) = TP_state[g11];
        out_h_phys(0, 1) = TP_state[g12];
        out_h_phys(0, 2) = TP_state[g13];
        out_h_phys(1, 1) = TP_state[g22];
        out_h_phys(1, 2) = TP_state[g23];
        out_h_phys(2, 2) = TP_state[g33];

        // extrinsic curvature
        out_K_tensor(0, 0) = TP_state[K11];
        out_K_tensor(0, 1) = TP_state[K12];
        out_K_tensor(0, 2) = TP_state[K13];
        out_K_tensor(1, 1) = TP_state[K22];
        out_K_tensor(1, 2) = TP_state[K23];
        out_K_tensor(2, 2) = TP_state[K33];

        // Z4 vector
        out_Z3(0) = TP_state[Z1];
        out_Z3(1) = TP_state[Z2];
        out_Z3(2) = TP_state[Z3];
        out_Theta = TP_state[Theta];

        // gauge
        out_lapse    = TP_state[lapse];
        out_shift(0) = TP_state[shift1];
        out_shift(1) = TP_state[shift2];
        out_shift(2) = TP_state[shift3];
    }
};

#endif /* TWOPUNCTURESINITIALDATA_HPP_ */
#endif /* USE_TWOPUNCTURES */
