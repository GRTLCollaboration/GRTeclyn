/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef SIMULATIONPARAMETERS_HPP
#define SIMULATIONPARAMETERS_HPP

// General includes
#include "BaseParameterChecker.hpp"
#include "GRParmParse.hpp"

// Problem specific includes:
#include "ArrayTools.hpp"
#include "BoostedBHInitialData.hpp"
#include "PunctureTracker.hpp"
#include "CCZ4RHS.hpp"
#ifdef USE_TWOPUNCTURES
#include "TP_Parameters.hpp"
#endif

class SimulationParameters : public BaseParameterChecker
{
  public:
    // NOLINTNEXTLINE(readability-identifier-length)
    SimulationParameters() : BaseParameterChecker()
    {
        GRParmParse pp;

        read_shared_params(pp);
#ifdef USE_TWOPUNCTURES
        read_tp_params(pp);
#else
        BoostedBHInitialData::params_t::check_params(1);
        BoostedBHInitialData::params_t::check_params(2);
#endif
        check_params();
    }

    /// Read shared parameters
    // NOLINTNEXTLINE(readability-identifier-length)
    void read_shared_params(GRParmParse &pp)
    {

        int formulation = CCZ4RHS<>::USE_CCZ4; // Whether to use BSSN or CCZ4
        pp.queryAdd("formulation", formulation);

        if (formulation != CCZ4RHS<>::USE_CCZ4 &&
            formulation != CCZ4RHS<>::USE_BSSN)
        {
            pp.error("formulation", "must be 0 or 1");
        }

        if (formulation == CCZ4RHS<>::USE_CCZ4)
        {
            CCZ4_params_t::check_params();
        }
        else if (formulation == CCZ4RHS<>::USE_BSSN)
        {
            if (pp.contains("ccz4.kappa1") || pp.contains("ccz4.kappa2") ||
                pp.contains("ccz4.kappa3"))
            {
                pp.warning("kappa1/2/3",
                           "should not be provided with BSSN formulation, "
                           "setting them all to zero");
            }
            pp.add("ccz4.kappa1", 0.0);
            pp.add("ccz4.kappa2", 0.0);
            pp.add("ccz4.kappa3", 0.0);
        }
        
        // Do we want puncture tracking and constraint norm calculation?
        bool puncture_tracking_enabled{false};
        pp.queryAdd("puncture_tracking.enabled", puncture_tracking_enabled);
        if (puncture_tracking_enabled)
        {
            puncture_tracker_params_t::check_params();
        }

        bool calculate_constraint_norms = false;
        pp.queryAdd("calculate_constraint_norms", calculate_constraint_norms);
    }

#ifdef USE_TWOPUNCTURES
    void read_tp_params(GRParmParse &pp)
    {
        pp.get("amrex.verbose", tp_params.verbose);
        // check whether to calculate the target ADM masses or use provided bare
        // masses
        bool calculate_target_masses;
        pp.load("TP_calculate_target_masses", calculate_target_masses, false);
        tp_params.give_bare_mass = !calculate_target_masses;

        // masses
        if (calculate_target_masses)
        {
            pp.load("TP_target_mass_plus", tp_params.target_M_plus);
            pp.load("TP_target_mass_minus", tp_params.target_M_minus);
            pp.load("TP_adm_tol", tp_params.adm_tol, 1e-10);
            amrex::Print() << "The black holes have target ADM masses of "
                           << tp_params.target_M_plus << " and "
                           << tp_params.target_M_minus << "\n";
            bh1_params.mass = tp_params.target_M_minus;
            bh2_params.mass = tp_params.target_M_plus;
        }
        else
        {
            pp.load("TP_mass_plus", tp_params.par_m_plus);
            pp.load("TP_mass_minus", tp_params.par_m_minus);
            bh1_params.mass = tp_params.par_m_plus;
            bh2_params.mass = tp_params.par_m_minus;
            amrex::Print() << "The black holes have bare masses of "
                           << std::setprecision(16) << tp_params.par_m_plus
                           << " and " << tp_params.par_m_minus << "\n";
            // reset precision
            amrex::Print() << std::setprecision(6);
        }

        // BH spin and momenta
        std::array<double, AMREX_SPACEDIM> spin_minus, spin_plus;
        pp.load("TP_momentum_minus", bh1_params.momentum);
        pp.load("TP_momentum_plus", bh2_params.momentum);
        pp.load("TP_spin_plus", spin_plus);
        pp.load("TP_spin_minus", spin_minus);
        FOR (i)
        {
            tp_params.par_P_minus[i] = bh1_params.momentum[i];
            tp_params.par_P_plus[i]  = bh2_params.momentum[i];
            tp_params.par_S_minus[i] = spin_minus[i];
            tp_params.par_S_plus[i]  = spin_plus[i];
        }

        amrex::Print() << "The corresponding momenta are:";
        amrex::Print() << "\nP_plus = ";
        FOR (i)
        {
            amrex::Print() << tp_params.par_P_plus[i] << " ";
        }
        amrex::Print() << "\nP_minus = ";
        FOR (i)
        {
            amrex::Print() << tp_params.par_P_minus[i] << " ";
        }

        amrex::Print() << "\nThe corresponding spins are:";
        amrex::Print() << "\nS_plus = ";
        FOR (i)
        {
            amrex::Print() << tp_params.par_S_plus[i] << " ";
        }
        amrex::Print() << "\nS_minus = ";
        FOR (i)
        {
            amrex::Print() << tp_params.par_S_minus[i] << " ";
        }
        amrex::Print() << "\n";

        // interpolation type
        bool use_spectral_interpolation;
        pp.load("TP_use_spectral_interpolation", use_spectral_interpolation,
                false);
        tp_params.grid_setup_method =
            (use_spectral_interpolation) ? "evaluation" : "Taylor expansion";

        // initial_lapse (default to psi^n)
        pp.load("TP_initial_lapse", tp_params.initial_lapse,
                std::string("psi^n"));
        if (tp_params.initial_lapse != "twopunctures-antisymmetric" &&
            tp_params.initial_lapse != "twopunctures-averaged" &&
            tp_params.initial_lapse != "psi^n" &&
            tp_params.initial_lapse != "brownsville")
        {
            std::string message  = "Parameter: TP_initial_lapse: ";
            message             += tp_params.initial_lapse;
            message             += " invalid";
            amrex::Abort(message.c_str());
        }
        if (tp_params.initial_lapse == "psi^n")
        {
            pp.load("TP_initial_lapse_psi_exponent",
                    tp_params.initial_lapse_psi_exponent, -2.0);
        }

        // Spectral grid parameters
        pp.load("TP_npoints_A", tp_params.npoints_A, 30);
        pp.load("TP_npoints_B", tp_params.npoints_B, 30);
        pp.load("TP_npoints_phi", tp_params.npoints_phi, 16);
        if (tp_params.npoints_phi % 4 != 0)
        {
            amrex::Abort("TP_npoints_phi must be a multiple of 4");
        }

        // Solver parameters and tolerances
        pp.load("TP_Newton_tol", tp_params.Newton_tol, 1e-10);
        pp.load("TP_Newton_maxit", tp_params.Newton_maxit, 5);
        pp.load("TP_epsilon", tp_params.TP_epsilon, 1e-6);
        pp.load("TP_Tiny", tp_params.TP_Tiny, 0.0);
        pp.load("TP_Extend_Radius", tp_params.TP_Extend_Radius, 0.0);

        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("amr.center", center);
        // BH positions
        pp.load("TP_offset_minus", tp_offset_minus);
        pp.load("TP_offset_plus", tp_offset_plus);
        bh1_params.center           = center;
        bh2_params.center           = center;
        bh1_params.center[0]       += tp_offset_minus;
        bh2_params.center[0]       += tp_offset_plus;
        double center_offset_x      = 0.5 * (tp_offset_plus + tp_offset_minus);
        tp_params.center_offset[0]  = center_offset_x;
        // par_b is half the distance between BH_minus and BH_plus
        tp_params.par_b = 0.5 * (tp_offset_plus - tp_offset_minus);
        pp.load("TP_swap_xz", tp_params.swap_xz, false);

        // Debug output
        pp.load("TP_do_residuum_debug_output",
                tp_params.do_residuum_debug_output, false);
        pp.load("TP_do_initial_debug_output", tp_params.do_initial_debug_output,
                false);

        // Irrelevant parameters set to default value
        tp_params.keep_u_around                   = false;
        tp_params.use_sources                     = false;
        tp_params.rescale_sources                 = true;
        tp_params.use_external_initial_guess      = false;
        tp_params.multiply_old_lapse              = false;
        tp_params.schedule_in_ADMBase_InitialData = true;
        tp_params.solve_momentum_constraint       = false;
        tp_params.metric_type                     = "something else";
        tp_params.conformal_storage               = "not conformal at all";
        tp_params.conformal_state                 = 0;
        tp_params.mp                              = 0;
        tp_params.mm                              = 0;
        tp_params.mp_adm                          = 0;
        tp_params.mm_adm                          = 0;
    }
#endif /* USE_TWOPUNCTURES */

    void check_params()
    {
#ifdef USE_TWOPUNCTURES
        // These checks are mostly taken from the Einstein Toolkit thorn
        // documentation:
        // https://einsteintoolkit.org/thornguide/EinsteinInitialData/TwoPunctures/documentation.html
        std::string mass_plus_name, mass_minus_name;
        if (tp_params.give_bare_mass)
        {
            mass_minus_name = "TP_mass_minus";
            mass_plus_name  = "TP_mass_plus";
        }
        else
        {
            mass_minus_name = "TP_target_mass_minus";
            mass_plus_name  = "TP_target_mass_plus";
            check_parameter("TP_adm_tol", tp_params.adm_tol,
                            tp_params.adm_tol > 0., "must be > 0.0");
        }
        check_parameter(mass_minus_name, bh1_params.mass, bh1_params.mass >= 0.,
                        "mustd be >= 0.0");
        check_parameter(mass_plus_name, bh2_params.mass, bh2_params.mass >= 0.,
                        "must be >= 0.0");

        int offset_dir = (!tp_params.swap_xz) ? 0 : 2;
        std::array<double, AMREX_SPACEDIM> center{};
        pp.get("amr.center", center);
        warn_parameter("TP_offset_minus", tp_offset_minus,
                       tp_offset_minus < (ivN[offset_dir] + 1) * coarsest_dx -
                                             center[offset_dir],
                       "should be within the computational domain");
        warn_parameter("TP_offset_plus", tp_offset_plus,
                       tp_offset_plus < (ivN[offset_dir] + 1) * coarsest_dx -
                                            center[offset_dir],
                       "should be within the computational domain");
        check_parameter("TP_npoints_A", tp_params.npoints_A,
                        tp_params.npoints_A >= 4, "must be >= 4");
        check_parameter("TP_npoints_B", tp_params.npoints_B,
                        tp_params.npoints_B >= 4, "must be >= 4");
        check_parameter("TP_npoints_phi", tp_params.npoints_phi,
                        tp_params.npoints_phi >= 4 &&
                            tp_params.npoints_phi % 2 == 0,
                        "must be >= 4 and divisible by 2");
        check_parameter("TP_Newton_maxit", tp_params.Newton_maxit,
                        tp_params.Newton_maxit >= 0, "must be >= 0");
        check_parameter("TP_Newton_tol", tp_params.Newton_tol,
                        tp_params.Newton_tol >= 0., "must be >= 0.0");
        check_parameter("TP_epsilon", tp_params.TP_epsilon,
                        tp_params.TP_epsilon >= 0., "must be >= 0.0");
        check_parameter("TP_Tiny", tp_params.TP_Tiny, tp_params.TP_Tiny >= 0.,
                        "must be >= 0.0");
        check_parameter("TP_Extend_Radius", tp_params.TP_Extend_Radius,
                        tp_params.TP_Extend_Radius >= 0., "must be >= 0.0");
#endif /* USE_TWOPUNCTURES */
    }

#ifdef USE_TWOPUNCTURES
    double tp_offset_plus, tp_offset_minus;
    TP::Parameters tp_params;
#endif
    SphericalExtraction::params_t
        extraction_params; // TODO: Remove once extraction is fixed
};

#endif /* SIMULATIONPARAMETERS_HPP */
