#include "SupportedWormholeLevel.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ChiTagger.hpp"
#include "ConstraintsWithMatter.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "RLMatterPumpParams.hpp"
#include "SmallDataIO.hpp"
#include "SpongeZone.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4WithMatter.hpp"
#include "WeylExtraction.hpp"
#include "EffectiveTeoMatter.hpp"
#include "ExternalGridInitialData.hpp"
#include "SupportedWormholeInitialData.hpp"
#include "ExoticScalarField.hpp"
#include "ComplexExoticScalarField.hpp"
#include "NoMatter.hpp"
#include "PhantomDecayPotential.hpp"
#include "ComplexScalarPotential.hpp"
#include "OscillonPotential.hpp"
#include "DustMatter.hpp"

#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <cmath>

//! Build the multi-site matter pump for the wormhole's complex_scalar branch.
//! Mirrors RadialRecipeMatterDispatch::build_rl_pump but uses the wormhole
//! SimulationParameters (no recipe_scalar_field_signs; single complex field).
//! When trajectory_mode == 0 every amplitude is 0, so the pump is a no-op.
inline RLMatterPumpParams
build_wormhole_pump(const SimulationParameters &params)
{
    RLMatterPumpParams pump;
    int n = params.trajectory_params.num_lumps;
    if (params.trajectory_mode != 1)
        n = 0;
    if (n > RL_MAX_LUMPS)
        n = RL_MAX_LUMPS;
    if (n < 0)
        n = 0;
    pump.num_sites       = n;
    pump.width           = params.rl_pump_width;
    pump.governor_center = params.rl_l2_ham_governor_center;
    pump.governor_width  = params.rl_l2_ham_governor_width;
    pump.governor        = RLRuntime::tanh_governor(
        RLRuntime::cached_L2_Ham(), pump.governor_center, pump.governor_width);
    pump.num_fields      = 1;
    pump.k_p             = params.rl_pump_kp;
    pump.k_d             = params.rl_pump_kd;
    pump.target_profile  = params.rl_pump_target_profile;
    pump.target_width    = params.rl_pump_target_width;
    pump.target_amp      = params.rl_pump_target_amp;
    for (int s = 0; s < n; ++s)
    {
        pump.sites[s].center_x  = RLRuntime::g_lump_state[s].x;
        pump.sites[s].center_y  = RLRuntime::g_lump_state[s].y;
        pump.sites[s].center_z  = RLRuntime::g_lump_state[s].z;
        pump.sites[s].amplitude = params.rl_pump_amplitude[s];
        pump.sites[s].frequency = params.rl_pump_frequency[s];
        pump.sites[s].phase     = params.rl_pump_phase[s];
        pump.sites[s].field_sign = 1; // single complex field
    }
    return pump;
}

void SupportedWormholeLevel::variableSetUp()
{
    BL_PROFILE("SupportedWormholeLevel::variableSetUp()");
    stateVariableSetUp();

    std::string matter_model = "exotic_scalar";
    {
        GRParmParse pp;
        pp.load("wormhole_matter_model", matter_model,
                std::string("exotic_scalar"));
    }

    if (matter_model == "no_matter")
    {
        ConstraintsWithMatter<NoMatter>::set_up(state_index);
        Weyl4WithMatter<NoMatter>::set_up(state_index);
    }
    else if (matter_model == "effective_teo")
    {
        ConstraintsWithMatter<EffectiveTeoMatter>::set_up(state_index);
        Weyl4WithMatter<EffectiveTeoMatter>::set_up(state_index);
    }
    else if (matter_model == "dust")
    {
        ConstraintsWithMatter<DustMatter>::set_up(state_index);
        Weyl4WithMatter<DustMatter>::set_up(state_index);
    }
    else if (matter_model == "oscillon_scalar")
    {
        ConstraintsWithMatter<ExoticScalarField<OscillonPotential>>::set_up(
            state_index);
        Weyl4WithMatter<ExoticScalarField<OscillonPotential>>::set_up(
            state_index);
    }
    else if (matter_model == "complex_scalar")
    {
        ConstraintsWithMatter<
            ComplexExoticScalarField<ComplexScalarPotential>>::set_up(
            state_index);
        Weyl4WithMatter<
            ComplexExoticScalarField<ComplexScalarPotential>>::set_up(
            state_index);
    }
    else
    {
        ConstraintsWithMatter<ExoticScalarField<PhantomDecayPotential>>::set_up(
            state_index);
        Weyl4WithMatter<ExoticScalarField<PhantomDecayPotential>>::set_up(
            state_index);
    }
}

void SupportedWormholeLevel::specificAdvance()
{
    amrex::MultiFab &S_new = get_new_data(state_index);
    const auto &arrs       = S_new.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);

    amrex::ParallelFor(S_new,
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, arrs[box_no]);
                           positive_chi_lapse(i, j, k, arrs[box_no]);
                       });
}

void SupportedWormholeLevel::initData()
{
    BL_PROFILE("SupportedWormholeLevel::initData");

    amrex::MultiFab &state = get_new_data(state_index);
    const auto &arrs       = state.arrays();

    if (!simParams().recipe_initial_data_file.empty())
    {
        ExternalGridInitialData ext_data(simParams().external_grid_params,
                                         Geom().CellSize(0));

        // Precollapsed-lapse override for the loaded (GRTresna) initial data.
        // The maximal-slicing solve writes lapse == 1, which relaxes violently
        // under the 1+log condition (a t~0.5 max_K spike that kicks the phantom
        // throat off its unstable equilibrium).  Seeding a precollapsed lapse
        // alpha = chi^p brings the gauge closer to its evolved value and damps
        // that transient.  lapse_type 0 keeps the loaded lapse (backward compat).
        const int lapse_type = simParams().wormhole_params.initial_lapse_type;

        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                for (int n = 0; n < NUM_VARS; ++n)
                {
                    arrs[box_no](i, j, k, n) = 0.;
                }
                ext_data.compute(i, j, k, arrs[box_no]);

                if (lapse_type != 0)
                {
                    const amrex::Real chi = amrex::max(
                        arrs[box_no](i, j, k, c_chi), amrex::Real(1.0e-10));
                    amrex::Real lapse = arrs[box_no](i, j, k, c_lapse);
                    if (lapse_type == 1)
                        lapse = std::sqrt(chi);
                    else if (lapse_type == 2)
                        lapse = amrex::Real(1.0) - amrex::Real(3.0) * std::log(chi);
                    else if (lapse_type == 3)
                        lapse = chi;
                    arrs[box_no](i, j, k, c_lapse) =
                        amrex::max(lapse, amrex::Real(1.0e-10));
                }
            });
    }
    else
    {
        SupportedWormholeInitialData wormhole(simParams().wormhole_params,
                                             Geom().CellSize(0));

        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                amrex::CellData<amrex::Real> cell =
                    arrs[box_no].cellData(i, j, k);
                for (int n = 0; n < cell.nComp(); ++n)
                {
                    cell[n] = 0.;
                }
                wormhole.compute(i, j, k, arrs[box_no]);
            });
    }

    amrex::Gpu::streamSynchronize();
}

void SupportedWormholeLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                    amrex::MultiFab &a_rhs,
                                    const double a_time)
{
    BL_PROFILE("SupportedWormholeLevel::specificEvalRHS()");
    const int soln_ghosts = a_soln.nGrowVect()[0];
    if (soln_ghosts > 0)
    {
        FillPatch(*this, a_soln, soln_ghosts, a_time, state_index, 0,
                  a_soln.nComp());
    }
    const auto &soln_arrs   = a_soln.arrays();
    const auto &soln_c_arrs = a_soln.const_arrays();
    const auto &rhs_arrs    = a_rhs.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);

    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    if (simParams().wormhole_matter_model == "no_matter")
    {
        NoMatter no_matter;
        CCZ4RHSWithMatter<NoMatter, MovingPunctureGaugeWithMatter,
                          FourthOrderDerivatives>
            ccz4rhs(no_matter, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
    else if (simParams().wormhole_matter_model == "effective_teo")
    {
        EffectiveTeoMatter teo_matter;
        CCZ4RHSWithMatter<EffectiveTeoMatter, MovingPunctureGaugeWithMatter,
                          FourthOrderDerivatives>
            ccz4rhs(teo_matter, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
    else if (simParams().wormhole_matter_model == "dust")
    {
        DustMatter dust_matter;
        CCZ4RHSWithMatter<DustMatter, MovingPunctureGaugeWithMatter,
                          FourthOrderDerivatives>
            ccz4rhs(dust_matter, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
    else if (simParams().wormhole_matter_model == "oscillon_scalar")
    {
        OscillonPotential potential;
        ExoticScalarField<OscillonPotential> oscillon_scalar(
            potential, simParams().wormhole_params.support_strength);
        CCZ4RHSWithMatter<ExoticScalarField<OscillonPotential>,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(oscillon_scalar, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
    else if (simParams().wormhole_matter_model == "complex_scalar")
    {
        ComplexScalarPotential potential(
            simParams().wormhole_params.phantom_mass,
            simParams().wormhole_params.phantom_lambda,
            simParams().wormhole_params.phantom_mu6);
        const RLMatterPumpParams pump = build_wormhole_pump(simParams());
        ComplexExoticScalarField<ComplexScalarPotential> complex_scalar(
            potential, simParams().wormhole_params.support_strength, pump);
        CCZ4RHSWithMatter<ComplexExoticScalarField<ComplexScalarPotential>,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(complex_scalar, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }
    else
    {
        PhantomDecayPotential potential(simParams().wormhole_params.phantom_mass);
        ExoticScalarField<PhantomDecayPotential> exotic_scalar(
            potential, simParams().wormhole_params.support_strength);
        CCZ4RHSWithMatter<ExoticScalarField<PhantomDecayPotential>,
                          MovingPunctureGaugeWithMatter, FourthOrderDerivatives>
            ccz4rhs(exotic_scalar, simParams().ccz4_params, Geom().CellSize(0),
                    simParams().sigma, simParams().formulation, 1.0,
                    simParams().wormhole_params.grid_center, a_time);

        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { ccz4rhs(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }

    // Sponge zone: extra radially-ramped KO dissipation in the outer shell to
    // absorb outgoing waves before they reach the Sommerfeld boundary (clean
    // GW extraction on a large box).  Applied once, after the matter-model
    // dispatch, as a second additive pass (matches RadialRecipe).  No-op when
    // sponge_enabled = false.
    if (simParams().sponge_params.enabled)
    {
        const SpongeZone sponge(simParams().sponge_params, Geom().CellSize(0));
        amrex::ParallelFor(
            a_rhs, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            { sponge.apply(i, j, k, rhs_arrs[box_no], soln_c_arrs[box_no]); });
    }

    amrex::Gpu::streamSynchronize();
}

void SupportedWormholeLevel::specificUpdateODE(amrex::MultiFab &a_soln)
{
    const auto &soln_arrs = a_soln.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);
    amrex::ParallelFor(a_soln, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    amrex::Gpu::streamSynchronize();
}

void SupportedWormholeLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto cur_time        = get_state_data(state_index).curTime();
    FillPatch(*this, state_new, 2, cur_time, state_index, c_chi, 1);
}

void SupportedWormholeLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                              amrex::Real a_regrid_threshold)
{
    BL_PROFILE("SupportedWormholeLevel::tag_cells()");
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto &tag_arrs       = a_tag_box_array.arrays();
    const auto &state_new_arrs = state_new.const_arrays();

    ChiTagger chi_tagger(Geom().CellSize(0), a_regrid_threshold);

    amrex::ParallelFor(state_new, amrex::IntVect(0),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           const auto &tags_arr  = tag_arrs[box_no];
                           const auto &state_arr = state_new_arrs[box_no];
                           chi_tagger(i, j, k, tags_arr, state_arr);
                       });
    amrex::Gpu::streamSynchronize();
}

BHAMR<SupportedWormholeLevel::num_punctures> *
SupportedWormholeLevel::get_bhamr_ptr()
{
    return dynamic_cast<BHAMR<num_punctures> *>(get_gramr_ptr());
}

void SupportedWormholeLevel::specificPostTimeStep()
{
    BL_PROFILE("SupportedWormholeLevel::specificPostTimeStep");

    if (simParams().calculate_constraint_norms && Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        amrex::MultiFab &state_new = get_new_data(state_index);
        FillPatch(*this, state_new, 2, time, state_index, 0, state_new.nComp());

        amrex::MultiFab cst(state_new.boxArray(), state_new.DistributionMap(), 4,
                            0);
        cst.setVal(0.0);
        const auto dx = Geom().CellSizeArray();
        if (simParams().wormhole_matter_model == "no_matter")
        {
            NoMatter no_matter;
            ConstraintsWithMatter<NoMatter> my_constraints(
                no_matter, dx[0], 1.0, 0, Interval(1, 3),
                simParams().wormhole_params.grid_center, time);
            for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto arr       = cst.array(mfi);
                const auto src_arr   = state_new.const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                    { my_constraints(ix, iy, iz, arr, src_arr); });
            }
        }
        else if (simParams().wormhole_matter_model == "effective_teo")
        {
            EffectiveTeoMatter teo_matter;
            ConstraintsWithMatter<EffectiveTeoMatter> my_constraints(
                teo_matter, dx[0], 1.0, 0, Interval(1, 3),
                simParams().wormhole_params.grid_center, time);
            for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto arr       = cst.array(mfi);
                const auto src_arr   = state_new.const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                    { my_constraints(ix, iy, iz, arr, src_arr); });
            }
        }
        else if (simParams().wormhole_matter_model == "complex_scalar")
        {
            ComplexScalarPotential potential(
                simParams().wormhole_params.phantom_mass,
                simParams().wormhole_params.phantom_lambda,
                simParams().wormhole_params.phantom_mu6);
            ComplexExoticScalarField<ComplexScalarPotential> complex_scalar(
                potential, simParams().wormhole_params.support_strength);
            ConstraintsWithMatter<ComplexExoticScalarField<ComplexScalarPotential>>
                my_constraints(complex_scalar, dx[0], 1.0, 0, Interval(1, 3),
                               simParams().wormhole_params.grid_center, time);
            for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto arr       = cst.array(mfi);
                const auto src_arr   = state_new.const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                { my_constraints(ix, iy, iz, arr, src_arr); });
            }
        }
        else
        {
            PhantomDecayPotential potential(
                simParams().wormhole_params.phantom_mass);
            ExoticScalarField<PhantomDecayPotential> exotic_scalar(
                potential, simParams().wormhole_params.support_strength);
            ConstraintsWithMatter<ExoticScalarField<PhantomDecayPotential>>
                my_constraints(exotic_scalar, dx[0], 1.0, 0, Interval(1, 3),
                               simParams().wormhole_params.grid_center, time);
            for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto arr       = cst.array(mfi);
                const auto src_arr   = state_new.const_array(mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                {
                    my_constraints(ix, iy, iz, arr, src_arr);
                });
            }
        }

        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(
            reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst.const_array(mfi);
            reduce_ops.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real ham = arr(i, j, k, 0);
                    const amrex::Real m1  = arr(i, j, k, 1);
                    const amrex::Real m2  = arr(i, j, k, 2);
                    const amrex::Real m3  = arr(i, j, k, 3);
                    const amrex::Real mom2 = (m1 * m1 + m2 * m2 + m3 * m3);
                    return {ham * ham * cell_vol, mom2 * cell_vol, cell_vol};
                });
        }

        auto [sum_ham2, sum_mom2, sum_vol] = reduce_data.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_ham2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mom2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_vol);

        const double L2_Ham =
            (sum_vol > 0.0) ? std::sqrt(sum_ham2 / sum_vol) : 0.0;
        const double L2_Mom =
            (sum_vol > 0.0) ? std::sqrt(sum_mom2 / sum_vol) : 0.0;

        // Publish L2_Ham for the pump governor (Rung 2 active support).
        RLRuntime::publish_cached_L2_Ham(L2_Ham);

        GRParmParse pp;
        std::string output_path = "./";
        pp.load("output_path", output_path, std::string("./"));
        std::string data_subpath;
        pp.load("data_subpath", data_subpath, std::string(""));

        if (!output_path.empty() && output_path.back() != '/')
            output_path += "/";
        if (!data_subpath.empty() && data_subpath.back() != '/')
            data_subpath += "/";

        const std::string out_dir = output_path + data_subpath;
        if (!out_dir.empty())
        {
            amrex::UtilCreateDirectory(out_dir, 0755, false);
        }

        const std::string prefix = out_dir + "constraint_norms";

        SmallDataIO constraints_file(prefix, dt, time, restart_time,
                                     SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            constraints_file.write_header_line({"L2_Ham", "L2_Mom"});
        }
        constraints_file.write_time_data_line({L2_Ham, L2_Mom});
    }

    if (Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        auto &diagnostic_level = parent->getLevel(0);
        amrex::MultiFab &state_diag = diagnostic_level.get_new_data(state_index);
        const auto &diag_geom       = parent->Geom(0);

        FillPatch(diagnostic_level, state_diag, 2, time, state_index, 0,
                  state_diag.nComp());

        {
            const auto &arrs = state_diag.arrays();
            TraceARemoval trace_A_removal;
            PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                                   simParams().min_lapse);
            amrex::ParallelFor(state_diag, amrex::IntVect(0),
                               [=] AMREX_GPU_DEVICE(int box_no, int i, int j,
                                                    int k)
                               {
                                   trace_A_removal(i, j, k, arrs[box_no]);
                                   positive_chi_lapse(i, j, k, arrs[box_no]);
                               });
            amrex::Gpu::streamSynchronize();
        }

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin,
                         amrex::ReduceOpMax, amrex::ReduceOpMax,
                         amrex::ReduceOpMin,
                         amrex::ReduceOpMin, amrex::ReduceOpMax,
                         amrex::ReduceOpMin, amrex::ReduceOpMax>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                          amrex::Real,
                          amrex::Real, amrex::Real,
                          amrex::Real, amrex::Real>
            reduce_data(reduce_ops);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        
        const auto prob_lo = diag_geom.ProbLoArray();
        const auto dx_arr = diag_geom.CellSizeArray();

        // The apparent-horizon / expansion proxy must be measured about the
        // physics center (where the initial data is centered), not the
        // coordinate origin at the domain corner.  Using the corner makes
        // r ~ |grid_center| at the actual center, which collapses the
        // 2*sqrt(chi)/r regularizing term and produces spurious theta_plus<0
        // (false trapped surfaces) offset to r ~ |grid_center|.
        const amrex::Real cx = simParams().wormhole_params.grid_center[0];
        const amrex::Real cy = simParams().wormhole_params.grid_center[1];
        const amrex::Real cz = simParams().wormhole_params.grid_center[2];

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const amrex::Real chi   = arr(i, j, k, c_chi);
                    const amrex::Real K     = arr(i, j, k, c_K);
                    
                    const amrex::Real x = prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] - cx;
                    const amrex::Real y = prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] - cy;
                    const amrex::Real z = prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] - cz;
                    const amrex::Real r2 = x*x + y*y + z*z;
                    const amrex::Real r = std::sqrt(r2);
                    
                    amrex::Real ah_radius = 0.0;
                    amrex::Real theta_plus_min_proxy = 1.0e30;
                    if (r > 1e-6)
                    {
                        const amrex::Real A11 = arr(i, j, k, c_A11);
                        const amrex::Real A22 = arr(i, j, k, c_A22);
                        const amrex::Real A33 = arr(i, j, k, c_A33);
                        const amrex::Real A12 = arr(i, j, k, c_A12);
                        const amrex::Real A13 = arr(i, j, k, c_A13);
                        const amrex::Real A23 = arr(i, j, k, c_A23);
                        
                        const amrex::Real Arr = (A11*x*x + A22*y*y + A33*z*z + 2.0*A12*x*y + 2.0*A13*x*z + 2.0*A23*y*z) / r2;
                        
                        const amrex::Real dx_chi =
                            (arr(i + 1, j, k, c_chi) - arr(i - 1, j, k, c_chi)) /
                            (2.0 * dx_arr[0]);
                        const amrex::Real dy_chi =
                            (arr(i, j + 1, k, c_chi) - arr(i, j - 1, k, c_chi)) /
                            (2.0 * dx_arr[1]);
                        const amrex::Real dz_chi =
                            (arr(i, j, k + 1, c_chi) - arr(i, j, k - 1, c_chi)) /
                            (2.0 * dx_arr[2]);
                        const amrex::Real dchi_dr =
                            (x * dx_chi + y * dy_chi + z * dz_chi) / r;
                        const amrex::Real sqrt_chi =
                            std::sqrt(amrex::max(chi, amrex::Real(1.0e-20)));

                        const amrex::Real theta_plus =
                            2.0 * sqrt_chi / r - dchi_dr / sqrt_chi + Arr -
                            (2.0 / 3.0) * K;
                        theta_plus_min_proxy = theta_plus;
                        
                        if (theta_plus <= 0.0)
                        {
                            ah_radius = r;
                        }
                    }

                    const amrex::Real sf_phi = arr(i, j, k, c_phi);
                    const amrex::Real sf_Pi  = arr(i, j, k, c_Pi);

                    return {lapse, chi, amrex::Math::abs(K), ah_radius, theta_plus_min_proxy,
                            sf_phi, sf_phi, sf_Pi, sf_Pi};
                });
        }

        const auto reduce_vals = reduce_data.value();
        amrex::Real min_lapse  = amrex::get<0>(reduce_vals);
        amrex::Real min_chi    = amrex::get<1>(reduce_vals);
        amrex::Real max_abs_K  = amrex::get<2>(reduce_vals);
        amrex::Real max_ah_r   = amrex::get<3>(reduce_vals);
        amrex::Real min_theta_plus = amrex::get<4>(reduce_vals);
        amrex::Real min_phi    = amrex::get<5>(reduce_vals);
        amrex::Real max_phi    = amrex::get<6>(reduce_vals);
        amrex::Real min_Pi     = amrex::get<7>(reduce_vals);
        amrex::Real max_Pi     = amrex::get<8>(reduce_vals);
        amrex::ParallelDescriptor::ReduceRealMin(min_lapse);
        amrex::ParallelDescriptor::ReduceRealMin(min_chi);
        amrex::ParallelDescriptor::ReduceRealMax(max_abs_K);
        amrex::ParallelDescriptor::ReduceRealMax(max_ah_r);
        amrex::ParallelDescriptor::ReduceRealMin(min_theta_plus);
        amrex::ParallelDescriptor::ReduceRealMin(min_phi);
        amrex::ParallelDescriptor::ReduceRealMax(max_phi);
        amrex::ParallelDescriptor::ReduceRealMin(min_Pi);
        amrex::ParallelDescriptor::ReduceRealMax(max_Pi);

        const amrex::Real tol =
            amrex::max(amrex::Real(1.0e-14), amrex::Real(1.0e-12) * amrex::Math::abs(min_lapse));

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_loc;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data_loc(reduce_ops_loc);
        using ReduceTupleLoc = typename decltype(reduce_data_loc)::Type;

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops_loc.eval(
                bx, reduce_data_loc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleLoc
                {
                    const amrex::Real lapse = arr(i, j, k, c_lapse);
                    const bool is_min       = (amrex::Math::abs(lapse - min_lapse) <= tol);
                    if (!is_min)
                    {
                        return {0.0, 0.0, 0.0, 0.0};
                    }
                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    return {x, y, z, 1.0};
                });
        }

        auto [sum_x, sum_y, sum_z, count] = reduce_data_loc.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_x);
        amrex::ParallelDescriptor::ReduceRealSum(sum_y);
        amrex::ParallelDescriptor::ReduceRealSum(sum_z);
        amrex::ParallelDescriptor::ReduceRealSum(count);

        const amrex::Real min_lapse_x = (count > 0.0) ? (sum_x / count) : 0.0;
        const amrex::Real min_lapse_y = (count > 0.0) ? (sum_y / count) : 0.0;
        const amrex::Real min_lapse_z = (count > 0.0) ? (sum_z / count) : 0.0;

        const amrex::Real tol_theta = amrex::max(
            amrex::Real(1.0e-12),
            amrex::Real(1.0e-8) * amrex::Math::abs(min_theta_plus));

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_theta_loc;
        amrex::ReduceData<amrex::Real, amrex::Real> reduce_data_theta_loc(
            reduce_ops_theta_loc);
        using ReduceTupleThetaLoc = typename decltype(reduce_data_theta_loc)::Type;

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops_theta_loc.eval(
                bx, reduce_data_theta_loc,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleThetaLoc
                {
                    const amrex::Real chi = arr(i, j, k, c_chi);
                    const amrex::Real K   = arr(i, j, k, c_K);

                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] - cx;
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] - cy;
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] - cz;
                    const amrex::Real r2 = x * x + y * y + z * z;
                    const amrex::Real r  = std::sqrt(r2);

                    if (r <= 1e-6)
                    {
                        return {0.0, 0.0};
                    }

                    const amrex::Real A11 = arr(i, j, k, c_A11);
                    const amrex::Real A22 = arr(i, j, k, c_A22);
                    const amrex::Real A33 = arr(i, j, k, c_A33);
                    const amrex::Real A12 = arr(i, j, k, c_A12);
                    const amrex::Real A13 = arr(i, j, k, c_A13);
                    const amrex::Real A23 = arr(i, j, k, c_A23);

                    const amrex::Real Arr =
                        (A11 * x * x + A22 * y * y + A33 * z * z +
                         2.0 * A12 * x * y + 2.0 * A13 * x * z +
                         2.0 * A23 * y * z) /
                        r2;

                    const amrex::Real dx_chi =
                        (arr(i + 1, j, k, c_chi) - arr(i - 1, j, k, c_chi)) /
                        (2.0 * dx_arr[0]);
                    const amrex::Real dy_chi =
                        (arr(i, j + 1, k, c_chi) - arr(i, j - 1, k, c_chi)) /
                        (2.0 * dx_arr[1]);
                    const amrex::Real dz_chi =
                        (arr(i, j, k + 1, c_chi) - arr(i, j, k - 1, c_chi)) /
                        (2.0 * dx_arr[2]);
                    const amrex::Real dchi_dr =
                        (x * dx_chi + y * dy_chi + z * dz_chi) / r;
                    const amrex::Real sqrt_chi =
                        std::sqrt(amrex::max(chi, amrex::Real(1.0e-20)));

                    const amrex::Real theta_plus =
                        2.0 * sqrt_chi / r - dchi_dr / sqrt_chi + Arr -
                        (2.0 / 3.0) * K;

                    const bool is_min =
                        (amrex::Math::abs(theta_plus - min_theta_plus) <=
                         tol_theta);
                    if (!is_min)
                    {
                        return {0.0, 0.0};
                    }
                    return {r, 1.0};
                });
        }

        auto [sum_r_theta, count_r_theta] = reduce_data_theta_loc.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_r_theta);
        amrex::ParallelDescriptor::ReduceRealSum(count_r_theta);
        const amrex::Real r_at_min_theta_plus =
            (count_r_theta > 0.0) ? (sum_r_theta / count_r_theta) : 0.0;

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            reduce_ops_bary;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            reduce_data_bary(reduce_ops_bary);
        using ReduceTupleBary = typename decltype(reduce_data_bary)::Type;

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops_bary.eval(
                bx, reduce_data_bary,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleBary
                {
                    const amrex::Real x =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                    const amrex::Real y =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                    const amrex::Real z =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                    const amrex::Real phi = arr(i, j, k, c_phi);
                    const amrex::Real pi  = arr(i, j, k, c_Pi);
                    const amrex::Real rho_scalar =
                        pi * pi + amrex::Real(0.5) * phi * phi;
                    const amrex::Real rho = amrex::max(
                        rho_scalar,
                        amrex::max(arr(i, j, k, c_teo_rho),
                                   arr(i, j, k, c_dust_rho)));
                    return {rho * x, rho * y, rho * z, rho};
                });
        }

        auto [sum_rx, sum_ry, sum_rz, sum_rho] = reduce_data_bary.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_rx);
        amrex::ParallelDescriptor::ReduceRealSum(sum_ry);
        amrex::ParallelDescriptor::ReduceRealSum(sum_rz);
        amrex::ParallelDescriptor::ReduceRealSum(sum_rho);
        const amrex::Real bary_x =
            (sum_rho > 0.0) ? (sum_rx / sum_rho) : 0.0;
        const amrex::Real bary_y =
            (sum_rho > 0.0) ? (sum_ry / sum_rho) : 0.0;
        const amrex::Real bary_z =
            (sum_rho > 0.0) ? (sum_rz / sum_rho) : 0.0;

        // Coordinate angular momentum about the z-axis through the throat,
        // J_z = sum (x p_y - y p_x) dV, with scalar momentum density
        // p_i = Pi1 d_i phi1 + Pi2 d_i phi2. Tracks the spin carried by the
        // complex phantom field (zero for the non-rotating real-scalar runs).
        const amrex::Real center_x = simParams().wormhole_params.grid_center[0];
        const amrex::Real center_y = simParams().wormhole_params.grid_center[1];
        const amrex::Real cell_vol_fine = dx_arr[0] * dx_arr[1] * dx_arr[2];

        amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops_jz;
        amrex::ReduceData<amrex::Real> reduce_data_jz(reduce_ops_jz);
        using ReduceTupleJz = typename decltype(reduce_data_jz)::Type;

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops_jz.eval(
                bx, reduce_data_jz,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleJz
                {
                    const amrex::Real Pi1 = arr(i, j, k, c_Pi);
                    const amrex::Real Pi2 = arr(i, j, k, c_Pi2);

                    const amrex::Real dx_phi1 =
                        (arr(i + 1, j, k, c_phi) - arr(i - 1, j, k, c_phi)) /
                        (2.0 * dx_arr[0]);
                    const amrex::Real dy_phi1 =
                        (arr(i, j + 1, k, c_phi) - arr(i, j - 1, k, c_phi)) /
                        (2.0 * dx_arr[1]);
                    const amrex::Real dx_phi2 =
                        (arr(i + 1, j, k, c_phi2) - arr(i - 1, j, k, c_phi2)) /
                        (2.0 * dx_arr[0]);
                    const amrex::Real dy_phi2 =
                        (arr(i, j + 1, k, c_phi2) - arr(i, j - 1, k, c_phi2)) /
                        (2.0 * dx_arr[1]);

                    const amrex::Real p_x = Pi1 * dx_phi1 + Pi2 * dx_phi2;
                    const amrex::Real p_y = Pi1 * dy_phi1 + Pi2 * dy_phi2;

                    const amrex::Real xr =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] -
                        center_x;
                    const amrex::Real yr =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] -
                        center_y;

                    return {(xr * p_y - yr * p_x) * cell_vol_fine};
                });
        }
        amrex::Real J_z = amrex::get<0>(reduce_data_jz.value());
        amrex::ParallelDescriptor::ReduceRealSum(J_z);

        // Conserved U(1) Noether charge of the complex phantom field,
        //   Q = integral (phi1 Pi2 - phi2 Pi1) sqrt(gamma) dV,  sqrt(gamma) =
        //   chi^{-3/2}.  Q is exactly conserved for a genuine bound (Q-ball)
        //   eigenstate and is the principled confinement diagnostic: unlike
        //   rho_sum (an AMR-refined-region sum that decays merely because the
        //   grid de-refines), Q_sphere is integrated over a FIXED coordinate
        //   sphere r < diag_sphere_radius about the throat, so a drop means the
        //   charge physically leaves the throat region rather than the grid
        //   coarsening.  Q_total is the whole-domain charge (radiation-robust).
        amrex::Real diag_sphere_radius = 10.0;
        {
            GRParmParse pp_r;
            pp_r.load("diag_sphere_radius", diag_sphere_radius, 10.0);
        }
        const amrex::Real center_z = simParams().wormhole_params.grid_center[2];
        const amrex::Real r_sphere2 = diag_sphere_radius * diag_sphere_radius;

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum>
            reduce_ops_q;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real>
            reduce_data_q(reduce_ops_q);
        using ReduceTupleQ = typename decltype(reduce_data_q)::Type;

        for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_diag.const_array(mfi);
            reduce_ops_q.eval(
                bx, reduce_data_q,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTupleQ
                {
                    const amrex::Real phi1 = arr(i, j, k, c_phi);
                    const amrex::Real phi2 = arr(i, j, k, c_phi2);
                    const amrex::Real Pi1  = arr(i, j, k, c_Pi);
                    const amrex::Real Pi2  = arr(i, j, k, c_Pi2);
                    const amrex::Real chi  = amrex::max(
                        arr(i, j, k, c_chi), amrex::Real(1.0e-10));
                    const amrex::Real sqrt_gamma =
                        amrex::Real(1.0) / (chi * std::sqrt(chi));
                    const amrex::Real q_density =
                        (phi1 * Pi2 - phi2 * Pi1) * sqrt_gamma * cell_vol_fine;

                    const amrex::Real xr =
                        prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] - center_x;
                    const amrex::Real yr =
                        prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] - center_y;
                    const amrex::Real zr =
                        prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] - center_z;
                    const bool inside =
                        (xr * xr + yr * yr + zr * zr) <= r_sphere2;

                    const amrex::Real rho_density =
                        (Pi1 * Pi1 + Pi2 * Pi2 +
                         amrex::Real(0.5) * (phi1 * phi1 + phi2 * phi2)) *
                        cell_vol_fine;

                    return {q_density,
                            inside ? q_density : amrex::Real(0.0),
                            inside ? rho_density : amrex::Real(0.0)};
                });
        }
        auto [Q_total, Q_sphere, rho_sphere] = reduce_data_q.value();
        amrex::ParallelDescriptor::ReduceRealSum(Q_total);
        amrex::ParallelDescriptor::ReduceRealSum(Q_sphere);
        amrex::ParallelDescriptor::ReduceRealSum(rho_sphere);

        // Pump control-effort diagnostic (Rung 2 active support): instantaneous
        // power the matter pump injects into the scalar field,
        //   P_pump = integral alpha * (Pi1 * S_Pi1 + Pi2 * S_Pi2) sqrt(gamma) dV,
        // where S_Pi{1,2} is the pump source added to d_t Pi{1,2} (the same term
        // ComplexExoticScalarField::add_matter_rhs applies).  This is the
        // Bianchi-violating energy injection; MAP-Elites penalises it so the
        // search favours minimal-intervention support.  Zero when the pump is
        // off (num_sites < 1), so passive runs report pump_work = 0.
        amrex::Real pump_work = 0.0;
        {
            const RLMatterPumpParams pump = build_wormhole_pump(simParams());
            if (pump.num_sites >= 1)
            {
                amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops_pw;
                amrex::ReduceData<amrex::Real> reduce_data_pw(reduce_ops_pw);
                using ReduceTuplePW = typename decltype(reduce_data_pw)::Type;
                const amrex::Real pw_time = time;
                for (amrex::MFIter mfi(state_diag, amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();
                    const auto arr       = state_diag.const_array(mfi);
                    reduce_ops_pw.eval(
                        bx, reduce_data_pw,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuplePW
                        {
                            const amrex::Real x =
                                prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0];
                            const amrex::Real y =
                                prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1];
                            const amrex::Real z =
                                prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2];
                            const amrex::Real lapse = arr(i, j, k, c_lapse);
                            const amrex::Real chi = amrex::max(
                                arr(i, j, k, c_chi), amrex::Real(1.0e-10));
                            const amrex::Real sqrt_gamma =
                                amrex::Real(1.0) / (chi * std::sqrt(chi));
                            const amrex::Real phi1 = arr(i, j, k, c_phi);
                            const amrex::Real phi2 = arr(i, j, k, c_phi2);
                            const amrex::Real Pi1  = arr(i, j, k, c_Pi);
                            const amrex::Real Pi2  = arr(i, j, k, c_Pi2);

                            // Reproduce the pump source (mirror
                            // ComplexExoticScalarField::add_matter_rhs).
                            amrex::Real s1 = 0.0, s2 = 0.0;
                            const amrex::Real gov = pump.governor;
                            if (pump.k_p > 0.0)
                            {
                                const amrex::Real kp = pump.k_p;
                                const amrex::Real kd = pump.k_d;
                                const amrex::Real inv_a = 1.0 / lapse;
                                const amrex::Real tw =
                                    (pump.target_width > 0.0) ? pump.target_width
                                                              : pump.width;
                                for (int s = 0; s < pump.num_sites; ++s)
                                {
                                    const auto &site = pump.sites[s];
                                    if (site.amplitude <= 0.0)
                                        continue;
                                    const amrex::Real env =
                                        RLRuntime::compute_site_envelope(
                                            x, y, z, site, tw,
                                            pump.target_profile);
                                    if (env < 1.0e-8)
                                        continue;
                                    const amrex::Real amp_t =
                                        (pump.target_amp > 0.0) ? pump.target_amp
                                                                : site.amplitude;
                                    const amrex::Real g   = amp_t * env;
                                    const amrex::Real arg =
                                        -site.frequency * pw_time + site.phase;
                                    const amrex::Real tphi1 = g * std::cos(arg);
                                    const amrex::Real tphi2 = g * std::sin(arg);
                                    const amrex::Real tPi1 =
                                        site.frequency * tphi2 * inv_a;
                                    const amrex::Real tPi2 =
                                        -site.frequency * tphi1 * inv_a;
                                    const amrex::Real w = gov * env;
                                    s1 += w * (-kp * (phi1 - tphi1) -
                                               kd * (Pi1 - tPi1));
                                    s2 += w * (-kp * (phi2 - tphi2) -
                                               kd * (Pi2 - tPi2));
                                }
                            }
                            else
                            {
                                for (int s = 0; s < pump.num_sites; ++s)
                                {
                                    const amrex::Real base =
                                        RLRuntime::compute_site_base(
                                            x, y, z, pump.sites[s], pump.width,
                                            gov);
                                    if (base <= 0.0)
                                        continue;
                                    const amrex::Real arg =
                                        pump.sites[s].frequency * pw_time +
                                        pump.sites[s].phase;
                                    s1 += base * std::cos(arg);
                                    s2 += base * std::sin(arg);
                                }
                            }
                            const amrex::Real power =
                                lapse * (Pi1 * s1 + Pi2 * s2) * sqrt_gamma *
                                cell_vol_fine;
                            return {power};
                        });
                }
                pump_work = amrex::get<0>(reduce_data_pw.value());
                amrex::ParallelDescriptor::ReduceRealSum(pump_work);
            }
        }

        GRParmParse pp;
        std::string output_path = "./";
        pp.load("output_path", output_path, std::string("./"));
        std::string data_subpath;
        pp.load("data_subpath", data_subpath, std::string(""));

        if (!output_path.empty() && output_path.back() != '/')
            output_path += "/";
        if (!data_subpath.empty() && data_subpath.back() != '/')
            data_subpath += "/";

        const std::string out_dir = output_path + data_subpath;
        if (!out_dir.empty())
        {
            amrex::UtilCreateDirectory(out_dir, 0755, false);
        }

        const std::string prefix = out_dir + "collapse_diagnostics";
        SmallDataIO diag_file(prefix, dt, time, restart_time, SmallDataIO::APPEND,
                              first_step);
        diag_file.remove_duplicate_time_data();
        if (first_step)
        {
            diag_file.write_header_line(
                {"min_lapse", "min_chi", "max_abs_K", "min_lapse_x",
                 "min_lapse_y", "min_lapse_z", "max_ah_r", "min_theta_plus",
                 "r_at_min_theta_plus",
                 "min_phi", "max_phi", "min_Pi", "max_Pi",
                 "barycenter_x", "barycenter_y", "barycenter_z", "rho_sum",
                 "J_z", "Q_total", "Q_sphere", "rho_sphere", "pump_work"});
        }
        diag_file.write_time_data_line({static_cast<double>(min_lapse),
                                        static_cast<double>(min_chi),
                                        static_cast<double>(max_abs_K),
                                        static_cast<double>(min_lapse_x),
                                        static_cast<double>(min_lapse_y),
                                        static_cast<double>(min_lapse_z),
                                        static_cast<double>(max_ah_r),
                                        static_cast<double>(min_theta_plus),
                                        static_cast<double>(r_at_min_theta_plus),
                                        static_cast<double>(min_phi),
                                        static_cast<double>(max_phi),
                                        static_cast<double>(min_Pi),
                                        static_cast<double>(max_Pi),
                                        static_cast<double>(bary_x),
                                        static_cast<double>(bary_y),
                                        static_cast<double>(bary_z),
                                        static_cast<double>(sum_rho),
                                        static_cast<double>(J_z),
                                        static_cast<double>(Q_total),
                                        static_cast<double>(Q_sphere),
                                        static_cast<double>(rho_sphere),
                                        static_cast<double>(pump_work)});
    }

    // ---- In-code Weyl4 / Psi4 spherical-harmonic extraction --------------
    // Proper GRTeclyn-collaboration SphericalExtraction: interpolates the
    // "Weyl4" derived variable onto spheres at the configured extraction_radii
    // and writes the requested (l,m) mode time series to the extraction_subpath
    // ("data/") every coarse step -- densely sampled and completely decoupled
    // from the plotfile / frame cadence, for clean GW waveform analysis.
    if (simParams().activate_extraction)
    {
        const int min_level =
            simParams().extraction_params.min_extraction_level();
        const bool calculate_weyl = at_level_timestep_multiple(min_level);

        if (calculate_weyl && Level() == min_level)
        {
            const amrex::Real m_time  = get_state_data(state_index).curTime();
            const amrex::Real m_dt    = get_gramr_ptr()->dtLevel(Level());
            const amrex::Real restart_time =
                get_gramr_ptr()->get_restart_time();
            const bool first_step = (m_time <= m_dt);

            WeylExtraction my_extraction(simParams().extraction_params, m_dt,
                                         m_time, first_step, restart_time);
            my_extraction.execute_query(&get_bhamr_ptr()->m_weyl_interpolator);
        }
    }
}