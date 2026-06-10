#include "ExternalGridInitialData.hpp"
#include "RadialRecipeInitialData.hpp"
#include "RadialRecipeLevel.hpp"
#include "RadialRecipeMatterDispatch.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ChiTagger.hpp"
#include "ConstraintsWithMatter.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "ScalarField.hpp"
#include "ExoticScalarField.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "Weyl4WithMatter.hpp"

#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <array>
#include <cmath>

namespace
{
// Observer-sampled (Warp Factory style) pointwise energy-condition margins.
//
// Given the 3+1 matter decomposition (rho, j_i, S_ij) measured by the Eulerian
// observer and the physical spatial metric gamma_ij, we build an orthonormal
// triad, transform (j, S) into that frame, and minimise the NEC/WEC/SEC/DEC
// contractions over a fixed sampling of spatial directions and observer speeds.
// Each margin is >= 0 where the corresponding condition is satisfied and < 0
// where it is violated; the rigorous Hawking-Ellis eigenvalue refinement runs
// post-hoc in the Python warpfactory module.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
ec_point_margins(amrex::Real rho, const amrex::Real jin[3],
                 const amrex::Real Sin[3][3], const amrex::Real g[3][3],
                 amrex::Real &nec, amrex::Real &wec, amrex::Real &sec,
                 amrex::Real &dec)
{
    constexpr amrex::Real tiny = 1.0e-30;

    // Gram-Schmidt orthonormalisation of the coordinate basis w.r.t. gamma_ij,
    // producing triad vectors e[a][i] (upper coordinate index i) satisfying
    // gamma_ij e[a]^i e[b]^j = delta_ab.
    auto gdot = [&](const amrex::Real a[3], const amrex::Real b[3])
    {
        amrex::Real s = 0.0;
        for (int p = 0; p < 3; ++p)
            for (int q = 0; q < 3; ++q)
                s += a[p] * g[p][q] * b[q];
        return s;
    };

    amrex::Real e[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    for (int a = 0; a < 3; ++a)
    {
        for (int b = 0; b < a; ++b)
        {
            const amrex::Real proj = gdot(e[a], e[b]);
            for (int p = 0; p < 3; ++p)
                e[a][p] -= proj * e[b][p];
        }
        const amrex::Real nrm = std::sqrt(amrex::max(gdot(e[a], e[a]), tiny));
        for (int p = 0; p < 3; ++p)
            e[a][p] /= nrm;
    }

    // Transform momentum density and spatial stress into the orthonormal frame.
    amrex::Real jhat[3] = {0.0, 0.0, 0.0};
    amrex::Real Shat[3][3] = {{0.0}};
    for (int a = 0; a < 3; ++a)
    {
        for (int p = 0; p < 3; ++p)
            jhat[a] += e[a][p] * jin[p];
        for (int b = 0; b < 3; ++b)
        {
            amrex::Real s = 0.0;
            for (int p = 0; p < 3; ++p)
                for (int q = 0; q < 3; ++q)
                    s += e[a][p] * e[b][q] * Sin[p][q];
            Shat[a][b] = s;
        }
    }
    const amrex::Real trS = Shat[0][0] + Shat[1][1] + Shat[2][2];

    constexpr int NTHETA = 4;
    constexpr int NPHI   = 6;
    constexpr int NSPEED = 3;
    constexpr amrex::Real vmax = 0.92;

    nec = wec = sec = dec = 1.0e30;

    for (int it = 0; it < NTHETA; ++it)
    {
        const amrex::Real theta = M_PI * (amrex::Real(it) + 0.5) / NTHETA;
        const amrex::Real st = std::sin(theta);
        const amrex::Real ct = std::cos(theta);
        for (int ip = 0; ip < NPHI; ++ip)
        {
            const amrex::Real phi = 2.0 * M_PI * amrex::Real(ip) / NPHI;
            const amrex::Real n[3] = {st * std::cos(phi), st * std::sin(phi),
                                      ct};

            amrex::Real Snn = 0.0;
            amrex::Real jn  = 0.0;
            amrex::Real Sn[3] = {0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a)
            {
                jn += jhat[a] * n[a];
                for (int b = 0; b < 3; ++b)
                {
                    Snn += Shat[a][b] * n[a] * n[b];
                    Sn[a] += Shat[a][b] * n[b];
                }
            }
            const amrex::Real abs_jn = std::abs(jn);

            // NEC: null k = (1, n); worst over +-n.
            nec = amrex::min(nec, rho + Snn - 2.0 * abs_jn);

            for (int iv = 0; iv < NSPEED; ++iv)
            {
                const amrex::Real v =
                    vmax * (amrex::Real(iv) + 1.0) / NSPEED;
                const amrex::Real g2 = 1.0 / (1.0 - v * v);
                const amrex::Real gam = std::sqrt(g2);

                // Timelike observer u = gamma (n_t + v n); worst over +-n.
                const amrex::Real Tuu =
                    g2 * (rho + v * v * Snn - 2.0 * v * abs_jn);
                wec = amrex::min(wec, Tuu);
                sec = amrex::min(sec, Tuu + 0.5 * (-rho + trS));

                // DEC: energy flux q^a = -T^a_b u^b must be causal and
                // future-directed; report q^0 - |q^spatial| for both +-n.
                for (int s = -1; s <= 1; s += 2)
                {
                    const amrex::Real q0 = gam * (rho + v * (s * jn));
                    amrex::Real fmag = 0.0;
                    for (int a = 0; a < 3; ++a)
                    {
                        const amrex::Real comp = jhat[a] + v * (s * Sn[a]);
                        fmag += comp * comp;
                    }
                    dec = amrex::min(dec, q0 - gam * std::sqrt(fmag));
                }
            }
        }
    }
}

// Fill the Hamiltonian/momentum constraint MultiFab using a given matter model.
// Templated so the same reduction downstream can be fed either the canonical
// ScalarField or the phantom ExoticScalarField, ensuring the constraint that is
// measured uses the same matter that actually sources the evolution.
template <class matter_t>
void fill_matter_constraints(amrex::MultiFab &cst,
                             const amrex::MultiFab &state_new,
                             const matter_t &matter, amrex::Real dx0,
                             const std::array<double, AMREX_SPACEDIM> &center,
                             amrex::Real time)
{
    ConstraintsWithMatter<matter_t> my_constraints(
        matter, dx0, 1.0, 0, Interval(1, 3), center, time);
    for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto arr       = cst.array(mfi);
        const auto src_arr   = state_new.const_array(mfi);
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
            { my_constraints(ix, iy, iz, arr, src_arr); });
    }
}

// Reduce the pointwise observer-sampled energy-condition margins (NEC/WEC/SEC/
// DEC) and the integrated NEC violation over the finest level, using the
// supplied matter model for the 3+1 stress-energy decomposition.  Returns
// {min_NEC, min_WEC, min_SEC, min_DEC, integral_NEC_violation}.
template <class matter_t>
std::array<amrex::Real, 5>
reduce_ec_margins(const amrex::MultiFab &state_fine, const matter_t &matter,
                  amrex::Real ec_dx, amrex::Real ec_cell_vol)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMin,
                     amrex::ReduceOpMin, amrex::ReduceOpMin, amrex::ReduceOpSum>
        ec_ops;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                      amrex::Real>
        ec_red(ec_ops);
    using ECTuple = typename decltype(ec_red)::Type;

    for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi)
    {
        const amrex::Box &bx = mfi.validbox();
        const auto st        = state_fine.const_array(mfi);
        ec_ops.eval(
            bx, ec_red,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ECTuple
            {
                const FourthOrderDerivatives deriv(ec_dx);
                const amrex::CellData<const amrex::Real> &cell =
                    st.cellData(i, j, k);
                typename matter_t::Vars vars(cell);
                const typename matter_t::D1Vars d1(i, j, k, st, deriv);
                const auto h_UU = CCZ4Geometry::compute_inverse_metric(vars);
                const auto chris =
                    CCZ4Geometry::compute_christoffel(d1, h_UU);
                const auto emt =
                    matter.compute_emtensor(vars, d1, h_UU, chris.ULL);

                const amrex::Real rho    = emt.rho;
                const amrex::Real jin[3] = {emt.j[0], emt.j[1], emt.j[2]};
                amrex::Real Sin[3][3];
                Sin[0][0] = emt.S[0][0];
                Sin[0][1] = Sin[1][0] = emt.S[0][1];
                Sin[0][2] = Sin[2][0] = emt.S[0][2];
                Sin[1][1] = emt.S[1][1];
                Sin[1][2] = Sin[2][1] = emt.S[1][2];
                Sin[2][2] = emt.S[2][2];

                const amrex::Real chi = st(i, j, k, c_chi);
                const amrex::Real inv_chi =
                    1.0 / amrex::max(chi, amrex::Real(1.0e-12));
                amrex::Real g[3][3];
                g[0][0] = st(i, j, k, c_h11) * inv_chi;
                g[0][1] = g[1][0] = st(i, j, k, c_h12) * inv_chi;
                g[0][2] = g[2][0] = st(i, j, k, c_h13) * inv_chi;
                g[1][1] = st(i, j, k, c_h22) * inv_chi;
                g[1][2] = g[2][1] = st(i, j, k, c_h23) * inv_chi;
                g[2][2] = st(i, j, k, c_h33) * inv_chi;

                amrex::Real nec, wec, sec, dec;
                ec_point_margins(rho, jin, Sin, g, nec, wec, sec, dec);

                const amrex::Real neg_nec =
                    (nec < 0.0) ? (-nec * ec_cell_vol) : 0.0;
                return {nec, wec, sec, dec, neg_nec};
            });
    }

    auto ec_vals             = ec_red.value();
    amrex::Real min_nec      = amrex::get<0>(ec_vals);
    amrex::Real min_wec      = amrex::get<1>(ec_vals);
    amrex::Real min_sec      = amrex::get<2>(ec_vals);
    amrex::Real min_dec      = amrex::get<3>(ec_vals);
    amrex::Real integral_nec = amrex::get<4>(ec_vals);
    amrex::ParallelDescriptor::ReduceRealMin(min_nec);
    amrex::ParallelDescriptor::ReduceRealMin(min_wec);
    amrex::ParallelDescriptor::ReduceRealMin(min_sec);
    amrex::ParallelDescriptor::ReduceRealMin(min_dec);
    amrex::ParallelDescriptor::ReduceRealSum(integral_nec);
    return {min_nec, min_wec, min_sec, min_dec, integral_nec};
}

// Derived plot field: geometric required energy density
//   rho_req = Ham_vac / (16 pi) = (R + (2/3) K^2 - A_ij A^ij) / (16 pi)
// This is the Eulerian energy density that Einstein's equations demand to
// source the *evolved geometry*, independent of whatever matter is actually
// present.  Cells with rho_req < 0 are precisely the exotic-matter requirement
// whose minimum/integral are logged as min_rho_req / integral_neg_rho in
// constraint_norms.dat, so this field lets the same quantity be mapped and
// plotted (it matches the scalar diagnostics by construction).
void compute_rho_req_mf(amrex::MultiFab &out_mf, int dcomp, int /*ncomp*/,
                        const amrex::MultiFab &src_mf,
                        const amrex::Geometry &geomdata,
                        amrex::Real /*time*/, const int * /*bcrec*/,
                        int /*level*/)
{
    const auto &out_arrays         = out_mf.arrays();
    const auto &src_arrays         = src_mf.const_arrays();
    const amrex::Real dx           = geomdata.CellSize(0);
    constexpr amrex::Real inv_16pi = 1.0 / (16.0 * M_PI);

    // Vacuum (geometry-only) Hamiltonian written into component dcomp; the empty
    // momentum interval means store_vars only fills Ham.
    Constraints vacuum_constraints(dx, dcomp, Interval());

    amrex::ParallelFor(
        out_mf, out_mf.nGrowVect(),
        [=] AMREX_GPU_DEVICE(int box_no, int ix, int iy, int iz) noexcept
        {
            vacuum_constraints(ix, iy, iz, out_arrays[box_no],
                               src_arrays[box_no]);
            out_arrays[box_no](ix, iy, iz, dcomp) *= inv_16pi;
        });
}
} // namespace

void RadialRecipeLevel::variableSetUp()
{
    BL_PROFILE("RadialRecipeLevel::variableSetUp()");
    stateVariableSetUp();

    RadialRecipeMatter::setup_derived_quantities(state_index, simParams());

    // Geometry-only required energy density rho_req = Ham_vac / (16 pi), so the
    // exotic-matter requirement (rho_req < 0) can be plotted as a field, not
    // just integrated into constraint_norms.dat.
    {
        int num_ghosts       = 2;
        auto &derive_lst     = amrex::AmrLevel::get_derive_lst();
        const auto &desc_lst = amrex::AmrLevel::get_desc_lst();
        derive_lst.add(
            "rho_req", amrex::IndexType::TheCellType(), 1,
            amrex::Vector<std::string>{"rho_req"}, compute_rho_req_mf,
            [=](const amrex::Box &box) { return amrex::grow(box, num_ghosts); },
            &amrex::cell_quartic_interp);
        derive_lst.addComponent("rho_req", desc_lst, state_index, 0, NUM_VARS);
    }
}

void RadialRecipeLevel::specificAdvance()
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

void RadialRecipeLevel::initData()
{
    BL_PROFILE("RadialRecipeLevel::initData");

    amrex::MultiFab &state = get_new_data(state_index);
    const auto &arrs       = state.arrays();

    if (!simParams().recipe_initial_data_file.empty())
    {
        ExternalGridInitialData ext_data(simParams().external_grid_params,
                                         Geom().CellSize(0));

        amrex::ParallelFor(
            state, state.nGrowVect(),
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
            {
                for (int n = 0; n < NUM_VARS; ++n)
                {
                    arrs[box_no](i, j, k, n) = 0.;
                }
                ext_data.compute(i, j, k, arrs[box_no]);
            });
    }
    else
    {
        RadialRecipeInitialData recipe(simParams().recipe_params,
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
                recipe.compute(i, j, k, arrs[box_no]);
            });
    }

    amrex::Gpu::streamSynchronize();
}

void RadialRecipeLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                        amrex::MultiFab &a_rhs,
                                        const double a_time)
{
    BL_PROFILE("RadialRecipeLevel::specificEvalRHS()");
    const int soln_ghosts = a_soln.nGrowVect()[0];
    if (soln_ghosts > 0)
    {
        FillPatch(*this, a_soln, soln_ghosts, a_time, state_index, 0,
                  a_soln.nComp());
    }
    const auto &soln_arrs = a_soln.arrays();
    TraceARemoval trace_A_removal;
    PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                           simParams().min_lapse);

    amrex::ParallelFor(a_soln, a_soln.nGrowVect(),
                       [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k)
                       {
                           trace_A_removal(i, j, k, soln_arrs[box_no]);
                           positive_chi_lapse(i, j, k, soln_arrs[box_no]);
                       });

    RadialRecipeMatter::eval_rhs(
        a_soln, a_rhs, simParams(), Geom().CellSize(0),
        simParams().recipe_params.grid_center, a_time);

    amrex::Gpu::streamSynchronize();
}

void RadialRecipeLevel::specificUpdateODE(amrex::MultiFab &a_soln)
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

void RadialRecipeLevel::pre_tag_cells()
{
    amrex::MultiFab &state_new = get_new_data(state_index);
    const auto cur_time        = get_state_data(state_index).curTime();
    FillPatch(*this, state_new, 2, cur_time, state_index, c_chi, 1);
}

void RadialRecipeLevel::tag_cells(amrex::TagBoxArray &a_tag_box_array,
                                  amrex::Real a_regrid_threshold)
{
    BL_PROFILE("RadialRecipeLevel::tag_cells()");
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

void RadialRecipeLevel::specificPostTimeStep()
{
    BL_PROFILE("RadialRecipeLevel::specificPostTimeStep");

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
        if (RadialRecipeMatter::uses_independent_scalars(simParams()))
        {
            GRTresnaIndependentScalars matter(
                simParams().recipe_num_scalar_fields,
                simParams().recipe_scalar_field_signs,
                simParams().recipe_scalar_mass);
            fill_matter_constraints(cst, state_new, matter, dx[0],
                                    simParams().recipe_params.grid_center,
                                    time);
        }
        else if (simParams().recipe_exotic_matter)
        {
            ExoticScalarField<DefaultPotential> exotic_scalar(
                DefaultPotential(), simParams().recipe_support_strength);
            fill_matter_constraints(cst, state_new, exotic_scalar, dx[0],
                                    simParams().recipe_params.grid_center,
                                    time);
        }
        else
        {
            ScalarField<DefaultPotential> scalar_field;
            fill_matter_constraints(cst, state_new, scalar_field, dx[0],
                                    simParams().recipe_params.grid_center,
                                    time);
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

        amrex::MultiFab cst_vac(state_new.boxArray(),
                                state_new.DistributionMap(), 4, 0);
        cst_vac.setVal(0.0);
        Constraints vacuum_constraints(dx[0], 0, Interval(1, 3));
        for (amrex::MFIter mfi(cst_vac, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst_vac.array(mfi);
            const auto src_arr   = state_new.const_array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int ix, int iy, int iz) noexcept
                { vacuum_constraints(ix, iy, iz, arr, src_arr); });
        }

        constexpr amrex::Real inv_16pi = 1.0 / (16.0 * M_PI);

        amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax,
                         amrex::ReduceOpSum, amrex::ReduceOpSum>
            rho_reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real>
            rho_reduce_data(rho_reduce_ops);
        using RhoReduceTuple = typename decltype(rho_reduce_data)::Type;

        for (amrex::MFIter mfi(cst_vac, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst_vac.const_array(mfi);
            rho_reduce_ops.eval(
                bx, rho_reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> RhoReduceTuple
                {
                    const amrex::Real ham_vac = arr(i, j, k, 0);
                    const amrex::Real rho_req = ham_vac * inv_16pi;
                    const amrex::Real neg_rho =
                        (rho_req < 0.0) ? (-rho_req * cell_vol) : 0.0;
                    return {rho_req, rho_req, neg_rho, cell_vol};
                });
        }

        auto rho_vals = rho_reduce_data.value();
        amrex::Real min_rho_req     = amrex::get<0>(rho_vals);
        amrex::Real max_rho_req     = amrex::get<1>(rho_vals);
        amrex::Real sum_neg_rho     = amrex::get<2>(rho_vals);
        amrex::Real sum_rho_vol     = amrex::get<3>(rho_vals);
        amrex::ParallelDescriptor::ReduceRealMin(min_rho_req);
        amrex::ParallelDescriptor::ReduceRealMax(max_rho_req);
        amrex::ParallelDescriptor::ReduceRealSum(sum_neg_rho);
        amrex::ParallelDescriptor::ReduceRealSum(sum_rho_vol);

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
            constraints_file.write_header_line(
                {"L2_Ham", "L2_Mom", "min_rho_req", "max_rho_req",
                 "integral_neg_rho"});
        }
        constraints_file.write_time_data_line(
            {L2_Ham, L2_Mom, static_cast<double>(min_rho_req),
             static_cast<double>(max_rho_req),
             static_cast<double>(sum_neg_rho)});
    }

    if (Level() == 0)
    {
        const amrex::Real time         = get_state_data(state_index).curTime();
        const amrex::Real dt           = parent->dtLevel(0);
        const amrex::Real restart_time = get_gramr_ptr()->get_restart_time();
        const bool first_step          = (time == 0.0);

        const int finest_lev = parent->finestLevel();
        auto &fine_level     = parent->getLevel(finest_lev);
        amrex::MultiFab &state_fine = fine_level.get_new_data(state_index);
        const auto &fine_geom       = parent->Geom(finest_lev);

        FillPatch(fine_level, state_fine, 2, time, state_index, 0,
                  state_fine.nComp());

        {
            const auto &arrs = state_fine.arrays();
            TraceARemoval trace_A_removal;
            PositiveChiAndLapse positive_chi_lapse(simParams().min_chi,
                                                   simParams().min_lapse);
            amrex::ParallelFor(state_fine, amrex::IntVect(0),
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

        const auto prob_lo = fine_geom.ProbLoArray();
        const auto dx_arr = fine_geom.CellSizeArray();

        // The apparent-horizon / expansion proxy must be measured about the
        // physics center (where the initial data is centered), not the
        // coordinate origin at the domain corner.  Using the corner makes
        // r ~ |grid_center| at the actual center, which collapses the
        // 2*sqrt(chi)/r regularizing term and produces spurious theta_plus<0
        // (false trapped surfaces) offset to r ~ |grid_center|.
        const amrex::Real cx = simParams().recipe_params.grid_center[0];
        const amrex::Real cy = simParams().recipe_params.grid_center[1];
        const amrex::Real cz = simParams().recipe_params.grid_center[2];

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
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

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
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

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
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
                 "min_phi", "max_phi", "min_Pi", "max_Pi"});
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
                                        static_cast<double>(max_Pi)});

        if (simParams().calculate_energy_conditions)
        {
            // In-situ observer-sampled energy conditions on the finest level.
            // We reconstruct the full 3+1 matter decomposition (rho, j_i, S_ij)
            // inline from the evolved state via ScalarField::compute_emtensor.
            // The spatial stresses S_ij require the evolved fields and are
            // unavailable to the t=0 Python proxies, which is why this lives in
            // C++.  We reduce the pointwise NEC/WEC/SEC/DEC margins to their
            // grid minima and the volume-integrated NEC violation.
            const amrex::Real ec_cell_vol =
                dx_arr[0] * dx_arr[1] * dx_arr[2];
            const amrex::Real ec_dx = dx_arr[0];

            // Use the same matter model as the evolution so that, when the
            // recipe requires exotic matter, the reported NEC/WEC genuinely
            // reflect the negative-energy phantom sector (rho <= 0) rather
            // than a canonical field that is never evolved.
            std::array<amrex::Real, 5> ec_res;
            if (RadialRecipeMatter::uses_independent_scalars(simParams()))
            {
                GRTresnaIndependentScalars matter(
                    simParams().recipe_num_scalar_fields,
                    simParams().recipe_scalar_field_signs,
                    simParams().recipe_scalar_mass);
                ec_res = reduce_ec_margins(state_fine, matter, ec_dx,
                                           ec_cell_vol);
            }
            else if (simParams().recipe_exotic_matter)
            {
                ExoticScalarField<DefaultPotential> exotic_scalar(
                    DefaultPotential(), simParams().recipe_support_strength);
                ec_res = reduce_ec_margins(state_fine, exotic_scalar, ec_dx,
                                           ec_cell_vol);
            }
            else
            {
                ScalarField<DefaultPotential> scalar_field;
                ec_res = reduce_ec_margins(state_fine, scalar_field, ec_dx,
                                           ec_cell_vol);
            }
            const amrex::Real min_nec      = ec_res[0];
            const amrex::Real min_wec      = ec_res[1];
            const amrex::Real min_sec      = ec_res[2];
            const amrex::Real min_dec      = ec_res[3];
            const amrex::Real integral_nec = ec_res[4];

            const std::string ec_prefix = out_dir + "energy_conditions";
            SmallDataIO ec_file(ec_prefix, dt, time, restart_time,
                                SmallDataIO::APPEND, first_step);
            ec_file.remove_duplicate_time_data();
            if (first_step)
            {
                // "matter_" prefix: these are the energy conditions of the
                // evolved matter sector (T_munu from compute_emtensor).  With
                // recipe_exotic_matter the matter is the phantom
                // ExoticScalarField and these margins go negative (NEC/WEC
                // violated) exactly where the geometry demands it; with the
                // canonical field they stay >= 0 by construction.  The
                // geometry-sourced effective stress energy T^eff = G/8pi is
                // evaluated post-hoc from plotfiles in the Python warpfactory
                // module (its spatial stresses require the time evolution).
                ec_file.write_header_line({"matter_min_NEC", "matter_min_WEC",
                                           "matter_min_SEC", "matter_min_DEC",
                                           "matter_integral_NEC_violation"});
            }
            ec_file.write_time_data_line(
                {static_cast<double>(min_nec), static_cast<double>(min_wec),
                 static_cast<double>(min_sec), static_cast<double>(min_dec),
                 static_cast<double>(integral_nec)});
        }

        if (simParams().calculate_curvature_invariants)
        {
            // Coordinate-invariant geometry diagnostics on the finest level:
            // the physical spatial Ricci scalar R^(3), the Ricci-tensor square
            // R_ij R^ij, and the extrinsic-curvature square K_ij K^ij.  Unlike
            // the matter energy conditions, these are sourced purely by the
            // evolved geometry and so directly probe the exotic tidal/curvature
            // structure of warp/wormhole candidates.  (Psi_4 is already emitted
            // to plotfiles via the Weyl4 derive.)
            const amrex::Real ci_dx       = dx_arr[0];
            const amrex::Real ci_cell_vol = dx_arr[0] * dx_arr[1] * dx_arr[2];

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax,
                             amrex::ReduceOpMax, amrex::ReduceOpSum>
                ci_ops;
            amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real,
                              amrex::Real>
                ci_red(ci_ops);
            using CITuple = typename decltype(ci_red)::Type;

            for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto st        = state_fine.const_array(mfi);
                ci_ops.eval(
                    bx, ci_red,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) -> CITuple
                    {
                        const FourthOrderDerivatives deriv(ci_dx);
                        const amrex::CellData<const amrex::Real> &cell =
                            st.cellData(i, j, k);
                        ScalarField<DefaultPotential>::Vars vars(cell);
                        const ScalarField<DefaultPotential>::D1Vars d1(
                            i, j, k, st, deriv);
                        const Tensor<2, amrex::Real> d2_chi =
                            deriv.diff2(i, j, k, st, c_chi);
                        const Tensor<4, amrex::Real> d2_h =
                            deriv.diff2_tensor(i, j, k, st, c_h11);

                        const auto h_UU =
                            CCZ4Geometry::compute_inverse_metric(vars);
                        const auto chris =
                            CCZ4Geometry::compute_christoffel(d1, h_UU);
                        const auto ricci = CCZ4Geometry::compute_ricci(
                            vars, d1, d2_chi, d2_h, h_UU, chris);

                        const amrex::Real R3 = ricci.scalar;

                        // Physical inverse spatial metric gamma^ij = chi h^ij.
                        const amrex::Real chi = vars.chi();
                        amrex::Real ricci_sq  = 0.0;
                        for (int a = 0; a < 3; ++a)
                            for (int b = 0; b < 3; ++b)
                                for (int c = 0; c < 3; ++c)
                                    for (int d = 0; d < 3; ++d)
                                        ricci_sq +=
                                            (chi * h_UU[a][c]) *
                                            (chi * h_UU[b][d]) *
                                            ricci.LL[a][b] * ricci.LL[c][d];

                        const amrex::Real Aij_sq =
                            CCZ4Geometry::compute_Aij_squared(vars, h_UU);
                        // K_ij K^ij = A_ij A^ij + K^2 / 3 (traceless split).
                        const amrex::Real K   = vars.K();
                        const amrex::Real KijKij =
                            Aij_sq + K * K / 3.0;

                        return {std::abs(R3), ricci_sq, KijKij,
                                R3 * R3 * ci_cell_vol};
                    });
            }

            auto ci_vals             = ci_red.value();
            amrex::Real max_abs_R3   = amrex::get<0>(ci_vals);
            amrex::Real max_ricci_sq = amrex::get<1>(ci_vals);
            amrex::Real max_KijKij   = amrex::get<2>(ci_vals);
            amrex::Real sum_R3sq_vol = amrex::get<3>(ci_vals);
            amrex::ParallelDescriptor::ReduceRealMax(max_abs_R3);
            amrex::ParallelDescriptor::ReduceRealMax(max_ricci_sq);
            amrex::ParallelDescriptor::ReduceRealMax(max_KijKij);
            amrex::ParallelDescriptor::ReduceRealSum(sum_R3sq_vol);

            const std::string ci_prefix = out_dir + "curvature_invariants";
            SmallDataIO ci_file(ci_prefix, dt, time, restart_time,
                                SmallDataIO::APPEND, first_step);
            ci_file.remove_duplicate_time_data();
            if (first_step)
            {
                ci_file.write_header_line({"max_abs_ricci_scalar",
                                           "max_ricci_tensor_sq",
                                           "max_Kij_sq", "L2_ricci_scalar"});
            }
            ci_file.write_time_data_line(
                {static_cast<double>(max_abs_R3),
                 static_cast<double>(max_ricci_sq),
                 static_cast<double>(max_KijKij),
                 static_cast<double>(std::sqrt(sum_R3sq_vol))});
        }
    }
}
