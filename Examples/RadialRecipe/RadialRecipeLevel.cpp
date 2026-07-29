#include "ExternalGridInitialData.hpp"
#include "RadialRecipeInitialData.hpp"
#include "RadialRecipeLevel.hpp"
#include "RadialRecipeMatterConstraints.hpp"
#include "RadialRecipeMatterDispatch.hpp"
#include "GRTresnaIndependentScalars.hpp"
#include "ComplexScalarField.hpp"
#include "CCZ4RHSWithMatter.hpp"
#include "ChiTagger.hpp"
#include "ConstraintsWithMatter.hpp"
#include "GRParmParse.hpp"
#include "PositiveChiAndLapse.hpp"
#include "ScalarField.hpp"
#include "ExoticScalarField.hpp"
#include "SmallDataIO.hpp"
#include "TraceARemoval.hpp"
#include "RLMatterPumpParams.hpp"
#include "RLPumpForce.hpp"
#include "Weyl4WithMatter.hpp"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Reduce.H>
#include <AMReX_Utility.H>
#include <AMReX_iMultiFab.H>
#include <algorithm>
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
// (Definition in RadialRecipeMatterConstraints.hpp)

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

namespace
{
// ---------------------------------------------------------------------------
// Composite (whole-hierarchy) constraint norms.
//
// The long-standing L2_Ham / L2_Mom in constraint_norms.dat are computed on
// LEVEL 0 ONLY, over every level-0 cell -- including the cells sitting
// underneath the refinement -- and averaged over the whole domain.  That is the
// right number for its actual consumer: it feeds RLRuntime::publish_cached_L2_Ham
// and hence the pump governor, where cheapness and run-to-run consistency are
// what matter.  It is NOT an accuracy measure, for three compounding reasons:
//
//   1. the refined levels never contribute, so violation living where the
//      matter actually is is invisible (finest dx is 2^max_level smaller);
//   2. covered coarse cells are still summed, so the active region is counted
//      at a resolution it is not evolved on, and its real residual never
//      appears at all;
//   3. the domain is overwhelmingly near-vacuum, so the volume-weighted mean is
//      diluted toward zero by empty space.
//
// The quantities below are the standard AMR composite: evaluate the constraints
// on every level, drop each coarse cell that has finer cells beneath it, weight
// by that level's own cell volume, accumulate.  Reported in NEW columns --
// cols 2/3 and the governor's input are deliberately left byte-for-byte alone,
// so this stays a diagnostics-only change and every pre-existing run remains
// comparable.  See research/neuralspacetime/Debug.md PART B fix item 0.
// ---------------------------------------------------------------------------
struct AmrConstraintNorms
{
    amrex::Real L2_Ham     = 0.0; // composite, covered cells + boundary dropped
    amrex::Real L2_Mom     = 0.0;
    amrex::Real L2_Ham_rel = 0.0; // / composite norm of the individual terms
    amrex::Real L2_Mom_rel = 0.0;
    amrex::Real Linf_Ham   = 0.0; // max |Ham| over the same set -- undiluted
    amrex::Real L2_Ham_ref = 0.0; // refined region (levels >= 1) only
};

// boundary_skip counts LEVEL-0 cells to drop at the outer domain boundary
// (their violation is boundary-condition, not physics); scaled per level.
// cst_level0, if given, is the caller's already-filled level-0 constraint
// MultiFab (same 8 components, same with_abs_terms).  Reusing it matters: level
// 0 holds ~99.9% of the cells, so recomputing it here would roughly double the
// per-step diagnostic cost for nothing.
AmrConstraintNorms
compute_amr_constraint_norms(amrex::Amr *parent,
                             const SimulationParameters &params,
                             amrex::Real time, int state_idx, int boundary_skip,
                             const amrex::MultiFab *cst_level0 = nullptr)
{
    AmrConstraintNorms out;

    amrex::Real sum_ham2 = 0.0, sum_mom2 = 0.0, sum_vol = 0.0;
    amrex::Real sum_ham_abs2 = 0.0, sum_mom_abs2 = 0.0;
    amrex::Real sum_ham2_ref = 0.0, sum_vol_ref = 0.0;
    amrex::Real max_ham = 0.0;

    const int finest = parent->finestLevel();
    const int n0     = parent->Geom(0).Domain().length(0);

    for (int lev = 0; lev <= finest; ++lev)
    {
        amrex::AmrLevel &level = parent->getLevel(lev);
        const auto dx          = level.Geom().CellSizeArray();

        amrex::MultiFab cst_local;
        const amrex::MultiFab *cstp = nullptr;
        if (lev == 0 && cst_level0 != nullptr)
        {
            cstp = cst_level0;
        }
        else
        {
            amrex::MultiFab &state = level.get_new_data(state_idx);
            // Safe here: Amr calls post_timestep on level 0 only after every
            // finer level has been advanced to the same time and averaged down.
            amrex::AmrLevel::FillPatch(level, state, 2, time, state_idx, 0,
                                       state.nComp());
            cst_local.define(state.boxArray(), state.DistributionMap(), 8, 0);
            cst_local.setVal(0.0);
            RadialRecipeMatter::fill_active_constraints(
                cst_local, state, params, dx[0],
                params.recipe_params.grid_center, time,
                /*with_abs_terms=*/true);
            cstp = &cst_local;
        }
        const amrex::MultiFab &cst = *cstp;

        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        // 1 exactly where this level's cell is covered by the next finer level.
        amrex::iMultiFab covered(cst.boxArray(), cst.DistributionMap(), 1, 0);
        if (lev < finest)
        {
            covered = amrex::makeFineMask(cst,
                                          parent->getLevel(lev + 1).boxArray(),
                                          parent->refRatio(lev), 0, 1);
        }
        else
        {
            covered.setVal(0);
        }

        // Outer boundary layer in this level's index space. Fine levels never
        // touch the domain edge, so in practice this only bites on level 0.
        const int scale     = level.Geom().Domain().length(0) / n0;
        amrex::Box interior = level.Geom().Domain();
        interior.grow(-boundary_skip * scale);
        const amrex::IntVect ilo = interior.smallEnd();
        const amrex::IntVect ihi = interior.bigEnd();

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpMax>
            ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                          amrex::Real, amrex::Real>
            data(ops);
        using Tuple = typename decltype(data)::Type;

        for (amrex::MFIter mfi(cst, amrex::TilingIfNotGPU()); mfi.isValid();
             ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = cst.const_array(mfi);
            const auto cov       = covered.const_array(mfi);
            ops.eval(
                bx, data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> Tuple
                {
                    if (cov(i, j, k) != 0 || i < ilo[0] || i > ihi[0] ||
                        j < ilo[1] || j > ihi[1] || k < ilo[2] || k > ihi[2])
                    {
                        return {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
                    }
                    const amrex::Real ham     = arr(i, j, k, 0);
                    const amrex::Real m1      = arr(i, j, k, 1);
                    const amrex::Real m2      = arr(i, j, k, 2);
                    const amrex::Real m3      = arr(i, j, k, 3);
                    const amrex::Real ham_abs = arr(i, j, k, 4);
                    const amrex::Real ma1     = arr(i, j, k, 5);
                    const amrex::Real ma2     = arr(i, j, k, 6);
                    const amrex::Real ma3     = arr(i, j, k, 7);
                    const amrex::Real mom2    = m1 * m1 + m2 * m2 + m3 * m3;
                    const amrex::Real mom_abs2 =
                        ma1 * ma1 + ma2 * ma2 + ma3 * ma3;
                    return {ham * ham * cell_vol, mom2 * cell_vol, cell_vol,
                            ham_abs * ham_abs * cell_vol, mom_abs2 * cell_vol,
                            std::abs(ham)};
                });
        }

        const auto v             = data.value();
        const amrex::Real lv_ham2 = amrex::get<0>(v);
        const amrex::Real lv_vol  = amrex::get<2>(v);
        sum_ham2     += lv_ham2;
        sum_mom2     += amrex::get<1>(v);
        sum_vol      += lv_vol;
        sum_ham_abs2 += amrex::get<3>(v);
        sum_mom_abs2 += amrex::get<4>(v);
        max_ham = std::max(max_ham, amrex::get<5>(v));
        if (lev >= 1)
        {
            sum_ham2_ref += lv_ham2;
            sum_vol_ref  += lv_vol;
        }
    }

    amrex::ParallelDescriptor::ReduceRealSum(sum_ham2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_mom2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_vol);
    amrex::ParallelDescriptor::ReduceRealSum(sum_ham_abs2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_mom_abs2);
    amrex::ParallelDescriptor::ReduceRealSum(sum_ham2_ref);
    amrex::ParallelDescriptor::ReduceRealSum(sum_vol_ref);
    amrex::ParallelDescriptor::ReduceRealMax(max_ham);

    if (sum_vol > 0.0)
    {
        out.L2_Ham = std::sqrt(sum_ham2 / sum_vol);
        out.L2_Mom = std::sqrt(sum_mom2 / sum_vol);
        const amrex::Real ham_abs = std::sqrt(sum_ham_abs2 / sum_vol);
        const amrex::Real mom_abs = std::sqrt(sum_mom_abs2 / sum_vol);
        out.L2_Ham_rel = (ham_abs > 0.0) ? out.L2_Ham / ham_abs : 0.0;
        out.L2_Mom_rel = (mom_abs > 0.0) ? out.L2_Mom / mom_abs : 0.0;
    }
    out.Linf_Ham = max_ham;
    out.L2_Ham_ref =
        (sum_vol_ref > 0.0) ? std::sqrt(sum_ham2_ref / sum_vol_ref) : 0.0;
    return out;
}
} // namespace

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

        // 8 components: Ham, Mom1-3, Ham_abs, Mom_abs1-3 for relative norms.
        amrex::MultiFab cst(state_new.boxArray(), state_new.DistributionMap(),
                            8, 0);
        cst.setVal(0.0);
        const auto dx = Geom().CellSizeArray();
        RadialRecipeMatter::fill_active_constraints(
            cst, state_new, simParams(), dx[0],
            simParams().recipe_params.grid_center, time,
            /*with_abs_terms=*/true);

        const amrex::Real cell_vol = dx[0] * dx[1] * dx[2];

        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum, amrex::ReduceOpSum,
                         amrex::ReduceOpSum>
            reduce_ops;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                          amrex::Real>
            reduce_data(reduce_ops);
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
                    const amrex::Real ham_abs = arr(i, j, k, 4);
                    const amrex::Real ma1 = arr(i, j, k, 5);
                    const amrex::Real ma2 = arr(i, j, k, 6);
                    const amrex::Real ma3 = arr(i, j, k, 7);
                    const amrex::Real mom2 = (m1 * m1 + m2 * m2 + m3 * m3);
                    const amrex::Real mom_abs2 =
                        (ma1 * ma1 + ma2 * ma2 + ma3 * ma3);
                    return {ham * ham * cell_vol, mom2 * cell_vol, cell_vol,
                            ham_abs * ham_abs * cell_vol,
                            mom_abs2 * cell_vol};
                });
        }

        auto [sum_ham2, sum_mom2, sum_vol, sum_ham_abs2, sum_mom_abs2] =
            reduce_data.value();
        amrex::ParallelDescriptor::ReduceRealSum(sum_ham2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mom2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_vol);
        amrex::ParallelDescriptor::ReduceRealSum(sum_ham_abs2);
        amrex::ParallelDescriptor::ReduceRealSum(sum_mom_abs2);

        const double L2_Ham =
            (sum_vol > 0.0) ? std::sqrt(sum_ham2 / sum_vol) : 0.0;
        const double L2_Mom =
            (sum_vol > 0.0) ? std::sqrt(sum_mom2 / sum_vol) : 0.0;
        const double L2_Ham_abs =
            (sum_vol > 0.0) ? std::sqrt(sum_ham_abs2 / sum_vol) : 0.0;
        const double L2_Mom_abs =
            (sum_vol > 0.0) ? std::sqrt(sum_mom_abs2 / sum_vol) : 0.0;
        const double L2_Ham_rel =
            (L2_Ham_abs > 0.0) ? (L2_Ham / L2_Ham_abs) : 0.0;
        const double L2_Mom_rel =
            (L2_Mom_abs > 0.0) ? (L2_Mom / L2_Mom_abs) : 0.0;

        RLRuntime::publish_cached_L2_Ham(L2_Ham);

        // Pump force L2 and governor (appended columns; see header below).
        // NOTE on the lapse: the pump adds S_A to ∂_t Π_A, a coordinate-time
        // rate, so the force density that actually enters ∂_t ρ and ∂_t j_i is
        // f = Σ_A s_A S_A (...) with NO lapse factor (the α in the "α f_⊥" of
        // the Duhamel bound is already absorbed converting the KG source J_A
        // to S_A = α J_A).  Matches ControllerReservoirMatter exactly.
        amrex::Real pump_force_L2 = 0.0;
        amrex::Real pump_fi_L2    = 0.0;
        const int n_fields =
            RadialRecipeMatter::uses_bicomplex_scalar(simParams()) ? 2 : 1;
        const bool bicomplex =
            RadialRecipeMatter::uses_bicomplex_scalar(simParams());
        const RLMatterPumpParams pump_diag =
            RadialRecipeMatter::build_rl_pump(simParams(), n_fields, time);
        const double governor_val = pump_diag.governor;
        if (pump_diag.num_sites >= 1)
        {
            amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum,
                             amrex::ReduceOpSum>
                pf_ops;
            amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> pf_data(
                pf_ops);
            const auto prob_lo = Geom().ProbLoArray();
            const auto dx_arr  = Geom().CellSizeArray();
            // Same centre-relative frame as Coordinates / RLLumpTracker / RHS.
            const amrex::Real cgx =
                static_cast<amrex::Real>(simParams().recipe_params.grid_center[0]);
            const amrex::Real cgy =
                static_cast<amrex::Real>(simParams().recipe_params.grid_center[1]);
            const amrex::Real cgz =
                static_cast<amrex::Real>(simParams().recipe_params.grid_center[2]);
            for (amrex::MFIter mfi(state_new, amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.validbox();
                const auto arr       = state_new.const_array(mfi);
                pf_ops.eval(
                    bx, pf_data,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        -> amrex::GpuTuple<amrex::Real, amrex::Real,
                                           amrex::Real>
                    {
                        const amrex::Real x = prob_lo[0] +
                                             (amrex::Real(i) + 0.5) * dx_arr[0] -
                                             cgx;
                        const amrex::Real y = prob_lo[1] +
                                             (amrex::Real(j) + 0.5) * dx_arr[1] -
                                             cgy;
                        const amrex::Real z = prob_lo[2] +
                                             (amrex::Real(k) + 0.5) * dx_arr[2] -
                                             cgz;
                        const amrex::Real lapse = arr(i, j, k, c_lapse);
                        const amrex::Real phi1p = arr(i, j, k, c_phi);
                        const amrex::Real phi2p = arr(i, j, k, c_phi2);
                        const amrex::Real Pi1p  = arr(i, j, k, c_Pi);
                        const amrex::Real Pi2p  = arr(i, j, k, c_Pi2);
                        const amrex::Real phi1m =
                            bicomplex ? arr(i, j, k, c_phi_m) : 0.0;
                        const amrex::Real phi2m =
                            bicomplex ? arr(i, j, k, c_phi2_m) : 0.0;
                        const amrex::Real Pi1m =
                            bicomplex ? arr(i, j, k, c_Pi_m) : 0.0;
                        const amrex::Real Pi2m =
                            bicomplex ? arr(i, j, k, c_Pi2_m) : 0.0;
                        RLPumpSources src;
                        if (bicomplex)
                        {
                            src = RLPumpForce::compute_bicomplex_sources(
                                pump_diag, x, y, z, time, lapse, phi1p, phi2p,
                                Pi1p, Pi2p, phi1m, phi2m, Pi1m, Pi2m);
                        }
                        else
                        {
                            src = RLPumpForce::compute_single_field_sources(
                                pump_diag, x, y, z, time, lapse, phi1p, phi2p,
                                Pi1p, Pi2p);
                        }
                        // Gravitational signs: +1 canonical, -1 phantom.
                        const amrex::Real f_perp =
                            (Pi1p * src.s1p + Pi2p * src.s2p) -
                            (Pi1m * src.s1m + Pi2m * src.s2m);

                        // f_i = Sum_A s_A S_A d_i phi_A.  state_new is
                        // FillPatch'd with 2 ghost cells above, so the same
                        // 4th-order centred stencil the RHS uses is available
                        // on the valid box.
                        auto d1c = [&](int comp, int dir) -> amrex::Real
                        {
                            const int di = (dir == 0) ? 1 : 0;
                            const int dj = (dir == 1) ? 1 : 0;
                            const int dk = (dir == 2) ? 1 : 0;
                            const amrex::Real fm2 =
                                arr(i - 2 * di, j - 2 * dj, k - 2 * dk, comp);
                            const amrex::Real fm1 =
                                arr(i - di, j - dj, k - dk, comp);
                            const amrex::Real fp1 =
                                arr(i + di, j + dj, k + dk, comp);
                            const amrex::Real fp2 =
                                arr(i + 2 * di, j + 2 * dj, k + 2 * dk, comp);
                            return (fm2 - 8.0 * fm1 + 8.0 * fp1 - fp2) /
                                   (12.0 * dx_arr[dir]);
                        };
                        amrex::Real fi2 = 0.0;
                        for (int dir = 0; dir < 3; ++dir)
                        {
                            amrex::Real f_dir =
                                src.s1p * d1c(c_phi, dir) +
                                src.s2p * d1c(c_phi2, dir);
                            if (bicomplex)
                            {
                                f_dir -= src.s1m * d1c(c_phi_m, dir) +
                                         src.s2m * d1c(c_phi2_m, dir);
                            }
                            fi2 += f_dir * f_dir;
                        }
                        return {f_perp * f_perp * cell_vol, fi2 * cell_vol,
                                cell_vol};
                    });
            }
            auto [sum_af2, sum_afi2, sum_pf_vol] = pf_data.value();
            amrex::ParallelDescriptor::ReduceRealSum(sum_af2);
            amrex::ParallelDescriptor::ReduceRealSum(sum_afi2);
            amrex::ParallelDescriptor::ReduceRealSum(sum_pf_vol);
            pump_force_L2 =
                (sum_pf_vol > 0.0) ? std::sqrt(sum_af2 / sum_pf_vol) : 0.0;
            pump_fi_L2 =
                (sum_pf_vol > 0.0) ? std::sqrt(sum_afi2 / sum_pf_vol) : 0.0;
        }

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

        // Whole-hierarchy norms (cols 12-17). These are the ones to quote as
        // an accuracy figure; cols 2-3 above are the governor's input and stay
        // exactly as they were. 4 level-0 cells of outer boundary are dropped.
        const AmrConstraintNorms amr_norms = compute_amr_constraint_norms(
            parent, simParams(), time, state_index, /*boundary_skip=*/4, &cst);

        const std::string prefix = out_dir + "constraint_norms";

        SmallDataIO constraints_file(prefix, dt, time, restart_time,
                                     SmallDataIO::APPEND, first_step);
        constraints_file.remove_duplicate_time_data();
        if (first_step)
        {
            // Append-only: existing readers use cols 1-3 (time, L2_Ham, L2_Mom)
            // and optionally 4-6. New columns follow.
            constraints_file.write_header_line(
                {"L2_Ham", "L2_Mom", "min_rho_req", "max_rho_req",
                 "integral_neg_rho", "L2_Ham_rel", "L2_Mom_rel",
                 "pump_force_L2", "governor", "pump_fi_L2", "L2_Ham_amr",
                 "L2_Mom_amr", "L2_Ham_amr_rel", "L2_Mom_amr_rel",
                 "Linf_Ham_amr", "L2_Ham_amr_ref"});
        }
        constraints_file.write_time_data_line(
            {L2_Ham, L2_Mom, static_cast<double>(min_rho_req),
             static_cast<double>(max_rho_req),
             static_cast<double>(sum_neg_rho), L2_Ham_rel, L2_Mom_rel,
             static_cast<double>(pump_force_L2), governor_val,
             static_cast<double>(pump_fi_L2),
             static_cast<double>(amr_norms.L2_Ham),
             static_cast<double>(amr_norms.L2_Mom),
             static_cast<double>(amr_norms.L2_Ham_rel),
             static_cast<double>(amr_norms.L2_Mom_rel),
             static_cast<double>(amr_norms.Linf_Ham),
             static_cast<double>(amr_norms.L2_Ham_ref)});
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

        // The dchi_dr stencil below reads i+-1 neighbours; on cells at the
        // finest level's outer boundary those are ghosts interpolated from
        // the coarse level, and the resulting theta_plus minimum parks at
        // the refinement edge (fixed r just outside recipe_basis_radius_max)
        // mimicking a growing horizon in runs whose interiors stay healthy.
        // Mask = 1 on this level's valid cells (faces shared with a sibling
        // box become 1 via FillBoundary), 0 on coarse-fine / domain ghosts;
        // the proxy is only evaluated where all six neighbours are valid.
        amrex::iMultiFab fine_mask(state_fine.boxArray(),
                                   state_fine.DistributionMap(), 1, 1);
        fine_mask.setVal(0);
        fine_mask.setVal(1, 0, 1, 0);
        fine_mask.FillBoundary(fine_geom.periodicity());

        for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.validbox();
            const auto arr       = state_fine.const_array(mfi);
            const auto mask_arr  = fine_mask.const_array(mfi);
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

                    const bool stencil_valid =
                        mask_arr(i - 1, j, k) != 0 && mask_arr(i + 1, j, k) != 0 &&
                        mask_arr(i, j - 1, k) != 0 && mask_arr(i, j + 1, k) != 0 &&
                        mask_arr(i, j, k - 1) != 0 && mask_arr(i, j, k + 1) != 0;

                    amrex::Real ah_radius = 0.0;
                    amrex::Real theta_plus_min_proxy = 1.0e30;
                    if (r > 1e-6 && stencil_valid)
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
            const auto mask_arr  = fine_mask.const_array(mfi);
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

                    const bool stencil_valid =
                        mask_arr(i - 1, j, k) != 0 && mask_arr(i + 1, j, k) != 0 &&
                        mask_arr(i, j - 1, k) != 0 && mask_arr(i, j + 1, k) != 0 &&
                        mask_arr(i, j, k - 1) != 0 && mask_arr(i, j, k + 1) != 0;

                    if (r <= 1e-6 || !stencil_valid)
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

        // Pump control-effort diagnostic: instantaneous power injected by the
        // PD/open-loop matter pump into the scalar momentum equations,
        //   P = integral alpha * f_perp * sqrt(gamma) dV,
        // with f_perp from RLPumpForce (same source law as the evolution RHS).
        // Zero when the pump is off (num_sites < 1) or past rl_pump_stop_time.
        amrex::Real pump_work = 0.0;
        {
            const bool bicomplex =
                RadialRecipeMatter::uses_bicomplex_scalar(simParams());
            const int n_fields = bicomplex ? 2 : 1;
            const RLMatterPumpParams pump =
                RadialRecipeMatter::build_rl_pump(simParams(), n_fields, time);
            if (pump.num_sites >= 1)
            {
                const amrex::Real cell_vol =
                    dx_arr[0] * dx_arr[1] * dx_arr[2];
                // Same centre-relative frame as Coordinates / RLLumpTracker / RHS.
                const amrex::Real cgx = static_cast<amrex::Real>(
                    simParams().recipe_params.grid_center[0]);
                const amrex::Real cgy = static_cast<amrex::Real>(
                    simParams().recipe_params.grid_center[1]);
                const amrex::Real cgz = static_cast<amrex::Real>(
                    simParams().recipe_params.grid_center[2]);
                amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops_pw;
                amrex::ReduceData<amrex::Real> reduce_data_pw(reduce_ops_pw);
                using ReduceTuplePW = typename decltype(reduce_data_pw)::Type;
                const amrex::Real pw_time = time;
                for (amrex::MFIter mfi(state_fine, amrex::TilingIfNotGPU());
                     mfi.isValid(); ++mfi)
                {
                    const amrex::Box &bx = mfi.validbox();
                    const auto arr       = state_fine.const_array(mfi);
                    reduce_ops_pw.eval(
                        bx, reduce_data_pw,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuplePW
                        {
                            const amrex::Real x =
                                prob_lo[0] + (amrex::Real(i) + 0.5) * dx_arr[0] -
                                cgx;
                            const amrex::Real y =
                                prob_lo[1] + (amrex::Real(j) + 0.5) * dx_arr[1] -
                                cgy;
                            const amrex::Real z =
                                prob_lo[2] + (amrex::Real(k) + 0.5) * dx_arr[2] -
                                cgz;
                            const amrex::Real lapse = arr(i, j, k, c_lapse);
                            const amrex::Real chi = amrex::max(
                                arr(i, j, k, c_chi), amrex::Real(1.0e-10));
                            const amrex::Real sqrt_gamma =
                                amrex::Real(1.0) / (chi * std::sqrt(chi));
                            const amrex::Real phi1p = arr(i, j, k, c_phi);
                            const amrex::Real phi2p = arr(i, j, k, c_phi2);
                            const amrex::Real Pi1p  = arr(i, j, k, c_Pi);
                            const amrex::Real Pi2p  = arr(i, j, k, c_Pi2);
                            const amrex::Real phi1m =
                                bicomplex ? arr(i, j, k, c_phi_m) : 0.0;
                            const amrex::Real phi2m =
                                bicomplex ? arr(i, j, k, c_phi2_m) : 0.0;
                            const amrex::Real Pi1m =
                                bicomplex ? arr(i, j, k, c_Pi_m) : 0.0;
                            const amrex::Real Pi2m =
                                bicomplex ? arr(i, j, k, c_Pi2_m) : 0.0;

                            RLPumpSources src;
                            if (bicomplex)
                            {
                                src = RLPumpForce::compute_bicomplex_sources(
                                    pump, x, y, z, pw_time, lapse, phi1p, phi2p,
                                    Pi1p, Pi2p, phi1m, phi2m, Pi1m, Pi2m);
                            }
                            else
                            {
                                src = RLPumpForce::compute_single_field_sources(
                                    pump, x, y, z, pw_time, lapse, phi1p, phi2p,
                                    Pi1p, Pi2p);
                            }
                            // f_perp with gravitational signs: +1 canonical,
                            // -1 phantom.  Power = f_perp * sqrt(gamma): the
                            // pump adds S_A to d_t Pi_A (a coordinate-time
                            // rate), so d_t rho|_pump = f_perp exactly, with
                            // NO lapse factor.  Same convention as
                            // ControllerReservoirMatter.
                            const amrex::Real f_perp =
                                (Pi1p * src.s1p + Pi2p * src.s2p) -
                                (Pi1m * src.s1m + Pi2m * src.s2m);
                            return {f_perp * sqrt_gamma * cell_vol};
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
                 "min_phi", "max_phi", "min_Pi", "max_Pi", "pump_work"});
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
                                        static_cast<double>(pump_work)});

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
                    simParams().recipe_scalar_mass,
                    simParams().recipe_scalar_lambda);
                ec_res = reduce_ec_margins(state_fine, matter, ec_dx,
                                           ec_cell_vol);
            }
            else if (RadialRecipeMatter::uses_complex_scalar(simParams()))
            {
                // NOTE: recipe_scalar_mu MUST be passed.  Until 2026-07-28 this
                // constructor omitted it, so every EC margin was evaluated with
                // the sextic term of the Q-ball potential silently zeroed --
                // a different V than the one the evolution integrates.
                ComplexScalarField matter(simParams().recipe_scalar_mass,
                                            simParams().recipe_scalar_lambda,
                                            simParams().recipe_scalar_sign,
                                            simParams().recipe_scalar_mu);
                ec_res = reduce_ec_margins(state_fine, matter, ec_dx,
                                           ec_cell_vol);
            }
            else if (RadialRecipeMatter::uses_bicomplex_scalar(simParams()))
            {
                // Same μ omission as above; at candidate-146 amplitudes the
                // sextic term is comparable to the mass/quartic terms.
                BiComplexScalarField matter(simParams().recipe_scalar_mass,
                                            simParams().recipe_scalar_lambda,
                                            simParams().recipe_scalar_mu);
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
