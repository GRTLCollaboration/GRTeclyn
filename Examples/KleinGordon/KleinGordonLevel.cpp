#include "KleinGordonLevel.hpp"
#include <AMReX_ParmParse.H>
#include <numeric>
// #include <Derive.H>

using namespace amrex;

amrex::Vector<std::string> KleinGordonLevel::plot_constraints;

namespace
{
struct WaveBCFill
{
    AMREX_GPU_DEVICE
    void operator()(const IntVect &iv, Array4<Real> const &dest,
                    const int /*dcomp*/, const int /*numcomp*/,
                    GeometryData const &geom, const Real /*time*/,
                    const BCRec * /*bcr*/, const int /*bcomp*/,
                    const int /*orig_comp*/) const
    {
        // removed because periodic

        // In this test, we only need to fill the x-direction bounary,
        // because it's periodic in other directions.  We also could
        // have used BCType::reflect_odd, and then we would not need to
        // do anything here.  However, this is an example of how to fill
        // external Dirichlet BC.
        // const int ilo = geom.Domain().smallEnd(0);
        // const int ihi = geom.Domain().bigEnd(0);
        // const auto [i,j,k] = iv.dim3();
        // if (i < ilo) {
        //     dest(i,j,k,0) = -dest(2*ilo-i-1,j,k,0);
        //     dest(i,j,k,1) = -dest(2*ilo-i-1,j,k,1);
        // }
        // if (i > ihi) {
        //     dest(i,j,k,0) = -dest(2*ihi-i+1,j,k,0);
        //     dest(i,j,k,1) = -dest(2*ihi-i+1,j,k,1);
        // }
    }
};

void wave_bcfill(Box const &bx, FArrayBox &data, const int dcomp,
                 const int numcomp, Geometry const &geom, const Real time,
                 const Vector<BCRec> &bcr, const int bcomp, const int scomp)
{
    GpuBndryFuncFab<WaveBCFill> gpu_bndry_func(WaveBCFill{});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}
} // namespace

// void
// KleinGordonLevel::variableSetUp ()
// {

//     BL_PROFILE("KleinGordonLevel::variableSetUp()");

//     amrex::Print() << "HERE in K-G variable setup \n";

//     const int nghost = simParams().num_ghosts;

//     desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
//                            StateDescriptor::Point, nghost, NUM_VARS,
//                            &cell_quartic_interp);

//     // int lo_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,    // external
//     Dirichlet
//     //                                        BCType::int_dir,    // periodic
//     //                                        BCType::int_dir) }; // periodic
//     // int hi_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::ext_dir,
//     //                                        BCType::int_dir,
//     //                                        BCType::int_dir) };
//     int lo_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::int_dir,    // periodic
//                                            BCType::int_dir,    // periodic
//                                            BCType::int_dir) }; // periodic
//     int hi_bc[BL_SPACEDIM] = {AMREX_D_DECL(BCType::int_dir,
//                                            BCType::int_dir,
//                                            BCType::int_dir) };

//     Vector<BCRec> bcs(NUM_VARS, BCRec(lo_bc, hi_bc));

//     StateDescriptor::BndryFunc bndryfunc(wave_bcfill);
//     bndryfunc.setRunOnGPU(true);

//     const int ncomp = simParams().nfields*2;
//     Vector<std::string> param_names(ncomp);// = {"u1", "v1", "u2", "v2"};

//     // amrex::Print() << ncomp <<  "\n";
//     // amrex::Print() << nfields <<  "\n";

//     for (int n = 0; n < simParams().nfields; n++)
//       {
//         char name[6];
// 	sprintf(name, "phi%d", n);
// 	param_names[2*n] = name;
// 	sprintf(name, "dphi%d", n);
// 	param_names[2*n+1] = name;
//       }

//     // for (int n = 0; n < ncomp; n++)
//     //   amrex::Print() << param_names[n] << " " << n <<  "\n";

//     amrex::Vector<std::string> name(NUM_VARS);
//     for (int i = 0; i < NUM_VARS; ++i)
//     {
//         name[i] = UserVariables::variable_names[i];
//     }

//     desc_lst.setComponent(State_Type, 0, name, bcs, bndryfunc);

//     //New diagnostic variable for testing interpolation between levels
//     (against analytic solution)

//     derive_lst.add(
//     		   "frac_error", amrex::IndexType::TheCellType(),
//        		   1, plot_constraints,
// 		   amrex::DeriveFuncFab(),		    // amrex null
// function
// 		   //		   derive_func_fab,
// 		   [=](const amrex::Box &box) { return amrex::grow(box,
// nghost);},
//     		   &amrex::cell_quartic_interp);

//     derive_lst.addComponent("frac_error", desc_lst, State_Type, 0, 1);

// }

// void KleinGordonLevel::variableCleanUp()
// {
//     desc_lst.clear();
//     derive_lst.clear();
// }

void KleinGordonLevel::initData()
{
    BL_PROFILE("KleinGordonLevel::initData()");

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    Real midpts[3];
    midpts[0] = 0.5 * (probhi[0] - problo[0]);
    midpts[1] = 0.5 * (probhi[1] - problo[1]);
    midpts[2] = 0.5 * (probhi[2] - problo[2]);

    MultiFab &S_new  = get_new_data(State_Type);
    auto const &snew = S_new.arrays();

    constexpr Real t0 = -5.4;

    amrex::Vector<amrex::Real> start_times{-5.4, 5.4};
    amrex::Vector<amrex::Real> start_pos{
        midpts[0], midpts[1], midpts[2] + 0.5 * midpts[2],
        midpts[0], midpts[1], midpts[2] - 0.5 * midpts[2]};

    // static_assert(std::is_trivially_copyable<BinaryBH>::value,
    //               "BinaryBH needs to be device copyable");

    InitialConditions SineGordon(simParams().alpha, simParams().k_r);

    int nfields = simParams().nfields;

    amrex::ParallelFor(
        S_new,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];
            Real z = problo[2] + (k + 0.5) * dx[2];

            //	Real rr2 = (x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)+(z-0.5)*(z-0.5);

            //        constexpr Real Pi = 3.1415926535897932384626;
            // constexpr Real mu = 0.7;
            // constexpr Real mu_coeff = mu/(std::sqrt(1-mu*mu));
            // constexpr Real t0 = 0;

            // constexpr Real k_r = 1;
            // constexpr Real omega = 1;
            for (int n = 0; n < nfields; n++)
            {
                //	    snew[bi](i,j,k,2*n) = 0.0;
                //	    snew[bi](i,j,k,2*n+1) = std::exp(-16.*rr2) *
                // std::pow(std::cos(Pi*rr2),6);

                // snew[bi](i,j,k,2*n) = 1+ampl[n]*std::exp(-width[n]*rr2);
                // snew[bi](i,j,k,2*n+1) = 0;

                snew[bi](i, j, k, 2 * n) =
                    SineGordon.breather_solution(x - midpts[0], 0);
                snew[bi](i, j, k, 2 * n + 1) =
                    SineGordon.breather_solution_deriv(x - midpts[0], 0);

                //	    snew[bi](i,j,k,2*n) =
                // SineGordon.breather_solution(x-start_pos[0], y-start_pos[1],
                // z-start_pos[2], start_times[0]) +
                // SineGordon.breather_solution(x-start_pos[3], y-start_pos[4],
                // z-start_pos[5], start_times[1]); snew[bi](i,j,k,2*n+1) = 0;
            }
        });
}

void KleinGordonLevel::specificAdvance()
{
    // amrex::MultiFab &S_new = get_new_data(State_Type);

    // // Check for nan's
    //   if (simParams().nan_check)
    //   {
    //       if (S_new.contains_nan(0, S_new.nComp(), amrex::IntVect(0), true))
    //       {
    //           amrex::Abort("NaN in specificAdvance");
    //       }
    //   }
}

void KleinGordonLevel::specificEvalRHS(amrex::MultiFab &a_soln,
                                       amrex::MultiFab &a_rhs,
                                       const double a_time)
{
    BL_PROFILE("KleinGordonLevel::specificEvalRHS()");

    const auto dxinv = Geom().InvCellSizeArray();
    AMREX_D_TERM(Real dx2inv = dxinv[0] * dxinv[0];
                 , Real dy2inv = dxinv[1] * dxinv[1];
                 , Real dz2inv = dxinv[2] * dxinv[2]);
    auto const &sa   = a_soln.arrays();
    auto const &sdot = a_rhs.arrays();

    Potential my_potential(simParams().scalar_mass);

    int nfields = simParams().nfields;

    amrex::ParallelFor(
        a_soln,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            auto const &s = sa[bi];
            auto const &f = sdot[bi];

            //      Real phi2 = std::pow(s(i,j,k,0),2)+std::pow(s(i,j,k,2),2);
            amrex::Real phi = 0;

            for (int n = 0; n < nfields; n++)
            {

                f(i, j, k, 2 * n) = s(i, j, k, 2 * n + 1);

                AMREX_D_TERM(
                    Real lapx = dx2inv * (-2.5 * s(i, j, k, 2 * n) +
                                          (4. / 3.) * (s(i - 1, j, k, 2 * n) +
                                                       s(i + 1, j, k, 2 * n)) -
                                          (1. / 12.) * (s(i - 2, j, k, 2 * n) +
                                                        s(i + 2, j, k, 2 * n)));
                    ,
                    Real lapy = dy2inv * (-2.5 * s(i, j, k, 2 * n) +
                                          (4. / 3.) * (s(i, j - 1, k, 2 * n) +
                                                       s(i, j + 1, k, 2 * n)) -
                                          (1. / 12.) * (s(i, j - 2, k, 2 * n) +
                                                        s(i, j + 2, k, 2 * n)));
                    , Real lapz =
                          dz2inv * (-2.5 * s(i, j, k, 2 * n) +
                                    (4. / 3.) * (s(i, j, k - 1, 2 * n) +
                                                 s(i, j, k + 1, 2 * n)) -
                                    (1. / 12.) * (s(i, j, k - 2, 2 * n) +
                                                  s(i, j, k + 2, 2 * n))));

                f(i, j, k, 2 * n + 1) = AMREX_D_TERM(lapx, +lapy, +lapz);

                f(i, j, k, 2 * n + 1) -= std::sin(s(i, j, k, 2 * n));
            }
        });
    Gpu::streamSynchronize();
}

void KleinGordonLevel::errorEst(TagBoxArray &tags, int /*clearval*/,
                                int /*tagval*/, Real /*time*/,
                                int /*n_error_buf*/, int /*ngrow*/)
{
    BL_PROFILE("KleinGordonLevel::errorEst()");

    auto const &S_new = get_new_data(State_Type);

    const char tagval = TagBox::SET;
    auto const &a     = tags.arrays();
    auto const &s     = S_new.const_arrays();
    amrex::ParallelFor(tags,
                       [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k)
                       {
                           // Just an example, not necessarily good choice.
                           if (amrex::Math::abs(s[bi](i, j, k, 1)) > 0.75)
                           {
                               a[bi](i, j, k) = tagval;
                           }
                       });
}

void KleinGordonLevel::derive(const std::string &name, amrex::Real time,
                              amrex::MultiFab &multifab, int dcomp)
{
    BL_PROFILE("KleinGordonLevel::derive()");

    InitialConditions SineGordon(simParams().alpha, simParams().k_r);

    const auto problo = geom.ProbLoArray();
    const auto probhi = geom.ProbHiArray();
    const auto dx     = geom.CellSizeArray();

    amrex::Real midpts[3];
    midpts[0] = simParams().center[0];
    midpts[1] = simParams().center[1];
    midpts[2] = simParams().center[2];

    auto const &sa = multifab.arrays();

    amrex::ParallelFor(
        multifab,
        [=] AMREX_GPU_DEVICE(int bi, int i, int j, int k) noexcept
        {
            auto const &s = sa[bi];

            amrex::Real x = problo[0] + (i + 0.5) * dx[0];
            amrex::Real y = problo[1] + (j + 0.5) * dx[1];
            amrex::Real z = problo[2] + (k + 0.5) * dx[2];

            amrex::Real exact_soln =
                SineGordon.breather_solution(x - midpts[0], time);

            s(i, j, k, dcomp) = exact_soln;
            // if (exact_soln==0)
            // 	 s_out(i,j,k,dcomp) = 0;
            // else
            //		       s_out(i,j,k,dcomp) =
            // std::log10(amrex::Math::abs((s(i,j,k,0)-exact_soln)/(exact_soln+1e-8)));//amrex::Math::abs(s(i,j,k,0)-exact_soln));
        });
}
