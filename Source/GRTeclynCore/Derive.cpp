/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */
#include "GRAMRLevel.hpp"
#include "NewConstraints.hpp"

// NOLINTBEGIN(bugprone-easily-swappable-parameters)
std::unique_ptr<amrex::MultiFab> GRAMRLevel::derive(const std::string &name,
                                                    amrex::Real time, int ngrow)
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    std::unique_ptr<amrex::MultiFab> multifab;
    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        multifab = std::make_unique<amrex::MultiFab>(
            this->boxArray(), this->DistributionMap(), rec->numState(), ngrow);
        derive(name, time, *multifab, 0);
    }
    else
    {
        amrex::Abort("UnKnown derived variable");
    }
    return multifab;
}

void GRAMRLevel::derive(const std::string &name, amrex::Real /*time*/,
                        amrex::MultiFab &multifab, int dcomp)
{
    const amrex::DeriveRec *rec = derive_lst.get(name);
    if (rec != nullptr)
    {
        AMREX_ALWAYS_ASSERT(
            m_is_writing_plotfile); // We can relax this if needed.
        // If we are in the middle of writing a plotfile, ghost cells have been
        // filled.
        const auto &state_mf = get_new_data(State_Type);
        const auto &src      = state_mf.const_arrays();
        if (name == "constraints")
        {
            const auto &dst = multifab.arrays();
            int iham        = -1;
            Interval imom;
            if (!plot_constraints.empty())
            {
                int inext                = dcomp;
                auto plot_constraints_it = std::find(
                    plot_constraints.begin(), plot_constraints.end(), "Ham");
                if (plot_constraints_it != std::end(plot_constraints))
                {
                    iham = inext++;
                }
                plot_constraints_it = std::find(plot_constraints.begin(),
                                                plot_constraints.end(), "Mom1");
                if (plot_constraints_it != std::end(plot_constraints))
                {
                    imom = Interval(inext, inext + AMREX_SPACEDIM - 1);
                    auto plot_constraints_it2 =
                        std::find(plot_constraints.begin(),
                                  plot_constraints.end(), "Mom2");
                    auto plot_constraints_it3 =
                        std::find(plot_constraints.begin(),
                                  plot_constraints.end(), "Mom3");
                    AMREX_ALWAYS_ASSERT(
                        plot_constraints_it2 != std::end(plot_constraints) &&
                        plot_constraints_it3 != std::end(plot_constraints));
                }
            }
            else
            {
                iham = dcomp;
                imom = Interval(dcomp + 1, dcomp + AMREX_SPACEDIM);
            }
            Constraints cst(Geom().CellSize(0), iham, imom);
            amrex::ParallelFor(
                multifab, amrex::IntVect(0),
                [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept
                { cst.compute(i, j, k, dst[box_no], src[box_no]); });
        }
        else
        {
            amrex::Abort("UnKnown derived variable");
        }
    }
    else
    {
        amrex::Abort("UnKnown derived variable");
    }
    amrex::Gpu::streamSynchronize();
}
