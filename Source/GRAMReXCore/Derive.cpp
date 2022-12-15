#include "GRAMRLevel.hpp"
#include "NewConstraints.hpp"

std::unique_ptr<amrex::MultiFab>
GRAMRLevel::derive (const std::string& name, amrex::Real time, int ngrow)
{
    std::unique_ptr<amrex::MultiFab> mf;
    const amrex::DeriveRec* rec = derive_lst.get(name);
    if (rec) {
        mf = std::make_unique<amrex::MultiFab>(this->boxArray(),
                                               this->DistributionMap(),
                                               rec->numState(), ngrow);
        derive(name, time, *mf, 0);
    } else {
        amrex::Abort("UnKnown derived variable");
    }
    return mf;
}

void GRAMRLevel::derive (const std::string& name, amrex::Real /*time*/,
                         amrex::MultiFab& mf, int dcomp)
{
    const amrex::DeriveRec* rec = derive_lst.get(name);
    if (rec) {
        AMREX_ALWAYS_ASSERT(m_is_writing_plotfile); // We can relax this if needed.
        // If we are in the middle of writing a plotfile, ghost cells have been filled.
        auto const& state_mf = get_new_data(State_Type);
        auto const& src = state_mf.const_arrays();
        if (name == "constraints") {
            auto const& dst = mf.arrays();
            int iham = -1;
            Interval imom;
            if ( ! plot_constraints.empty() ) {
                int inext = dcomp;
                auto r = std::find(plot_constraints.begin(), plot_constraints.end(),
                                   "Ham");
                if (r != std::end(plot_constraints)) {
                    iham = inext++;
                }
                r = std::find(plot_constraints.begin(), plot_constraints.end(),
                              "Mom");
                if (r != std::end(plot_constraints)) {
                    imom = Interval(inext, inext+AMREX_SPACEDIM-1);
                    inext += AMREX_SPACEDIM;
                }
            } else {
                iham = dcomp;
                imom = Interval(dcomp+1, dcomp+AMREX_SPACEDIM);
            }
            Constraints cst(Geom().CellSize(0), iham, imom);
            amrex::ParallelFor(mf, amrex::IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                cst.compute(i,j,k,dst[box_no], src[box_no]);
            });
        } else {
            amrex::Abort("UnKnown derived variable");
        }
    } else {
        amrex::Abort("UnKnown derived variable");
    }
    amrex::Gpu::streamSynchronize();
}
