/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#if !defined(AMRREDUCTIONS_HPP_)
#error "This file should only be included through AMRReductions.hpp"
#endif

#ifndef AMRREDUCTIONS_IMPL_HPP_
#define AMRREDUCTIONS_IMPL_HPP_

#include "AMRReductions.hpp"

#include <AMReX_AmrLevel.H>
#include <AMReX_Derive.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <type_traits>
#include <utility>

namespace AMRReductions
{
// AMREX_D_TERM (rather than dx[0]*dx[0]*dx[0]) is kept deliberately: it makes
// the cell measure correct for whatever AMREX_SPACEDIM the code is compiled
// with, e.g. a future Cartoon-formalism build running on a lower-dimensional
// grid. Note that Cartoon norms would still need an extra symmetry-reduction
// measure factor (e.g. 2*pi*r) on top of this to be physically meaningful.
namespace detail
{
// The identity of the reduce_op_t of a_reduction_t, i.e. the value that a
// cell which is skipped must contribute (0 for sums, the largest number for
// minima, ...). AMReX's reduce operations define this as init().
template <Reduction a_reduction_t> amrex::Real identity()
{
    amrex::Real value{};
    typename a_reduction_t::reduce_op_t{}.init(value);
    return value;
}

template <class> inline constexpr bool always_false_v = false;

// amrex::ReduceOps only reduces over the data on this MPI rank, so the results
// still need combining across ranks with the matching MPI reduction.
template <class ReduceOp> void all_reduce(amrex::Real &a_value)
{
    const auto communicator = amrex::ParallelContext::CommunicatorSub();
    if constexpr (std::is_same_v<ReduceOp, amrex::ReduceOpSum>)
    {
        amrex::ParallelAllReduce::Sum(a_value, communicator);
    }
    else if constexpr (std::is_same_v<ReduceOp, amrex::ReduceOpMin>)
    {
        amrex::ParallelAllReduce::Min(a_value, communicator);
    }
    else if constexpr (std::is_same_v<ReduceOp, amrex::ReduceOpMax>)
    {
        amrex::ParallelAllReduce::Max(a_value, communicator);
    }
    else
    {
        static_assert(always_false_v<ReduceOp>,
                      "Only ReduceOpSum, ReduceOpMin and ReduceOpMax are "
                      "supported as a reduce_op_t");
    }
}

/* The per-cell work of compute(). This is a functor rather than a lambda so
 * that operator() can be force-inlined into the reduction loop: the compiler
 * would otherwise be free to call it out of line for every cell, which makes
 * the reduction several times slower. */
template <Reduction... ReductionTypes> struct ReductionKernel
{
    using tuple_t = amrex::GpuTuple<typename ReductionTypes::real_t...>;

    amrex::GpuTuple<ReductionTypes...> reductions;
    //! What a skipped cell returns (see identity())
    tuple_t identities;
    amrex::MultiArray4<amrex::Real const> data_arrays;
    //! Nonzero where a cell is covered by a finer level; evaluates to false
    //! on the finest level
    amrex::MultiArray4<int const> mask_arrays;
    amrex::Real cell_volume{};

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE tuple_t
    operator()(int box_no, int ix, int iy, int iz) const noexcept
    {
        if (mask_arrays && mask_arrays[box_no](ix, iy, iz) != 0)
        {
            return identities;
        }
        return cell_values(data_arrays[box_no], ix, iy, iz,
                           std::index_sequence_for<ReductionTypes...>{});
    }

    // A member function rather than a lambda so that it is inlined too
    template <std::size_t... Is>
    [[nodiscard]] AMREX_GPU_DEVICE AMREX_FORCE_INLINE tuple_t
    cell_values(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
                int iz, std::index_sequence<Is...> /*unused*/) const noexcept
    {
        return {amrex::get<Is>(reductions)
                    .cell_value(a_cell, ix, iy, iz, cell_volume)...};
    }
};
} // namespace detail

inline DerivedHierarchy::DerivedHierarchy(
    amrex::Amr &a_amr, amrex::Real a_time,
    const amrex::Vector<std::string> &a_names)
{
    // What to derive, and the first component of the result it fills
    struct Source
    {
        std::string name;
        int dcomp;
    };
    amrex::Vector<Source> sources;
    const auto is_included = [&](const std::string &a_name)
    {
        return std::ranges::any_of(sources, [&](const Source &a_source)
                                   { return a_source.name == a_name; });
    };

    for (const std::string &name : a_names)
    {
        const int dcomp = static_cast<int>(m_component_names.size());
        int state_index{};
        int state_comp{};
        if (amrex::AmrLevel::isStateVariable(name, state_index, state_comp))
        {
            if (!is_included(name))
            {
                sources.push_back({.name = name, .dcomp = dcomp});
                m_component_names.push_back(name);
            }
            continue;
        }

        // A derived quantity, or one of its components: either way the whole
        // quantity is derived (once, even if several of its components are
        // requested)
        const amrex::DeriveRec *record =
            amrex::AmrLevel::get_derive_lst().get(name);
        if (record == nullptr)
        {
            amrex::Abort("DerivedHierarchy: " + name +
                         " is neither a state variable nor a derived quantity");
        }
        if (!is_included(record->name()))
        {
            sources.push_back({.name = record->name(), .dcomp = dcomp});
            for (int icomp = 0; icomp < record->numDerive(); ++icomp)
            {
                m_component_names.push_back(record->variableName(icomp));
            }
        }
    }
    AMREX_ALWAYS_ASSERT(!m_component_names.empty());

    const int num_components = static_cast<int>(m_component_names.size());
    m_data.reserve(a_amr.finestLevel() + 1);
    for (int level = 0; level <= a_amr.finestLevel(); ++level)
    {
        m_data.emplace_back(a_amr.boxArray(level), a_amr.DistributionMap(level),
                            num_components, 0);
    }
    const auto data_pointers = amrex::GetVecOfPtrs(m_data);
    for (const Source &source : sources)
    {
        a_amr.derive(source.name, a_time, data_pointers, source.dcomp);
    }
}

inline int DerivedHierarchy::comp(const std::string &a_name) const
{
    const auto found = std::ranges::find(m_component_names, a_name);
    if (found == m_component_names.end())
    {
        std::string available;
        for (const std::string &name : m_component_names)
        {
            available += " " + name;
        }
        amrex::Abort("DerivedHierarchy: no component called " + a_name +
                     "; the components are" + available);
    }
    return static_cast<int>(found - m_component_names.begin());
}

inline amrex::Vector<const amrex::MultiFab *> DerivedHierarchy::data() const
{
    return amrex::GetVecOfConstPtrs(m_data);
}

inline AMRReducer::AMRReducer(const amrex::Amr &a_amr)
    : AMRReducer(a_amr.Geom(), a_amr.refRatio())
{
}

inline AMRReducer::AMRReducer(amrex::Vector<amrex::Geometry> a_geometry,
                              amrex::Vector<amrex::IntVect> a_refinement_ratios)
    : m_geometry(std::move(a_geometry)),
      m_refinement_ratios(std::move(a_refinement_ratios))
{
    AMREX_ALWAYS_ASSERT(!m_geometry.empty());
}

inline amrex::Real AMRReducer::domain_volume() const
{
    const auto dx                 = m_geometry[0].CellSizeArray();
    const amrex::Real cell_volume = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    return cell_volume * m_geometry[0].Domain().d_numPts();
}

template <Reduction... ReductionTypes>
std::array<amrex::Real, sizeof...(ReductionTypes)>
AMRReducer::compute(const amrex::Vector<const amrex::MultiFab *> &a_data,
                    const ReductionTypes &...a_reductions) const
{
    static_assert(sizeof...(ReductionTypes) > 0,
                  "AMRReducer::compute needs at least one reduction");
    AMREX_ALWAYS_ASSERT(!a_data.empty());
    AMREX_ALWAYS_ASSERT(std::none_of(a_data.begin(), a_data.end(),
                                     [](const amrex::MultiFab *a_mf)
                                     { return a_mf == nullptr; }));
    AMREX_ALWAYS_ASSERT(a_data.size() <= m_geometry.size());
    // m_refinement_ratios may be sized to amr.max_level rather than the
    // current finest level, so it only needs to be at least as long as
    // required to cover the levels actually present in a_data.
    AMREX_ALWAYS_ASSERT(a_data.size() <= m_refinement_ratios.size() + 1);
    const auto check_components = [&](const auto &a_reduction)
    {
        AMREX_ALWAYS_ASSERT(
            a_reduction.num_comps >= 1 && a_reduction.comp >= 0 &&
            a_reduction.comp + a_reduction.num_comps <= a_data[0]->nComp());
    };
    (check_components(a_reductions), ...);

    using index_sequence = std::index_sequence_for<ReductionTypes...>;

    detail::ReductionKernel<ReductionTypes...> kernel{
        .reductions = amrex::makeTuple(a_reductions...),
        .identities = {detail::identity<ReductionTypes>()...},
    };

    // As in amrex::volumeWeightedSum(), a single ReduceData accumulates the
    // reductions over all levels
    amrex::ReduceOps<typename ReductionTypes::reduce_op_t...> reduce_ops;
    amrex::ReduceData<typename ReductionTypes::real_t...> reduce_data(
        reduce_ops);

    for (int level = 0; level < a_data.size(); ++level)
    {
        const auto &data   = *a_data[level];
        const auto dx      = m_geometry[level].CellSizeArray();
        kernel.cell_volume = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
        kernel.data_arrays = data.const_arrays();

        // Cells covered by the next finer level are masked out; the finest
        // level has no mask.
        // The mask uses the asynchronous arena, as in volumeWeightedSum(), as
        // on GPUs it may be freed before the reduction kernel has finished.
        amrex::iMultiFab mask;
        kernel.mask_arrays = {};
        if (level + 1 < a_data.size())
        {
            mask = amrex::makeFineMask(
                data, *a_data[level + 1], amrex::IntVect(0),
                m_refinement_ratios[level], amrex::Periodicity::NonPeriodic(),
                0, 1, amrex::MFInfo().SetArena(amrex::The_Async_Arena()));
            kernel.mask_arrays = mask.const_arrays();
        }

        reduce_ops.eval(data, amrex::IntVect(0), reduce_data, kernel);
    }

    const auto totals        = reduce_data.value(reduce_ops);
    const amrex::Real volume = domain_volume();
    std::array<amrex::Real, sizeof...(ReductionTypes)> results{};
    [&]<std::size_t... Is>(std::index_sequence<Is...>)
    {
        ((results[Is] = amrex::get<Is>(totals),
          detail::all_reduce<typename ReductionTypes::reduce_op_t>(results[Is]),
          results[Is] = a_reductions.finalize(results[Is], volume)),
         ...);
    }(index_sequence{});
    return results;
}
} // namespace AMRReductions

#endif /* AMRREDUCTIONS_IMPL_HPP_ */
