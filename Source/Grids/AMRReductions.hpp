/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef AMRREDUCTIONS_HPP_
#define AMRREDUCTIONS_HPP_

#include <AMReX_Amr.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>

#include <array>
#include <cmath>
#include <concepts>
#include <string>

namespace AMRReductions
{
namespace detail
{
// Square of f, or the sum of squares of a multi-component f
// NOLINTBEGIN(bugprone-easily-swappable-parameters)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
square(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy, int iz,
       int a_comp, int a_num_comps) noexcept
// NOLINTEND(bugprone-easily-swappable-parameters)
{
    amrex::Real sum_of_squares = 0.0;
    for (int comp = a_comp; comp < a_comp + a_num_comps; ++comp)
    {
        const amrex::Real value  = a_cell(ix, iy, iz, comp);
        sum_of_squares          += value * value;
    }
    return sum_of_squares;
}

// f, or the magnitude sqrt(sum_c f_c^2) of a multi-component f
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy, int iz,
      int a_comp, int a_num_comps) noexcept
{
    return (a_num_comps == 1)
               ? a_cell(ix, iy, iz, a_comp)
               : std::sqrt(square(a_cell, ix, iy, iz, a_comp, a_num_comps));
}
} // namespace detail

/* The reductions AMRReducer can compute. Each acts on the variable comp or,
 * if num_comps > 1, on the Euclidean magnitude sqrt(sum_c f_c^2) of the
 * components [comp, comp + num_comps), e.g. to combine the momentum
 * constraint components into a single quantity.
 *
 * A reduction provides
 * - reduce_op_t: how the values of different cells are combined
 *   (amrex::ReduceOpSum, amrex::ReduceOpMin or amrex::ReduceOpMax)
 * - cell_value(): the value of a single cell
 * - finalize(): the result, given the combined cell values and the volume of
 *   the domain
 *
 * Further reductions can be added by writing a similar struct. */

//! Minimum over all cells
struct Min
{
    int comp      = 0;
    int num_comps = 1;

    using reduce_op_t = amrex::ReduceOpMin;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    cell_value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
               int iz, amrex::Real /*a_cell_volume*/) const noexcept
    {
        return detail::value(a_cell, ix, iy, iz, comp, num_comps);
    }

    [[nodiscard]] static amrex::Real finalize(amrex::Real a_min,
                                              amrex::Real /*a_domain_volume*/)
    {
        return a_min;
    }
};

//! Maximum over all cells
struct Max
{
    int comp      = 0;
    int num_comps = 1;

    using reduce_op_t = amrex::ReduceOpMax;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    cell_value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
               int iz, amrex::Real /*a_cell_volume*/) const noexcept
    {
        return detail::value(a_cell, ix, iy, iz, comp, num_comps);
    }

    [[nodiscard]] static amrex::Real finalize(amrex::Real a_max,
                                              amrex::Real /*a_domain_volume*/)
    {
        return a_max;
    }
};

//! Volume integral
struct Integral
{
    int comp      = 0;
    int num_comps = 1;

    using reduce_op_t = amrex::ReduceOpSum;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    cell_value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
               int iz, amrex::Real a_cell_volume) const noexcept
    {
        return a_cell_volume *
               detail::value(a_cell, ix, iy, iz, comp, num_comps);
    }

    [[nodiscard]] static amrex::Real finalize(amrex::Real a_integral,
                                              amrex::Real /*a_domain_volume*/)
    {
        return a_integral;
    }
};

/** Volume-weighted p-norm, (int |f|^p dV [/ V])^(1/p). The exponent p is a
 *  template parameter so that the common cases (see L1Norm and L2Norm below)
 *  need no runtime branching or std::pow in every cell. */
template <int p> struct PNorm
{
    static_assert(p >= 1, "PNorm needs p >= 1");

    int comp                 = 0;
    int num_comps            = 1;
    bool normalize_by_volume = true;

    using reduce_op_t = amrex::ReduceOpSum;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
    cell_value(const amrex::Array4<amrex::Real const> &a_cell, int ix, int iy,
               int iz, amrex::Real a_cell_volume) const noexcept
    {
        if constexpr (p == 2)
        {
            return a_cell_volume *
                   detail::square(a_cell, ix, iy, iz, comp, num_comps);
        }
        else
        {
            const amrex::Real abs_value =
                std::abs(detail::value(a_cell, ix, iy, iz, comp, num_comps));
            if constexpr (p == 1)
            {
                return a_cell_volume * abs_value;
            }
            else
            {
                return a_cell_volume * std::pow(abs_value, p);
            }
        }
    }

    [[nodiscard]] amrex::Real finalize(amrex::Real a_integral,
                                       amrex::Real a_domain_volume) const
    {
        const amrex::Real integral =
            normalize_by_volume ? a_integral / a_domain_volume : a_integral;
        if constexpr (p == 1)
        {
            return integral;
        }
        else if constexpr (p == 2)
        {
            return std::sqrt(integral);
        }
        else
        {
            return std::pow(integral, 1.0 / p);
        }
    }
};

//! The most common norms. Being aliases for concrete types, these can be
//! written without template arguments, e.g. L2Norm{.comp = 0}.
using L1Norm = PNorm<1>;
using L2Norm = PNorm<2>;

//! The requirements on a reduction passed to AMRReducer::compute
template <class T>
concept Reduction =
    requires(const T &a_reduction,
             const amrex::Array4<amrex::Real const> &a_cell, amrex::Real a_real)
{
    typename T::reduce_op_t;
    {a_reduction.comp}->std::convertible_to<int>;
    {a_reduction.num_comps}->std::convertible_to<int>;
    {a_reduction.cell_value(a_cell, 0, 0, 0, a_real)}
        ->std::same_as<amrex::Real>;
    {a_reduction.finalize(a_real, a_real)}->std::same_as<amrex::Real>;
};

/** Derives state variables and/or derived quantities (as registered in the
 *  derive list, e.g. Constraints::name) at a_time into a single multi-level
 *  hierarchy, so that AMRReducer::compute can reduce all of them in one pass
 *  over the grid. A name may also be one component of a derived quantity
 *  (e.g. "Ham"), in which case the whole quantity is derived. Components are
 *  looked up by name with comp(). */
class DerivedHierarchy
{
  public:
    DerivedHierarchy(amrex::Amr &a_amr, amrex::Real a_time,
                     const amrex::Vector<std::string> &a_names);

    //! The index of the component called a_name; aborts if there is none
    [[nodiscard]] int comp(const std::string &a_name) const;

    //! One MultiFab per level, as AMRReducer::compute expects
    [[nodiscard]] amrex::Vector<const amrex::MultiFab *> data() const;

    [[nodiscard]] const amrex::Vector<std::string> &component_names() const
    {
        return m_component_names;
    }

  private:
    amrex::Vector<amrex::MultiFab> m_data;
    amrex::Vector<std::string> m_component_names;
};

/** TEMPORARY: How many cells to exclude next to the low and high face of the
 * domain in each direction, e.g. to exclude the cells whose evolution equation
 * is replaced by an approximate (Sommerfeld) outer boundary condition, but not
 *  those next to a reflective symmetry plane. The widths are in coarsest-level
 *  cells and are scaled by the refinement ratio on finer levels, so that the
 *  excluded layer has the same physical thickness on every level. Faces in
 *  periodic directions are never excluded. */
struct BoundaryExclusion
{
    amrex::IntVect lo{0};
    amrex::IntVect hi{0};
};

/** Computes reductions over a multi-level hierarchy of MultiFabs (e.g. a
 *  DerivedHierarchy or the result of amrex::Amr::derive). All the reductions
 *  passed to a single compute() call are computed in one pass over the grid,
 *  e.g.
 *
 *    const AMRReductions::DerivedHierarchy diagnostics(
 *        amr, time, {Constraints::name, "chi"});
 *    const auto [ham_l2_norm, chi_min] = reducer.compute(
 *        diagnostics.data(),
 *        AMRReductions::L2Norm{.comp = diagnostics.comp("Ham")},
 *        AMRReductions::Min{.comp = diagnostics.comp("chi")});
 *
 *  Cells covered by a finer level are excluded so that every point in space
 *  is counted once.
 *
 *  Layers of cells next to chosen faces of the domain can be excluded with
 *  a BoundaryExclusion (see below). The domain volume used by normalized
 *  norms excludes the same layers. */
class AMRReducer
{
  public:
    explicit AMRReducer(const amrex::Amr &a_amr,
                        const BoundaryExclusion &a_exclusion = {});

    //! a_geometry and a_refinement_ratios must cover every level of the data
    //! passed to compute() (sizing them to amr.max_level is fine)
    AMRReducer(amrex::Vector<amrex::Geometry> a_geometry,
               amrex::Vector<amrex::IntVect> a_refinement_ratios,
               const BoundaryExclusion &a_exclusion = {});

    //! Returns the result of each of a_reductions, in the same order.
    //! a_data holds one MultiFab per level, e.g. DerivedHierarchy::data(), or
    //! amrex::GetVecOfConstPtrs() of what Amr::derive returns.
    template <Reduction... ReductionTypes>
    [[nodiscard]] std::array<amrex::Real, sizeof...(ReductionTypes)>
    compute(const amrex::Vector<const amrex::MultiFab *> &a_data,
            const ReductionTypes &...a_reductions) const;

    //! Volume of the coarse-level domain (minus any excluded boundary
    //! layer), i.e. the volume the hierarchy would have if it were unrefined
    [[nodiscard]] amrex::Real domain_volume() const;

  private:
    amrex::Vector<amrex::Geometry> m_geometry;
    amrex::Vector<amrex::IntVect> m_refinement_ratios;
    BoundaryExclusion m_exclusion;
};
} // namespace AMRReductions

#include "AMRReductions.impl.hpp"

#endif /* AMRREDUCTIONS_HPP_ */
