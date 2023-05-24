/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRINTERVAL_FDF5A7A_HPP_
#define GRINTERVAL_FDF5A7A_HPP_

// clang-format off

// Namespace to avoid conflicts with current code
namespace Old
{
/// A templated version of Chombo's Interval - allows compile time checking.
/**Note: iend is included in the interval, i.e. the interval <1,3> has
 * values 1,2,3 and therefore size 3.
 */
template <int ibegin, int iend> struct GRInterval
{
    static constexpr int begin() { return ibegin; }

    /// The largest component contained in the interval
    static constexpr int end() { return iend; }

    /// Returns the size of the interval
    static constexpr int size() { return iend - ibegin + 1; }

    /// Checks whether a values is in the interval (Note: the end component is
    /// also contained)
    static constexpr bool contains(int i)
    {
        return ((i > ibegin) && (i <= iend));
    }
};
} // namespace Old

#endif /* GRINTERVAL_HPP_ */
