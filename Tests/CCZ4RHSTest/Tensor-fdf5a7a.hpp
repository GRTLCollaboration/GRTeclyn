/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TENSOR_FDF5A7A_HPP_
#define TENSOR_FDF5A7A_HPP_

// clang-format off
// NOLINTBEGIN

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"

// Namespace to avoid conflicts with current code
namespace Old
{
/// This class implements a Tensor with given rank, element data type, and
/// dimension.  By default the dimension is equal to DEFAULT_TENSOR_DIM.
template <int rank, class data_t, int size = DEFAULT_TENSOR_DIM> class Tensor
{
    template <int, class, int> friend class Tensor;
    typedef typename Tensor<rank - 1, data_t, size>::arr_t arr_t[size];
    arr_t arr;

  public:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE constexpr Tensor() {}

    //    ALWAYS_INLINE
    //    Tensor(std::initializer_list<data_t> list) :
    //        arr (list)
    //    {}

    template <typename... T>
    AMREX_GPU_HOST_DEVICE Tensor(T... data) : arr{data...}
    {
    }

    constexpr operator arr_t &() { return arr; }

    constexpr operator const arr_t &() const { return arr; }
};

template <class data_t, int size> class Tensor<0, data_t, size>
{
    template <int, class, int> friend class Tensor;
    typedef data_t arr_t;
    arr_t arr;

  public:
    constexpr Tensor() {}

    constexpr Tensor(data_t val) : arr(val) {}

    constexpr operator arr_t &() { return arr; }

    constexpr operator const arr_t &() const { return arr; }
};
} // namespace Old

// NOLINTEND

#endif /* TENSOR_FDF5A7A_HPP_ */
