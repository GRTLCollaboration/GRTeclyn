/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef GRUTILS_HPP_
#define GRUTILS_HPP_

// #ifndef GR_SPACEDIM
// #define GR_SPACEDIM 3
// #endif
constexpr int GR_SPACEDIM = 3;

#ifndef DEFAULT_TENSOR_DIM
#define DEFAULT_TENSOR_DIM AMREX_SPACEDIM
#endif

// A function to return the right index for the tensors based on the
// ordering below 0: T11, 1: T12, 2: T13, 3: T22, 4: T23, 5: T33
#define VAR_IDX(ivar, i, j) ((ivar) + (i) + (j) + (((i) * (j) != 0) ? 1 : 0))

// A version for where the base reference for the tensor is 0
#define VAR_IDX0(i, j) VAR_IDX(0, (i), (j))

// NOLINTBEGIN(cppcoreguidelines-macro-usage)
// Fancy 'for' loop macros to iterate through spatial tensors
// use as "FOR(i, j) { ... }"
#define FOR1(IDX) for (int(IDX) = 0; (IDX) < DEFAULT_TENSOR_DIM; ++(IDX))
#define FOR2(IDX1, IDX2)                                                       \
    FOR1 (IDX1)                                                                \
        FOR1 (IDX2)
#define FOR3(IDX1, IDX2, IDX3)                                                 \
    FOR2 (IDX1, IDX2)                                                          \
        FOR1 (IDX3)
#define FOR4(IDX1, IDX2, IDX3, IDX4)                                           \
    FOR2 (IDX1, IDX2)                                                          \
        FOR2 (IDX3, IDX4)
#define FOR5(IDX1, IDX2, IDX3, IDX4, IDX5)                                     \
    FOR4 (IDX1, IDX2, IDX3, IDX4)                                              \
        FOR1 (IDX5)
#define DUMMYFOR() // prevents warning that appeared in debug mode:
                   // 'ISO C++11 requires at least one argument for the "..." in
                   // a variadic macro'

#define FOR2_SYM(IDX1, IDX2)                                                   \
    for (int(IDX1) = 0; (IDX1) < DEFAULT_TENSOR_DIM; ++(IDX1))                 \
        for (int(IDX2) = IDX1; (IDX2) < DEFAULT_TENSOR_DIM; ++(IDX2))

#define GET_MACRO6(_1, _2, _3, _4, _5, NAME, ...) NAME
#define FOR(...)                                                               \
    GET_MACRO6(__VA_ARGS__, FOR5, FOR4, FOR3, FOR2, FOR1, DUMMYFOR)(__VA_ARGS__)
// NOLINTEND(cppcoreguidelines-macro-usage)

#endif /* GRUTILS_HPP_*/
