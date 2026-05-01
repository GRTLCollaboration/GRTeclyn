/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef STATE_TYPES_HPP_
#define STATE_TYPES_HPP_

// We use this StateType to index various "groups" of state data. In our
// applications however, we only have one type of state data for all evolution
// vars. In general, amrex allows for multiple types of state data to be
// defined. These data types can differ by e.g. numbers of components, ghost
// cells, boundary conditions, and so on.
enum StateType
{
    state_index = 0,
    NUM_STATE_TYPE
};

#endif /* STATE_TYPES_HPP_ */