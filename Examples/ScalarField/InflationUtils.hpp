/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef INFLATIONUTILS_HPP_
#define INFLATIONUTILS_HPP_

using namespace amrex;

namespace InflationUtils
{
    // Look-up table 
    // Used to construct polarisation basis tensors
    const Tensor<2, int> lut{0, 1, 2, 1, 3, 4, 2, 4, 5};
    const Real tolerance = 1.e-12;
}

#endif /* INFLATIONUTILS_HPP_ */