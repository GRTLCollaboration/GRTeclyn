/* GRTeclyn
 * Copyright 2022 The GRTL collaboration.
 * Please refer to LICENSE in GRTeclyn's root directory.
 */

#ifndef TPAMR_HPP_
#define TPAMR_HPP_

#include "BHAMR.hpp"

#ifdef USE_TWOPUNCTURES
#include "TwoPunctures.hpp"

class TPAMR : public BHAMR
{
  public:
    TP::TwoPunctures m_two_punctures;

    void set_two_punctures_parameters(const TP::Parameters &params)
    {
        m_two_punctures.Parameters::operator=(params);
    }
};

#endif /* USE_TWOPUNCTURES */

#endif /* TPAMR_HPP_ */
