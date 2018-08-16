#pragma once

#include "RationalNumber.hpp"
#include "FFInt.hpp"

namespace firefly{

  class Monomial{
  public:
    Monomial(RationalNumber &rn, uint deg_);
    Monomial(FFInt &ff, uint deg_);
    bool is_coef_rn = false;
    RationalNumber coef_rn;
    FFInt coef_ff;
    uint deg;
  };
}