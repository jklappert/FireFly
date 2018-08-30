#pragma once

#include "RationalNumber.hpp"
#include "FFInt.hpp"

namespace firefly {

  class Monomial {
  public:
    Monomial(const std::vector<uint>& powers_, const RationalNumber& coef_);
    bool operator<(const Monomial&);
    bool operator>(const Monomial&);
    Monomial operator*(const Monomial&);
    std::vector<uint> powers;
    RationalNumber coef;
  };
}
