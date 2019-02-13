#pragma once

#include "RationalNumber.hpp"
#include "FFInt.hpp"

namespace firefly {

  class Monomial {
  public:
    Monomial(const std::vector<uint32_t>& powers_, const RationalNumber& coef_);
    bool operator<(const Monomial&);
    bool operator>(const Monomial&);
    Monomial operator*(const Monomial&);
    Monomial operator-() const;
    std::vector<uint32_t> powers;
    RationalNumber coef;
  };
}
