#pragma once

#include <vector>
#include "RationalNumber.hpp"

namespace firefly {

  class Polynomial {
  public:
    /**
     *    A constructor for a polynomial with RationalNumber objects as
     *    coefficients
     *    @param coefs_ a vector of RationalNumber coefficients for the polynomial
     *    in ascending order (x^0, x^2,...)
     */
    Polynomial(std::vector<RationalNumber> coefs_);
    std::vector<RationalNumber> coefs;  /**< The vector which holds all coefficients*/
  };

  std::ostream &operator<< (std::ostream &out, const Polynomial &pol);

}
