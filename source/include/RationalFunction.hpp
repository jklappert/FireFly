#pragma once

#include "RationalNumber.hpp"
#include "Polynomial.hpp"

namespace firefly {

  class RationalFunction {
  public:
    /**
     *    Constructor of RationalFunction
     *    @param n the polynomial of the numerator
     *    @param d the polynomial of the denominator
     */
    RationalFunction(Polynomial n, Polynomial d);
    RationalFunction();
    Polynomial numerator;  /**< The coefficients of the numerator */
    Polynomial denominator; /**< The coefficients of the denominator */
  };

  std::ostream& operator<< (std::ostream& out, const RationalFunction& rf);
}
