// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include <string>
#include <vector>

#include "RationalNumber.hpp"
#include "Polynomial.hpp"

namespace firefly {
  /**
   * @class RationalFunction
   * @brief A container class representing rational functions
   */
  class RationalFunction {
  public:
    /**
     *    Constructor of RationalFunction
     *    @param n the polynomial of the numerator
     *    @param d the polynomial of the denominator
     */
    RationalFunction(Polynomial n, Polynomial d);
    RationalFunction();
    /**
    *  Transforms the Polynomial object to a string where each variable
    *  is replaced by the corresponding symbol in a given vector
    *  @param symbols a vector of symbols, e.g. {"x","y","z"}.
    */
    std::string to_string(const std::vector<std::string>& symbols) const;

    Polynomial numerator;  /**< The coefficients of the numerator */
    Polynomial denominator; /**< The coefficients of the denominator */
  };

  std::ostream& operator<<(std::ostream& out, const RationalFunction& rf);
}
