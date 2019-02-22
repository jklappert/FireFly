// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

#pragma once

#include "RationalNumber.hpp"
#include "FFInt.hpp"

namespace firefly {
  /**
   * @class Monomial
   * @brief A container class representing monomials
   */
  class Monomial {
  public:
    /**
     *  Constructor of the Monomial class
     *  @param powers_ the power of the monomial as a vector
     *  @param coef_ the coef of the monomial as a RationalNumber object
     */
    Monomial(const std::vector<uint32_t>& powers_, const RationalNumber& coef_);
    bool operator<(const Monomial&);
    bool operator>(const Monomial&);
    Monomial operator*(const Monomial&);
    Monomial operator-() const;
    std::vector<uint32_t> powers;
    RationalNumber coef;
  };
}
