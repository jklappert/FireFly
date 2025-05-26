//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

#pragma once

#include "firefly/config.hpp"
#include "firefly/RationalNumber.hpp"

#include <cstdint>
#include <vector>

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
    std::vector<uint32_t> powers; /**< The degree of the monomial */
    RationalNumber coef; /**< The coefficient of the monomial */
  };
}
