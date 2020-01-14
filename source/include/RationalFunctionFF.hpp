//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
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

#include "PolynomialFF.hpp"

namespace firefly {
  /**
   * @class RationalFunctionFF
   * @brief A container class representing rational functions over finite fields
   */
  class RationalFunctionFF {
  public:
    /**
     *    Constructor of RationalFunctionFF
     *    @param n the polynomial of the numerator
     *    @param d the polynomial of the denominator
     */
    RationalFunctionFF(const PolynomialFF& n, const PolynomialFF& d);
    RationalFunctionFF();
    /**
    *  Transforms the Polynomial object to a string where each variable
    *  is replaced by the corresponding symbol in a given vector
    *  @param vars a vector of variables, e.g. {"x","y","z"}.
    */
    std::string to_string(const std::vector<std::string>& vars) const;
    PolynomialFF numerator;  /**< The coefficients of the numerator */
    PolynomialFF denominator; /**< The coefficients of the denominator */
  };

  std::ostream& operator<<(std::ostream& out, const RationalFunctionFF& rf);
}
