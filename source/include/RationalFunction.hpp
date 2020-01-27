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
    RationalFunction(const Polynomial& n, const Polynomial& d);
    RationalFunction();
    /**
    *  Transforms the Polynomial object to a string where each variable
    *  is replaced by the corresponding symbol in a given vector
    *  @param vars a vector of variables, e.g. {"x","y","z"}.
    */
    std::string to_string(const std::vector<std::string>& vars) const;
    /**
     *  Generates a Horner form for this rational function
     *  @param vars A vector of the variables as strings
     *  @return The horner form for this rational function
     */
    std::string generate_horner(const std::vector<std::string>& vars) const;
    /**
     *  Adds a rational function as factor for this rational function
     *  @param factor the factor that will be added to the internal factors
     */
    void add_factor(const RationalFunction& factor);
    /**
     *  @return the factors of this rational function
     */
    std::vector<RationalFunction> get_factors() const;
    Polynomial numerator;  /**< The coefficients of the numerator */
    Polynomial denominator; /**< The coefficients of the denominator */
  private:
    std::vector<RationalFunction> factors {};
  };

  std::ostream& operator<<(std::ostream& out, const RationalFunction& rf);
}
