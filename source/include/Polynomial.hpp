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

#include <string>

#include "PolynomialFF.hpp"
#include "Monomial.hpp"

namespace firefly {
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;

  /**
   * @class Polynomial
   * @brief A container class representing polynomials
   */
  class Polynomial {
  public:
    /**
     *    A constructor for a polynomial with RationalNumber objects as
     *    coefficients
     *    @param coef an unordered map of coefficients and degrees
     */
    Polynomial(const rn_map& coef);
    /**
     *    A constructor for a polynomial by a Monomial object
     *    @param coef a Monomial object
     */
    Polynomial(const Monomial& coef);
    Polynomial();
    Polynomial& operator*=(const RationalNumber&);
    /**
     *  Converts the Polynomial object to a PolynomialFF object
     */
    PolynomialFF convert_to_PolynomialFF();
    /**
     *  Sorts the degress of the polynomial in lexographically order
     */
    void sort();
    /**
     *  Clears the polynomial
     */
    void clear();
    /**
     *  Transforms the Polynomial object to a string where each variable
     *  is replaced by the corresponding symbol in a given vector
     *  @param vars a vector of symbols, e.g. {"x","y","z"}.
     */
    std::string to_string(const std::vector<std::string>& vars) const;
    /**
     *  Generates a Horner form for this polynomial
     *  @param vars A vector of the variables as strings
     *  @return The horner form for this polynomial
     */
    std::string generate_horner(std::vector<std::string> vars) const;
    std::vector<Monomial> coefs;
  private:
    uint32_t n;
  };

  std::ostream& operator<< (std::ostream& out, const Polynomial& pol);

}
