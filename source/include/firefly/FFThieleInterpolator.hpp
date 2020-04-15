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

#include "PolynomialFF.hpp"

namespace firefly {

    /**
   * @class ThieleInterpolator
   * @brief A interpolator class for univariate rational function interpolation using Thiele's formula
   */
  class ThieleInterpolator {
  public:
    /**
     *  Constructor
     */
    ThieleInterpolator();
    /**
     *  Adds a point to the Thiele interpolation
     *  @param num the numerical value of the black box
     *  @param yi the numerical value at which the black box has been probed
     *  @return a bool which indicates if the interpolation terminated
     */
    bool add_point(const FFInt& num, const FFInt& yi);
    /**
     *  @return returns a pair of ff_maps which correspond to the numerator and denominator
     */
    std::pair<ff_map, ff_map> get_result();
    ThieleInterpolator(const ThieleInterpolator& other);
    ThieleInterpolator(ThieleInterpolator && other);
    ThieleInterpolator& operator=(const ThieleInterpolator& other);
    ThieleInterpolator& operator=(ThieleInterpolator && other);
  private:
    std::vector<FFInt> ai {}; /**< A vector which holds all coefficients a_i */
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
    /**
    *    Computes the coefficient a(i) = ai.at(i) recursively
    *    @param i The order of a(i)
    *    @param num f(y_i)
    *    @return a(i)
    */
    FFInt comp_ai(uint32_t i, const FFInt& num);
    /**
    *    Constructs the canonical form of the rational function recursivly
    *    @return the rational function in its canonical form
    */
    std::pair<ff_map, ff_map>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical();
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param i order of the highest coefficient a_i
     *    @param y y_i
     *    @returns f(y_i)
     */
    FFInt comp_fyi(uint32_t i, const FFInt& y);
  };
}
