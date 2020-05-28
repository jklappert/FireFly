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
#include "firefly/Polynomial.hpp"

namespace firefly {

  /**
   *  Generates a Horner form of a polynomial over a finite field
   *  @param monomials A ff_map representing a polynomial over a finite field
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_ff(const ff_map& monomials, const std::vector<std::string>& vars, uint32_t index = 0);
  /**
   *  Generates a Horner form of a polynomial over the rationals
   *  @param monomials A rn_map representing a polynomial over the rationals
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_rn(const rn_map& monomials, const std::vector<std::string>& vars, uint32_t index = 0);
  /**
   *  Generates a Horner form of a polynomial over the rationals
   *  @param monomials A vector of monomial objects representing a polynomial over the rationals
   *  @param vars A vector of the occurring variables
   *  @param index The index of the recursion. Should not be set when calling this functions.
   */
  std::string generate_horner_mon(const std::vector<Monomial> monomials, const std::vector<std::string>& vars, uint32_t index = 0);
}
