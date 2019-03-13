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

#include "FFInt.hpp"
#include <vector>

namespace firefly {
  /**
   *  example for n = 4 which uses the singular_solver
   */
  FFInt singular_solver(std::vector<FFInt> yis);
  /**
   *  example for n = 1
   */
  FFInt n_eq_1(FFInt z1);
  /**
   *  example for n = 4 and the usage of the Chinese Remainder Theorem
   */
  FFInt n_eq_4(std::vector<FFInt> yis);
  /**
   *  example for a large interpolation problem augmented with large coefficients
   */
  FFInt gghh(std::vector<FFInt> yis);
  /**
   *  example for the reconstruction of a polynomial
   */
  FFInt pol_n_eq_3(std::vector<FFInt> yis);
  /**
   *  example for a three loop gg -> h integral coefficient
   */
  FFInt ggh(std::vector<FFInt> yis);
  /**
   *  a benchmark function with 20 variables
   */
  FFInt bench_1(std::vector<FFInt> yis);
  /**
   *  a benchmark function with 5 variables and almost complete dense numerator
   */
  FFInt bench_2(std::vector<FFInt> yis);
  /**
   *  a benchmark function with 5 variables and almost complete sparse with high degrees
   */
  FFInt bench_3(std::vector<FFInt> yis);
}
