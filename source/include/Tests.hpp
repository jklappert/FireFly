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
  FFInt singular_solver(const std::vector<FFInt>& yis);
  /**
   *  example for n = 1
   */
  FFInt n_eq_1(const FFInt& z1);
  /**
   *  example for n = 4 and the usage of the Chinese Remainder Theorem
   */
  FFInt n_eq_4(const std::vector<FFInt>& yis);
  /**
   *  example for a large interpolation problem augmented with large coefficients
   */
  FFInt gghh(const std::vector<FFInt>& yis);
  /**
   *  examples for the reconstruction of a polynomial with 3 variables
   */
  FFInt pol_n_eq_3(const std::vector<FFInt>& yis);
  /**
   *  example for the reconstruction of a sparse polynomial with 20 variables
   */
  FFInt pol_20_20(const std::vector<FFInt>& yis);
  /**
   *  example for reconstruction of a polynomial with 9 variables
   */
  FFInt pol_1(const std::vector<FFInt>& yis);
  /**
   *  example for reconstruction of a more complex polynomial with 5 variables
   */
  FFInt pol_6(const std::vector<FFInt>& yis);
  /**
   *  example for reconstuction of a polynomial with 5 variables
   */
  FFInt pol_7(const std::vector<FFInt>& yis);
  /**
   *  example for a three loop gg -> h integral coefficient
   */
  FFInt ggh(const std::vector<FFInt>& yis);
  /**
   *  a benchmark function with 20 variables
   */
  FFInt bench_1(const std::vector<FFInt>& yis);
  /**
   *  a benchmark function with 5 variables and almost complete dense numerator
   */
  FFInt bench_2(const std::vector<FFInt>& yis);
  /**
   *  a benchmark function with 5 variables and almost complete sparse with high degrees
   */
  FFInt bench_3(const std::vector<FFInt>& yis);
}
