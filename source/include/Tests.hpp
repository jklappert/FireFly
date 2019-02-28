// ====================================================================
// This file is part of FireFly.
//
// FireFly is licenced under the GNU General Public License (GNU GPL)
// version 3.
// ====================================================================

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
