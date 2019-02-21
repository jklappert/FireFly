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
}
