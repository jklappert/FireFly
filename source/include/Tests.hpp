#pragma once

#include "FFInt.hpp"
#include <vector>

namespace firefly{
  void black_box(std::vector<FFInt>& result, std::vector<FFInt> values);
  FFInt singular_solver(std::vector<FFInt> yis);
  FFInt n_eq_1(FFInt z1);
  FFInt n_eq_4(std::vector<FFInt> yis);
  FFInt gghh(std::vector<FFInt> yis);
  FFInt pol_n_eq_3(std::vector<FFInt> yis);
  FFInt ggh(std::vector<FFInt> yis);
}
