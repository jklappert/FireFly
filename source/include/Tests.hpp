#pragma once

#include "FFInt.hpp"
#include <vector>

namespace firefly{
  FFInt singular_solver(std::vector<FFInt> yis);
  FFInt n_eq_1(FFInt z1);
  FFInt n_eq_4(std::vector<FFInt> yis);
  FFInt gghh(std::vector<FFInt> yis);
}