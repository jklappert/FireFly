#include <iostream>
#include <vector>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "FFInt.hpp"
#include "Polynomial.hpp"
#include "Logger.hpp"
#include "ReconstHelper.hpp"
#include <gmpxx.h>
#include <string>
#include "utils.hpp"
#include "RationalNumber.hpp"

int main() {
  firefly::RatReconst rec_rat(1);
  firefly::PolyReconst rec_pol(1);
  auto vec_rat = rec_rat.reconst();
  auto vec_pol = rec_pol.reconst();
  INFO_MSG("f(x) = (" << vec_rat.first.at(0) << " + " << vec_rat.first.at(1) << " * x)/(" << vec_rat.second.at(0) << " + " << vec_rat.second.at(100) << " * x^4)");
  INFO_MSG(vec_pol.at(0) << " + " << vec_pol.at(2));

  return 0;
}
