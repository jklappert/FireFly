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
  //firefly::RatReconst rec (1);
  firefly::PolyReconst rec (1);
  auto vec = rec.reconst();
  INFO_MSG("f(x) = " << vec.at(0) << " + " << vec.at(1) << " * x");
  

  return 0;
}
