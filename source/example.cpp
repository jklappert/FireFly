#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "Logger.hpp"

int main() {
  firefly::RatReconst rec_rat(1);
  firefly::PolyReconst rec_pol(1);

  try {
    auto rat_fun = rec_rat.reconst();
    auto pol_fun = rec_pol.reconst();
    std::cout << rat_fun;
    std::cout << "f(x) = " << pol_fun;
  } catch (std::exception &e) {
    ERROR_MSG(e.what());
  }

  return 0;
}
