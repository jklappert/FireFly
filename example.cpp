#include <iostream>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "Logger.hpp"

int main() {
  //firefly::RatReconst rec_rat(1);
  firefly::PolyReconst rec_pol(4);

  try {
    //auto rat_fun = rec_rat.reconst();
    auto pol_fun = rec_pol.reconst();
    //std::cout << rat_fun;
    std::cout << "f(x) = " << pol_fun;
  } catch (std::exception &e) {
    ERROR_MSG(e.what());
  }

  std::vector<std::pair<std::string, uint64_t>> replacements;
  replacements.push_back(std::make_pair("d", 3));

  firefly::FFInt d("d", 7, replacements);
  firefly::FFInt a("4", 7, replacements);
  firefly::FFInt b("497823478997845978789348974597878986789430893490308508340894568045608904568456", 7, replacements);

  INFO_MSG(d << " " << a << " " << b);

  try {
    firefly::FFInt c("sdf", 7, replacements);
  } catch (const std::exception &) {
    firefly::FFInt c("", 7, replacements);
  }

  return 0;
}
