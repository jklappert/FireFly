#include <iostream>
#include <vector>
#include "PolyReconst.hpp"
#include "RatReconst.hpp"
#include "FFInt.hpp"
#include "Polynomial.hpp"
#include "Logger.hpp"

int main() {
  /*   std::vector<firefly::FFInt> v1;
     v1.emplace_back(firefly::FFInt(4,7));
     v1.emplace_back(firefly::FFInt(2,7));
     v1.emplace_back(firefly::FFInt(3,7));
     firefly::Polynomial p1(v1);
     std::vector<firefly::FFInt> v2;
     v2.emplace_back(firefly::FFInt(1,7));
     v2.emplace_back(firefly::FFInt(5,7));
     firefly::Polynomial p2(v2);
     firefly::Polynomial p3 = p2*p1;
     std::cout << p3 << "\n";*/
  firefly::RatReconst rec (1);
  //firefly::PolyReconst rec (1);
  auto vec = rec.reconst();
  rec.constrCanonical();
  INFO_MSG ("Coefficient size: " << vec.size());
  INFO_MSG("f(x) = (" << rec.canonical.first << ")/(" << rec.canonical.second << ")");
  //rec.constrCanonical();
  //std::cout << rec.canonical << std::endl;

  return 0;
}
