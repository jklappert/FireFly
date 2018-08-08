#include <iostream>
#include <vector>
#include "PolyReconst.hpp"
#include "FFInt.hpp"
#include "Polynomial.hpp"

int main(int argc, char **argv) {
   std::vector<firefly::FFInt> v1;
   v1.emplace_back(firefly::FFInt(4,7));
   v1.emplace_back(firefly::FFInt(2,7));
   v1.emplace_back(firefly::FFInt(3,7));
   firefly::Polynomial p1(v1);
   std::vector<firefly::FFInt> v2;
   v2.emplace_back(firefly::FFInt(1,7));
   v2.emplace_back(firefly::FFInt(5,7));
   firefly::Polynomial p2(v2);
   firefly::Polynomial p3 = p2*p1;
   std::cout << p3 << "\n";
   //firefly::PolyReconst rec (1);
   //auto vec = rec.reconst();
   //std::cout << vec.at(0).n << " " << vec.at(1).n << " " << vec.at(2).n << std::endl;
   return 0;
}
