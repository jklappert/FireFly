#pragma once

#include <vector>
#include <iostream>
#include "FFInt.hpp"

namespace firefly{

   class Polynomial{
   public:
      Polynomial();
      Polynomial(std::vector<FFInt> coef_);
      Polynomial operator+(const Polynomial&);
      Polynomial operator-(const Polynomial&);
      Polynomial operator*(const Polynomial&);
      Polynomial& operator=(const Polynomial&);
      int deg;
      std::vector<FFInt> coef {};
   };

   std::ostream& operator<<(std::ostream& out, const Polynomial& a);
}