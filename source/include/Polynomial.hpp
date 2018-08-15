#pragma once

#include <vector>
#include <iostream>
#include "FFInt.hpp"

namespace firefly {

  class Polynomial {
  public:
    Polynomial();
    Polynomial(std::vector<FFInt> coef_);
    Polynomial operator+ (const Polynomial &);
    Polynomial operator- (const Polynomial &);
    Polynomial operator* (const Polynomial &);
    Polynomial &operator= (const Polynomial &);
    Polynomial operator* (const FFInt &);
    Polynomial operator/ (const FFInt &);
    int deg;
    std::vector<FFInt> coef {};
    FFInt calc(FFInt x);
  };

  std::ostream &operator<< (std::ostream &out, const Polynomial &a);
}
