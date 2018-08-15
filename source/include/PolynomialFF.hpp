#pragma once

#include <vector>
#include <iostream>
#include "FFInt.hpp"

namespace firefly {

  class PolynomialFF {
  public:
    PolynomialFF();
    PolynomialFF(std::vector<FFInt> coef_);
    PolynomialFF operator+ (const PolynomialFF &);
    PolynomialFF operator- (const PolynomialFF &);
    PolynomialFF operator* (const PolynomialFF &);
    PolynomialFF &operator= (const PolynomialFF &);
    PolynomialFF operator* (const FFInt &);
    PolynomialFF operator/ (const FFInt &);
    int deg;
    std::vector<FFInt> coef {};
    FFInt calc(FFInt x);
  };

  std::ostream &operator<< (std::ostream &out, const PolynomialFF &a);
}
