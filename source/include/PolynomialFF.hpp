#pragma once

#include <iostream>
#include <vector>
#include "MultiPolynomialFF.hpp"

namespace firefly {

  class PolynomialFF : public MultiPolynomialFF {
  public:
    PolynomialFF();
    PolynomialFF(std::vector<FFInt> coef_);
    PolynomialFF operator+(const PolynomialFF &);
    PolynomialFF operator-(const PolynomialFF &);
    PolynomialFF operator*(const PolynomialFF &);
    PolynomialFF &operator=(const PolynomialFF &);
    PolynomialFF operator*(const FFInt &);
    PolynomialFF operator/(const FFInt &);
    uint deg;
    FFInt calc(FFInt x);
    std::vector<FFInt> coef {};
  };

  std::ostream &operator<<(std::ostream &out, const PolynomialFF &a);
}
