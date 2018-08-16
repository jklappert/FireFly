#pragma once

#include <vector>
#include "FFInt.hpp"

namespace firefly {

  class MultiPolynomialFF {
  public:
    MultiPolynomialFF(int n, const std::vector<MultiPolynomialFF> &coef);
    MultiPolynomialFF();
    MultiPolynomialFF operator+(const MultiPolynomialFF &);
    int n;
    std::vector<MultiPolynomialFF> coef {};
    std::vector<int> monomials {};
    //sort function -> macht aus coef monomials
    // save in canonical form
  };
}
