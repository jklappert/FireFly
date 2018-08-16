#pragma once

#include <map>
#include "FFInt.hpp"

namespace firefly {

// Not tested yet!
  class MultiPolynomialFF {
  public:
    MultiPolynomialFF(uint n, const std::map<uint, MultiPolynomialFF> & coef_);
    MultiPolynomialFF();
    MultiPolynomialFF operator+(const MultiPolynomialFF &);
    MultiPolynomialFF operator-(const MultiPolynomialFF &);
    MultiPolynomialFF operator*(const FFInt &);
    uint n;
    uint deg;
    std::map<uint, MultiPolynomialFF> coef {};
    //std::vector<int> monomials {};
    //sort function -> macht aus coef monomials
    // save in canonical form
  };
}
