#pragma once

#include <iostream>
#include "UintHasher.hpp"
#include <unordered_map>
#include "FFInt.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint>, FFInt, UintHasher> ff_map;

  class PolynomialFF {
  public:
    PolynomialFF();
    PolynomialFF(uint n_, ff_map coef_);
    PolynomialFF operator+(const PolynomialFF&);
    PolynomialFF operator-(const PolynomialFF&);
    PolynomialFF& operator=(const PolynomialFF&) = default;
    PolynomialFF operator*(const FFInt&);
    PolynomialFF operator/(const FFInt&);
    uint n;
    FFInt calc(std::vector<FFInt> x);
    ff_map coef;
    bool zero();
    PolynomialFF mul(const uint zi);
    std::vector<uint> min_deg();
    std::vector<uint> max_deg();
  private:
    std::vector<uint> min_degree {};
    std::vector<uint> max_degree {};
  };

  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a);
}
