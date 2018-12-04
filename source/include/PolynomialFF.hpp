#pragma once

#include <iostream>
#include <unordered_map>
#include "RationalNumber.hpp"
#include "UintHasher.hpp"
#include "FFInt.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint>, FFInt, UintHasher> ff_map;
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;

  class PolynomialFF {
  public:
    PolynomialFF();
    PolynomialFF(uint n_, ff_map coefs_);
    PolynomialFF operator+(const PolynomialFF&);
    PolynomialFF operator-(const PolynomialFF&);
    PolynomialFF& operator=(const PolynomialFF&) = default;
    PolynomialFF& operator-=(const PolynomialFF&);
    PolynomialFF& operator+=(const PolynomialFF&);
    PolynomialFF operator*(const PolynomialFF&);
    PolynomialFF operator*(const FFInt&);
    PolynomialFF operator/(const FFInt&);
    uint n;
    FFInt calc(std::vector<FFInt> x);
    ff_map coefs;
    bool zero();
    PolynomialFF homogenize(uint degree);
    PolynomialFF mul(const uint zi);
    std::vector<uint> min_deg();
    std::vector<uint> max_deg();
    PolynomialFF add_shift(std::vector<FFInt> &shift);
  private:
    std::vector<uint> min_degree {};
    std::vector<uint> max_degree {};
  };

  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a);
}
