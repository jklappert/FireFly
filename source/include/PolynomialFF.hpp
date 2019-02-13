#pragma once

#include <iostream>
#include <unordered_map>
#include "RationalNumber.hpp"
#include "UintHasher.hpp"
#include "FFInt.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint32_t>, FFInt, UintHasher> ff_map;
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;

  class PolynomialFF {
  public:
    PolynomialFF();
    PolynomialFF(uint32_t n_, ff_map coefs_);
    PolynomialFF& operator=(const PolynomialFF&) = default;
    PolynomialFF& operator-=(const PolynomialFF&);
    PolynomialFF& operator+=(const PolynomialFF&);
    PolynomialFF operator*(const PolynomialFF&);
    PolynomialFF operator*(const FFInt&);
    PolynomialFF operator/(const FFInt&);
    uint32_t n = 0;
    FFInt calc(std::vector<FFInt> x);
    ff_map coefs {};
    bool zero();
    PolynomialFF homogenize(uint32_t degree);
    PolynomialFF mul(const uint32_t zi);
    std::vector<uint32_t> min_deg();
    std::vector<uint32_t> max_deg();
    PolynomialFF add_shift(const std::vector<FFInt>& shift);
  private:
    std::vector<uint32_t> min_degree {};
    std::vector<uint32_t> max_degree {};
    FFInt bin_coef(uint32_t n, uint32_t k);
  };
    PolynomialFF operator+(const PolynomialFF& a, const PolynomialFF& b);
    PolynomialFF operator-(const PolynomialFF& a, const PolynomialFF& b);

  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a);
}
