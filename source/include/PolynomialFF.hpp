//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2019  Jonas Klappert and Fabian Lange
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//==================================================================================

#pragma once

#include <iostream>
#include <unordered_map>
#include "RationalNumber.hpp"
#include "UintHasher.hpp"
#include "FFInt.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint32_t>, FFInt, UintHasher> ff_map;
  typedef std::unordered_map<std::vector<uint32_t>, RationalNumber, UintHasher> rn_map;

  /**
   * @class PolynomialFF
   * @brief A container class representing polynomials over finite fields
   */
  class PolynomialFF {
  public:
    /**
     *  Default constructor
     */
    PolynomialFF();
    /**
     *  Constructor with a given number of variables and coefficients
     *  @param n_ the number of variables
     *  @param coefs_ an ff_map type which holds all degrees and coefficients
     */
    PolynomialFF(uint32_t n_, ff_map coefs_);
    PolynomialFF& operator=(const PolynomialFF&) = default;
    PolynomialFF& operator-=(const PolynomialFF&);
    PolynomialFF& operator+=(const PolynomialFF&);
    PolynomialFF& operator*=(const FFInt&);
    PolynomialFF& operator/=(const FFInt&);
    uint32_t n = 0; /**< An integer indicating the number of variables */
    /**
     *  Evaluates the polynomial at a given parameter point
     *  @param x the parameter point at which the polynomial should be evaluated
     *  @return f(x)
     */
    FFInt calc(const std::vector<FFInt>& x) const;
    /**
     *  Evaluates the polynomial at a given parameter point omitting the first
     *  variable
     *  @param x the parameter point which is of length n - 1
     *  @return f(x)
     */
    FFInt calc_n_m_1(const std::vector<FFInt>& x) const;
    ff_map coefs {};
    /**
     *  @return true if the PolynomialFF object has no coefficients or only one which is zero
     */
    bool zero();
    /**
     *  Adds an additional functional dependence of a new variable with respect to
     *  a given degree. The sum of each variable degree is thus filled up to a
     *  given upper bound.
     *  @param degree the upper bound up to which the monomials should be homogenized
     */
    PolynomialFF homogenize(const uint32_t degree);
    /**
     *  Multiplies the PolynomialFF with a monomial of form zi^1
     *  @param zi the zi which should be multiplied with the PolynomialFF
     */
    PolynomialFF mul(const uint32_t zi);
    /**
     *  @return the minimal degree of the PolynomialFF
     */
    std::vector<uint32_t> min_deg();
    /**
     *  @return the maximal degree of the PolynomialFF
     */
    std::vector<uint32_t> max_deg();
    /**
     *  Calculates a new PolynomialFF if one shifts its variables
     *  @param shift the shift of each variable
     */
    PolynomialFF add_shift(const std::vector<FFInt>& shift);
  private:
    std::vector<uint32_t> min_degree {};
    std::vector<uint32_t> max_degree {};
    /**
     *  Calculates a binomial coefficient n over k
     *  @param n the n
     *  @param k the k
     */
    FFInt bin_coef(uint32_t n, uint32_t k);
    /**
     *  Multiplies two polynomials with no overlap of degrees
     *  @param a first polynimial
     *  @param b second polynomial
     *  @curr_deg degree number of b which is absent in a
     */
    PolynomialFF mul_shift(const ff_map& a, const ff_map& b, uint32_t curr_deg);
  };
  PolynomialFF operator*(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator+(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator-(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator*(const PolynomialFF& a, const FFInt&);
  PolynomialFF operator/(const PolynomialFF& a, const FFInt&);
  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a);
}
