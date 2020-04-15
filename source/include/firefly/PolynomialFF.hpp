//==================================================================================
//    FireFly - Reconstructing rational functions and polynomial over finite fields.
//    Copyright (C) 2020  Jonas Klappert and Fabian Lange
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

#include "RationalNumber.hpp"
#include "UintHasher.hpp"
#include "ShuntingYardParser.hpp"

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
    FFInt calc_n_m_1(const std::vector<FFInt>& x);
    /**
     *  @return true if the PolynomialFF object has no coefficients or only one which is zero
     */
    bool zero() const;
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
    PolynomialFF mul(const uint32_t zi) const;
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
    PolynomialFF add_shift(const std::vector<FFInt>& shift) const;
    /**
     *  Removes any zero coefficient
     */
    void remove_zero_coefs();
    uint32_t n = 0; /**< An integer indicating the number of variables */
    ff_map coefs {};
    /**
     *  Transforms the Polynomial object to a string where each variable
     *  is replaced by the corresponding symbol in a given vector
     *  @param vars a vector of symbols, e.g. {"x","y","z"}.
     */
    std::string to_string(const std::vector<std::string>& vars) const;
  private:
    /**
     *  Calculates a binomial coefficient n over k
     *  @param n the n
     *  @param k the k
     */
    FFInt bin_coef(uint32_t n, uint32_t k) const;
    /**
     *  Multiplies two polynomials with no overlap of degrees
     *  @param a first polynimial
     *  @param b second polynomial
     *  @param curr_deg degree number of b which is absent in a
     */
    PolynomialFF mul_shift(const ff_map& a, const ff_map& b, uint32_t curr_deg) const;
    /**
     *  Executes prerequired steps to generate a Horner form
     */
    void generate_hornerff();
    /**
     *  Generates Horner scheme coefficients of a polynomial recursively
     *  @param var An integer representing the current variable
     *  @param monomials A map of monomials that build a polynomial
     *  @param vars the list of variables
     *  @return A Horner form of a polynomial
     */
    std::string generate_horner_coefs(int var, const ff_map& monomials, const std::vector<std::string>& vars);
    ShuntingYardParser s_y_fun; /**< A ShuntingYardParser object which is used to evaluate the polynomial in Horner form */
    std::vector<uint32_t> min_degree {}; /**< The minimal degree of the polynomial. Only available after min_deg() is called */
    std::vector<uint32_t> max_degree {}; /**< The maximal degree of the polynomial. Only available after max_deg() is called */
    std::vector<std::string> vars {}; /**< A vector holding the variables of this polynomial as strings */
    bool generate_new_horner = true; /**< Indicates whether one needs to generate a new Horner form */
    bool eval_horner = false; /**< Indicates whether one should evaluate this polynomial in Horner form */
  };

  PolynomialFF operator*(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator+(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator-(const PolynomialFF& a, const PolynomialFF& b);
  PolynomialFF operator*(const PolynomialFF& a, const FFInt&);
  PolynomialFF operator/(const PolynomialFF& a, const FFInt&);
  std::ostream& operator<<(std::ostream& out, const PolynomialFF& a);
}
