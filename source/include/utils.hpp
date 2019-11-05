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

#include "PolynomialFF.hpp"
#include "RationalNumber.hpp"

namespace firefly {
  /**
   *    Applies the cinese remainder theorem
   *    @param p1 a pair of a coefficient a and a prime p
   *    @param p2 a pair of a coefficient a and a prime p
   *    @returns the combination of p1 and p2 corresponding to the chinese
   *    remainder theorem
   */
  std::pair<mpz_class, mpz_class> run_chinese_remainder(
    const std::pair<mpz_class, mpz_class>& p1,
    const std::pair<mpz_class, mpz_class>& p2);

  /**
   *    Applies the rational reconstruction algorithm
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef(const mpz_class& a, const mpz_class& p);

  /**
   *    Applies the rational reconstruction algorithm MQRR from
   *    Maximal Quotient Rational Reconstruction: An Almost Optimal Algorithm for Rational Reconstruction
   *    by M. Monagan
   *    @param a a number over a finite field
   *    @param p a prime number defining the finite field
   *    @return a RationalNumber which has been reconstruction using the
   *    rational reconstruction algorithm
   */
  std::pair<bool, RationalNumber> get_rational_coef_mqrr(const mpz_class& a, const mpz_class& p);

  /**
  *  Solves the given modified Vandermonde system
  *  @param degs the contributing degrees
  *  @param nums the evaluated numerical values
  *  @param val the anchor points
  *  @return the polynomial
  */
  PolynomialFF solve_vandermonde_system(std::vector<std::vector<uint32_t>>& degs,
                                        const std::vector<FFInt>& nums,
                                        const std::vector<FFInt>& val);
  /**
   *  Compares two vetors colexographically, i.e. (1,0,0) < (0,1,0), and returns
   *  true if the first arguement is greater than the second
   *  @param a first vector which should be probed if it its greater
   *  @param b the reference vector for the comparison
   */
  bool a_grt_b(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);

  /**
   *  Compares two vetors colexographically, i.e. (1,0,0) < (0,1,0), and returns
   *  true if the first arguement is greater than the second. This function
   *  is particulary written for tuples of 0 and 1
   *  @param a first vector which should be probed if it its greater
   *  @param b the reference vector for the comparison
   */
  bool a_grt_b_s(const std::vector<uint32_t>& a, const std::vector<uint32_t>& b);

  /**
   *  Generates a vector of possible tuples of 1 and 0 for a given length r
   *  @param r the length of the vector
   */
  std::vector<std::vector<uint32_t>> generate_possible_shifts(uint32_t r);

  // TODO
  uint32_t compute_bunch_size(const uint32_t queue_length, const uint32_t thr_n, const uint32_t bunch_size);
  uint32_t compute_job_number(const uint32_t queue_length, const uint32_t threads, const uint32_t total_threads, const uint32_t bunch_size);

#ifdef DEFAULT
  uint64_t mod_mul(uint64_t a, uint64_t b, uint64_t m);
  /**
   *  Performs a exponentiation modulo m
   *  @param base the base
   *  @param exp the exponent
   *  @param m the modulus
   *  @return (base^exp) mod m
   */
  uint64_t mod_pow(uint64_t base, uint64_t exp, uint64_t m);
  /**
   *  Calculates the multiplicative inverse using the Extended Euclidean Algorithm
   *  @param a the integer of which the multiplicative inverse should be calculated
   *  @param m the modulus
   */
  uint64_t mod_inv(uint64_t a, uint64_t m);
#endif
}
