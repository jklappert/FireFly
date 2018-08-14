#pragma once

#include <cstdint>
#include <vector>
#include <gmpxx.h>
#include "FFInt.hpp"
#include "Polynomial.hpp"
#include "RationalNumber.hpp"

namespace firefly {

  class PolyReconst {
  public:
    /**
     *    A constructor
     *    @param n_ The number of parameters as an integer
     */
    PolyReconst(int n_);
    /**
     *    Calls the reconstruction algorithm
     *    @return A vector with a pair of integers corresponding to the
     *    coefficients of the polynomial. The vector is ordered in an ascending
     *    way such that the first coefficient corresponds to x^0,...
     */
    std::vector<RationalNumber> reconst();

  private:
    /**
     *    The reconstruction algorithm for polynomials with coefficients which
     *    are members of a finite field
     *    @param prime The corresponding prime number
     *    @return A vector which contains the coefficients in ascending order
     */
    std::vector<mpz_class> reconst_ff(const uint64_t prime);
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param ai The vector of previously computed ai
     *    @param num f(y_i)
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @return a(i)
     */
    FFInt comp_ai(const std::vector<FFInt> &ai, const FFInt &num, int i, int ip);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @return The vector of coefficients of the canonical form
     */
    std::vector<FFInt> constr_canonical(const std::vector<FFInt> &ai, const uint64_t prime) const;
    /**
     *    Iterative construction of the canonical form
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @param i The iteration step; stops at ai.size()
     *    @return One iteration step of the canonical polynomial
     */
    Polynomial iterate_canonical(const std::vector<FFInt> &ai, const uint64_t prime, uint i) const;
    /**
     *    A numerical black box function which provides the reconstruction
     *    algorithm with the finite field member f(y)
     *    @param p The prime number which defines the finite field
     *    @param y A member of the finite field which should be evaluated in f(y)
     *    @return The numerical result in the finite field of prime
     */
    FFInt num(uint64_t p, const FFInt &y) const;
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param p The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const uint64_t prime);
    std::vector<mpz_class> convert_to_mpz(const std::vector<FFInt>& ai) const;
    std::vector<FFInt> convert_to_ffint(const uint64_t prime) const;
    const uint breakCondition = 3;
    int n; /**< The number of parameters */
    std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
    mpz_class combined_prime {}; /**< The combination of the used prime numbers
                                  with the chinese remained theorem */
    std::vector<mpz_class> combined_ci {}; /**< The combination of the finite field
                                  results with the chinese remained theorem */
    std::vector<RationalNumber> guessi {}; /**< The guesses of the rational coefficients */
  };
}
