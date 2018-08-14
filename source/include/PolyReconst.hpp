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
     *    @param n_ the number of parameters as an integer
     */
    PolyReconst (int n_);
    /**
     *    Calls the reconstruction algorithm
     *    @returns a vector with a pair of integers corrisponding to the
     *    coefficients of the polynom. The vector is ordered in an anscending
     *    way such that the first coefficient corresponds to to x^0,...
     */
    std::vector<RationalNumber> reconst();

  private:
    /**
     *    The reconstruction algorithm for polynimials with coefficients which
     *    are members of a finite field
     *    @param prime the corresponding prime number
     *    @return A vector which contains the coefficients in ascending order
     */
    std::vector<mpz_class> reconst_ff(const uint64_t prime);
    /**
     *    A numerical black box function which provides the reconstruction
     *    algorithm with the finite field member f(y)
     *    @param p the prime number which defines the finite field
     *    @param y a member of the finite field which should be evaluated in f(y)
     */
    FFInt num (uint64_t p, const FFInt &y) const;
    /**
     *    Computes the coefficient a_i recursevly using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param i the order of a_i
     *    @param ip recursion order
     *    @param num f(y_i)
     *    @returns a_i
     */
    const uint breakCondition = 3;
    FFInt comp_ai (const std::vector<FFInt> &ai, int i, int ip, const FFInt &num);
    Polynomial iterate_canonical (const std::vector<FFInt> &ai, uint i, const uint64_t prime) const;
    std::vector<FFInt> constr_canonical(const std::vector<FFInt> &ai, const uint64_t prime) const;
    std::vector<mpz_class> convert_to_mpz(const std::vector<FFInt>& ai) const;
    std::vector<FFInt> convert_to_ffint(const uint64_t prime) const;
    bool test_guess(const uint64_t prime);
    int n; /**< The number of parameters */
    std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
    mpz_class combined_prime {};
    std::vector<mpz_class> combined_ci {};
    std::vector<RationalNumber> guessi {};
  };
}
