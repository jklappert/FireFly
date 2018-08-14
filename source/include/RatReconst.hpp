#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"
#include "Polynomial.hpp"
#include "RationalNumber.hpp"
#include "gmpxx.h"

namespace firefly {

  class RatReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst (int n_);
    /**
     *    Calls the reconstruction algorithm
     *    @returns a vector with a pair of integers corrisponding to the
     *    coefficient of the polynom. The vector is ordered in an anscending
     *    way such that the first coefficient corresponds to to z^0,...
     */
    std::pair<std::vector<RationalNumber>, std::vector<RationalNumber>> reconst();
  private:
    std::pair<std::vector<mpz_class>, std::vector<mpz_class>> reconst_ff(const uint64_t prime);
    /**
     *    Constructs the canonical form of the rational function recursevly
     */
    std::pair<Polynomial, Polynomial>  construct_canonical(std::vector<FFInt> &ai, const uint64_t prime) const;
    int n; /**< The number of parameters */
    /**
     *    Computes the coefficient a_i recursevly using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param i the order of a_i
     *    @param ip recursion order
     *    @param num f(y_i)
     *    @returns a_i
     */
    FFInt comp_ai (std::vector<FFInt> &ai, int i, int ip, const FFInt &num);
    std::pair<Polynomial, Polynomial> normalize (std::pair<Polynomial, Polynomial> &ratFun, const uint64_t prime) const;
    std::pair<Polynomial, Polynomial> iterate_canonical (std::vector<FFInt> &ai, uint i, const uint64_t prime) const;
    /**
     *    A numerical implementation of Thiele's interpolation formula
     *    from arXiv:1608.01902 eq. (3.10)
     *    @param i the order of the highest a_i
     *    @param ip recusrion order
     *    @param num the argument of f(z) which is a finite field member
     *    @returns a finite field member which corresponds to one recusion
     *    step
     */
    FFInt iterateCanonicalNum (std::vector<FFInt> &ai, uint i, uint ip, const FFInt &num, const uint64_t prime) const;
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param i order of the highest coefficient a_i
     *    @param y y_i
     *    @returns f(y_i)
     */
    FFInt compFyi (std::vector<FFInt> &ai, int i, const FFInt &y, const uint64_t prime) const;
    std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
    /**
     *    A numerical black box function which provides the reconstruction
     *    algorithm with the finite field member f(y)
     *    @param p the prime number which defines the finite field
     *    @param y a member of the finite field which should be evaluated in f(y)
     */
    FFInt num (uint64_t p, const FFInt &y) const;
    mpz_class combined_prime {};
    const uint breakCondition = 3;
    std::pair<std::vector<mpz_class>, std::vector<mpz_class>> convert_to_mpz(const std::pair<Polynomial, Polynomial> &p) const;
    std::vector<FFInt> convert_to_ffint(const uint64_t prime, const std::vector<RationalNumber> &guess) const;
    bool test_guess(const uint64_t prime);
    std::vector<RationalNumber> g_ni {};
    std::vector<RationalNumber> g_di {};
    std::vector<mpz_class> combined_ni {};
    std::vector<mpz_class> combined_di {};
  };
}
