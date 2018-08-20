#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"
#include "PolynomialFF.hpp"
#include "RationalNumber.hpp"
#include "RationalFunction.hpp"
#include "gmpxx.h"

namespace firefly {

  class RatReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(int n_);
    /**
     *    Calls the reconstruction algorithm
     *    @returns a RationalFunction object
     *    @throw runtimeerror if the prime numbers are not sufficient to reconstruct
     *    rational coefficients
     */
    //RationalFunction reconst();
  private:
    //std::pair<std::vector<mpz_class>, std::vector<mpz_class>> reconst_ff(const uint64_t prime);
    /**
     *    Computes the coefficient a_i recursivly using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param i the order of a_i
     *    @param ip recursion order
     *    @param num f(y_i)
     *    @returns a_i
     */
    //FFInt comp_ai(std::vector<FFInt> &ai, int i, int ip, const FFInt &num);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param ratFun the rational function
     *    @param prime a prime number defining the finite field
     *    @return A normalized version of ratFun
     */
    //std::pair<PolynomialFF, PolynomialFF> normalize(std::pair<PolynomialFF, PolynomialFF> &ratFun, const uint64_t prime) const;
    /**
     *    Constructs the canonical form of the rational function recursivly
     *    @param ai a vector of the coefficients ai as FFints
     *    @param prime a prime number defining the current finite field
     *    @return the rational function in its canonical form
     */
    //std::pair<PolynomialFF, PolynomialFF>  construct_canonical(std::vector<FFInt> &ai, const uint64_t prime) const;
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param ai a vector of the coefficients ai as FFInts
     *    @param i an integer telling the current degree of the rational function
     *    @param prime a prime number defining the current finite field
     *    @return the recursivly iterated rational function in its canonical form
     */
    //std::pair<PolynomialFF, PolynomialFF> iterate_canonical(std::vector<FFInt> &ai, uint i, const uint64_t prime) const;
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param ai a vector of coefficients ai
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @param prime a prime number defining the finite field
     *    @returns f(y_i)
     */
    //FFInt comp_fyi(std::vector<FFInt> &ai, uint i, uint ip, const FFInt &y, const uint64_t prime) const;
    /**
     *    A numerical black box function which provides the reconstruction
     *    algorithm with the finite field member f(y)
     *    @param prime the prime number which defines the finite field
     *    @param y a member of the finite field which should be evaluated in f(y)
     */
    //FFInt num(uint64_t prime, const FFInt &y) const;
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    //bool test_guess(const uint64_t prime);
    /**
     *    Converts the coefficients of a rational function from FFInts to mpz_class
     *    objects
     *    @param rf a rational function
     *    @return the coefficients of the given rational function converted to
     *    mpz_class objects
     */
    //std::pair<std::vector<mpz_class>, std::vector<mpz_class>> convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF> &rf) const;
    /**
     *    Converts the elements of a vector of RationalNumber objects to FFInts
     *    @param ri the vector of RationalNumber objects
     *    @param prime the prime number defining the corresponding finite field
     *    @return elements of ri converted to FFInts
     */
    //std::vector<FFInt> convert_to_ffint(const std::vector<RationalNumber> &ri, const uint64_t prime) const;
    int n; /**< The number of parameters */
    //mpz_class combined_prime {};  /**< The combination of the used prime numbers with the chinese remained theorem */
    //const uint breakCondition = 3;  /**< The number of additional evaluations of the black box function to validate the termination criterion */
    //std::vector<FFInt> yi {}; /**< A vector which holds all arguments y_i */
    //std::vector<RationalNumber> g_ni {}; /**< rational coefficient guesses for the numerator*/
    //std::vector<RationalNumber> g_di {}; /**< rational coefficient guesses for the denominator*/
    //std::vector<mpz_class> combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    //std::vector<mpz_class> combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
  };
}
