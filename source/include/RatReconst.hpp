#pragma once

#include <cstdint>
#include <vector>
#include "FFInt.hpp"
#include "PolyReconst.hpp"
#include "PolynomialFF.hpp"
#include "RationalNumber.hpp"
#include "RationalFunction.hpp"
#include "gmpxx.h"

namespace firefly {
  typedef std::unordered_map<std::vector<uint>, mpz_class, UintHasher> mpz_map;
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;

  class RatReconst {
  public:
    /**
     *    A constructor
     *    @param n_ the number of parameters
     */
    RatReconst(uint n_);
    /**
     *
     */
    void feed(const FFInt& new_ti, const std::vector<FFInt>& yis, const FFInt& num);
    /**
     *
     */
    RationalFunction get_result();
    std::vector<FFInt> shift {};
    bool done = false;
    bool new_prime = false;
    uint zi = 1;
  private:
    FFInt comp_ai(int i, int ip, const FFInt& num);
    /**
     *    Normalize the rational function such that the first non-zero coefficient
     *    of the denominator is normalized to 1
     *    @param ratFun the rational function
     *    @param prime a prime number defining the finite field
     *    @return A normalized version of ratFun
     */
    void normalize();
    /**
     *    Constructs the canonical form of the rational function recursivly
     *    @param ai a vector of the coefficients ai as FFints
     *    @param prime a prime number defining the current finite field
     *    @return the rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF>  construct_canonical();
    /**
     *    Iterates Thiele's interpolation formula to get the canonical form
     *    of the rational function
     *    @param ai a vector of the coefficients ai as FFInts
     *    @param i an integer telling the current degree of the rational function
     *    @param prime a prime number defining the current finite field
     *    @return the recursivly iterated rational function in its canonical form
     */
    std::pair<PolynomialFF, PolynomialFF> iterate_canonical(uint i);
    /**
     *    Calculates f(y_i) using  Thiele's interpolation formula
     *    @param ai a vector of coefficients ai
     *    @param i order of the highest coefficient a_i
     *    @param ip order of sub coefficient a_ip
     *    @param y y_i
     *    @param prime a prime number defining the finite field
     *    @returns f(y_i)
     */
    FFInt comp_fyi(uint i, uint ip, const FFInt& y);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const FFInt& num);
    /**
     *    Converts the coefficients of a rational function from FFInts to mpz_class
     *    objects
     *    @param rf a rational function
     *    @return the coefficients of the given rational function converted to
     *    mpz_class objects
     */
    std::pair<mpz_map, mpz_map> convert_to_mpz(const std::pair<PolynomialFF, PolynomialFF>& rf) const;
    /**
     *    Converts the elements of a vector of RationalNumber objects to FFInts
     *    @param ri the vector of RationalNumber objects
     *    @param prime the prime number defining the corresponding finite field
     *    @return elements of ri converted to FFInts
     */
    ff_map convert_to_ffint(const rn_map& ri) const;
    /**
     *
     */
    void remove_shift();
    bool rec_rat_coef();
    uint n; /**< The number of parameters */
    bool check = false;
    bool use_chinese_remainder = false;
    bool poly_new_prime = false;
    bool shifted = false;
    uint curr_zi = 0;
    std::vector<FFInt> ai {};
    std::unordered_map<uint, PolyReconst> coef_n {};
    std::unordered_map<uint, PolyReconst> coef_d {};
    RationalFunction result;
    mpz_class combined_prime {};  /**< The combination of the used prime numbers with the chinese remained theorem */
    std::vector<FFInt> ti {}; /**< A vector which holds all arguments t_i */
    rn_map g_ni {}; /**< rational coefficient guesses for the numerator*/
    rn_map g_di {}; /**< rational coefficient guesses for the denominator*/
    mpz_map combined_ni {};  /**< The combination of the coefficients of the numerator over finite field with the chinese remained theorem */
    mpz_map combined_di {};  /**< The combination of the coefficients of the denominator over finite field with the chinese remained theorem */
  };
}
