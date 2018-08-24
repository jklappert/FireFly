#pragma once

#include <cstdint>
#include <vector>
#include <unordered_map>
#include <gmpxx.h>
#include "FFInt.hpp"
#include "Polynomial.hpp"
#include "PolynomialFF.hpp"
#include "RationalNumber.hpp"

namespace firefly {

  typedef std::unordered_map<std::vector<uint>, mpz_class, UintHasher> mpz_map;
  typedef std::unordered_map<std::vector<uint>, RationalNumber, UintHasher> rn_map;

  class PolyReconst {
  public:
    /**
     *    A constructor
     *    @param n_ The number of parameters as an integer
     */
    PolyReconst(uint n_, uint64_t prime);
    /**
     *    Calls the reconstruction algorithm
     *    @return the reconstructed Polynomial
     *    @throw runtimeerror if the prime numbers are not sufficient to reconstruct
     *    rational coefficients
     */
    void feed(uint64_t prime, const std::vector<FFInt> &yis, FFInt &num);
    bool done = false;
    bool new_prime = false;
    uint next_zi = 1;
    Polynomial get_result();
  private:
    /**
     *    Computes the coefficient a(i) = ai.at(i) recursively using eq. (3.11) of
     *    arXiv:1608.01902
     *    @param zi the integer i to a zi
     *    @param ai The vector of previously computed ai
     *    @param num f(y_i)
     *    @param i The order of a(i)
     *    @param ip Recursion order
     *    @return a(i)
     */
    PolynomialFF comp_ai(const uint zi, int i, int ip, const PolynomialFF &num, std::vector<PolynomialFF>& ai);
    /**
     *    Convert the reconstructed polynomial to the canonical form
     *    @param zi the integer i to a zi
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @return The vector of coefficients of the canonical form
     */
    PolynomialFF construct_canonical(const uint zi, const uint64_t prime, std::vector<PolynomialFF>& ai);
    /**
     *    Iterative construction of the canonical form
     *    @param zi the integer i to a z
     *    @param ai The computed ai
     *    @param prime The prime number of the finite field
     *    @param i The iteration step; stops at ai.size()
     *    @return One iteration step of the canonical polynomial
     */
    PolynomialFF iterate_canonical(const uint zi, const uint64_t prime, uint i, std::vector<PolynomialFF>& ai);
    /**
     *    Test if the guess yields the same answer for the function in the finite
     *    field of prime
     *    @param prime The prime number which defines the finite field
     *    @return true or false
     */
    bool test_guess(const uint64_t prime, const FFInt& num);
    /**
     *    Converts a vector of FFInts to a vector of mpz_class
     *    @param ai a vector of FFInts
     *    @return The vector ai converted to mpz_class objects
     */
    mpz_map convert_to_mpz(const PolynomialFF &poly) const;
    /**
     *    Convert a vector of RationalNumber objects to FFInts
     *    @param ri a vector of RationalNumber objects
     *    @param prime a prime defining the current finite field
     *    @return the vector ri converted to FFInt objects
     */
    ff_map convert_to_ffint(const rn_map &ri, const uint64_t prime) const;
    uint n; /**< The number of parameters */
    bool use_chinese_remainder = false;
    bool check = false;
    uint curr_zi = 1;
    mpz_class combined_prime; /**< The combination of the used prime numbers with the chinese remained theorem */
    mpz_map combined_ci; /**< The combination of the finite field results with the chinese remained theorem */
    rn_map gi {}; /**< The guesses of the rational coefficients */
    std::unordered_map<uint, std::vector<FFInt>> yis {};
    std::unordered_map<uint, std::vector<PolynomialFF>> ais{};
    std::unordered_map<uint, int> max_deg{};
  };
}
